"""Cell-count adjustments that let an exact MPI decomposition balance better.

SNaC gives every rank a box of whole cells, and blocks that touch must share
their tangential cell counts and partitions. With awkward block sizes, even the
best exact decomposition for a rank count can balance poorly; moving a few cell
counts by a small fraction often fixes that. This module finds the trade-off
between the largest relative change of a cell count and the balance it buys.

Per axis, blocks that must share a grid form a component with one cell count n
and one partition count d. A block's largest local box is the product of
ceil(n/d) over its axes (SNaC's initmpi gives remainder cells to the first
ranks), so the imbalance is N * max(largest box) / sum(block cells).

Changes are measured from each component's original count: the count before
earlier proposals changed it (``Block.origin_ng``). Applying proposals one
after another therefore never drifts beyond the tolerance.

For fixed partitions, only the largest n with a given local extent ceil(n/d)
can balance best: it keeps the largest box and adds cells. Each component thus
has few candidate counts, and a branch and bound over components finds the
best counts within a change cap for each of a few promising partition sets. A
second search then looks for the same balance, give or take a little, with
fewer changed counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, floor, prod
from time import monotonic
from typing import Any, Callable

import numpy as np

from .axis_format import EXPLICIT, MONOTONE, SYMMETRIC
from .grid import AXIS_NAMES, axis_grid_arrays
from .model import Block, Project
from .topology import FACE_INFO, FaceConnection, axis_components
from .validation import GRID_TOL, SPACING_JUMP_WARNING_RATIO, copy_project

# Relative tolerance when comparing imbalances.
_TIE = 1.0e-12
# Alternatives must cut the remaining excess load by this fraction ...
_EXCESS_CUT = 0.25
# ... and lower the imbalance by at least this much.
_MIN_GAIN = 0.01
_MAX_ALTERNATIVES = 4
# To change fewer counts, a proposal may give up this much of its imbalance
# (relative), but never more than this fraction of its remaining excess load.
_SIMPLER_SLACK = 0.005
_SIMPLER_EXCESS = 0.1
# Boundary cell widths may change as on a uniform grid, plus this slack.
_WIDTH_SLACK = 1.05
# Partition sets explored per rank count, best bound first.
MAX_PARTITION_SETS = 48


@dataclass(frozen=True)
class CountComponent:
    """Blocks that share one cell count along an axis."""

    axis: int
    block_ids: tuple[int, ...]
    current: int
    # The count changes are measured from, and the window it allows.
    origin: int
    low: int
    high: int
    reason: str | None = None

    @property
    def free(self) -> bool:
        return self.reason is None


@dataclass
class BalancePoint:
    """Cell counts per component and the partitions that balance them."""

    imbalance: float
    # Largest relative difference from an original count, and the number of
    # components whose count differs from the current one.
    change: float
    changed: int
    counts: tuple[int, ...]
    dims: dict[int, tuple[int, int, int]]


@dataclass
class _Partition:
    dims: dict[int, tuple[int, int, int]]
    splits: tuple[int, ...]


@dataclass
class BalanceSearch:
    """Cell-count alternatives found for one rank count."""

    points: list[BalancePoint] = field(default_factory=list)
    optimal: bool = True
    fixed_components: list[dict[str, Any]] = field(default_factory=list)
    components: list[CountComponent] = field(default_factory=list)


class _Balancer:
    def __init__(
        self,
        project: Project,
        connections: list[FaceConnection],
        target: int,
        *,
        deadline: float,
        node_limit: int,
    ) -> None:
        self.project = project
        self.by_id = {block.id: block for block in project.blocks}
        self.block_ids = [block.id for block in project.blocks]
        self.target = target
        self.tolerance = project.decomposition.ng_tolerance
        self.min_local = project.decomposition.min_local_cells
        self.max_aspect = project.decomposition.max_local_aspect
        self.deadline = deadline
        self.nodes_left = node_limit
        self.truncated = False
        self.components: list[CountComponent] = []
        self.component_of: dict[int, list[int]] = {block_id: [0, 0, 0] for block_id in self.block_ids}
        for axis in range(3):
            for group in axis_components(project.blocks, connections, axis):
                ids = tuple(sorted(group))
                for block_id in ids:
                    self.component_of[block_id][axis] = len(self.components)
                self.components.append(self._component(axis, ids))
        self._widths_cache: dict[tuple[int, int], dict[int, tuple[float, float]] | None] = {}
        self._classes_cache: dict[tuple[int, int, float | None], list[tuple[int, list[int]]]] = {}
        self._base_widths = {
            index: self._grid_widths(component, component.origin)
            for index, component in enumerate(self.components)
            if component.free
        }
        self.pairs: dict[int, list[tuple[int, int, int, int, int, int, float]]] = {}
        for connection in connections:
            self._add_pair(connection)
        self.partners = {
            index: {pair[0] for pair in pairs} | {pair[3] for pair in pairs} for index, pairs in self.pairs.items()
        }

    # Components -------------------------------------------------------------

    def _component(self, axis: int, ids: tuple[int, ...]) -> CountComponent:
        blocks = [self.by_id[block_id] for block_id in ids]
        counts = {block.ng[axis] for block in blocks}
        current = min(counts)
        origins = {block.origin_ng()[axis] for block in blocks}
        origin = origins.pop() if len(origins) == 1 else current
        tolerance = self.tolerance
        low = max(1, ceil(origin * (1.0 - tolerance) - 1.0e-9))
        high = max(low, floor(origin * (1.0 + tolerance) + 1.0e-9))
        reason = None
        if axis == self.project.fft_axis_index:
            reason = "FFT axis"
        elif len(counts) > 1:
            reason = "blocks disagree on the cell count"
        elif any(_count_follows_grading(block, axis) for block in blocks):
            reason = "the grading sets the cell count"
        elif any(block.axis_locks[axis] for block in blocks):
            reason = "locked"
        elif low == high == current:
            reason = f"±{tolerance * 100:.3g}% allows no other count"
        if reason is not None:
            origin = low = high = current
        return CountComponent(axis, ids, current, origin, low, high, reason)

    def fixed_components(self) -> list[dict[str, Any]]:
        return [
            {
                "axis": AXIS_NAMES[component.axis],
                "blockIds": list(component.block_ids),
                "ng": component.current,
                "reason": component.reason,
            }
            for component in self.components
            if not component.free
        ]

    def windows(self) -> dict[int, tuple[tuple[int, int], ...]]:
        """The lowest and highest count per block axis, current counts included."""

        return {
            block_id: tuple(
                (
                    min(self.components[index].low, self.components[index].current),
                    max(self.components[index].high, self.components[index].current),
                )
                for index in self.component_of[block_id]
            )
            for block_id in self.block_ids
        }

    def caps(self) -> list[float]:
        values = {0.0}
        for component in self.components:
            if component.free:
                values.update(
                    abs(n - component.origin) / component.origin for n in range(component.low, component.high + 1)
                )
        return sorted(values)

    # Grid checks ------------------------------------------------------------

    def _grid_widths(self, component: CountComponent, n: int) -> dict[int, tuple[float, float]] | None:
        """First and last cell widths per block, or None if the grids break."""

        axis = component.axis
        widths: dict[int, tuple[float, float]] = {}
        reference: np.ndarray | None = None
        try:
            for block_id in component.block_ids:
                block = self.by_id[block_id]
                arrays = axis_grid_arrays(
                    block.axes[AXIS_NAMES[axis]].to_dict(), block.lmin[axis], block.lmax[axis], n
                )
                spacing = arrays.face_spacing[1:-1]
                if spacing.size != n or not np.all(np.isfinite(arrays.faces)) or np.any(spacing <= 0.0):
                    return None
                # Blocks of a component share their grid lines along the axis.
                if reference is not None:
                    scale = block.lmax[axis] - block.lmin[axis]
                    if reference.shape != arrays.faces.shape or np.max(np.abs(reference - arrays.faces)) > max(
                        GRID_TOL * scale, 1.0e-12 * scale
                    ):
                        return None
                reference = arrays.faces
                widths[block_id] = (float(spacing[0]), float(spacing[-1]))
        except (ValueError, ArithmeticError, FloatingPointError):
            return None
        return widths

    def widths(self, index: int, n: int) -> dict[int, tuple[float, float]] | None:
        key = (index, n)
        if key not in self._widths_cache:
            component = self.components[index]
            widths = self._grid_widths(component, n)
            base = self._base_widths.get(index)
            if widths is not None and base is not None and n not in (component.origin, component.current):
                # Boundary cells may scale as on a uniform grid; a grading whose
                # end cells move much more than that changes character.
                lowest = 1.0 / ((1.0 + self.tolerance) * _WIDTH_SLACK)
                highest = _WIDTH_SLACK / (1.0 - min(self.tolerance, 0.9))
                for block_id, pair in widths.items():
                    if not all(lowest <= new / old <= highest for new, old in zip(pair, base[block_id])):
                        widths = None
                        break
            self._widths_cache[key] = widths
        return self._widths_cache[key]

    def _add_pair(self, connection: FaceConnection) -> None:
        a_axis, a_side = FACE_INFO[connection.a_face]
        b_axis, b_side = FACE_INFO[connection.b_face]
        a_index = self.component_of[connection.a_id][a_axis]
        b_index = self.component_of[connection.b_id][b_axis]
        if not (self.components[a_index].free or self.components[b_index].free):
            return
        a_widths = self.widths(a_index, self.components[a_index].current)
        b_widths = self.widths(b_index, self.components[b_index].current)
        if a_widths is None or b_widths is None:
            return
        a_width = a_widths[connection.a_id][a_side]
        b_width = b_widths[connection.b_id][b_side]
        # Never add a spacing-jump warning, nor worsen an existing jump.
        limit = max(SPACING_JUMP_WARNING_RATIO, max(a_width, b_width) / min(a_width, b_width)) * (1.0 + 1.0e-9)
        pair = (a_index, connection.a_id, a_side, b_index, connection.b_id, b_side, limit)
        self.pairs.setdefault(a_index, []).append(pair)
        if b_index != a_index:
            self.pairs.setdefault(b_index, []).append(pair)

    def pairs_ok(self, index: int, counts: list[int | None]) -> bool:
        for a_index, a_id, a_side, b_index, b_id, b_side, limit in self.pairs.get(index, ()):
            if counts[a_index] is None or counts[b_index] is None:
                continue
            a_widths = self.widths(a_index, counts[a_index])
            b_widths = self.widths(b_index, counts[b_index])
            if a_widths is None or b_widths is None:
                return False
            a_width = a_widths[a_id][a_side]
            b_width = b_widths[b_id][b_side]
            if max(a_width, b_width) / min(a_width, b_width) > limit:
                return False
        return True

    # Partitions -------------------------------------------------------------

    def partition(self, dims: dict[int, tuple[int, int, int]]) -> _Partition | None:
        splits = [0] * len(self.components)
        for index, component in enumerate(self.components):
            values = {dims[block_id][component.axis] for block_id in component.block_ids}
            if len(values) != 1:
                return None
            splits[index] = values.pop()
            if not component.free and splits[index] > 1 and component.current // splits[index] < self.min_local:
                return None
        return _Partition(dims, tuple(splits))

    def classes(self, index: int, split: int, cap: float | None) -> list[tuple[int, list[int]]]:
        """Candidate counts by local extent, largest count first in each.

        *cap* bounds the relative difference from the original count; None
        allows the current count only.
        """

        component = self.components[index]
        if not component.free:
            return [(ceil(component.current / split), [component.current])]
        key = (index, split, cap)
        if key not in self._classes_cache:
            if cap is None:
                counts = [component.current]
            else:
                counts = [
                    n
                    for n in range(component.high, component.low - 1, -1)
                    if abs(n - component.origin) <= cap * component.origin * (1.0 + 1.0e-12)
                ]
            by_extent: dict[int, list[int]] = {}
            for n in counts:
                if split == 1 or n // split >= self.min_local:
                    by_extent.setdefault(-(-n // split), []).append(n)
            self._classes_cache[key] = sorted(by_extent.items())
        return self._classes_cache[key]

    # Search -----------------------------------------------------------------

    def solve(
        self,
        partitions: list[_Partition],
        cap: float | None,
        *,
        limit: float | None = None,
        fewest: bool = False,
    ) -> BalancePoint | None:
        """The best counts with every change at most *cap*.

        Best means the lowest imbalance, then the smallest change, then the
        fewest changed counts. With *fewest*, the fewest changed counts come
        first. *limit* rejects imbalances above it.
        """

        prepared = []
        for partition in partitions:
            classes = [self.classes(index, partition.splits[index], cap) for index in range(len(self.components))]
            if any(not options for options in classes):
                continue
            # Smallest extents and largest counts bound the imbalance from below.
            largest_box = max(
                prod(classes[index][0][0] for index in self.component_of[block_id]) for block_id in self.block_ids
            )
            most_cells = sum(
                prod(classes[index][-1][1][0] for index in self.component_of[block_id]) for block_id in self.block_ids
            )
            bound = self.target * largest_box / most_cells
            if limit is None or bound <= limit * (1.0 + _TIE):
                prepared.append((bound, partition, classes))
        prepared.sort(key=lambda item: item[0])
        best: BalancePoint | None = None
        for bound, partition, classes in prepared:
            if not fewest and best is not None and bound > best.imbalance * (1.0 + _TIE):
                break
            if self._out_of_budget():
                break
            found = self._search(partition, classes, best, limit=limit, fewest=fewest)
            if found is not None:
                best = found
        return best

    def _out_of_budget(self) -> bool:
        if self.nodes_left <= 0 or monotonic() >= self.deadline:
            self.truncated = True
            return True
        return False

    def _search(
        self,
        partition: _Partition,
        classes: list[list[tuple[int, list[int]]]],
        incumbent: BalancePoint | None,
        *,
        limit: float | None,
        fewest: bool,
    ) -> BalancePoint | None:
        components = self.components
        target = self.target
        free = [index for index, component in enumerate(components) if component.free]
        weight = {
            index: sum(prod(self.by_id[block_id].ng) for block_id in components[index].block_ids) for index in free
        }
        free.sort(key=lambda index: (-weight[index], index))
        # Optimistic values for unassigned components: smallest extent, largest count.
        extent = [options[0][0] for options in classes]
        count = [options[-1][1][0] for options in classes]
        keeps_current = {
            index: any(components[index].current in values for _, values in classes[index]) for index in free
        }
        assigned: list[int | None] = [
            None if components[index].free else components[index].current for index in range(len(components))
        ]
        box = {block_id: prod(extent[index] for index in self.component_of[block_id]) for block_id in self.block_ids}
        cells = {block_id: prod(count[index] for index in self.component_of[block_id]) for block_id in self.block_ids}
        total = sum(cells.values())
        best = incumbent
        found: BalancePoint | None = None

        def better(imbalance: float, change: float, changed: int) -> bool:
            if best is None:
                return True
            if fewest and changed != best.changed:
                return changed < best.changed
            if abs(imbalance - best.imbalance) > _TIE * best.imbalance:
                return imbalance < best.imbalance
            if abs(change - best.change) > _TIE:
                return change < best.change
            return changed < best.changed

        def hopeless(bound: float, change: float, changed: int) -> bool:
            """Whether no completion can beat the best solution so far."""

            if limit is not None and bound > limit * (1.0 + _TIE):
                return True
            if best is None:
                return False
            if fewest and changed != best.changed:
                return changed > best.changed
            if bound > best.imbalance * (1.0 + _TIE):
                return True
            return bound >= best.imbalance * (1.0 - _TIE) and change > best.change + _TIE

        def leaf(change: float, changed: int) -> None:
            nonlocal best, found
            imbalance = target * max(box.values()) / total
            if limit is not None and imbalance > limit * (1.0 + _TIE):
                return
            if not better(imbalance, change, changed) or not self._aspect_ok(partition, assigned):
                return
            best = found = BalancePoint(
                imbalance=imbalance,
                change=change,
                changed=changed,
                counts=tuple(int(value) for value in assigned),
                dims=partition.dims,
            )

        def visit(position: int, change: float, changed: int) -> bool:
            nonlocal total
            if self.nodes_left <= 0 or monotonic() >= self.deadline:
                self.truncated = True
                return False
            self.nodes_left -= 1
            if position == len(free):
                leaf(change, changed)
                return True
            index = free[position]
            component = components[index]
            # The largest count of an extent wins unless a constraint still to
            # be checked might reject it while accepting a smaller count.
            dominated = not self.max_aspect and all(
                assigned[partner] is not None for partner in self.partners.get(index, ()) if partner != index
            )
            members = component.block_ids
            member_set = set(members)
            others = max((value for block_id, value in box.items() if block_id not in member_set), default=0)
            old_extent, old_count = extent[index], count[index]

            def bound_for(local_extent: int, n: int) -> float:
                trial_total = total + sum(cells[block_id] // old_count * (n - old_count) for block_id in members)
                trial_max = max(others, max(box[block_id] // old_extent * local_extent for block_id in members))
                return target * trial_max / trial_total

            options = []
            for local_extent, values in classes[index]:
                # When fewer changes matter, keeping the current count is its own option.
                values = [n for n in values if not (fewest and n == component.current)]
                if values:
                    step = abs(values[0] - component.origin) / component.origin
                    options.append((bound_for(local_extent, values[0]), step, local_extent, values))
            options.sort(key=lambda item: (item[0], item[1]))
            unsorted = 0
            if fewest and keeps_current[index]:
                local_extent = -(-component.current // partition.splits[index])
                options.insert(0, (bound_for(local_extent, component.current), 0.0, local_extent, [component.current]))
                unsorted = 1
            for order, (bound, _, local_extent, values) in enumerate(options):
                if order >= unsorted and hopeless(bound, 0.0, changed):
                    # Later options have larger bounds.
                    break
                for n in values:
                    step = abs(n - component.origin) / component.origin
                    is_change = int(n != component.current)
                    if hopeless(bound, max(change, step), changed + is_change) or self.widths(index, n) is None:
                        continue
                    assigned[index] = n
                    if not self.pairs_ok(index, assigned):
                        assigned[index] = None
                        continue
                    saved = {block_id: (box[block_id], cells[block_id]) for block_id in members}
                    for block_id in members:
                        box[block_id] = box[block_id] // old_extent * local_extent
                        cells[block_id] = cells[block_id] // old_count * n
                    previous_total = total
                    total += sum(cells[block_id] - saved[block_id][1] for block_id in members)
                    extent[index], count[index] = local_extent, n
                    keep_going = visit(position + 1, max(change, step), changed + is_change)
                    extent[index], count[index] = old_extent, old_count
                    total = previous_total
                    for block_id, (old_box, old_cells) in saved.items():
                        box[block_id], cells[block_id] = old_box, old_cells
                    assigned[index] = None
                    if not keep_going:
                        return False
                    if dominated:
                        # Smaller counts with the same extent only lose cells.
                        break
            return True

        visit(0, 0.0, 0)
        return found

    def _aspect_ok(self, partition: _Partition, counts: list[int | None]) -> bool:
        if not self.max_aspect:
            return True
        for block_id in self.block_ids:
            lengths = [counts[index] / partition.splits[index] for index in self.component_of[block_id]]
            if max(lengths) / min(lengths) > self.max_aspect * (1.0 + 1.0e-12):
                return False
        return True


def _count_follows_grading(block: Block, axis: int) -> bool:
    spec = block.axes[AXIS_NAMES[axis]]
    if spec.kind in {EXPLICIT, SYMMETRIC}:
        return True
    return spec.kind == MONOTONE and "n" not in spec.controls


def _simpler_limit(imbalance: float) -> float:
    """The imbalance a proposal may accept in exchange for fewer changed counts."""

    return min(imbalance * (1.0 + _SIMPLER_SLACK), imbalance + _SIMPLER_EXCESS * (imbalance - 1.0))


def search_balance_alternatives(
    project: Project,
    connections: list[FaceConnection],
    target: int,
    reference_imbalance: float,
    collect_partitions: Callable[..., tuple[list[tuple[float, dict[int, tuple[int, int, int]]]], bool]],
    *,
    deadline: float,
    node_limit: int = 400_000,
) -> BalanceSearch:
    """Points on the trade-off between cell-count change and load balance.

    *reference_imbalance* is the balance reached without changes.
    *collect_partitions(windows, incumbent, deadline)* returns promising
    partition sets as (bound, dims per block) and whether it had to leave some
    out. Walking by increasing change, a point is kept when it cuts the
    remaining excess load by at least a quarter; the best point overall is
    always kept. Each point then changes as few counts as it can while giving
    up at most a little of its balance (``_simpler_limit``).
    """

    search = BalanceSearch()
    if project.decomposition.ng_tolerance <= 0.0 or reference_imbalance <= 1.0 + 1.0e-9:
        return search
    balancer = _Balancer(project, connections, target, deadline=deadline, node_limit=node_limit)
    search.components = balancer.components
    search.fixed_components = balancer.fixed_components()
    if not any(component.free for component in balancer.components):
        return search

    collected, dropped = collect_partitions(balancer.windows(), reference_imbalance, deadline)
    partitions = [
        partition
        for partition in (balancer.partition(dims) for _, dims in collected[:MAX_PARTITION_SETS])
        if partition is not None
    ]
    search.optimal = not dropped and len(collected) <= MAX_PARTITION_SETS
    if not partitions:
        search.optimal = search.optimal and not balancer.truncated
        return search

    caps = balancer.caps()
    cache: dict[float | None, BalancePoint | None] = {}

    def solve(cap: float | None) -> BalancePoint | None:
        if cap not in cache:
            cache[cap] = balancer.solve(partitions, cap)
        return cache[cap]

    def smallest_cap(level: float) -> tuple[float, BalancePoint] | None:
        """The smallest cap whose best imbalance reaches *level*, with that solution."""

        low, high = 0, len(caps) - 1
        top = solve(caps[high])
        if top is None or top.imbalance > level * (1.0 + _TIE):
            return None
        while low < high:
            middle = (low + high) // 2
            found = solve(caps[middle])
            if found is not None and found.imbalance <= level * (1.0 + _TIE):
                high = middle
            else:
                low = middle + 1
        return caps[high], solve(caps[high])

    def simplest(cap: float, point: BalancePoint, limit: float) -> BalancePoint:
        simpler = balancer.solve(partitions, cap, limit=limit, fewest=True)
        return simpler if simpler is not None and simpler.changed < point.changed else point

    current = solve(None)
    start = min(reference_imbalance, current.imbalance if current else reference_imbalance)
    best = solve(caps[-1])
    points: list[BalancePoint] = []
    reached_best = False
    if best is not None and best.imbalance <= start - _MIN_GAIN:
        level = start
        while len(points) < _MAX_ALTERNATIVES:
            goal = min(1.0 + (1.0 - _EXCESS_CUT) * (level - 1.0), level - _MIN_GAIN)
            found = smallest_cap(goal)
            if found is None or found[1].changed == 0:
                break
            cap, point = found
            reached_best = reached_best or point.imbalance <= best.imbalance * (1.0 + _TIE)
            point = simplest(cap, point, min(goal, _simpler_limit(point.imbalance)))
            points.append(point)
            level = point.imbalance
            if level <= 1.0 + _TIE:
                break
        if not reached_best:
            cap, overall = smallest_cap(best.imbalance) or (caps[-1], best)
            overall = simplest(cap, overall, _simpler_limit(overall.imbalance))
            points = points[: _MAX_ALTERNATIVES - 1] + [overall]
    search.points = [point for point in points if point.changed]
    search.optimal = search.optimal and not balancer.truncated
    return search


def apply_balance_alternative(project: Project | dict[str, Any], alternative: dict[str, Any]) -> Project:
    """Set the cell counts and partitions of one balance alternative.

    Each block remembers its counts from before the first proposal it took
    (``ngOriginal``) and the counts set now (``ngBalanced``), so later
    proposals stay within the tolerance of the original counts.
    """

    project = copy_project(project)
    by_id = {block.id: block for block in project.blocks}
    entries = alternative.get("blocks") if isinstance(alternative, dict) else None
    if not isinstance(entries, list) or not entries:
        raise ValueError("balance alternative must list the blocks it changes")
    origins: dict[int, list[int]] = {}
    for entry in entries:
        block_id = int(entry.get("blockId", 0))
        if block_id not in by_id:
            raise ValueError(f"balance alternative refers to missing block {block_id}")
        ng = [int(value) for value in entry.get("ng", [])]
        dims = [int(value) for value in entry.get("dims", [])]
        if len(ng) != 3 or len(dims) != 3 or min(ng + dims) < 1:
            raise ValueError(f"balance alternative for block {block_id} needs three positive ng and dims")
        origins[block_id] = by_id[block_id].origin_ng()
        by_id[block_id].ng = ng
        by_id[block_id].dims = dims
    project.validate()
    for entry in entries:
        block = by_id[int(entry["blockId"])]
        if block.ng != [int(value) for value in entry["ng"]]:
            raise ValueError(f"block {block.id}: its grading does not allow the proposed cell counts")
        origin = origins[block.id]
        if block.ng == origin:
            block.ng_original = block.ng_balanced = None
        else:
            block.ng_original, block.ng_balanced = origin, list(block.ng)
    return project
