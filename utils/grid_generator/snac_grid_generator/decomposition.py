"""Topology-aware MPI decomposition for SNaC grid projects."""

from __future__ import annotations

from dataclasses import dataclass, field
import heapq
from itertools import combinations
from math import ceil, floor, gcd, log, prod
from time import monotonic
from typing import Any

from .balance import MAX_PARTITION_SETS, BalancePoint, BalanceSearch, search_balance_alternatives
from .diagnostics import error, warning
from .grid import AXIS_NAMES
from .model import DEFAULT_MIN_LOCAL_CELLS, Block, Project
from .topology import FaceConnection, Topology, build_topology
from .validation import (
    CheckResult,
    check_project,
    connected_component_errors,
    copy_project,
)

MIN_LOCAL_CELLS = DEFAULT_MIN_LOCAL_CELLS  # public alias of the project default
_SHAPE_WEIGHT = 0.08
_COMMUNICATION_WEIGHT = 0.02
_SPREAD_WEIGHT = 0.02
# Above this load imbalance, nearby rank totals that balance clearly better are suggested.
_SUGGEST_BETTER_ABOVE = 1.05
# Relative slack when comparing imbalance bounds, against rounding.
_BOUND_SLACK = 1.0e-12


@dataclass(frozen=True)
class BlockDecomposition:
    """Achieved MPI decomposition metrics for one block."""

    block_id: int
    dims: tuple[int, int, int]
    ranks: int
    min_cells: int
    max_cells: int
    average_cells: float
    min_local_extent: int
    max_local_aspect: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "blockId": self.block_id,
            "dims": list(self.dims),
            "ranks": self.ranks,
            "minCells": self.min_cells,
            "maxCells": self.max_cells,
            "averageCells": self.average_cells,
            "minLocalExtent": self.min_local_extent,
            "maxLocalAspect": self.max_local_aspect,
        }


@dataclass
class DecompositionResult(CheckResult):
    """Result and diagnostics from exact-rank decomposition optimization."""

    target_ranks: int = 0
    total_ranks: int = 0
    active_axes: tuple[str, ...] = ()
    min_cells: int = 0
    max_cells: int = 0
    average_cells: float = 0.0
    imbalance: float = 0.0
    communication: int = 0
    min_local_extent: int = 0
    max_local_aspect: float = 0.0
    quality_score: float = 0.0
    optimal: bool = True
    searched_nodes: int = 0
    nearby_ranks: list[int] = field(default_factory=list)
    better_ranks: list[dict[str, float | int]] = field(default_factory=list)
    blocks: list[BlockDecomposition] = field(default_factory=list)
    # Cell-count changes that balance better, by increasing change; see balance.py.
    ng_tolerance: float = 0.0
    alternatives: list[dict[str, Any]] = field(default_factory=list)
    alternatives_optimal: bool = True
    fixed_components: list[dict[str, Any]] = field(default_factory=list)
    # Counts that earlier proposals changed; the tolerance is measured from these.
    original_counts: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            **super().to_dict(),
            "targetRanks": self.target_ranks,
            "totalRanks": self.total_ranks,
            "activeAxes": list(self.active_axes),
            "minCells": self.min_cells,
            "maxCells": self.max_cells,
            "averageCells": self.average_cells,
            "imbalance": self.imbalance,
            "communication": self.communication,
            "minLocalExtent": self.min_local_extent,
            "maxLocalAspect": self.max_local_aspect,
            "qualityScore": self.quality_score,
            "optimal": self.optimal,
            "searchedNodes": self.searched_nodes,
            "nearbyRanks": self.nearby_ranks,
            "betterRanks": self.better_ranks,
            "blocks": [block.to_dict() for block in self.blocks],
            "ngTolerance": self.ng_tolerance,
            "alternatives": self.alternatives,
            "alternativesOptimal": self.alternatives_optimal,
            "fixedComponents": self.fixed_components,
            "originalCounts": self.original_counts,
        }


@dataclass(frozen=True)
class _DecompositionCandidate:
    dims: tuple[int, int, int]
    ranks: int
    min_cells: int
    max_cells: int
    communication: int
    min_local_extent: int
    max_local_aspect: float


@dataclass
class _SearchOutcome:
    candidates: dict[int, _DecompositionCandidate] | None = None
    score: tuple[float, float, float, int, float] | None = None
    nodes: int = 0
    truncated: bool = False
    # Collect mode: (imbalance bound, dims per block), best bound first, and
    # whether assignments that might have balanced well were left out.
    collected: list[tuple[float, dict[int, tuple[int, int, int]]]] = field(default_factory=list)
    dropped: bool = False


def optimize_project_decomposition(
    project: Project | dict[str, Any],
    *,
    node_limit: int = 500_000,
    time_limit: float = 4.0,
    alternatives: bool = False,
    alternatives_time_limit: float = 4.0,
) -> tuple[Project, DecompositionResult]:
    """Choose an exact-rank, topology-compatible MPI decomposition.

    The returned project keeps its cell counts. With *alternatives*, the result
    also lists cell-count changes within the project's ``ng_tolerance`` that
    balance better; ``apply_balance_alternative`` applies one.
    """

    project = copy_project(project)
    target = project.decomposition.target_ranks
    min_local_cells = project.decomposition.min_local_cells
    max_local_aspect = project.decomposition.max_local_aspect
    result = DecompositionResult(target_ranks=target)
    topology = build_topology(project.blocks, project.periodic_axes)
    result.add(*topology.issues)
    result.add(*connected_component_errors(project.blocks, topology.connections))
    if not result.ok:
        return project, result
    if not project.blocks:
        result.add(error("decomposition", "cannot decompose a project with no blocks"))
        return project, result
    if target < len(project.blocks):
        result.add(
            error("decomposition", f"target MPI ranks must be at least the number of blocks ({len(project.blocks)})")
        )
        return project, result

    active_sets = _active_axis_sets(project, target)
    if not active_sets:
        result.add(error("decomposition", "the decomposition mode and allowed axes do not define a valid search"))
        return project, result

    deadline = monotonic() + max(0.05, time_limit)
    best: _SearchOutcome | None = None
    best_axes: tuple[int, ...] = ()
    nodes_remaining = max(1, node_limit)
    any_truncated = False
    for set_index, active_axes in enumerate(active_sets):
        sets_left = len(active_sets) - set_index
        allowance = max(1, nodes_remaining // sets_left)
        outcome = _search_decomposition(
            project.blocks,
            topology.connections,
            target,
            active_axes,
            min_local_cells=min_local_cells,
            max_local_aspect=max_local_aspect,
            node_limit=allowance,
            deadline=deadline,
        )
        result.searched_nodes += outcome.nodes
        nodes_remaining = max(0, nodes_remaining - outcome.nodes)
        any_truncated = any_truncated or outcome.truncated
        if outcome.candidates is not None and (best is None or outcome.score < best.score):
            best = outcome
            best_axes = active_axes
        if set_index < len(active_sets) - 1 and (monotonic() >= deadline or nodes_remaining <= 0):
            any_truncated = True
            break

    if best is None or best.candidates is None:
        result.optimal = not any_truncated
        if any_truncated:
            result.add(error("decomposition", "decomposition search limit reached before finding an exact solution"))
        else:
            result.add(
                error("decomposition", f"no exact {project.decomposition.mode} decomposition uses {target} MPI ranks")
            )
            result.nearby_ranks = _nearby_decomposition_counts(
                project,
                topology.connections,
                target,
                deadline=monotonic() + min(1.5, max(0.2, time_limit / 2.0)),
            )
        return project, result

    result.optimal = not any_truncated
    if any_truncated:
        result.add(
            warning(
                "decomposition",
                "decomposition search was bounded; the result is valid but may not be globally optimal",
            )
        )
    by_id = {block.id: block for block in project.blocks}
    for block_id, candidate in best.candidates.items():
        by_id[block_id].dims = list(candidate.dims)
    project.validate()
    _fill_decomposition_result(result, project.blocks, topology.connections, best.candidates, best_axes)
    if result.imbalance > _SUGGEST_BETTER_ABOVE:
        result.better_ranks = _better_balanced_counts(
            project,
            topology.connections,
            target,
            result.imbalance,
            deadline=monotonic() + min(1.5, max(0.2, time_limit / 2.0)),
        )
    result.ng_tolerance = project.decomposition.ng_tolerance
    if alternatives:
        _add_balance_alternatives(project, topology, result, deadline=monotonic() + max(0.05, alternatives_time_limit))
    return project, result


def describe_project_decomposition(project: Project | dict[str, Any]) -> DecompositionResult:
    """Report the balance of the project's current block dims."""

    project = copy_project(project)
    target = project.decomposition.target_ranks
    result = DecompositionResult(target_ranks=target, ng_tolerance=project.decomposition.ng_tolerance)
    topology = build_topology(project.blocks, project.periodic_axes)
    result.add(*topology.issues)
    if not project.blocks or not result.ok:
        result.add(error("decomposition", "cannot describe the decomposition of this project"))
        return result
    candidates = {block.id: _decomposition_candidate(block, tuple(block.dims)) for block in project.blocks}
    split_axes = tuple(axis for axis in range(3) if any(block.dims[axis] > 1 for block in project.blocks))
    result.target_ranks = target or sum(candidate.ranks for candidate in candidates.values())
    _fill_decomposition_result(result, project.blocks, topology.connections, candidates, split_axes)
    if target and result.total_ranks != target:
        result.add(error("decomposition", f"decomposition uses {result.total_ranks} MPI ranks, expected {target}"))
    return result


def _add_balance_alternatives(
    project: Project,
    topology: Topology,
    result: DecompositionResult,
    *,
    deadline: float,
) -> None:
    target = result.target_ranks
    decomposition = project.decomposition

    def collect(
        windows: dict[int, tuple[tuple[int, int], ...]],
        incumbent: float,
        until: float,
    ) -> tuple[list[tuple[float, dict[int, tuple[int, int, int]]]], bool]:
        collected: list[tuple[float, dict[int, tuple[int, int, int]]]] = []
        dropped = False
        sets = _active_axis_sets(project, target)
        for index, active_axes in enumerate(sets):
            outcome = _search_decomposition(
                project.blocks,
                topology.connections,
                target,
                active_axes,
                min_local_cells=decomposition.min_local_cells,
                max_local_aspect=decomposition.max_local_aspect,
                node_limit=200_000,
                deadline=monotonic() + max(0.02, (until - monotonic()) / (2 * (len(sets) - index))),
                collect=MAX_PARTITION_SETS + 1,
                ng_windows=windows,
                incumbent=incumbent,
            )
            collected.extend(outcome.collected)
            dropped = dropped or outcome.dropped or outcome.truncated
        collected.sort(key=lambda item: item[0])
        return collected, dropped

    search = search_balance_alternatives(
        project,
        topology.connections,
        target,
        result.imbalance,
        collect,
        deadline=deadline,
    )
    result.fixed_components = search.fixed_components
    result.alternatives_optimal = search.optimal
    result.original_counts = [
        {
            "axis": AXIS_NAMES[component.axis],
            "blockIds": list(component.block_ids),
            "original": component.origin,
            "current": component.current,
        }
        for component in search.components
        if component.origin != component.current
    ]
    if not search.points:
        return
    base_errors = set(check_project(project, topology=topology).errors)
    seen: set[str] = set()
    for point in search.points:
        alternative = _balance_alternative(project, topology, search, point, base_errors, deadline)
        if alternative is None or alternative["imbalance"] > result.imbalance - 1.0e-9:
            continue
        key = repr(alternative["blocks"])
        if key not in seen:
            seen.add(key)
            result.alternatives.append(alternative)
    result.alternatives.sort(key=lambda item: (item["maxChange"], item["imbalance"]))


def _balance_alternative(
    project: Project,
    topology: Topology,
    search: BalanceSearch,
    point: BalancePoint,
    base_errors: set[str],
    deadline: float,
) -> dict[str, Any] | None:
    """Describe one balance point as the project Balance would then produce."""

    candidate = copy_project(project)
    by_id = {block.id: block for block in candidate.blocks}
    changes = []
    for component, count in zip(search.components, point.counts):
        for block_id in component.block_ids:
            by_id[block_id].ng[component.axis] = count
        if count != component.current:
            change = {
                "axis": AXIS_NAMES[component.axis],
                "blockIds": list(component.block_ids),
                "from": component.current,
                "to": count,
            }
            if component.origin != component.current:
                change["original"] = component.origin
            changes.append(change)
    proposed = {block.id: list(block.ng) for block in candidate.blocks}
    try:
        candidate.validate()
    except ValueError:
        return None
    if any(block.ng != proposed[block.id] for block in candidate.blocks):
        return None

    # Prefer what Balance itself picks for the new counts, unless it balances worse.
    target = candidate.decomposition.target_ranks
    chosen: dict[int, _DecompositionCandidate] | None = None
    chosen_axes: tuple[int, ...] = ()
    best_score = None
    for active_axes in _active_axis_sets(candidate, target):
        outcome = _search_decomposition(
            candidate.blocks,
            topology.connections,
            target,
            active_axes,
            min_local_cells=candidate.decomposition.min_local_cells,
            max_local_aspect=candidate.decomposition.max_local_aspect,
            node_limit=50_000,
            deadline=min(deadline, monotonic() + 0.5),
        )
        if outcome.candidates is not None and (best_score is None or outcome.score < best_score):
            best_score, chosen, chosen_axes = outcome.score, outcome.candidates, active_axes
    if chosen is None or best_score[1] > point.imbalance * (1.0 + 1.0e-9):
        chosen = {block_id: _decomposition_candidate(by_id[block_id], dims) for block_id, dims in point.dims.items()}
        chosen_axes = tuple(axis for axis in range(3) if any(dims[axis] > 1 for dims in point.dims.values()))
        if any(
            not _candidate_meets_constraints(
                by_id[block_id],
                item,
                candidate.decomposition.min_local_cells,
                candidate.decomposition.max_local_aspect,
            )
            for block_id, item in chosen.items()
        ):
            return None
    for block_id, item in chosen.items():
        by_id[block_id].dims = list(item.dims)

    check = check_project(candidate, topology=topology)
    if not set(check.errors) <= base_errors:
        return None
    metrics = DecompositionResult(target_ranks=target)
    _fill_decomposition_result(metrics, candidate.blocks, topology.connections, chosen, chosen_axes)
    old_cells = sum(prod(block.ng) for block in project.blocks)
    new_cells = sum(prod(block.ng) for block in candidate.blocks)
    return {
        "imbalance": metrics.imbalance,
        # Measured from the original counts, as the tolerance is.
        "maxChange": point.change,
        "exact": all(block.ng[axis] % block.dims[axis] == 0 for block in candidate.blocks for axis in range(3)),
        "totalCells": new_cells,
        "cellChange": new_cells / old_cells - 1.0,
        "maxInterfaceRatio": max((item["ratio"] for item in check.interfaces), default=1.0),
        "maxLocalAspect": metrics.max_local_aspect,
        "minLocalExtent": metrics.min_local_extent,
        "activeAxes": list(metrics.active_axes),
        "changes": changes,
        "blocks": [
            {"blockId": block.id, "ng": list(block.ng), "dims": list(block.dims)}
            for block in sorted(candidate.blocks, key=lambda item: item.id)
        ],
    }


def _active_axis_sets(project: Project, target: int) -> list[tuple[int, ...]]:
    if target == len(project.blocks):
        return [()]
    allowed = tuple(
        index
        for index, enabled in enumerate(project.decomposition.axes)
        if enabled and index != project.fft_axis_index
    )
    mode = project.decomposition.mode
    sizes = range(len(allowed), 0, -1) if mode == "auto" else (int(mode[0]),)
    return [axes for size in sizes for axes in combinations(allowed, size)]


def _search_decomposition(
    blocks: list[Block],
    connections: list[FaceConnection],
    target: int,
    active_axes: tuple[int, ...],
    *,
    min_local_cells: int,
    max_local_aspect: float,
    node_limit: int,
    deadline: float,
    first_only: bool = False,
    collect: int = 0,
    ng_windows: dict[int, tuple[tuple[int, int], ...]] | None = None,
    incumbent: float = float("inf"),
) -> _SearchOutcome:
    """Branch and bound over per-block dims for exactly *target* ranks.

    By default this finds the best assignment for the current cell counts.

    With *collect* and *ng_windows* (the lowest and highest cell count each
    block axis may take), it instead gathers up to *collect* assignments that
    could balance best once cell counts move within their windows, ranked by a
    lower bound on that imbalance. Assignments whose bound exceeds an imbalance
    already reachable without changes (*incumbent*, or a better one found on
    the way) cannot win and are pruned. The local aspect limit is left to the
    caller, since it depends on the final cell counts.
    """

    outcome = _SearchOutcome()
    nblocks = len(blocks)
    collecting = collect > 0 and ng_windows is not None
    heap: list[tuple[float, int, dict[int, tuple[int, int, int]]]] = []
    reachable = incumbent

    def limit(block: Block, axis: int) -> int:
        return ng_windows[block.id][axis][1] if collecting else block.ng[axis]

    def load(block_id: int, candidate: _DecompositionCandidate) -> float:
        if not collecting:
            return candidate.max_cells
        # The smallest largest-local-box this split can reach inside the windows.
        return prod(
            max(-(-low // dims), min_local_cells) if dims > 1 else low
            for (low, _), dims in zip(ng_windows[block_id], candidate.dims)
        )

    if not active_axes and not collecting:
        if target == nblocks:
            candidates = {block.id: _decomposition_candidate(block, (1, 1, 1)) for block in blocks}
            if not all(
                _candidate_meets_constraints(block, candidates[block.id], min_local_cells, max_local_aspect)
                for block in blocks
            ):
                return outcome
            outcome.candidates = candidates
            outcome.score = _decomposition_score(blocks, connections, outcome.candidates, target)
        return outcome

    capacity = sum(
        prod(max(1, limit(block, axis) // min_local_cells) for axis in active_axes)
        for block in blocks
    )
    if target > capacity:
        return outcome

    max_block_ranks = target - nblocks + 1
    total_cells = sum(prod(block.ng) for block in blocks)
    mean_cells = total_cells / target
    load_scale = (
        sum(prod(high for _, high in ng_windows[block.id]) for block in blocks) / target if collecting else mean_cells
    )
    candidates: dict[int, list[_DecompositionCandidate]] = {}
    for block in blocks:
        ideal = target * prod(block.ng) / total_cells
        block_candidates = _block_decomposition_candidates(
            block,
            active_axes,
            max_block_ranks,
            min_local_cells=min_local_cells,
            max_local_aspect=max_local_aspect,
            ng_limit=tuple(high for _, high in ng_windows[block.id]) if collecting else None,
        )
        block_candidates.sort(key=lambda candidate: _candidate_priority(block, candidate, ideal, mean_cells))
        candidates[block.id] = block_candidates

    neighbor_rules: dict[int, list[tuple[int, int]]] = {block.id: [] for block in blocks}
    for connection in connections:
        for axis_index in range(3):
            if axis_index == connection.axis_index:
                continue
            neighbor_rules[connection.a_id].append((connection.b_id, axis_index))
            neighbor_rules[connection.b_id].append((connection.a_id, axis_index))
    by_id = {block.id: block for block in blocks}
    order = sorted(
        by_id,
        key=lambda block_id: (-len(neighbor_rules[block_id]), -prod(by_id[block_id].ng), block_id),
    )
    assigned: dict[int, _DecompositionCandidate] = {}
    stop = False

    def compatible_with_assigned(block_id: int, candidate: _DecompositionCandidate) -> bool:
        return not any(
            neighbor_id in assigned and candidate.dims[axis] != assigned[neighbor_id].dims[axis]
            for neighbor_id, axis in neighbor_rules[block_id]
        )

    def remaining_rank_bounds(position: int) -> tuple[int, int, int] | None:
        minimum = 0
        maximum = 0
        rank_step = 0
        for future_id in order[position:]:
            fixed: dict[int, int] = {}
            for neighbor_id, axis in neighbor_rules[future_id]:
                if neighbor_id not in assigned:
                    continue
                value = assigned[neighbor_id].dims[axis]
                if axis in fixed and fixed[axis] != value:
                    return None
                fixed[axis] = value
            minimum_ranks = prod(fixed.values())
            maximum_ranks = minimum_ranks * prod(
                max(1, limit(by_id[future_id], axis) // min_local_cells)
                for axis in active_axes
                if axis not in fixed
            )
            minimum += minimum_ranks
            maximum += min(max_block_ranks, maximum_ranks)
            rank_step = gcd(rank_step, minimum_ranks)
        return minimum, maximum, rank_step

    def search(position: int, ranks_used: int, current_max_cells: float) -> None:
        nonlocal stop, reachable
        if stop:
            return
        if outcome.nodes >= node_limit or monotonic() >= deadline:
            outcome.truncated = True
            stop = True
            return
        outcome.nodes += 1
        bounds = remaining_rank_bounds(position)
        if bounds is None or ranks_used + bounds[0] > target or ranks_used + bounds[1] < target:
            return
        if bounds[2] and (target - ranks_used - bounds[0]) % bounds[2]:
            return
        if position == len(order):
            if ranks_used != target:
                return
            if any(not any(candidate.dims[axis] > 1 for candidate in assigned.values()) for axis in active_axes):
                return
            if collecting:
                if all(
                    _candidate_meets_constraints(by_id[block_id], candidate, min_local_cells, max_local_aspect)
                    for block_id, candidate in assigned.items()
                ):
                    score = _decomposition_score(blocks, connections, assigned, target)
                    reachable = min(reachable, score[1])
                    if outcome.score is None or score < outcome.score:
                        outcome.score = score
                        outcome.candidates = dict(assigned)
                bound = current_max_cells / load_scale
                if bound <= reachable * (1.0 + _BOUND_SLACK):
                    entry = (-bound, -outcome.nodes, {block_id: c.dims for block_id, c in assigned.items()})
                    if len(heap) < collect:
                        heapq.heappush(heap, entry)
                    else:
                        outcome.dropped = True
                        if bound < -heap[0][0]:
                            heapq.heapreplace(heap, entry)
                return
            score = _decomposition_score(blocks, connections, assigned, target)
            if outcome.score is None or score < outcome.score:
                outcome.score = score
                outcome.candidates = dict(assigned)
            if first_only:
                stop = True
            return
        if collecting:
            bound = current_max_cells / load_scale
            if bound > reachable * (1.0 + _BOUND_SLACK):
                return
            if len(heap) >= collect and bound >= -heap[0][0]:
                outcome.dropped = True
                return
        elif outcome.score is not None and current_max_cells / mean_cells > outcome.score[0]:
            return

        block_id = order[position]
        blocks_after = len(order) - position - 1
        max_ranks = target - ranks_used - blocks_after
        for candidate in candidates[block_id]:
            if candidate.ranks > max_ranks or not compatible_with_assigned(block_id, candidate):
                continue
            assigned[block_id] = candidate
            search(
                position + 1,
                ranks_used + candidate.ranks,
                max(current_max_cells, load(block_id, candidate)),
            )
            assigned.pop(block_id, None)
            if stop and (first_only or outcome.truncated):
                return

    search(0, 0, 0)
    kept = sorted(
        (-negative, -order, dims) for negative, order, dims in heap if -negative <= reachable * (1.0 + _BOUND_SLACK)
    )
    outcome.collected = [(bound, dims) for bound, _, dims in kept]
    return outcome


def _block_decomposition_candidates(
    block: Block,
    active_axes: tuple[int, ...],
    max_ranks: int,
    *,
    min_local_cells: int,
    max_local_aspect: float,
    ng_limit: tuple[int, int, int] | None = None,
) -> list[_DecompositionCandidate]:
    dims = [1, 1, 1]
    result: list[_DecompositionCandidate] = []
    cells = block.ng if ng_limit is None else ng_limit

    def admissible(candidate: _DecompositionCandidate) -> bool:
        if ng_limit is None:
            return _candidate_meets_constraints(block, candidate, min_local_cells, max_local_aspect)
        # With adjustable cell counts, only the largest allowed counts matter here.
        return all(
            cells[axis] // candidate.dims[axis] >= min_local_cells for axis in range(3) if candidate.dims[axis] > 1
        )

    def visit(position: int, ranks: int) -> None:
        if position == len(active_axes):
            candidate = _decomposition_candidate(block, tuple(dims))
            if admissible(candidate):
                result.append(candidate)
            return
        axis = active_axes[position]
        partition_limit = max(1, cells[axis] // min_local_cells)
        upper = min(partition_limit, max_ranks // ranks)
        for value in range(1, upper + 1):
            dims[axis] = value
            visit(position + 1, ranks * value)
        dims[axis] = 1

    visit(0, 1)
    return result


def _candidate_meets_constraints(
    block: Block,
    candidate: _DecompositionCandidate,
    min_local_cells: int,
    max_local_aspect: float,
) -> bool:
    split_extents = [
        block.ng[index] // candidate.dims[index]
        for index in range(3)
        if candidate.dims[index] > 1
    ]
    enough_cells = not split_extents or min(split_extents) >= min_local_cells
    return enough_cells and (
        max_local_aspect == 0.0 or candidate.max_local_aspect <= max_local_aspect
    )


def _decomposition_candidate(block: Block, dims: tuple[int, int, int]) -> _DecompositionCandidate:
    minimum_extents = tuple(floor(block.ng[index] / dims[index]) for index in range(3))
    maximum_extents = tuple(ceil(block.ng[index] / dims[index]) for index in range(3))
    local_lengths = tuple(block.ng[index] / dims[index] for index in range(3))
    communication = (
        (dims[0] - 1) * block.ng[1] * block.ng[2]
        + (dims[1] - 1) * block.ng[0] * block.ng[2]
        + (dims[2] - 1) * block.ng[0] * block.ng[1]
    )
    return _DecompositionCandidate(
        dims=dims,
        ranks=prod(dims),
        min_cells=prod(minimum_extents),
        max_cells=prod(maximum_extents),
        communication=communication,
        min_local_extent=min(minimum_extents),
        max_local_aspect=max(local_lengths) / min(local_lengths),
    )


def _candidate_priority(
    block: Block,
    candidate: _DecompositionCandidate,
    ideal_ranks: float,
    mean_cells: float,
) -> tuple[float, float, int, tuple[int, int, int]]:
    rank_error = abs(log(candidate.ranks / ideal_ranks))
    shape_penalty = _SHAPE_WEIGHT * log(max(1.0, candidate.max_local_aspect))
    communication_penalty = _COMMUNICATION_WEIGHT * candidate.communication / prod(block.ng)
    load_penalty = abs(candidate.max_cells / mean_cells - 1.0)
    return (
        load_penalty + rank_error + shape_penalty + communication_penalty,
        candidate.max_local_aspect,
        candidate.communication,
        candidate.dims,
    )


def _decomposition_score(
    blocks: list[Block],
    connections: list[FaceConnection],
    candidates: dict[int, _DecompositionCandidate],
    target: int,
) -> tuple[float, float, float, int, float]:
    total_cells = sum(prod(block.ng) for block in blocks)
    mean_cells = total_cells / target
    maximum = max(candidate.max_cells for candidate in candidates.values())
    minimum = min(candidate.min_cells for candidate in candidates.values())
    imbalance = maximum / mean_cells
    spread = maximum / minimum
    aspect = max(candidate.max_local_aspect for candidate in candidates.values())
    communication = _communication_cost(connections, candidates)
    quality = (
        imbalance
        + _SHAPE_WEIGHT * log(max(1.0, aspect))
        + _COMMUNICATION_WEIGHT * communication / total_cells
        + _SPREAD_WEIGHT * log(max(1.0, spread))
    )
    return quality, imbalance, aspect, communication, spread


def _communication_cost(
    connections: list[FaceConnection],
    candidates: dict[int, _DecompositionCandidate],
) -> int:
    communication = sum(candidate.communication for candidate in candidates.values())
    for connection in connections:
        dims = candidates[connection.a_id].dims
        communication += prod(dims[axis] for axis in range(3) if axis != connection.axis_index)
    return communication


def _fill_decomposition_result(
    result: DecompositionResult,
    blocks: list[Block],
    connections: list[FaceConnection],
    candidates: dict[int, _DecompositionCandidate],
    active_axes: tuple[int, ...],
) -> None:
    total_cells = sum(prod(block.ng) for block in blocks)
    score = _decomposition_score(blocks, connections, candidates, result.target_ranks)
    result.total_ranks = sum(candidate.ranks for candidate in candidates.values())
    result.active_axes = tuple(AXIS_NAMES[index] for index in active_axes)
    result.min_cells = min(candidate.min_cells for candidate in candidates.values())
    result.max_cells = max(candidate.max_cells for candidate in candidates.values())
    result.average_cells = total_cells / result.total_ranks
    result.imbalance = score[1]
    result.max_local_aspect = score[2]
    result.communication = score[3]
    result.min_local_extent = min(candidate.min_local_extent for candidate in candidates.values())
    result.quality_score = score[0]
    by_id = {block.id: block for block in blocks}
    result.blocks = [
        BlockDecomposition(
            block_id=block_id,
            dims=candidate.dims,
            ranks=candidate.ranks,
            min_cells=candidate.min_cells,
            max_cells=candidate.max_cells,
            average_cells=prod(by_id[block_id].ng) / candidate.ranks,
            min_local_extent=candidate.min_local_extent,
            max_local_aspect=candidate.max_local_aspect,
        )
        for block_id, candidate in sorted(candidates.items())
    ]


def _better_balanced_counts(
    project: Project,
    connections: list[FaceConnection],
    target: int,
    imbalance: float,
    *,
    deadline: float,
) -> list[dict[str, float | int]]:
    """Nearby rank totals whose best decomposition at least halves the excess load.

    The interface constraints can make an exact rank count balance poorly while
    a total one or two ranks away balances well; these are the suggestions.
    """

    excess = imbalance - 1.0
    found: list[dict[str, float | int]] = []
    maximum_offset = min(16, max(4, target // 10))
    for offset in range(1, maximum_offset + 1):
        for trial in (target - offset, target + offset):
            if trial < len(project.blocks) or monotonic() >= deadline:
                continue
            best: float | None = None
            for active_axes in _active_axis_sets(project, trial):
                outcome = _search_decomposition(
                    project.blocks,
                    connections,
                    trial,
                    active_axes,
                    min_local_cells=project.decomposition.min_local_cells,
                    max_local_aspect=project.decomposition.max_local_aspect,
                    node_limit=20_000,
                    deadline=deadline,
                )
                if outcome.score is not None and (best is None or outcome.score[1] < best):
                    best = outcome.score[1]
            if best is not None and best - 1.0 <= 0.5 * excess and imbalance - best >= 0.05:
                found.append({"ranks": trial, "imbalance": best})
        if len(found) >= 3:
            break
    return sorted(found, key=lambda item: (abs(int(item["ranks"]) - target), int(item["ranks"])))[:3]


def _nearby_decomposition_counts(
    project: Project,
    connections: list[FaceConnection],
    target: int,
    *,
    deadline: float,
) -> list[int]:
    nearby: list[int] = []
    maximum_offset = min(16, max(4, target // 10))
    for offset in range(1, maximum_offset + 1):
        for trial in (target - offset, target + offset):
            if trial < len(project.blocks) or trial in nearby:
                continue
            active_sets = _active_axis_sets(project, trial)
            found = False
            for active_axes in active_sets:
                outcome = _search_decomposition(
                    project.blocks,
                    connections,
                    trial,
                    active_axes,
                    min_local_cells=project.decomposition.min_local_cells,
                    max_local_aspect=project.decomposition.max_local_aspect,
                    node_limit=40_000,
                    deadline=deadline,
                    first_only=True,
                )
                if outcome.candidates is not None:
                    nearby.append(trial)
                    found = True
                    break
                if monotonic() >= deadline:
                    return sorted(nearby, key=lambda value: (abs(value - target), value))[:4]
            if found and len(nearby) >= 4:
                return sorted(nearby, key=lambda value: (abs(value - target), value))[:4]
    return sorted(nearby, key=lambda value: (abs(value - target), value))[:4]
