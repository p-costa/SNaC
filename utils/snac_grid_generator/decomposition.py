"""Topology-aware MPI decomposition for SNaC grid projects."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from math import ceil, floor, gcd, log, prod
from time import monotonic
from typing import Any

from .grid import AXIS_NAMES
from .model import Block, Project
from .validation import (
    CheckResult,
    FaceConnection,
    _build_topology,
    _check_connected_components,
    _project_copy,
)

MIN_LOCAL_CELLS = 4
_SHAPE_WEIGHT = 0.08
_COMMUNICATION_WEIGHT = 0.02
_SPREAD_WEIGHT = 0.02


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
    blocks: list[BlockDecomposition] = field(default_factory=list)

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
            "blocks": [block.to_dict() for block in self.blocks],
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


def optimize_project_decomposition(
    project: Project | dict[str, Any],
    *,
    node_limit: int = 500_000,
    time_limit: float = 4.0,
) -> tuple[Project, DecompositionResult]:
    """Choose an exact-rank, topology-compatible MPI decomposition."""

    project = _project_copy(project)
    target = project.decomposition.target_ranks
    result = DecompositionResult(target_ranks=target)
    topology = _build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    result.errors.extend(_check_connected_components(project.blocks, topology.connections))
    if result.errors:
        return project, result
    if not project.blocks:
        result.errors.append("cannot decompose a project with no blocks")
        return project, result
    if target < len(project.blocks):
        result.errors.append(
            f"target MPI ranks must be at least the number of blocks ({len(project.blocks)})"
        )
        return project, result

    active_sets = _active_axis_sets(project, target)
    if not active_sets:
        result.errors.append("the decomposition mode and allowed axes do not define a valid search")
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
            result.errors.append("decomposition search limit reached before finding an exact solution")
        else:
            result.errors.append(f"no exact {project.decomposition.mode} decomposition uses {target} MPI ranks")
            result.nearby_ranks = _nearby_decomposition_counts(
                project,
                topology.connections,
                target,
                deadline=monotonic() + min(1.5, max(0.2, time_limit / 2.0)),
            )
        return project, result

    result.optimal = not any_truncated
    if any_truncated:
        result.warnings.append("decomposition search was bounded; the result is valid but may not be globally optimal")
    by_id = {block.id: block for block in project.blocks}
    for block_id, candidate in best.candidates.items():
        by_id[block_id].dims = list(candidate.dims)
    project.validate()
    _fill_decomposition_result(result, project.blocks, topology.connections, best.candidates, best_axes)
    return project, result


def _active_axis_sets(project: Project, target: int) -> list[tuple[int, ...]]:
    if target == len(project.blocks):
        return [()]
    allowed = tuple(index for index, enabled in enumerate(project.decomposition.axes) if enabled)
    mode = project.decomposition.mode
    sizes = range(len(allowed), 0, -1) if mode == "auto" else (int(mode[0]),)
    return [axes for size in sizes for axes in combinations(allowed, size)]


def _search_decomposition(
    blocks: list[Block],
    connections: list[FaceConnection],
    target: int,
    active_axes: tuple[int, ...],
    *,
    node_limit: int,
    deadline: float,
    first_only: bool = False,
) -> _SearchOutcome:
    outcome = _SearchOutcome()
    nblocks = len(blocks)
    if not active_axes:
        if target == nblocks:
            outcome.candidates = {
                block.id: _decomposition_candidate(block, (1, 1, 1)) for block in blocks
            }
            outcome.score = _decomposition_score(blocks, connections, outcome.candidates, target)
        return outcome

    capacity = sum(
        prod(max(1, block.ng[axis] // MIN_LOCAL_CELLS) for axis in active_axes)
        for block in blocks
    )
    if target > capacity:
        return outcome

    max_block_ranks = target - nblocks + 1
    total_cells = sum(prod(block.ng) for block in blocks)
    mean_cells = total_cells / target
    candidates: dict[int, list[_DecompositionCandidate]] = {}
    for block in blocks:
        ideal = target * prod(block.ng) / total_cells
        block_candidates = _block_decomposition_candidates(block, active_axes, max_block_ranks)
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
                max(1, by_id[future_id].ng[axis] // MIN_LOCAL_CELLS)
                for axis in active_axes
                if axis not in fixed
            )
            minimum += minimum_ranks
            maximum += min(max_block_ranks, maximum_ranks)
            rank_step = gcd(rank_step, minimum_ranks)
        return minimum, maximum, rank_step

    def search(position: int, ranks_used: int, current_max_cells: int) -> None:
        nonlocal stop
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
            score = _decomposition_score(blocks, connections, assigned, target)
            if outcome.score is None or score < outcome.score:
                outcome.score = score
                outcome.candidates = dict(assigned)
            if first_only:
                stop = True
            return
        if outcome.score is not None and current_max_cells / mean_cells > outcome.score[0]:
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
                max(current_max_cells, candidate.max_cells),
            )
            assigned.pop(block_id, None)
            if stop and (first_only or outcome.truncated):
                return

    search(0, 0, 0)
    return outcome


def _block_decomposition_candidates(
    block: Block,
    active_axes: tuple[int, ...],
    max_ranks: int,
) -> list[_DecompositionCandidate]:
    dims = [1, 1, 1]
    result: list[_DecompositionCandidate] = []

    def visit(position: int, ranks: int) -> None:
        if position == len(active_axes):
            result.append(_decomposition_candidate(block, tuple(dims)))
            return
        axis = active_axes[position]
        partition_limit = max(1, block.ng[axis] // MIN_LOCAL_CELLS)
        upper = min(partition_limit, max_ranks // ranks)
        for value in range(1, upper + 1):
            dims[axis] = value
            visit(position + 1, ranks * value)
        dims[axis] = 1

    visit(0, 1)
    return result


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
