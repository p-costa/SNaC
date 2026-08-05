"""Grid-congruence and interface-spacing repair helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isclose, log, sqrt
from typing import Any

import numpy as np

from .grid import AXIS_NAMES, axis_grid_arrays
from .model import AxisSpec, Block, Project
from .numeric import coordinates_close, geometry_tolerance
from .topology import (
    FACE_INFO,
    FACE_ORDER,
    FaceConnection,
    axis_extent,
    build_topology,
)
from .validation import (
    GRID_TOL,
    SPACING_JUMP_WARNING_RATIO,
    CheckResult,
    _boundary_width,
    _interior_faces,
    _project_copy,
)


@dataclass(frozen=True)
class GridRepair:
    """One proposed axis-grid replacement."""

    block_id: int
    axis: str
    source_block_id: int
    old_n: int
    new_n: int
    kind: str = "congruence"
    face: str | None = None
    old_spacing: float | None = None
    new_spacing: float | None = None

    def to_dict(self) -> dict[str, Any]:
        data = {
            "blockId": self.block_id,
            "axis": self.axis,
            "sourceBlockId": self.source_block_id,
            "oldN": self.old_n,
            "newN": self.new_n,
        }
        if self.kind != "congruence":
            data.update(
                {
                    "kind": self.kind,
                    "face": self.face,
                    "oldSpacing": self.old_spacing,
                    "newSpacing": self.new_spacing,
                }
            )
        return data


@dataclass
class RepairResult(CheckResult):
    """Validation messages and replacements proposed by grid repair."""

    changes: list[GridRepair] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {**super().to_dict(), "changes": [change.to_dict() for change in self.changes]}


def repair_project_grids(
    project: Project | dict[str, Any],
    source_block_id: int | None = None,
) -> tuple[Project, RepairResult]:
    """Propose congruent grids and matched normal interface spacing."""

    original = _project_copy(project)
    project = _project_copy(original)
    result = RepairResult()
    topology = build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    if result.errors:
        return original, result

    by_id = {block.id: block for block in project.blocks}
    selected = source_block_id if source_block_id in by_id else None
    for axis_index, axis in enumerate(AXIS_NAMES):
        for component in _axis_components(project.blocks, topology.connections, axis_index):
            if len(component) < 2:
                continue
            locked = sorted(block_id for block_id in component if by_id[block_id].axis_locks[axis_index])
            if len(locked) > 1:
                reference = by_id[locked[0]]
                conflicts = [
                    block_id
                    for block_id in locked[1:]
                    if not _same_axis_grid(reference, by_id[block_id], axis_index)
                ]
                if conflicts:
                    result.errors.append(
                        f"locked {axis} grids conflict between blocks {locked[0]} and {conflicts}"
                    )
                    continue
            owner_id = locked[0] if locked else selected if selected in component else min(component)
            owner = by_id[owner_id]
            for block_id in sorted(component):
                if block_id == owner_id:
                    continue
                block = by_id[block_id]
                if block.axis_locks[axis_index] and not _same_axis_grid(owner, block, axis_index):
                    result.errors.append(
                        f"block {block_id}: locked {axis} grid conflicts with authoritative block {owner_id}"
                    )
                    continue
                if _same_axis_grid(owner, block, axis_index):
                    continue
                old_n = block.ng[axis_index]
                block.axes[axis] = AxisSpec.from_dict(owner.axes[axis].to_dict())
                block.ng[axis_index] = owner.ng[axis_index]
                result.changes.append(GridRepair(block_id, axis, owner_id, old_n, block.ng[axis_index]))

    if result.errors:
        result.changes.clear()
        return original, result
    _repair_interface_spacing(project, topology.connections, selected, result)
    if result.errors:
        result.changes.clear()
        return original, result
    project.validate()
    return project, result


def _repair_interface_spacing(
    project: Project,
    connections: list[FaceConnection],
    selected: int | None,
    result: RepairResult,
) -> None:
    by_id = {block.id: block for block in project.blocks}
    constraints: dict[tuple[int, int, int], tuple[float, int]] = {}

    def constrain(block_id: int, axis: int, side: int, spacing: float, source_id: int) -> None:
        key = (block_id, axis, side)
        previous = constraints.get(key)
        if previous is not None and not isclose(previous[0], spacing, rel_tol=GRID_TOL, abs_tol=0.0):
            result.errors.append(
                f"block {block_id}: conflicting {AXIS_NAMES[axis]} spacing targets on {FACE_ORDER[2 * axis + side]}"
            )
        else:
            constraints[key] = (spacing, source_id)

    for connection in connections:
        a = by_id[connection.a_id]
        b = by_id[connection.b_id]
        axis = connection.axis_index
        a_width = _boundary_width(a, axis, connection.a_face)
        b_width = _boundary_width(b, axis, connection.b_face)
        ratio = max(a_width, b_width) / min(a_width, b_width)
        if ratio <= SPACING_JUMP_WARNING_RATIO:
            continue
        _, a_side = FACE_INFO[connection.a_face]
        _, b_side = FACE_INFO[connection.b_face]
        a_locked = a.axis_locks[axis]
        b_locked = b.axis_locks[axis]
        if a.id == b.id:
            if a_locked:
                result.errors.append(
                    f"block {a.id}: locked periodic {AXIS_NAMES[axis]} grid has an interface spacing jump of {ratio:.3g}"
                )
                continue
            target = sqrt(a_width * b_width)
            constrain(a.id, axis, a_side, target, a.id)
            constrain(a.id, axis, b_side, target, a.id)
            continue
        if a_locked and b_locked:
            result.errors.append(
                f"blocks {a.id} and {b.id}: locked {AXIS_NAMES[axis]} grids have an interface spacing jump of {ratio:.3g}"
            )
            continue
        if a_locked:
            constrain(b.id, axis, b_side, a_width, a.id)
        elif b_locked:
            constrain(a.id, axis, a_side, b_width, b.id)
        elif selected == a.id or (selected != b.id and a.id < b.id):
            constrain(b.id, axis, b_side, a_width, a.id)
        else:
            constrain(a.id, axis, a_side, b_width, b.id)

    grouped: dict[tuple[int, int], dict[int, tuple[float, int]]] = {}
    for (block_id, axis, side), target in constraints.items():
        grouped.setdefault((block_id, axis), {})[side] = target
    if result.errors:
        return

    spacing_changes = False
    for (block_id, axis_index), targets in sorted(grouped.items()):
        block = by_id[block_id]
        axis = AXIS_NAMES[axis_index]
        arrays = axis_grid_arrays(
            block.axes[axis].to_dict(),
            block.lmin[axis_index],
            block.lmax[axis_index],
            block.ng[axis_index],
        )
        widths = arrays.face_spacing[1:-1]
        old_endpoints = (float(widths[0]), float(widths[-1]))
        lower_target = targets.get(0, (None, 0))[0]
        upper_target = targets.get(1, (None, 0))[0]
        try:
            matched_axis, repaired = _matched_axis_grid(
                block,
                axis_index,
                lower_target,
                upper_target,
            )
        except ValueError as exc:
            result.errors.append(f"block {block_id}: cannot match {axis} interface spacing ({exc})")
            continue
        block.axes[axis] = matched_axis
        if matched_axis.kind == "explicit":
            result.warnings.append(f"block {block_id} {axis}: converted to explicit coordinates")
        for side, (target, source_id) in sorted(targets.items()):
            result.changes.append(
                GridRepair(
                    block_id=block_id,
                    axis=axis,
                    source_block_id=source_id,
                    old_n=block.ng[axis_index],
                    new_n=block.ng[axis_index],
                    kind="spacing",
                    face=FACE_ORDER[2 * axis_index + side],
                    old_spacing=old_endpoints[side],
                    new_spacing=float(repaired[0 if side == 0 else -1]),
                )
            )
        spacing_changes = True
    if spacing_changes and not project.write_external_grid:
        project.write_external_grid = True
        result.warnings.append("external grid writing was enabled for matched interface spacing")


def _matched_axis_grid(
    block: Block,
    axis_index: int,
    target_lower: float | None,
    target_upper: float | None,
) -> tuple[AxisSpec, np.ndarray]:
    axis_name = AXIS_NAMES[axis_index]
    current = block.axes[axis_name]
    if current.kind == "simple_ratio" and (target_lower is None) != (target_upper is None):
        candidate = AxisSpec.from_dict(current.to_dict())
        control = "width_start" if target_lower is not None else "width_end"
        target = float(target_lower if target_lower is not None else target_upper)
        candidate.controls = ["n", control]
        setattr(candidate, control, target)
        try:
            arrays = axis_grid_arrays(
                candidate.to_dict(),
                block.lmin[axis_index],
                block.lmax[axis_index],
                block.ng[axis_index],
            )
            widths = arrays.face_spacing[1:-1]
            achieved = float(widths[0 if target_lower is not None else -1])
            if isclose(achieved, target, rel_tol=1.0e-8, abs_tol=1.0e-12):
                return candidate, widths
        except ValueError:
            pass

    arrays = axis_grid_arrays(
        current.to_dict(),
        block.lmin[axis_index],
        block.lmax[axis_index],
        block.ng[axis_index],
    )
    repaired = _smooth_endpoint_match(
        arrays.face_spacing[1:-1],
        float(block.lmax[axis_index] - block.lmin[axis_index]),
        target_lower,
        target_upper,
    )
    faces = np.concatenate(
        ([block.lmin[axis_index]], block.lmin[axis_index] + np.cumsum(repaired))
    )
    faces[-1] = block.lmax[axis_index]
    return AxisSpec(kind="explicit", faces=faces.tolist()), repaired


def _smooth_endpoint_match(
    widths: np.ndarray,
    length: float,
    target_lower: float | None,
    target_upper: float | None,
) -> np.ndarray:
    widths = np.asarray(widths, dtype=float)
    if widths.size == 1:
        target = target_lower if target_lower is not None else target_upper
        if target is not None and not isclose(target, length, rel_tol=GRID_TOL, abs_tol=0.0):
            raise ValueError("a one-cell grid cannot match a different interface spacing")
        return np.asarray([length])

    lower = float(widths[0] if target_lower is None else target_lower)
    upper = float(widths[-1] if target_upper is None else target_upper)
    if lower <= 0.0 or upper <= 0.0 or lower + upper >= length:
        raise ValueError("interface spacing targets leave no positive interior grid")
    if widths.size == 2:
        if target_lower is not None and target_upper is not None and not isclose(
            lower + upper, length, rel_tol=GRID_TOL, abs_tol=0.0
        ):
            raise ValueError("two-cell grid interface targets do not span the block")
        if target_lower is None:
            lower = length - upper
        elif target_upper is None:
            upper = length - lower
        return np.asarray([lower, upper])

    position = np.linspace(0.0, 1.0, widths.size)
    endpoint_scale = np.exp(
        log(lower / widths[0]) * (1.0 - position) + log(upper / widths[-1]) * position
    )
    base = widths * endpoint_scale
    shape = 4.0 * position * (1.0 - position)

    def total(curvature: float) -> float:
        return float(np.sum(base * np.exp(curvature * shape)))

    low = -64.0
    high = 64.0
    if total(low) > length or total(high) < length:
        raise ValueError("could not construct a positive grid for the interface spacing targets")
    for _ in range(100):
        middle = 0.5 * (low + high)
        if total(middle) < length:
            low = middle
        else:
            high = middle
    repaired = base * np.exp(0.5 * (low + high) * shape)
    repaired *= length / float(repaired.sum())
    repaired[0] = lower
    repaired[-1] = upper
    repaired[1:-1] *= (length - lower - upper) / float(repaired[1:-1].sum())
    return repaired


def _axis_components(
    blocks: list[Block],
    connections: list[FaceConnection],
    axis_index: int,
) -> list[set[int]]:
    graph: dict[int, set[int]] = {block.id: set() for block in blocks}
    for connection in connections:
        if connection.axis_index == axis_index:
            continue
        graph[connection.a_id].add(connection.b_id)
        graph[connection.b_id].add(connection.a_id)

    result: list[set[int]] = []
    unseen = set(graph)
    while unseen:
        seed = min(unseen)
        component = {seed}
        queue = [seed]
        unseen.remove(seed)
        while queue:
            current = queue.pop(0)
            for neighbor in graph[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        result.append(component)
    return result


def _same_axis_grid(a: Block, b: Block, axis_index: int) -> bool:
    if not (
        coordinates_close(
            a.lmin[axis_index],
            b.lmin[axis_index],
            scale=max(axis_extent(a, axis_index), axis_extent(b, axis_index)),
        )
        and coordinates_close(
            a.lmax[axis_index],
            b.lmax[axis_index],
            scale=max(axis_extent(a, axis_index), axis_extent(b, axis_index)),
        )
    ):
        return False
    try:
        a_faces = _interior_faces(a, axis_index)
        b_faces = _interior_faces(b, axis_index)
    except Exception:
        return False
    if a_faces.shape != b_faces.shape:
        return False
    scale = max(axis_extent(a, axis_index), axis_extent(b, axis_index))
    tolerance = max(
        GRID_TOL * scale,
        geometry_tolerance(
            float(a_faces[0]),
            float(a_faces[-1]),
            float(b_faces[0]),
            float(b_faces[-1]),
            scale=scale,
        ),
    )
    return float(np.max(np.abs(a_faces - b_faces))) <= tolerance
