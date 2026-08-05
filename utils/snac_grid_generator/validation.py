"""Structured-grid checks for the grid generator."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import prod
from typing import Any

import numpy as np

from .grid import AXIS_NAMES, axis_grid_arrays
from .model import AxisSpec, Block, Project
from .numeric import coordinates_close, geometry_tolerance
from .snac_bc import clear_friend_boundaries, connect_friend_boundaries
from .snac_validation import validate_snac_project
from .topology import (
    FACE_INFO,
    FACE_ORDER,
    FaceConnection,
    axis_extent,
    build_topology,
)

GRID_TOL = 1.0e-9
SPACING_JUMP_WARNING_RATIO = 3.0


@dataclass
class CheckResult:
    """Validation result returned by the GUI API."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    interfaces: list[dict[str, Any]] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def extend(self, other: "CheckResult") -> None:
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        self.interfaces.extend(other.interfaces)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ok": self.ok,
            "errors": self.errors,
            "warnings": self.warnings,
            "interfaces": self.interfaces,
        }


def check_project(project: Project | dict[str, Any]) -> CheckResult:
    """Check that a project can be exported as a structured SNaC grid."""

    project = _project_copy(project)
    result = CheckResult()
    topology = build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    if not project.blocks:
        result.errors.append("project has no blocks and cannot be exported as a SNaC case")
    result.errors.extend(_check_connected_components(project.blocks, topology.connections))
    by_id = {block.id: block for block in project.blocks}

    result.errors.extend(_check_block_basics(project.blocks))
    result.errors.extend(validate_snac_project(project, topology.connections))
    result.errors.extend(_check_connected_block_grids(by_id, topology.connections))
    result.errors.extend(_check_decomposition_request(project))
    result.interfaces, spacing_warnings = _interface_spacing_metrics(by_id, topology.connections)
    result.warnings.extend(spacing_warnings)
    return result


def infer_project_connectivity(project: Project | dict[str, Any]) -> tuple[Project, CheckResult]:
    """Return a copy with full-face ``F`` boundary connectivity inferred."""

    project = _project_copy(project)
    result = CheckResult()
    topology = build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    if result.errors:
        return project, result

    clear_friend_boundaries(project.blocks)
    by_id = {block.id: block for block in project.blocks}
    for connection in topology.connections:
        connect_friend_boundaries(
            by_id[connection.a_id],
            connection.a_face,
            by_id[connection.b_id],
            connection.b_face,
        )
    return project, result


def update_project_structure(
    project: Project | dict[str, Any],
    source_block_id: int | None = None,
) -> tuple[Project, CheckResult]:
    """Infer friend BCs and propagate MPI partitions across structured faces."""

    project, result = infer_project_connectivity(project)
    if result.errors:
        return project, result

    topology = build_topology(project.blocks, project.periodic_axes)
    result.warnings.extend(_propagate_mpi_partitions(project.blocks, topology.connections, source_block_id))
    result.extend(check_project(project))
    return project, result


def apply_axis_to_aligned_blocks(
    project: Project | dict[str, Any],
    source_block_id: int,
    axis: str,
) -> tuple[Project, list[int]]:
    """Copy one axis grid to face-connected blocks with the same axis extent."""

    project = _project_copy(project)
    if axis not in AXIS_NAMES:
        raise ValueError(f"unknown axis {axis!r}")
    source = next((block for block in project.blocks if block.id == source_block_id), None)
    if source is None:
        raise ValueError(f"block {source_block_id} does not exist")

    axis_index = AXIS_NAMES.index(axis)
    topology = build_topology(project.blocks, project.periodic_axes)
    if topology.errors:
        raise ValueError(topology.errors[0])
    by_id = {block.id: block for block in project.blocks}
    adjacency: dict[int, list[tuple[int, int]]] = {block.id: [] for block in project.blocks}
    for connection in topology.connections:
        adjacency[connection.a_id].append((connection.b_id, connection.axis_index))
        adjacency[connection.b_id].append((connection.a_id, connection.axis_index))

    aligned = {source.id}
    pending = [source.id]
    while pending:
        current = pending.pop()
        for neighbor_id, normal_axis in adjacency[current]:
            neighbor = by_id[neighbor_id]
            if normal_axis == axis_index or neighbor_id in aligned:
                continue
            if not (
                coordinates_close(
                    neighbor.lmin[axis_index],
                    source.lmin[axis_index],
                    scale=max(axis_extent(neighbor, axis_index), axis_extent(source, axis_index)),
                )
                and coordinates_close(
                    neighbor.lmax[axis_index],
                    source.lmax[axis_index],
                    scale=max(axis_extent(neighbor, axis_index), axis_extent(source, axis_index)),
                )
            ):
                continue
            aligned.add(neighbor_id)
            pending.append(neighbor_id)

    changed: list[int] = []
    for block_id in sorted(aligned - {source.id}):
        block = by_id[block_id]
        block.axes[axis] = AxisSpec.from_dict(source.axes[axis].to_dict())
        block.ng[axis_index] = source.ng[axis_index]
        changed.append(block_id)
    project.validate()
    return project, changed


def _project_copy(project: Project | dict[str, Any]) -> Project:
    if isinstance(project, Project):
        return Project.from_dict(project.to_dict())
    return Project.from_dict(project)


def _check_block_basics(blocks: list[Block]) -> list[str]:
    errors: list[str] = []
    for block in blocks:
        for idx, axis in enumerate(AXIS_NAMES):
            if block.dims[idx] > block.ng[idx]:
                errors.append(f"block {block.id}: mpi partitions along {axis} exceed ng")
            try:
                axis_grid_arrays(block.axes[axis].to_dict(), block.lmin[idx], block.lmax[idx], block.ng[idx])
            except Exception as exc:
                errors.append(f"block {block.id}: invalid {axis} grid ({exc})")
    return errors


def _check_decomposition_request(project: Project) -> list[str]:
    target = project.decomposition.target_ranks
    if target == 0:
        return []
    errors: list[str] = []
    if target < len(project.blocks):
        errors.append(f"target MPI ranks must be at least the number of blocks ({len(project.blocks)})")
    total = sum(prod(block.dims) for block in project.blocks)
    if total != target:
        errors.append(f"MPI decomposition uses {total} ranks, expected {target}")

    active = set()
    min_local_cells = project.decomposition.min_local_cells
    max_local_aspect = project.decomposition.max_local_aspect
    for block in project.blocks:
        local_extents = [block.ng[index] / block.dims[index] for index in range(3)]
        split_extents = [
            block.ng[index] // block.dims[index]
            for index in range(3)
            if block.dims[index] > 1
        ]
        if split_extents and min(split_extents) < min_local_cells:
            errors.append(
                f"block {block.id}: decomposition leaves fewer than {min_local_cells} cells on a rank"
            )
        aspect = max(local_extents) / min(local_extents)
        if max_local_aspect and aspect > max_local_aspect:
            errors.append(
                f"block {block.id}: local partition aspect {aspect:.3g} exceeds {max_local_aspect:.3g}"
            )
    for axis_index, axis in enumerate(AXIS_NAMES):
        if any(block.dims[axis_index] > 1 for block in project.blocks):
            active.add(axis_index)
        if not project.decomposition.axes[axis_index] and any(
            block.dims[axis_index] > 1 for block in project.blocks
        ):
            errors.append(f"MPI decomposition uses disabled {axis} partitions")
    if project.decomposition.mode != "auto" and target > len(project.blocks):
        expected = int(project.decomposition.mode[0])
        if len(active) != expected:
            errors.append(
                f"{project.decomposition.mode} decomposition requires {expected} active axes, found {len(active)}"
            )
    return errors


def _check_connected_components(blocks: list[Block], connections: list[FaceConnection]) -> list[str]:
    if len(blocks) < 2:
        return []

    graph: dict[int, set[int]] = {block.id: set() for block in blocks}
    for connection in connections:
        graph[connection.a_id].add(connection.b_id)
        graph[connection.b_id].add(connection.a_id)

    seed = min(graph)
    seen = {seed}
    queue = [seed]
    while queue:
        current = queue.pop(0)
        for neighbor in graph[current]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    missing = sorted(set(graph) - seen)
    if missing:
        return [f"blocks {missing} are disconnected from block {seed}"]
    return []


def _check_connected_block_grids(by_id: dict[int, Block], connections: list[FaceConnection]) -> list[str]:
    errors: list[str] = []
    for connection in connections:
        a = by_id[connection.a_id]
        b = by_id[connection.b_id]
        normal = connection.axis_index
        for axis_index, axis in enumerate(AXIS_NAMES):
            if axis_index == normal:
                continue
            if a.ng[axis_index] != b.ng[axis_index]:
                errors.append(f"blocks {a.id} and {b.id}: touching grids have different ng along {axis}")
            if a.dims[axis_index] != b.dims[axis_index]:
                errors.append(f"blocks {a.id} and {b.id}: touching MPI partitions have different dims along {axis}")
            try:
                a_faces = _interior_faces(a, axis_index)
                b_faces = _interior_faces(b, axis_index)
            except Exception as exc:
                errors.append(f"blocks {a.id} and {b.id}: could not compare {axis} grid ({exc})")
                continue
            if a_faces.shape != b_faces.shape:
                continue
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
            if float(np.max(np.abs(a_faces - b_faces))) > tolerance:
                errors.append(f"blocks {a.id} and {b.id}: touching grid lines do not match along {axis}")
    return errors


def interface_spacing_metrics(project: Project | dict[str, Any]) -> list[dict[str, Any]]:
    """Return normal-spacing ratios for every valid full-face connection."""

    project = _project_copy(project)
    topology = build_topology(project.blocks, project.periodic_axes)
    metrics, _ = _interface_spacing_metrics(
        {block.id: block for block in project.blocks},
        topology.connections,
    )
    return metrics


def _interface_spacing_metrics(
    by_id: dict[int, Block],
    connections: list[FaceConnection],
) -> tuple[list[dict[str, Any]], list[str]]:
    metrics: list[dict[str, Any]] = []
    warnings: list[str] = []
    for connection in connections:
        a = by_id[connection.a_id]
        b = by_id[connection.b_id]
        axis = connection.axis_index
        try:
            a_width = _boundary_width(a, axis, connection.a_face)
            b_width = _boundary_width(b, axis, connection.b_face)
        except Exception as exc:
            warnings.append(f"blocks {a.id} and {b.id}: could not compare interface spacing ({exc})")
            continue
        ratio = max(a_width, b_width) / min(a_width, b_width)
        metrics.append(
            {
                "blockId": a.id,
                "face": FACE_ORDER[connection.a_face],
                "neighborId": b.id,
                "neighborFace": FACE_ORDER[connection.b_face],
                "axis": AXIS_NAMES[axis],
                "ratio": ratio,
                "spacing": a_width,
                "neighborSpacing": b_width,
            }
        )
        if ratio > SPACING_JUMP_WARNING_RATIO:
            warnings.append(
                f"blocks {a.id} and {b.id}: {AXIS_NAMES[axis]} interface spacing jumps by {ratio:.3g}"
            )
    return metrics, warnings


def _propagate_mpi_partitions(
    blocks: list[Block],
    connections: list[FaceConnection],
    source_block_id: int | None,
) -> list[str]:
    warnings: list[str] = []
    by_id = {block.id: block for block in blocks}
    source_id = source_block_id if source_block_id in by_id else None

    for axis_index, axis in enumerate(AXIS_NAMES):
        graph: dict[int, set[int]] = {block.id: set() for block in blocks}
        for connection in connections:
            if connection.axis_index == axis_index:
                continue
            graph[connection.a_id].add(connection.b_id)
            graph[connection.b_id].add(connection.a_id)

        seen: set[int] = set()
        for block in sorted(blocks, key=lambda item: item.id):
            if block.id in seen:
                continue
            component = []
            queue = [block.id]
            seen.add(block.id)
            while queue:
                current = queue.pop(0)
                component.append(current)
                for neighbor in graph[current]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        queue.append(neighbor)

            if len(component) < 2:
                continue
            if source_id in component:
                target = by_id[source_id].dims[axis_index]
                owner = source_id
            else:
                owner = min(component)
                target = by_id[owner].dims[axis_index]
            for block_id in component:
                block = by_id[block_id]
                if block.dims[axis_index] != target:
                    warnings.append(
                        f"block {block_id}: mpi {axis} partition changed from {block.dims[axis_index]} "
                        f"to {target} to match block {owner}"
                    )
                    block.dims[axis_index] = target
    return warnings


def _interior_faces(block: Block, axis_index: int) -> np.ndarray:
    axis = AXIS_NAMES[axis_index]
    arrays = axis_grid_arrays(
        block.axes[axis].to_dict(),
        block.lmin[axis_index],
        block.lmax[axis_index],
        block.ng[axis_index],
    )
    return arrays.faces[:-1]


def _boundary_width(block: Block, axis_index: int, face: int) -> float:
    axis = AXIS_NAMES[axis_index]
    arrays = axis_grid_arrays(
        block.axes[axis].to_dict(),
        block.lmin[axis_index],
        block.lmax[axis_index],
        block.ng[axis_index],
    )
    widths = arrays.face_spacing[1:-1]
    _, side = FACE_INFO[face]
    return float(widths[0] if side == 0 else widths[-1])
