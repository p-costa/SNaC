"""Structured-grid checks and repair helpers for the grid generator."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isclose
from typing import Any

import numpy as np

from .grid import AXIS_NAMES, axis_grid_arrays
from .model import Block, Project

FACE_ORDER = ("x-", "x+", "y-", "y+", "z-", "z+")
FACE_INDEX = {
    (0, 0): 0,
    (0, 1): 1,
    (1, 0): 2,
    (1, 1): 3,
    (2, 0): 4,
    (2, 1): 5,
}
FACE_INFO = {
    0: (0, 0),
    1: (0, 1),
    2: (1, 0),
    3: (1, 1),
    4: (2, 0),
    5: (2, 1),
}

GEOM_TOL = 1.0e-10
GRID_TOL = 1.0e-9
SPACING_JUMP_WARNING_RATIO = 3.0

_BC_PAIRS = {"ND", "DN", "NN", "DD", "FD", "DF", "FF", "FN", "NF"}
_NORMAL_BC_PAIRS = {
    ("FF", "FF"),
    ("ND", "DN"),
    ("DN", "ND"),
    ("DD", "NN"),
    ("FD", "FN"),
    ("DF", "NF"),
    ("FN", "FD"),
    ("NF", "DF"),
    ("FN", "FN"),
    ("NF", "NF"),
    ("NN", "NN"),
    ("DN", "NN"),
    ("ND", "NN"),
    ("NN", "DD"),
}


@dataclass(frozen=True)
class FaceConnection:
    """Full-face connection between two blocks."""

    a_id: int
    a_face: int
    b_id: int
    b_face: int
    axis_index: int


@dataclass
class CheckResult:
    """Validation result returned by the GUI API."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def extend(self, other: "CheckResult") -> None:
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)

    def to_dict(self) -> dict[str, Any]:
        return {"ok": self.ok, "errors": self.errors, "warnings": self.warnings}


@dataclass
class _Topology:
    connections: list[FaceConnection] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def check_project(project: Project | dict[str, Any]) -> CheckResult:
    """Check that a project can be exported as a structured SNaC grid."""

    project = _project_copy(project)
    result = CheckResult()
    result.errors.extend(_check_periodic_geometry(project.blocks, project.periodic_axes))
    topology = _build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    by_id = {block.id: block for block in project.blocks}

    result.errors.extend(_check_block_basics(project.blocks))
    result.errors.extend(_check_boundary_pairs(project.blocks))
    result.errors.extend(_check_friend_boundaries(project.blocks, topology.connections))
    result.errors.extend(_check_connected_block_grids(by_id, topology.connections))
    result.warnings.extend(_check_spacing_jumps(by_id, topology.connections))
    return result


def infer_project_connectivity(project: Project | dict[str, Any]) -> tuple[Project, CheckResult]:
    """Return a copy with full-face ``F`` boundary connectivity inferred."""

    project = _project_copy(project)
    result = CheckResult()
    result.errors.extend(_check_periodic_geometry(project.blocks, project.periodic_axes))
    topology = _build_topology(project.blocks, project.periodic_axes)
    result.errors.extend(topology.errors)
    result.warnings.extend(topology.warnings)
    if result.errors:
        return project, result

    _clear_friend_boundaries(project.blocks)
    by_id = {block.id: block for block in project.blocks}
    for connection in topology.connections:
        _connect(by_id[connection.a_id], connection.a_face, by_id[connection.b_id], connection.b_face)
    return project, result


def update_project_structure(
    project: Project | dict[str, Any],
    source_block_id: int | None = None,
) -> tuple[Project, CheckResult]:
    """Infer friend BCs and propagate MPI partitions across structured faces."""

    project, result = infer_project_connectivity(project)
    if result.errors:
        return project, result

    topology = _build_topology(project.blocks, project.periodic_axes)
    result.warnings.extend(_propagate_mpi_partitions(project.blocks, topology.connections, source_block_id))
    result.extend(check_project(project))
    return project, result


def _project_copy(project: Project | dict[str, Any]) -> Project:
    if isinstance(project, Project):
        return Project.from_dict(project.to_dict())
    return Project.from_dict(project)


def _build_topology(blocks: list[Block], periodic_axes: list[bool] | None = None) -> _Topology:
    topology = _Topology()
    face_owner: dict[tuple[int, int], FaceConnection] = {}
    periodic_axes = periodic_axes or [False, False, False]

    for i, left in enumerate(blocks):
        for right in blocks[i + 1 :]:
            if _volume_overlap(left, right):
                topology.errors.append(f"blocks {left.id} and {right.id} overlap")
                continue
            for axis_index in range(3):
                if _touches(left.lmax[axis_index], right.lmin[axis_index]):
                    _add_face_contact(topology, face_owner, left, 1, right, 0, axis_index)
                if _touches(right.lmax[axis_index], left.lmin[axis_index]):
                    _add_face_contact(topology, face_owner, right, 1, left, 0, axis_index)

    for block in blocks:
        for axis_index, is_periodic in enumerate(periodic_axes[:3]):
            if is_periodic:
                lo_face = FACE_INDEX[(axis_index, 0)]
                hi_face = FACE_INDEX[(axis_index, 1)]
                topology.connections.append(FaceConnection(block.id, lo_face, block.id, hi_face, axis_index))

    return topology


def _check_periodic_geometry(blocks: list[Block], periodic_axes: list[bool]) -> list[str]:
    errors: list[str] = []
    if not blocks:
        return errors
    reference = blocks[0]
    for axis_index, is_periodic in enumerate(periodic_axes[:3]):
        if not is_periodic:
            continue
        axis = AXIS_NAMES[axis_index]
        for block in blocks[1:]:
            if not (
                _touches(block.lmin[axis_index], reference.lmin[axis_index])
                and _touches(block.lmax[axis_index], reference.lmax[axis_index])
            ):
                errors.append(
                    f"periodic {axis} requires every block to share the same {axis} min/max extent"
                )
                break
    return errors


def _add_face_contact(
    topology: _Topology,
    face_owner: dict[tuple[int, int], FaceConnection],
    lower_block: Block,
    lower_side: int,
    upper_block: Block,
    upper_side: int,
    axis_index: int,
) -> None:
    if _same_cross_section(lower_block, upper_block, axis_index):
        lower_face = FACE_INDEX[(axis_index, lower_side)]
        upper_face = FACE_INDEX[(axis_index, upper_side)]
        connection = FaceConnection(lower_block.id, lower_face, upper_block.id, upper_face, axis_index)
        for block_id, face in ((connection.a_id, connection.a_face), (connection.b_id, connection.b_face)):
            old = face_owner.get((block_id, face))
            if old is not None:
                topology.errors.append(f"block {block_id} face {FACE_ORDER[face]} has multiple full-face neighbors")
            face_owner[(block_id, face)] = connection
        topology.connections.append(connection)
    elif _cross_section_overlap(lower_block, upper_block, axis_index):
        topology.errors.append(
            f"blocks {lower_block.id} and {upper_block.id} touch on a partial {AXIS_NAMES[axis_index]} face"
        )


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


def _check_boundary_pairs(blocks: list[Block]) -> list[str]:
    errors: list[str] = []
    for block in blocks:
        for axis_index, axis in enumerate(AXIS_NAMES):
            lo = FACE_INDEX[(axis_index, 0)]
            hi = FACE_INDEX[(axis_index, 1)]
            pressure_pair = f"{block.cbcpre[lo]}{block.cbcpre[hi]}"
            if pressure_pair not in _BC_PAIRS:
                errors.append(f"block {block.id}: pressure BC pair {pressure_pair} along {axis} is invalid")
            for component in range(3):
                velocity_pair = f"{block.cbcvel[component][lo]}{block.cbcvel[component][hi]}"
                if velocity_pair not in _BC_PAIRS:
                    errors.append(
                        f"block {block.id}: velocity component {component + 1} BC pair "
                        f"{velocity_pair} along {axis} is invalid"
                    )
            for iscal, cbcscal in enumerate(block.cbcscal, start=1):
                scalar_pair = f"{cbcscal[lo]}{cbcscal[hi]}"
                if scalar_pair not in _BC_PAIRS:
                    errors.append(f"block {block.id}: scalar {iscal} BC pair {scalar_pair} along {axis} is invalid")
            normal_velocity_pair = f"{block.cbcvel[axis_index][lo]}{block.cbcvel[axis_index][hi]}"
            if (normal_velocity_pair, pressure_pair) not in _NORMAL_BC_PAIRS:
                errors.append(
                    f"block {block.id}: normal velocity BC pair {normal_velocity_pair} is incompatible with "
                    f"pressure pair {pressure_pair} along {axis}"
                )
    return errors


def _check_friend_boundaries(blocks: list[Block], connections: list[FaceConnection]) -> list[str]:
    errors: list[str] = []
    expected: dict[tuple[int, int], int] = {}
    existing_ids = {block.id for block in blocks}
    for connection in connections:
        expected[(connection.a_id, connection.a_face)] = connection.b_id
        expected[(connection.b_id, connection.b_face)] = connection.a_id

    for block in blocks:
        for face, label in enumerate(FACE_ORDER):
            expected_friend = expected.get((block.id, face))
            friend_codes = [
                block.cbcpre[face],
                *(block.cbcvel[component][face] for component in range(3)),
                *(cbcscal[face] for cbcscal in block.cbcscal),
            ]
            friend_values = [
                block.bcpre[face],
                *(block.bcvel[component][face] for component in range(3)),
                *(bcscal[face] for bcscal in block.bcscal),
            ]
            actual_friends = [int(round(value)) for code, value in zip(friend_codes, friend_values) if code == "F"]

            if expected_friend is None:
                if actual_friends:
                    errors.append(f"block {block.id} face {label} is marked F but has no full-face neighbor")
                continue

            if any(code != "F" for code in friend_codes):
                errors.append(
                    f"block {block.id} face {label} touches block {expected_friend} but is not fully marked F"
                )
                continue

            for friend in actual_friends:
                if friend not in existing_ids:
                    errors.append(f"block {block.id} face {label} references missing block {friend}")
                elif friend != expected_friend:
                    errors.append(
                        f"block {block.id} face {label} references block {friend}, expected {expected_friend}"
                    )
    return errors


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
            scale = max(1.0, abs(a.lmax[axis_index] - a.lmin[axis_index]))
            if float(np.max(np.abs(a_faces - b_faces))) > GRID_TOL * scale:
                errors.append(f"blocks {a.id} and {b.id}: touching grid lines do not match along {axis}")
    return errors


def _check_spacing_jumps(by_id: dict[int, Block], connections: list[FaceConnection]) -> list[str]:
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
        if ratio > SPACING_JUMP_WARNING_RATIO:
            warnings.append(
                f"blocks {a.id} and {b.id}: {AXIS_NAMES[axis]} interface spacing jumps by {ratio:.3g}"
            )
    return warnings


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


def _clear_friend_boundaries(blocks: list[Block]) -> None:
    for block in blocks:
        for face in range(6):
            if block.cbcpre[face] == "F":
                block.cbcpre[face] = "N"
                block.bcpre[face] = 0.0
            for component in range(3):
                if block.cbcvel[component][face] == "F":
                    block.cbcvel[component][face] = "D"
                    block.bcvel[component][face] = 0.0
            for iscal in range(len(block.cbcscal)):
                if block.cbcscal[iscal][face] == "F":
                    block.cbcscal[iscal][face] = "N"
                    block.bcscal[iscal][face] = 0.0


def _connect(a: Block, a_face: int, b: Block, b_face: int) -> None:
    for component in range(3):
        a.cbcvel[component][a_face] = "F"
        b.cbcvel[component][b_face] = "F"
        a.bcvel[component][a_face] = float(b.id)
        b.bcvel[component][b_face] = float(a.id)
    for iscal in range(len(a.cbcscal)):
        a.cbcscal[iscal][a_face] = "F"
        b.cbcscal[iscal][b_face] = "F"
        a.bcscal[iscal][a_face] = float(b.id)
        b.bcscal[iscal][b_face] = float(a.id)
    a.cbcpre[a_face] = "F"
    b.cbcpre[b_face] = "F"
    a.bcpre[a_face] = float(b.id)
    b.bcpre[b_face] = float(a.id)


def _same_cross_section(a: Block, b: Block, axis_index: int) -> bool:
    for idx in range(3):
        if idx == axis_index:
            continue
        if not (_touches(a.lmin[idx], b.lmin[idx]) and _touches(a.lmax[idx], b.lmax[idx])):
            return False
    return True


def _cross_section_overlap(a: Block, b: Block, axis_index: int) -> bool:
    return all(
        _interval_overlap(a.lmin[idx], a.lmax[idx], b.lmin[idx], b.lmax[idx])
        for idx in range(3)
        if idx != axis_index
    )


def _volume_overlap(a: Block, b: Block) -> bool:
    return all(_interval_overlap(a.lmin[idx], a.lmax[idx], b.lmin[idx], b.lmax[idx]) for idx in range(3))


def _interval_overlap(a0: float, a1: float, b0: float, b1: float) -> bool:
    return min(a1, b1) - max(a0, b0) > GEOM_TOL


def _touches(a: float, b: float) -> bool:
    return isclose(a, b, rel_tol=0.0, abs_tol=GEOM_TOL)


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
