"""SNaC-specific project checks mirrored from the solver input contract."""

from __future__ import annotations

from math import isfinite

from .grid import AXIS_NAMES
from .model import Block, Project
from .snac_grid import SNAC_GRID_FUNCTIONS
from .topology import FACE_INDEX, FACE_ORDER, FaceConnection

SNAC_INITIAL_VELOCITY_FIELDS = (
    "zer",
    "uni",
    "cou",
    "poi",
    "log",
    "hcp",
    "hcl",
    "tgv",
    "kov",
    "pdc",
)

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


def validate_snac_project(project: Project, connections: list[FaceConnection]) -> list[str]:
    """Return errors for SNaC fields represented by the grid project."""

    errors = _check_boundary_pairs(project.blocks)
    errors.extend(_check_friend_boundaries(project.blocks, connections))
    errors.extend(_check_native_grids(project.blocks))
    errors.extend(_check_initial_velocity_fields(project.blocks))
    errors.extend(_check_periodic_inflow(project))
    if not project.write_external_grid and any(
        block.axes[axis].kind != "snac" for block in project.blocks for axis in AXIS_NAMES
    ):
        errors.append("non-native grids require external grid writing")
    return errors


def _check_native_grids(blocks: list[Block]) -> list[str]:
    errors: list[str] = []
    for block in blocks:
        for axis in AXIS_NAMES:
            specification = block.axes[axis]
            if specification.kind != "snac":
                continue
            if specification.gt not in SNAC_GRID_FUNCTIONS:
                errors.append(f"block {block.id}: unsupported SNaC {axis} mapping function {specification.gt}")
            if not isfinite(specification.gr) or specification.gr < 0.0:
                errors.append(f"block {block.id}: SNaC {axis} grid growth parameter must be finite and non-negative")
    return errors


def _check_initial_velocity_fields(blocks: list[Block]) -> list[str]:
    errors: list[str] = []
    supported = set(SNAC_INITIAL_VELOCITY_FIELDS)
    for block in blocks:
        if block.inivel not in supported:
            errors.append(f"block {block.id}: unsupported SNaC initial velocity field {block.inivel!r}")
    return errors


def _check_periodic_inflow(project: Project) -> list[str]:
    errors: list[str] = []
    for axis_index, periodic in enumerate(project.periodic_axes):
        if not periodic:
            continue
        for block in project.blocks:
            for side in (0, 1):
                face = FACE_INDEX[(axis_index, side)]
                if block.inflow[face] > 0:
                    errors.append(
                        f"block {block.id}: inflow type on periodic face {FACE_ORDER[face]} is invalid"
                    )
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
