"""SNaC-specific project checks mirrored from the solver input contract."""

from __future__ import annotations

from math import isfinite, sqrt

import numpy as np

from .axis_format import NATIVE
from .diagnostics import Diagnostic, error
from .grid import AXIS_NAMES, axis_grid_arrays
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


def validate_snac_project(project: Project, connections: list[FaceConnection]) -> list[Diagnostic]:
    """Return errors for SNaC fields represented by the grid project."""

    errors = _check_boundary_pairs(project.blocks)
    errors.extend(_check_friend_boundaries(project.blocks, connections))
    errors.extend(_check_native_grids(project.blocks))
    errors.extend(_check_initial_velocity_fields(project.blocks))
    errors.extend(_check_periodic_inflow(project))
    errors.extend(_check_fft_axis(project))
    non_native = sorted(
        {block.id for block in project.blocks for axis in AXIS_NAMES if block.axes[axis].kind != NATIVE}
    )
    if not project.write_external_grid and non_native:
        errors.append(
            error("external_grid_required", "non-native grids require external grid writing", *non_native)
        )
    return errors


def fft_spacing_tolerance(lmin: float, lmax: float, n: int, precision: str) -> float:
    """Absolute spacing tolerance for a uniform FFT-direction grid.

    SNaC rejects a loaded FFT-direction grid whose spacing deviates from
    (lmax - lmin)/n by more than sqrt(eps)*dl + 16*eps*max(|lmin|, |lmax|) in its
    working precision. The generator uses a ten times tighter relative bound.
    """

    eps = float(np.finfo(np.float32 if precision == "single" else np.float64).eps)
    spacing = (lmax - lmin) / n
    return 0.1 * sqrt(eps) * spacing + 16.0 * eps * max(abs(lmin), abs(lmax))


def _check_fft_axis(project: Project) -> list[Diagnostic]:
    """Mirror SNaC's FFT-direction sanity checks and require a uniform grid."""

    index = project.fft_axis_index
    if index is None or not project.blocks:
        return []
    axis = AXIS_NAMES[index]
    errors: list[Diagnostic] = []
    blocks = sorted(project.blocks, key=lambda item: item.id)
    for block in blocks:
        if block.dims[index] != 1:
            errors.append(
                error(
                    "fft_axis",
                    f"block {block.id}: SNaC allows no MPI partitions along the FFT axis {axis}",
                    block.id,
                    axis=axis,
                )
            )

    # SNaC compares these values exactly across blocks.
    reference = blocks[0]
    key = (reference.ng[index], reference.lmin[index], reference.lmax[index])
    differing = [
        block.id for block in blocks[1:] if (block.ng[index], block.lmin[index], block.lmax[index]) != key
    ]
    if differing:
        errors.append(
            error(
                "fft_axis",
                f"blocks {differing} differ from block {reference.id} in ng, lmin, or lmax "
                f"along the FFT axis {axis}",
                *differing,
                reference.id,
                axis=axis,
            )
        )

    for block in blocks:
        specification = block.axes[axis]
        if specification.kind == NATIVE:
            if specification.gt not in range(4) or specification.gr != 0.0:
                errors.append(
                    error(
                        "fft_axis",
                        f"block {block.id}: the FFT axis {axis} needs a uniform native grid, "
                        "so use mapping function 0 to 3 with gr = 0",
                        block.id,
                        axis=axis,
                    )
                )
            continue
        try:
            arrays = axis_grid_arrays(
                specification.to_dict(), block.lmin[index], block.lmax[index], block.ng[index]
            )
        except ValueError:
            continue  # reported as an invalid grid by the structural checks
        widths = arrays.face_spacing[1:-1]
        uniform = (block.lmax[index] - block.lmin[index]) / block.ng[index]
        tolerance = fft_spacing_tolerance(
            block.lmin[index], block.lmax[index], block.ng[index], project.external_grid_precision
        )
        if float(np.max(np.abs(widths - uniform))) > tolerance:
            errors.append(
                error(
                    "fft_axis",
                    f"block {block.id}: the grid along the FFT axis {axis} must be uniform, "
                    f"but its spacing varies by a factor of {float(widths.max() / widths.min()):.6g}",
                    block.id,
                    axis=axis,
                )
            )
    return errors


def _check_native_grids(blocks: list[Block]) -> list[Diagnostic]:
    errors: list[Diagnostic] = []
    for block in blocks:
        for axis in AXIS_NAMES:
            specification = block.axes[axis]
            if specification.kind != NATIVE:
                continue
            if specification.gt not in SNAC_GRID_FUNCTIONS:
                errors.append(
                    error(
                        "grid_mapping",
                        f"block {block.id}: unsupported SNaC {axis} mapping function {specification.gt}",
                        block.id,
                        axis=axis,
                    )
                )
            if not isfinite(specification.gr) or specification.gr < 0.0:
                errors.append(
                    error(
                        "grid_mapping",
                        f"block {block.id}: SNaC {axis} grid growth parameter must be finite and non-negative",
                        block.id,
                        axis=axis,
                    )
                )
    return errors


def _check_initial_velocity_fields(blocks: list[Block]) -> list[Diagnostic]:
    errors: list[Diagnostic] = []
    supported = set(SNAC_INITIAL_VELOCITY_FIELDS)
    for block in blocks:
        if block.inivel not in supported:
            errors.append(
                error(
                    "initial_velocity",
                    f"block {block.id}: unsupported SNaC initial velocity field {block.inivel!r}",
                    block.id,
                )
            )
    return errors


def _check_periodic_inflow(project: Project) -> list[Diagnostic]:
    errors: list[Diagnostic] = []
    for axis_index, periodic in enumerate(project.periodic_axes):
        if not periodic:
            continue
        for block in project.blocks:
            for side in (0, 1):
                face = FACE_INDEX[(axis_index, side)]
                if block.inflow[face] > 0:
                    errors.append(
                        error(
                            "boundary_condition",
                            f"block {block.id}: inflow type on periodic face {FACE_ORDER[face]} is invalid",
                            (block.id, FACE_ORDER[face]),
                            axis=AXIS_NAMES[axis_index],
                        )
                    )
    return errors


def _check_boundary_pairs(blocks: list[Block]) -> list[Diagnostic]:
    errors: list[Diagnostic] = []
    for block in blocks:
        for axis_index, axis in enumerate(AXIS_NAMES):
            lo = FACE_INDEX[(axis_index, 0)]
            hi = FACE_INDEX[(axis_index, 1)]
            faces = ((block.id, FACE_ORDER[lo]), (block.id, FACE_ORDER[hi]))

            def invalid(message: str) -> None:
                errors.append(error("boundary_condition", f"block {block.id}: {message}", *faces, axis=axis))

            pressure_pair = f"{block.cbcpre[lo]}{block.cbcpre[hi]}"
            if pressure_pair not in _BC_PAIRS:
                invalid(f"pressure BC pair {pressure_pair} along {axis} is invalid")
            for component in range(3):
                velocity_pair = f"{block.cbcvel[component][lo]}{block.cbcvel[component][hi]}"
                if velocity_pair not in _BC_PAIRS:
                    invalid(f"velocity component {component + 1} BC pair {velocity_pair} along {axis} is invalid")
            for iscal, cbcscal in enumerate(block.cbcscal, start=1):
                scalar_pair = f"{cbcscal[lo]}{cbcscal[hi]}"
                if scalar_pair not in _BC_PAIRS:
                    invalid(f"scalar {iscal} BC pair {scalar_pair} along {axis} is invalid")
            normal_velocity_pair = f"{block.cbcvel[axis_index][lo]}{block.cbcvel[axis_index][hi]}"
            if (normal_velocity_pair, pressure_pair) not in _NORMAL_BC_PAIRS:
                invalid(
                    f"normal velocity BC pair {normal_velocity_pair} is incompatible with "
                    f"pressure pair {pressure_pair} along {axis}"
                )
    return errors


def _check_friend_boundaries(blocks: list[Block], connections: list[FaceConnection]) -> list[Diagnostic]:
    errors: list[Diagnostic] = []
    expected: dict[tuple[int, int], tuple[int, int]] = {}
    existing_ids = {block.id for block in blocks}
    for connection in connections:
        expected[(connection.a_id, connection.a_face)] = (connection.b_id, connection.b_face)
        expected[(connection.b_id, connection.b_face)] = (connection.a_id, connection.a_face)

    for block in blocks:
        for face, label in enumerate(FACE_ORDER):
            neighbor = expected.get((block.id, face))
            expected_friend = None if neighbor is None else neighbor[0]
            neighbor_face = None if neighbor is None else (neighbor[0], FACE_ORDER[neighbor[1]])
            here = (block.id, label)
            axis = label[0]
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
                    errors.append(
                        error(
                            "friend_boundary",
                            f"block {block.id} face {label} is marked F but has no full-face neighbor",
                            here,
                            axis=axis,
                        )
                    )
                continue

            if any(code != "F" for code in friend_codes):
                errors.append(
                    error(
                        "friend_boundary",
                        f"block {block.id} face {label} touches block {expected_friend} but is not fully marked F",
                        here,
                        neighbor_face,
                        axis=axis,
                    )
                )
                continue

            for friend in dict.fromkeys(actual_friends):  # one message per referenced block
                if friend not in existing_ids:
                    errors.append(
                        error(
                            "friend_boundary",
                            f"block {block.id} face {label} references missing block {friend}",
                            here,
                            neighbor_face,
                            axis=axis,
                        )
                    )
                elif friend != expected_friend:
                    errors.append(
                        error(
                            "friend_boundary",
                            f"block {block.id} face {label} references block {friend}, expected {expected_friend}",
                            here,
                            neighbor_face,
                            friend,
                            axis=axis,
                        )
                    )
    return errors
