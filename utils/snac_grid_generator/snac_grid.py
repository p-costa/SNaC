"""SNaC-native grid mappings and binary serialization."""

from __future__ import annotations

from dataclasses import dataclass
from math import tanh
from pathlib import Path
from typing import Callable, Protocol

import numpy as np

SNAC_GRID_FUNCTIONS = {
    0: "Cluster both ends",
    1: "Cluster lower end",
    2: "Cluster middle",
    3: "Cluster upper end",
    4: "Geometric lower end",
    5: "Geometric upper end",
    6: "Geometric both ends",
    7: "Geometric middle",
}


class GridArraysLike(Protocol):
    faces: np.ndarray
    centers: np.ndarray
    face_spacing: np.ndarray
    center_spacing: np.ndarray


@dataclass(frozen=True)
class BinaryGrid:
    """Validated arrays read from one SNaC grid binary."""

    faces: np.ndarray
    centers: np.ndarray
    face_spacing: np.ndarray
    center_spacing: np.ndarray
    dtype: str


def native_faces(n: int, lmin: float, lmax: float, gt: int, gr: float) -> np.ndarray:
    """Evaluate a mapping from ``src/initgrid.f90`` at SNaC face positions."""

    mapper = _SNAC_MAPPERS.get(gt, _gridpoint_cluster_two_end)
    r0 = np.arange(1, n + 1, dtype=float) / float(n)
    mapped = np.array([mapper(gr, value) for value in r0], dtype=float)
    faces = np.concatenate(([lmin], lmin + mapped * (lmax - lmin)))
    faces[-1] = lmax
    return faces


def binary_payload(arrays: GridArraysLike, precision: str = "double") -> np.ndarray:
    """Return the four-array payload expected by SNaC's ``load_grid``."""

    dtype = {"single": "<f4", "double": "<f8"}.get(precision)
    if dtype is None:
        raise ValueError("grid binary precision must be 'single' or 'double'")
    physical_faces = np.asarray(arrays.faces[:-1], dtype=dtype)
    if not np.all(np.isfinite(physical_faces)) or np.any(np.diff(physical_faces) <= 0.0):
        raise ValueError(
            f"{precision}-precision grid files cannot represent all face coordinates distinctly"
        )
    payload = np.concatenate(
        (
            arrays.faces[1:-1],
            arrays.centers[1:-1],
            arrays.face_spacing[1:-1],
            arrays.center_spacing[1:-1],
        )
    ).astype(dtype, copy=False)
    if not np.all(np.isfinite(payload)):
        raise ValueError(f"{precision}-precision grid files contain non-finite values")
    return payload


def read_grid_binary(
    path: str | Path,
    n: int,
    lmin: float,
    lmax: float,
) -> BinaryGrid:
    """Read and validate the four arrays written by SNaC ``save_grid``."""

    if n < 1:
        raise ValueError("grid binary cell count must be positive")
    if not lmax > lmin:
        raise ValueError("grid binary extent must be positive")
    path = Path(path)
    payload = path.read_bytes()
    expected_sizes = {4 * n * itemsize for itemsize in (4, 8)}
    if len(payload) not in expected_sizes:
        sizes = " or ".join(str(value) for value in sorted(expected_sizes))
        raise ValueError(
            f"{path} has {len(payload)} bytes; expected {sizes} bytes for {n} cells"
        )

    itemsize = len(payload) // (4 * n)
    failures: list[str] = []
    for endian in ("<", ">"):
        dtype = np.dtype(f"{endian}f{itemsize}")
        values = np.frombuffer(payload, dtype=dtype).reshape(4, n).astype(float)
        try:
            return _validated_binary_grid(values, dtype, lmin, lmax)
        except ValueError as exc:
            failures.append(str(exc))
    raise ValueError(
        f"{path} is not a consistent SNaC grid binary for {n} cells "
        f"({'; '.join(dict.fromkeys(failures))})"
    )


def _validated_binary_grid(
    values: np.ndarray,
    dtype: np.dtype,
    lmin: float,
    lmax: float,
) -> BinaryGrid:
    if not np.all(np.isfinite(values)):
        raise ValueError("contains non-finite values")

    stored_faces, centers, face_spacing, center_spacing = values
    faces = np.concatenate(([float(lmin)], stored_faces))
    length = abs(float(lmax) - float(lmin))
    relative_tolerance = 5.0e-5 if dtype.itemsize == 4 else 5.0e-11
    coordinate_tolerance = max(
        relative_tolerance * length,
        np.finfo(dtype).eps * max(abs(lmin), abs(lmax), length) * 16.0,
    )
    spacing_tolerance = max(
        relative_tolerance * length,
        np.finfo(dtype).eps * length * 16.0,
    )
    if abs(float(faces[-1]) - float(lmax)) > coordinate_tolerance:
        raise ValueError("upper face does not match the block extent")
    faces[-1] = float(lmax)
    widths = np.diff(faces)
    if np.any(widths <= 0.0):
        raise ValueError("face coordinates are not strictly increasing")

    expected_centers = 0.5 * (faces[:-1] + faces[1:])
    expected_center_spacing = np.concatenate(
        (0.5 * (widths[:-1] + widths[1:]), [widths[-1]])
    )
    checks = (
        ("cell centers", centers, expected_centers, 0.0, coordinate_tolerance),
        ("face spacings", face_spacing, widths, relative_tolerance, spacing_tolerance),
        (
            "center spacings",
            center_spacing,
            expected_center_spacing,
            relative_tolerance,
            spacing_tolerance,
        ),
    )
    for label, actual, expected, rtol, atol in checks:
        if not np.allclose(
            actual,
            expected,
            rtol=rtol,
            atol=atol,
        ):
            raise ValueError(f"stored {label} disagree with the face coordinates")

    return BinaryGrid(
        faces=faces,
        centers=centers.copy(),
        face_spacing=face_spacing.copy(),
        center_spacing=center_spacing.copy(),
        dtype=dtype.str,
    )


def _gridpoint_cluster_two_end(alpha: float, r0: float) -> float:
    if alpha != 0.0:
        return 0.5 * (1.0 + tanh((r0 - 0.5) * alpha) / tanh(alpha / 2.0))
    return r0


def _gridpoint_cluster_one_end(alpha: float, r0: float) -> float:
    if alpha != 0.0:
        return 1.0 + tanh((r0 - 1.0) * alpha) / tanh(alpha)
    return r0


def _gridpoint_cluster_one_end_r(alpha: float, r0: float) -> float:
    if alpha != 0.0:
        return 1.0 - (1.0 + tanh(((1.0 - r0) - 1.0) * alpha) / tanh(alpha))
    return r0


def _gridpoint_cluster_middle(alpha: float, r0: float) -> float:
    if alpha == 0.0:
        return r0
    if r0 <= 0.5:
        return 0.5 * tanh(2.0 * alpha * r0) / tanh(alpha)
    return 0.5 * (2.0 + tanh(2.0 * alpha * (r0 - 1.0)) / tanh(alpha))


def _gridpoint_cluster_geometric_one_end(alpha: float, r0: float) -> float:
    power_value = 1.0 / 3.0
    if r0 == 1.0:
        return r0
    if alpha == 0.0:
        return r0
    return r0 * (((1.0 - r0**alpha) / (1.0 - r0) / alpha) ** power_value)


def _gridpoint_cluster_geometric_one_end_r(alpha: float, r0: float) -> float:
    return 1.0 - _gridpoint_cluster_geometric_one_end(alpha, 1.0 - r0)


def _gridpoint_cluster_geometric_two_ends(alpha: float, r0: float) -> float:
    if r0 <= 0.5:
        return 0.5 * _gridpoint_cluster_geometric_one_end(alpha / 2.0, 2.0 * r0)
    return 1.0 - 0.5 * _gridpoint_cluster_geometric_one_end(alpha / 2.0, 2.0 * (1.0 - r0))


def _gridpoint_cluster_geometric_middle(alpha: float, r0: float) -> float:
    if r0 <= 0.5:
        return 0.5 * (1.0 - _gridpoint_cluster_geometric_one_end(alpha / 2.0, (0.5 - r0) * 2.0))
    return 0.5 * (1.0 + _gridpoint_cluster_geometric_one_end(alpha / 2.0, (r0 - 0.5) * 2.0))


_SNAC_MAPPERS: dict[int, Callable[[float, float], float]] = {
    0: _gridpoint_cluster_two_end,
    1: _gridpoint_cluster_one_end,
    2: _gridpoint_cluster_middle,
    3: _gridpoint_cluster_one_end_r,
    4: _gridpoint_cluster_geometric_one_end,
    5: _gridpoint_cluster_geometric_one_end_r,
    6: _gridpoint_cluster_geometric_two_ends,
    7: _gridpoint_cluster_geometric_middle,
}
