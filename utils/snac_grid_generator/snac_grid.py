"""SNaC-native grid mappings and binary serialization."""

from __future__ import annotations

from math import tanh
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


def native_faces(n: int, lmin: float, lmax: float, gt: int, gr: float) -> np.ndarray:
    """Evaluate a mapping from ``src/initgrid.f90`` at SNaC face positions."""

    mapper = _SNAC_MAPPERS.get(gt, _gridpoint_cluster_two_end)
    r0 = np.arange(1, n + 1, dtype=float) / float(n)
    mapped = np.array([mapper(gr, value) for value in r0], dtype=float)
    faces = np.concatenate(([lmin], lmin + mapped * (lmax - lmin)))
    faces[-1] = lmax
    return faces


def binary_payload(arrays: GridArraysLike) -> np.ndarray:
    """Return the four-array payload expected by SNaC's ``load_grid``."""

    return np.concatenate(
        (
            arrays.faces[1:-1],
            arrays.centers[1:-1],
            arrays.face_spacing[1:-1],
            arrays.center_spacing[1:-1],
        )
    ).astype("<f8", copy=False)


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
