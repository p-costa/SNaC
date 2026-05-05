"""Grid spacing and SNaC binary-grid helpers.

The functions in this module intentionally mirror ``src/initgrid.f90`` for the
native SNaC mapping functions, and add OpenFOAM-style multi-grading as an
external-grid option.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, isclose, pow, tanh
from typing import Callable, Sequence

import numpy as np

AXIS_NAMES = ("x", "y", "z")

GRID_FUNCTIONS = {
    0: "Cluster both ends",
    1: "Cluster lower end",
    2: "Cluster middle",
    3: "Cluster upper end",
    4: "Geometric lower end",
    5: "Geometric upper end",
    6: "Geometric both ends",
    7: "Geometric middle",
}

_REL_TOL = 1.0e-11
_RMAX = 1.0e5


@dataclass(frozen=True)
class GridArrays:
    """SNaC-compatible one-dimensional grid arrays for one block axis."""

    faces: np.ndarray
    centers: np.ndarray
    face_spacing: np.ndarray
    center_spacing: np.ndarray

    @property
    def binary_payload(self) -> np.ndarray:
        """Return the array order written by ``save_grid`` in SNaC."""

        return np.concatenate(
            (
                self.faces[1:-1],
                self.centers[1:-1],
                self.face_spacing[1:-1],
                self.center_spacing[1:-1],
            )
        ).astype("<f8", copy=False)


@dataclass(frozen=True)
class SpacingSolution:
    """Derived values for a geometric spacing sequence."""

    n: int
    length: float
    ratio: float
    cell_ratio: float
    width_start: float
    width_end: float


def solve_spacing(
    *,
    length: float | None = None,
    n: int | None = None,
    ratio: float | None = None,
    cell_ratio: float | None = None,
    width_start: float | None = None,
    width_end: float | None = None,
) -> SpacingSolution:
    """Solve the common blockMesh grading relations from any useful subset.

    This is the backend equivalent of the blockMesh grading calculator: once
    three independent quantities are supplied, the remaining quantities are
    derived where possible.
    """

    values: dict[str, float | int] = {}
    if length is not None:
        values["length"] = _positive(length, "length")
    if n is not None:
        if int(n) < 1:
            raise ValueError("n must be >= 1")
        values["n"] = int(n)
    if ratio is not None:
        values["ratio"] = _positive(ratio, "ratio")
    if cell_ratio is not None:
        values["cell_ratio"] = _positive(cell_ratio, "cell_ratio")
    if width_start is not None:
        values["width_start"] = _positive(width_start, "width_start")
    if width_end is not None:
        values["width_end"] = _positive(width_end, "width_end")

    for _ in range(24):
        old_count = len(values)
        _derive(values)
        if len(values) == old_count:
            break

    required = {"length", "n", "ratio", "cell_ratio", "width_start", "width_end"}
    missing = required - values.keys()
    if missing:
        raise ValueError(f"not enough independent spacing parameters: missing {sorted(missing)}")

    return SpacingSolution(
        n=int(values["n"]),
        length=float(values["length"]),
        ratio=float(values["ratio"]),
        cell_ratio=float(values["cell_ratio"]),
        width_start=float(values["width_start"]),
        width_end=float(values["width_end"]),
    )


def axis_grid_arrays(axis: dict, lmin: float, lmax: float, n: int) -> GridArrays:
    """Build SNaC grid arrays for one axis specification."""

    if n < 1:
        raise ValueError("axis grid must contain at least one cell")
    if lmax <= lmin:
        raise ValueError("axis lmax must be greater than lmin")

    kind = axis.get("kind", "snac")
    length = float(lmax) - float(lmin)

    if kind == "snac":
        faces = _snac_faces(
            n=n,
            lmin=float(lmin),
            lmax=float(lmax),
            gt=int(axis.get("gt", 0)),
            gr=float(axis.get("gr", 0.0)),
        )
    elif kind in {"multi", "simple_ratio"}:
        segments = _segments_from_axis(axis, length)
        widths = _multi_grading_widths(length, n, segments, profile=str(axis.get("profile", "geometric")))
        faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths)))
    elif kind == "max_min":
        segments = _segments_from_axis(axis, length)
        widths = _multi_grading_widths(length, n, segments, profile=str(axis.get("profile", "geometric")))
        faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths)))
    else:
        raise ValueError(f"unknown axis grid kind {kind!r}")

    return _arrays_from_faces(faces)


def spacing_from_widths(widths: Sequence[float], lmin: float) -> GridArrays:
    """Build SNaC grid arrays directly from cell widths."""

    if not widths:
        raise ValueError("widths may not be empty")
    widths_np = np.asarray(widths, dtype=float)
    if np.any(widths_np <= 0.0):
        raise ValueError("all widths must be positive")
    faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths_np)))
    return _arrays_from_faces(faces)


def _derive(values: dict[str, float | int]) -> None:
    n = int(values["n"]) if "n" in values else None
    ratio = float(values["ratio"]) if "ratio" in values else None
    cell_ratio = float(values["cell_ratio"]) if "cell_ratio" in values else None
    width_start = float(values["width_start"]) if "width_start" in values else None
    width_end = float(values["width_end"]) if "width_end" in values else None
    length = float(values["length"]) if "length" in values else None

    if n is not None and n > 1 and ratio is not None and "cell_ratio" not in values:
        values["cell_ratio"] = pow(ratio, 1.0 / (n - 1))
    if n is not None and n > 1 and cell_ratio is not None and "ratio" not in values:
        values["ratio"] = pow(cell_ratio, n - 1)
    if width_start is not None and width_end is not None and "ratio" not in values:
        values["ratio"] = width_end / width_start
    if width_start is not None and ratio is not None and "width_end" not in values:
        values["width_end"] = width_start * ratio
    if width_end is not None and ratio is not None and "width_start" not in values:
        values["width_start"] = width_end / ratio

    n = int(values["n"]) if "n" in values else None
    cell_ratio = float(values["cell_ratio"]) if "cell_ratio" in values else None
    width_start = float(values["width_start"]) if "width_start" in values else None
    width_end = float(values["width_end"]) if "width_end" in values else None
    length = float(values["length"]) if "length" in values else None

    if n is not None and cell_ratio is not None and width_start is not None and "length" not in values:
        values["length"] = _sum_geometric(width_start, cell_ratio, n)
    if n is not None and cell_ratio is not None and width_end is not None and "length" not in values:
        values["length"] = _sum_geometric(width_end, 1.0 / cell_ratio, n)
    if n is not None and cell_ratio is not None and length is not None and "width_start" not in values:
        values["width_start"] = _first_width_from_cell_ratio(length, n, cell_ratio)
    if n is not None and cell_ratio is not None and length is not None and "width_end" not in values:
        values["width_end"] = _first_width_from_cell_ratio(length, n, cell_ratio) * pow(cell_ratio, max(n - 1, 0))
    if n is not None and width_start is not None and length is not None and "cell_ratio" not in values:
        values["cell_ratio"] = _solve_cell_ratio_from_start(n, width_start, length)
    if n is not None and width_end is not None and length is not None and "cell_ratio" not in values:
        values["cell_ratio"] = _solve_cell_ratio_from_end(n, width_end, length)

    ratio = float(values["ratio"]) if "ratio" in values else None
    width_start = float(values["width_start"]) if "width_start" in values else None
    length = float(values["length"]) if "length" in values else None
    if ratio is not None and width_start is not None and length is not None and "n" not in values:
        values["n"] = _solve_n_from_ratio_start(ratio, width_start, length)


def _segments_from_axis(axis: dict, length: float) -> list[dict[str, float]]:
    kind = axis.get("kind", "snac")
    if kind == "multi":
        segments = axis.get("segments") or []
        if not segments:
            return [{"length": 1.0, "cells": 1.0, "ratio": 1.0}]
        return [
            {
                "length": _positive(seg.get("length", 1.0), "segment length weight"),
                "cells": _positive(seg.get("cells", 1.0), "segment cell weight"),
                "ratio": _positive(seg.get("ratio", 1.0), "segment ratio"),
            }
            for seg in segments
        ]
    if kind == "simple_ratio":
        ratio = _positive(axis.get("ratio", 1.0), "ratio")
        if axis.get("side", "end") == "start":
            ratio = 1.0 / ratio
        return [{"length": 1.0, "cells": 1.0, "ratio": ratio}]
    if kind == "max_min":
        min_width = _positive(axis.get("min", length), "min width")
        max_width = _positive(axis.get("max", length), "max width")
        if max_width < min_width:
            min_width, max_width = max_width, min_width
        ratio = max_width / min_width
        side = axis.get("side", "end")
        if side == "start":
            return [{"length": 1.0, "cells": 1.0, "ratio": 1.0 / ratio}]
        if side == "both":
            return [
                {"length": 0.5, "cells": 0.5, "ratio": ratio},
                {"length": 0.5, "cells": 0.5, "ratio": 1.0 / ratio},
            ]
        if side == "middle":
            return [
                {"length": 0.5, "cells": 0.5, "ratio": 1.0 / ratio},
                {"length": 0.5, "cells": 0.5, "ratio": ratio},
            ]
        return [{"length": 1.0, "cells": 1.0, "ratio": ratio}]
    raise ValueError(f"unsupported axis kind {kind!r}")


def _multi_grading_widths(length: float, n: int, segments: Sequence[dict[str, float]], *, profile: str = "geometric") -> np.ndarray:
    if n < len(segments):
        raise ValueError(f"axis has {n} cells but {len(segments)} grading segments")

    length_weights = np.asarray([seg["length"] for seg in segments], dtype=float)
    cell_weights = np.asarray([seg["cells"] for seg in segments], dtype=float)
    ratios = np.asarray([seg["ratio"] for seg in segments], dtype=float)

    segment_lengths = length * length_weights / length_weights.sum()
    segment_cells = _allocate_cells(n, cell_weights)

    widths: list[np.ndarray] = []
    for seg_length, seg_n, ratio in zip(segment_lengths, segment_cells, ratios, strict=True):
        widths.append(_profile_widths(float(seg_length), int(seg_n), float(ratio), profile))
    result = np.concatenate(widths)
    result[-1] += length - float(result.sum())
    return result


def _allocate_cells(n: int, weights: np.ndarray) -> np.ndarray:
    scaled = weights / weights.sum() * n
    cells = np.maximum(1, np.floor(scaled).astype(int))
    while int(cells.sum()) > n:
        candidates = np.where(cells > 1)[0]
        if candidates.size == 0:
            raise ValueError("cannot allocate at least one cell to each grading segment")
        idx = candidates[np.argmin(scaled[candidates] - np.floor(scaled[candidates]))]
        cells[idx] -= 1
    while int(cells.sum()) < n:
        remainder = scaled - np.floor(scaled)
        idx = int(np.argmax(remainder))
        cells[idx] += 1
        remainder[idx] = 0.0
    return cells


def _geometric_widths(length: float, n: int, ratio: float) -> np.ndarray:
    if n == 1:
        return np.array([length], dtype=float)
    cell_ratio = pow(ratio, 1.0 / (n - 1))
    first = _first_width_from_cell_ratio(length, n, cell_ratio)
    powers = np.power(cell_ratio, np.arange(n, dtype=float))
    widths = first * powers
    widths[-1] += length - float(widths.sum())
    return widths


def _profile_widths(length: float, n: int, ratio: float, profile: str) -> np.ndarray:
    profile_key = profile.lower()
    if profile_key == "geometric":
        return _geometric_widths(length, n, ratio)
    if profile_key in {"tanh", "erf"}:
        return _mapped_ratio_widths(length, n, ratio, profile_key)
    raise ValueError(f"unknown grid profile {profile!r}")


def _mapped_ratio_widths(length: float, n: int, ratio: float, profile: str) -> np.ndarray:
    if n == 1 or isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        return np.full(n, length / n, dtype=float)
    target = max(ratio, 1.0 / ratio)
    alpha = _solve_mapping_alpha(n, target, profile)
    r0 = np.arange(0, n + 1, dtype=float) / float(n)
    faces = np.array([_map_one_end(alpha, value, profile) for value in r0], dtype=float)
    widths = np.diff(faces) * length
    if ratio < 1.0:
        widths = widths[::-1]
    widths[-1] += length - float(widths.sum())
    return widths


def _solve_mapping_alpha(n: int, ratio: float, profile: str) -> float:
    if ratio <= 1.0:
        return 0.0

    def ratio_for(alpha: float) -> float:
        r0 = np.arange(0, n + 1, dtype=float) / float(n)
        faces = np.array([_map_one_end(alpha, value, profile) for value in r0], dtype=float)
        widths = np.diff(faces)
        return float(widths[-1] / widths[0])

    high = 1.0
    while ratio_for(high) < ratio and high < 80.0:
        high *= 2.0
    if high >= 80.0 and ratio_for(high) < ratio:
        raise ValueError(f"{profile} profile cannot reach spacing ratio {ratio:g} with {n} cells")
    return _bisect(lambda alpha: ratio_for(alpha) - ratio, 0.0, high)


def _map_one_end(alpha: float, r0: float, profile: str) -> float:
    if alpha == 0.0:
        return r0
    if profile == "tanh":
        return 1.0 + tanh((r0 - 1.0) * alpha) / tanh(alpha)
    return 1.0 + erf((r0 - 1.0) * alpha) / erf(alpha)


def _snac_faces(n: int, lmin: float, lmax: float, gt: int, gr: float) -> np.ndarray:
    mapper = _SNAC_MAPPERS.get(gt, _gridpoint_cluster_two_end)
    r0 = np.arange(1, n + 1, dtype=float) / float(n)
    mapped = np.array([mapper(gr, value) for value in r0], dtype=float)
    faces = np.concatenate(([lmin], lmin + mapped * (lmax - lmin)))
    faces[-1] = lmax
    return faces


def _arrays_from_faces(faces: np.ndarray) -> GridArrays:
    if faces.ndim != 1 or faces.size < 2:
        raise ValueError("faces must be a one-dimensional array with at least two entries")
    widths = np.diff(faces)
    if np.any(widths <= 0.0):
        raise ValueError("grid faces must be strictly increasing")

    n = widths.size
    face_spacing = np.empty(n + 2, dtype=float)
    center_spacing = np.empty(n + 2, dtype=float)
    faces_out = np.empty(n + 2, dtype=float)
    centers_out = np.empty(n + 2, dtype=float)

    face_spacing[0] = widths[0]
    face_spacing[1 : n + 1] = widths
    face_spacing[n + 1] = widths[-1]

    center_spacing[0] = widths[0]
    if n > 1:
        center_spacing[1:n] = 0.5 * (widths[:-1] + widths[1:])
    center_spacing[n] = widths[-1]
    center_spacing[n + 1] = widths[-1]

    faces_out[0] = faces[0]
    faces_out[1 : n + 1] = faces[1:]
    faces_out[n + 1] = faces[-1] + widths[-1]

    centers_out[0] = faces[0] - 0.5 * widths[0]
    centers_out[1 : n + 1] = 0.5 * (faces[:-1] + faces[1:])
    centers_out[n + 1] = centers_out[n] + widths[-1]

    return GridArrays(
        faces=faces_out,
        centers=centers_out,
        face_spacing=face_spacing,
        center_spacing=center_spacing,
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


def _sum_geometric(width_start: float, cell_ratio: float, n: int) -> float:
    if n < 1:
        return 0.0
    if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        return width_start * n
    return width_start * (1.0 - pow(cell_ratio, n)) / (1.0 - cell_ratio)


def _first_width_from_cell_ratio(length: float, n: int, cell_ratio: float) -> float:
    if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        return length / n
    return length * (1.0 - cell_ratio) / (1.0 - pow(cell_ratio, n))


def _solve_cell_ratio_from_start(n: int, width_start: float, length: float) -> float:
    if isclose(n * width_start, length, rel_tol=_REL_TOL, abs_tol=_REL_TOL * length):
        return 1.0
    if n * width_start < length:
        c_min = pow(1.0 + _REL_TOL, 1.0 / max(n - 1, 1))
        c_max = pow(_RMAX, 1.0 / max(n - 1, 1))
    else:
        c_min = pow(1.0 / _RMAX, 1.0 / max(n - 1, 1))
        c_max = pow(1.0 - _REL_TOL, 1.0 / max(n - 1, 1))
    return _bisect(lambda c: _sum_geometric(width_start, c, n) - length, c_min, c_max)


def _solve_cell_ratio_from_end(n: int, width_end: float, length: float) -> float:
    if isclose(n * width_end, length, rel_tol=_REL_TOL, abs_tol=_REL_TOL * length):
        return 1.0
    if n * width_end > length:
        c_min = pow(1.0 + _REL_TOL, 1.0 / max(n - 1, 1))
        c_max = pow(_RMAX, 1.0 / max(n - 1, 1))
    else:
        c_min = pow(1.0 / _RMAX, 1.0 / max(n - 1, 1))
        c_max = pow(1.0 - _REL_TOL, 1.0 / max(n - 1, 1))
    return _bisect(lambda c: _sum_geometric(width_end, 1.0 / c, n) - length, c_min, c_max)


def _solve_n_from_ratio_start(ratio: float, width_start: float, length: float) -> int:
    if isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        return max(1, round(length / width_start))

    def f(n_float: float) -> float:
        if n_float <= 1.0:
            return width_start - length
        c = pow(ratio, 1.0 / (n_float - 1.0))
        return width_start * (1.0 - pow(c, n_float)) / (1.0 - c) - length

    d_min = min(width_start, width_start * ratio)
    high = max(2.0, length / d_min + 1.0)
    n_float = _bisect(f, 1.000001, high)
    return max(1, round(n_float))


def _bisect(f: Callable[[float], float], left: float, right: float, *, max_steps: int = 500) -> float:
    f_left = f(left)
    f_right = f(right)
    if f_left == 0.0:
        return left
    if f_right == 0.0:
        return right
    if f_left * f_right > 0.0:
        raise ValueError("root is not bracketed by the supplied interval")
    for _ in range(max_steps):
        mid = 0.5 * (left + right)
        if abs(right - left) < _REL_TOL:
            return mid
        f_mid = f(mid)
        if f_left * f_mid <= 0.0:
            right = mid
            f_right = f_mid
        else:
            left = mid
            f_left = f_mid
    raise ValueError("root finding did not converge")


def _positive(value: float | int | str, name: str) -> float:
    value_f = float(value)
    if value_f <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value_f
