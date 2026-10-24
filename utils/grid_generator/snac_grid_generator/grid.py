"""Reusable one-dimensional grid spacing and grading helpers."""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, exp, isclose, log, pow, tanh
from typing import Callable, Sequence

import numpy as np

from .numeric import geometry_tolerance
from .snac_grid import SNAC_GRID_FUNCTIONS, native_faces

AXIS_NAMES = ("x", "y", "z")

GRID_FUNCTIONS = SNAC_GRID_FUNCTIONS

_REL_TOL = 1.0e-11
_RMAX = 1.0e5


@dataclass(frozen=True)
class GridArrays:
    """One-dimensional face, center, and spacing arrays."""

    faces: np.ndarray
    centers: np.ndarray
    face_spacing: np.ndarray
    center_spacing: np.ndarray


@dataclass(frozen=True)
class SpacingSolution:
    """Derived values for a geometric spacing sequence."""

    n: int
    length: float
    ratio: float
    cell_ratio: float
    width_start: float
    width_end: float


@dataclass(frozen=True)
class MinMaxFit:
    """Integer cell count and achieved spacings for min/max targets."""

    n: int
    width_min: float
    width_max: float


@dataclass(frozen=True)
class MonotoneFit:
    """Resolved one-way grading controls and their achieved endpoint widths."""

    n: int
    ratio: float
    cell_ratio: float | None
    width_start: float
    width_end: float
    ideal_n: float | None


MONOTONE_CONTROLS = ("n", "ratio", "cell_ratio", "width_start", "width_end")


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

    n_value = int(values["n"])
    length_value = float(values["length"])
    cell_ratio_value = float(values["cell_ratio"])
    width_start_value = _first_width_from_cell_ratio(length_value, n_value, cell_ratio_value)
    width_end_value = width_start_value * pow(cell_ratio_value, max(n_value - 1, 0))
    ratio_value = width_end_value / width_start_value

    return SpacingSolution(
        n=n_value,
        length=length_value,
        ratio=ratio_value,
        cell_ratio=cell_ratio_value,
        width_start=width_start_value,
        width_end=width_end_value,
    )


def fit_monotone_spacing(axis: dict, length: float, n: int) -> MonotoneFit:
    """Resolve two user-selected controls for a monotone spacing profile."""

    length = _positive(length, "axis length")
    profile = str(axis.get("profile", "geometric")).lower()
    raw_controls = axis.get("controls")
    legacy = not isinstance(raw_controls, (list, tuple))
    controls = tuple(dict.fromkeys(str(value) for value in (raw_controls or ("n", "ratio"))))
    if len(controls) != 2 or any(value not in MONOTONE_CONTROLS for value in controls):
        raise ValueError("monotone grading requires two distinct supported controls")
    if profile not in {"geometric", "tanh", "erf"}:
        raise ValueError(f"unknown grid profile {profile!r}")
    if profile != "geometric" and "cell_ratio" in controls:
        raise ValueError("cell-to-cell ratio is available only for geometric grading")

    ratio_value = _positive(axis.get("ratio", 1.0), "ratio")
    if legacy and str(axis.get("side", "end")) == "start":
        ratio_value = 1.0 / ratio_value
    supplied = {
        "n": int(n),
        "ratio": ratio_value,
        "cell_ratio": _positive(axis.get("cell_ratio", 1.0), "cell ratio"),
        "width_start": _positive(axis.get("width_start", length / max(n, 1)), "lower spacing"),
        "width_end": _positive(axis.get("width_end", length / max(n, 1)), "upper spacing"),
    }

    if profile == "geometric":
        kwargs: dict[str, float | int] = {"length": length}
        for control in controls:
            kwargs[control] = supplied[control]
        solution = solve_spacing(**kwargs)
        ideal_n = _ideal_geometric_count(length, controls, supplied, solution)
        return MonotoneFit(
            n=solution.n,
            ratio=solution.ratio,
            cell_ratio=solution.cell_ratio,
            width_start=solution.width_start,
            width_end=solution.width_end,
            ideal_n=ideal_n,
        )

    if "n" in controls:
        n_value = int(supplied["n"])
        other = next(value for value in controls if value != "n")
        if other == "ratio":
            ratio = float(supplied["ratio"])
        else:
            ratio = _solve_profile_ratio_from_endpoint(
                length,
                n_value,
                float(supplied[other]),
                profile,
                endpoint=other,
            )
        ideal_n: float | None = float(n_value)
    else:
        if set(controls) == {"width_start", "width_end"}:
            ratio = float(supplied["width_end"]) / float(supplied["width_start"])
        elif "ratio" in controls:
            ratio = float(supplied["ratio"])
        else:
            raise ValueError("this profile requires cells, endpoint ratio, or both endpoint spacings")
        targets = {control: float(supplied[control]) for control in controls if control.startswith("width_")}
        n_value, ideal_n = _fit_profile_cell_count(length, ratio, profile, targets)

    widths = _profile_widths(length, n_value, ratio, profile)
    return MonotoneFit(
        n=n_value,
        ratio=float(widths[-1] / widths[0]),
        cell_ratio=None,
        width_start=float(widths[0]),
        width_end=float(widths[-1]),
        ideal_n=ideal_n,
    )


def fit_min_max_cell_count(axis: dict, length: float, *, max_cells: int = 1_000_000) -> MinMaxFit:
    """Fit an integer cell count to absolute min/max spacing targets.

    The requested spacing ratio and clustering profile are retained exactly.
    Since an integer number of cells generally cannot satisfy both absolute
    targets exactly, the closest common scaling of the two targets is used.
    """

    length = _positive(length, "axis length")
    target_min = _positive(axis.get("min", length), "min width")
    target_max = _positive(axis.get("max", length), "max width")
    if target_max < target_min:
        raise ValueError("max width must be greater than or equal to min width")

    estimated_max = int(length / target_min) + 2
    upper = max(2, min(max_cells, estimated_max))
    lower = 1

    while lower <= upper:
        mid = (lower + upper) // 2
        widths = _max_min_widths(length, mid, axis)
        if float(widths.min()) > target_min:
            lower = mid + 1
        else:
            upper = mid - 1

    if lower > max_cells:
        raise ValueError(f"min/max spacing requires more than {max_cells} cells")

    candidates = {
        max(1, min(max_cells, value))
        for value in range(max(1, lower - 4), min(max_cells, lower + 4) + 1)
    }
    best: tuple[float, int, float, float] | None = None
    for candidate in sorted(candidates):
        widths = _max_min_widths(length, candidate, axis)
        width_min = float(widths.min())
        width_max = float(widths.max())
        error = max(abs(_log_ratio(width_min, target_min)), abs(_log_ratio(width_max, target_max)))
        score = (error, candidate, width_min, width_max)
        if best is None or score < best:
            best = score

    if best is None:
        raise ValueError("could not fit min/max spacing targets")
    return MinMaxFit(n=best[1], width_min=best[2], width_max=best[3])


def axis_grid_arrays(axis: dict, lmin: float, lmax: float, n: int) -> GridArrays:
    """Build grid arrays for one axis specification."""

    if n < 1:
        raise ValueError("axis grid must contain at least one cell")
    if lmax <= lmin:
        raise ValueError("axis lmax must be greater than lmin")

    kind = axis.get("kind", "snac")
    length = float(lmax) - float(lmin)

    if kind == "snac":
        faces = native_faces(
            n=n,
            lmin=float(lmin),
            lmax=float(lmax),
            gt=int(axis.get("gt", 0)),
            gr=float(axis.get("gr", 0.0)),
        )
    elif kind == "multi":
        segments = _segments_from_axis(axis, length)
        widths = _multi_grading_widths(length, n, segments, profile=str(axis.get("profile", "geometric")))
        faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths)))
    elif kind == "simple_ratio":
        fit = fit_monotone_spacing(axis, length, n)
        widths = _profile_widths(length, fit.n, fit.ratio, str(axis.get("profile", "geometric")))
        faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths)))
    elif kind == "max_min":
        widths = _max_min_widths(length, n, axis)
        faces = np.concatenate(([float(lmin)], float(lmin) + np.cumsum(widths)))
    elif kind == "explicit":
        faces = np.asarray(axis.get("faces", []), dtype=float)
        if faces.size != n + 1:
            raise ValueError(f"explicit grid has {faces.size} faces, expected {n + 1}")
        if not np.all(np.isfinite(faces)) or np.any(np.diff(faces) <= 0.0):
            raise ValueError("explicit grid faces must be finite and strictly increasing")
        tolerance = geometry_tolerance(float(faces[0]), float(faces[-1]), lmin, lmax, scale=length)
        if abs(float(faces[0]) - float(lmin)) > tolerance or abs(float(faces[-1]) - float(lmax)) > tolerance:
            raise ValueError("explicit grid faces must span the axis extent")
        faces = faces.copy()
    else:
        raise ValueError(f"unknown axis grid kind {kind!r}")

    faces[-1] = float(lmax)
    return _arrays_from_faces(faces)


def axis_grid_diagnostics(axis: dict, lmin: float, lmax: float, n: int) -> dict:
    """Return achieved spacing metrics used by the GUI and validation."""

    arrays = axis_grid_arrays(axis, lmin, lmax, n)
    widths = arrays.face_spacing[1:-1]
    effective_n = int(widths.size)
    adjacent = widths[1:] / widths[:-1] if effective_n > 1 else np.ones(1, dtype=float)
    growth = np.maximum(adjacent, 1.0 / adjacent)
    result: dict[str, object] = {
        "n": effective_n,
        "lower": float(widths[0]),
        "upper": float(widths[-1]),
        "min": float(widths.min()),
        "max": float(widths.max()),
        "expansion": float(widths[-1] / widths[0]),
        "ratio": float(widths.max() / widths.min()),
        "maxGrowth": float(growth.max()),
        "cellRatio": None,
        "idealN": None,
        "residuals": {},
        "segments": [],
    }

    kind = str(axis.get("kind", "snac"))
    length = float(lmax) - float(lmin)
    if kind == "simple_ratio":
        fit = fit_monotone_spacing(axis, length, n)
        result["cellRatio"] = fit.cell_ratio
        result["idealN"] = fit.ideal_n
        achieved = {
            "n": float(fit.n),
            "ratio": fit.ratio,
            "cell_ratio": fit.cell_ratio,
            "width_start": fit.width_start,
            "width_end": fit.width_end,
        }
        controls = axis.get("controls") or ("n", "ratio")
        residuals: dict[str, float] = {}
        for control in controls:
            requested = float(n) if control == "n" else float(axis.get(control, achieved.get(control) or 1.0))
            value = achieved.get(control)
            if value is not None and requested > 0.0:
                residuals[str(control)] = float(value / requested - 1.0)
        result["residuals"] = residuals
    elif kind == "multi":
        segments = _segments_from_axis(axis, length)
        _, segment_data = _multi_grading_data(
            length,
            effective_n,
            segments,
            profile=str(axis.get("profile", "geometric")),
        )
        result["segments"] = segment_data

    return result


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

    if n == 1 and length is not None:
        if ratio is not None and not isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
            raise ValueError("n=1 requires a unit expansion ratio")
        if cell_ratio is not None and not isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
            raise ValueError("n=1 requires a unit cell ratio")
        if width_start is not None and not isclose(width_start, length, rel_tol=_REL_TOL, abs_tol=_REL_TOL):
            raise ValueError("n=1 requires width_start to equal length")
        if width_end is not None and not isclose(width_end, length, rel_tol=_REL_TOL, abs_tol=_REL_TOL):
            raise ValueError("n=1 requires width_end to equal length")
        values.update(ratio=1.0, cell_ratio=1.0, width_start=length, width_end=length)
        return

    if n is not None and n > 1 and ratio is not None and "cell_ratio" not in values:
        values["cell_ratio"] = pow(ratio, 1.0 / (n - 1))
    if n is not None and n > 1 and cell_ratio is not None and "ratio" not in values:
        values["ratio"] = pow(cell_ratio, n - 1)
    if n is None and ratio is not None and cell_ratio is not None:
        values["n"] = _solve_n_from_ratio_cell_ratio(ratio, cell_ratio)
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
    if n is None and cell_ratio is not None and width_start is not None and length is not None:
        values["n"] = _solve_n_from_cell_ratio_start(cell_ratio, width_start, length)
    if n is None and cell_ratio is not None and width_end is not None and length is not None:
        values["n"] = _solve_n_from_cell_ratio_end(cell_ratio, width_end, length)

    ratio = float(values["ratio"]) if "ratio" in values else None
    width_start = float(values["width_start"]) if "width_start" in values else None
    length = float(values["length"]) if "length" in values else None
    if ratio is not None and width_start is not None and length is not None and "n" not in values:
        values["n"] = _solve_n_from_ratio_start(ratio, width_start, length)


def _segments_from_axis(axis: dict, length: float) -> list[dict[str, float | bool]]:
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
                "continuous": bool(seg.get("continuous", False)),
            }
            for seg in segments
        ]
    if kind == "max_min":
        min_width = _positive(axis.get("min", length), "min width")
        max_width = _positive(axis.get("max", length), "max width")
        if max_width < min_width:
            raise ValueError("max width must be greater than or equal to min width")
        ratio = max_width / min_width
        side = axis.get("side", "end")
        if side == "start":
            return [{"length": 1.0, "cells": 1.0, "ratio": 1.0 / ratio}]
        return [{"length": 1.0, "cells": 1.0, "ratio": ratio}]
    raise ValueError(f"unsupported axis kind {kind!r}")


def _multi_grading_widths(
    length: float,
    n: int,
    segments: Sequence[dict[str, float | bool]],
    *,
    profile: str = "geometric",
) -> np.ndarray:
    result, _ = _multi_grading_data(length, n, segments, profile=profile)
    return result


def _multi_grading_data(
    length: float,
    n: int,
    segments: Sequence[dict[str, float | bool]],
    *,
    profile: str = "geometric",
) -> tuple[np.ndarray, list[dict[str, float | int | bool]]]:
    if n < len(segments):
        raise ValueError(f"axis has {n} cells but {len(segments)} grading segments")

    length_weights = np.asarray([seg["length"] for seg in segments], dtype=float)
    cell_weights = np.asarray([seg["cells"] for seg in segments], dtype=float)
    ratios = np.asarray([seg["ratio"] for seg in segments], dtype=float)

    segment_lengths = length * length_weights / length_weights.sum()
    segment_cells = _allocate_cells(n, cell_weights)

    widths: list[np.ndarray] = []
    diagnostics: list[dict[str, float | int | bool]] = []
    previous_end: float | None = None
    for index, (seg_length, seg_n, ratio) in enumerate(zip(segment_lengths, segment_cells, ratios)):
        linked = index > 0 and bool(segments[index].get("continuous", False))
        if linked and previous_end is not None:
            ratio = _solve_profile_ratio_from_endpoint(
                float(seg_length),
                int(seg_n),
                previous_end,
                profile,
                endpoint="width_start",
            )
        segment_widths = _profile_widths(float(seg_length), int(seg_n), float(ratio), profile)
        start = float(segment_widths[0])
        end = float(segment_widths[-1])
        jump = 1.0 if previous_end is None else max(start / previous_end, previous_end / start)
        widths.append(segment_widths)
        diagnostics.append(
            {
                "index": index,
                "length": float(seg_length),
                "cells": int(seg_n),
                "ratio": end / start,
                "widthStart": start,
                "widthEnd": end,
                "jump": jump,
                "continuous": linked,
            }
        )
        previous_end = end
    result = np.concatenate(widths)
    result[-1] += length - float(result.sum())
    if diagnostics:
        diagnostics[-1]["widthEnd"] = float(result[-1])
        diagnostics[-1]["ratio"] = float(result[-1]) / float(diagnostics[-1]["widthStart"])
    return result, diagnostics


def _allocate_cells(n: int, weights: np.ndarray) -> np.ndarray:
    scaled = weights / weights.sum() * n
    cells = np.maximum(1, np.floor(scaled).astype(int))
    while int(cells.sum()) > n:
        candidates = np.where(cells > 1)[0]
        if candidates.size == 0:
            raise ValueError("cannot allocate at least one cell to each grading segment")
        idx = candidates[np.argmax(cells[candidates] - scaled[candidates])]
        cells[idx] -= 1
    while int(cells.sum()) < n:
        idx = int(np.argmax(scaled - cells))
        cells[idx] += 1
    return cells


def _max_min_widths(length: float, n: int, axis: dict) -> np.ndarray:
    min_width = _positive(axis.get("min", length), "min width")
    max_width = _positive(axis.get("max", length), "max width")
    if max_width < min_width:
        raise ValueError("max width must be greater than or equal to min width")
    ratio = max_width / min_width
    side = str(axis.get("side", "end"))
    profile = str(axis.get("profile", "geometric"))

    if side not in {"end", "start", "both", "middle"}:
        raise ValueError(f"unknown min/max clustering side {side!r}")

    if side == "start":
        return _profile_widths(length, n, 1.0 / ratio, profile)
    if side == "both":
        return _symmetric_profile_widths(length, n, ratio, profile, cluster_middle=False)
    if side == "middle":
        return _symmetric_profile_widths(length, n, ratio, profile, cluster_middle=True)
    return _profile_widths(length, n, ratio, profile)


def _symmetric_profile_widths(
    length: float,
    n: int,
    ratio: float,
    profile: str,
    *,
    cluster_middle: bool,
) -> np.ndarray:
    half_n = (n + 1) // 2
    half = _profile_widths(1.0, half_n, ratio, profile)
    if cluster_middle:
        half = half[::-1]
    if n % 2:
        widths = np.concatenate((half, half[-2::-1]))
    else:
        widths = np.concatenate((half, half[::-1]))
    widths *= length / float(widths.sum())
    widths[-1] += length - float(widths.sum())
    return widths


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


def _solve_profile_ratio_from_endpoint(
    length: float,
    n: int,
    target: float,
    profile: str,
    *,
    endpoint: str,
) -> float:
    target = _positive(target, "endpoint spacing")
    if endpoint not in {"width_start", "width_end"}:
        raise ValueError(f"unknown spacing endpoint {endpoint!r}")
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        if not isclose(target, length, rel_tol=_REL_TOL, abs_tol=_REL_TOL * length):
            raise ValueError("one cell requires its spacing to equal the axis length")
        return 1.0

    index = 0 if endpoint == "width_start" else -1

    def error(log_ratio: float) -> float:
        widths = _profile_widths(length, n, exp(log_ratio), profile)
        return float(widths[index] - target)

    limit = log(_RMAX)
    try:
        root = _bisect(error, -limit, limit)
    except ValueError as exc:
        uniform = length / n
        raise ValueError(
            f"{endpoint.replace('_', ' ')} {target:g} cannot be reached with {n} cells "
            f"and {profile} grading (uniform spacing is {uniform:g})"
        ) from exc
    return exp(root)


def _fit_profile_cell_count(
    length: float,
    ratio: float,
    profile: str,
    targets: dict[str, float],
    *,
    max_cells: int = 1_000_000,
) -> tuple[int, float | None]:
    if not targets:
        raise ValueError("a spacing target is required when cell count is derived")
    ratio = _positive(ratio, "ratio")
    for key, value in targets.items():
        targets[key] = _positive(value, key.replace("_", " "))

    primary = "width_start" if "width_start" in targets else "width_end"
    index = 0 if primary == "width_start" else -1
    target = targets[primary]
    minimum_n = 1 if isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL) else 2

    def width_for(n_value: int) -> float:
        return float(_profile_widths(length, n_value, ratio, profile)[index])

    if width_for(minimum_n) < target:
        crossing = minimum_n
    else:
        high = max(minimum_n + 1, int(length / target) + 2)
        high = min(high, max_cells)
        while width_for(high) > target and high < max_cells:
            high = min(max_cells, high * 2)
        if high >= max_cells and width_for(high) > target:
            raise ValueError(f"spacing controls require more than {max_cells} cells")
        low = minimum_n
        while low < high:
            mid = (low + high) // 2
            if width_for(mid) > target:
                low = mid + 1
            else:
                high = mid
        crossing = low

    candidates = range(max(minimum_n, crossing - 4), min(max_cells, crossing + 4) + 1)

    def score(n_value: int) -> float:
        widths = _profile_widths(length, n_value, ratio, profile)
        errors = []
        if "width_start" in targets:
            errors.append(abs(_log_ratio(float(widths[0]), targets["width_start"])))
        if "width_end" in targets:
            errors.append(abs(_log_ratio(float(widths[-1]), targets["width_end"])))
        return max(errors)

    best = min(candidates, key=lambda value: (score(value), value))
    lower = max(minimum_n, crossing - 1)
    upper = max(lower + 1, crossing)
    lower_width = width_for(lower)
    upper_width = width_for(upper)
    if lower_width > 0.0 and upper_width > 0.0 and not isclose(lower_width, upper_width):
        fraction = log(target / lower_width) / log(upper_width / lower_width)
        ideal_n = float(lower + min(1.0, max(0.0, fraction)))
    else:
        ideal_n = float(best)
    return best, ideal_n


def _ideal_geometric_count(
    length: float,
    controls: tuple[str, ...],
    supplied: dict[str, float | int],
    solution: SpacingSolution,
) -> float | None:
    if "n" in controls:
        return float(solution.n)
    if set(controls) == {"ratio", "cell_ratio"}:
        ratio = float(supplied["ratio"])
        cell_ratio = float(supplied["cell_ratio"])
        if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
            return None
        return 1.0 + log(ratio) / log(cell_ratio)
    if "cell_ratio" in controls:
        cell_ratio = float(supplied["cell_ratio"])
        if "width_start" in controls:
            width = float(supplied["width_start"])
            if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
                return length / width
            argument = 1.0 - length / width * (1.0 - cell_ratio)
            return log(argument) / log(cell_ratio) if argument > 0.0 else None
        if "width_end" in controls:
            width = float(supplied["width_end"])
            if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
                return length / width
            argument = 1.0 / (1.0 + length / width * (1.0 - cell_ratio) / cell_ratio)
            return log(argument) / log(cell_ratio) if argument > 0.0 else None

    ratio = solution.ratio
    if "width_start" in controls:
        width_start = float(supplied["width_start"])
    elif "width_end" in controls:
        width_start = float(supplied["width_end"]) / ratio
    else:
        return None
    return _continuous_n_from_ratio_start(length, ratio, width_start)


def _continuous_n_from_ratio_start(length: float, ratio: float, width_start: float) -> float | None:
    if isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        return length / width_start

    def implied_length(n_value: float) -> float:
        cell_ratio = pow(ratio, 1.0 / (n_value - 1.0))
        if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
            return width_start * n_value
        return width_start * (1.0 - pow(cell_ratio, n_value)) / (1.0 - cell_ratio)

    low = 2.0
    if implied_length(low) >= length:
        return low
    high = max(4.0, length / min(width_start, width_start * ratio) + 2.0)
    while implied_length(high) < length and high < 1_000_000.0:
        high *= 2.0
    if high >= 1_000_000.0 and implied_length(high) < length:
        return None
    return _bisect(lambda value: implied_length(value) - length, low, high)


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
        if widths[0] <= 0.0:
            return float("inf")
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


def _arrays_from_faces(faces: np.ndarray) -> GridArrays:
    if faces.ndim != 1 or faces.size < 2:
        raise ValueError("faces must be a one-dimensional array with at least two entries")
    if not np.all(np.isfinite(faces)):
        raise ValueError("grid faces must be finite")
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

    def implied_length(n_value: int) -> float:
        if n_value <= 1:
            return width_start
        cell_ratio_value = pow(ratio, 1.0 / (n_value - 1))
        return _sum_geometric(width_start, cell_ratio_value, n_value)

    estimate = max(2, int(length / min(width_start, width_start * ratio)) + 2)
    return _closest_integer_count(length, implied_length, estimate)


def _solve_n_from_ratio_cell_ratio(ratio: float, cell_ratio: float) -> int:
    if isclose(ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL) and isclose(
        cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL
    ):
        raise ValueError("ratio and cell_ratio do not determine n for a uniform grid")
    if isclose(cell_ratio, 1.0, rel_tol=0.0, abs_tol=_REL_TOL):
        raise ValueError("cell_ratio=1 is incompatible with a non-unit expansion ratio")
    n_float = 1.0 + np.log(ratio) / np.log(cell_ratio)
    if not np.isfinite(n_float) or n_float < 1.0:
        raise ValueError("ratio and cell_ratio imply an invalid cell count")
    return max(1, round(float(n_float)))


def _solve_n_from_cell_ratio_start(cell_ratio: float, width_start: float, length: float) -> int:
    def implied_length(n_value: int) -> float:
        return _sum_geometric(width_start, cell_ratio, n_value)

    estimate = max(2, int(length / min(width_start, width_start * cell_ratio)) + 2)
    return _closest_integer_count(length, implied_length, estimate)


def _solve_n_from_cell_ratio_end(cell_ratio: float, width_end: float, length: float) -> int:
    def implied_length(n_value: int) -> float:
        return _sum_geometric(width_end, 1.0 / cell_ratio, n_value)

    estimate = max(2, int(length / min(width_end, width_end / cell_ratio)) + 2)
    return _closest_integer_count(length, implied_length, estimate)


def _closest_integer_count(length: float, implied_length: Callable[[int], float], estimate: int) -> int:
    def evaluate(n_value: int) -> float:
        try:
            value = float(implied_length(n_value))
        except OverflowError:
            return float("inf")
        return value if np.isfinite(value) else float("inf")

    high = max(2, estimate)
    while evaluate(high) < length and high < 1_000_000:
        high = min(1_000_000, high * 2)
    if high >= 1_000_000 and evaluate(high) < length:
        raise ValueError("spacing parameters require more than 1000000 cells")

    low = 1
    while low < high:
        mid = (low + high) // 2
        if evaluate(mid) < length:
            low = mid + 1
        else:
            high = mid
    candidates = range(max(1, low - 2), min(1_000_000, low + 2) + 1)
    return min(candidates, key=lambda n_value: abs(_log_ratio(evaluate(n_value), length)))


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
    if not np.isfinite(value_f) or value_f <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value_f


def _log_ratio(value: float, reference: float) -> float:
    return float(np.log(value / reference))
