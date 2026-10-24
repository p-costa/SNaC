"""Scale-aware floating-point comparisons shared by grid modules."""

from __future__ import annotations

from math import isfinite
from sys import float_info

GEOMETRY_REL_TOL = 1.0e-10
_ROUND_OFF_FACTOR = 8.0


def geometry_tolerance(*coordinates: float, scale: float) -> float:
    """Return a tolerance based on local extent and coordinate round-off."""

    local_scale = abs(float(scale))
    if not isfinite(local_scale) or local_scale <= 0.0:
        raise ValueError("geometry comparison scale must be finite and positive")
    magnitude = max(
        (abs(float(value)) for value in coordinates),
        default=local_scale,
    )
    round_off = (
        _ROUND_OFF_FACTOR * float_info.epsilon * max(magnitude, local_scale)
    )
    return max(GEOMETRY_REL_TOL * local_scale, round_off)


def coordinates_close(a: float, b: float, *, scale: float) -> bool:
    """Compare coordinates using a local geometric scale."""

    return abs(float(a) - float(b)) <= geometry_tolerance(a, b, scale=scale)
