"""Canonical names for axis grid specifications, and their legacy aliases.

Project files since schema 4 use these grading kinds and camelCase keys, like
the rest of the project JSON. Older spellings are still accepted everywhere an
axis specification is read, so existing files, scripts, and tests keep working.
"""

from __future__ import annotations

from typing import Any

NATIVE = "native"
MONOTONE = "monotone"
SYMMETRIC = "symmetric"
MULTI = "multi"
EXPLICIT = "explicit"
AXIS_KINDS = (NATIVE, MONOTONE, SYMMETRIC, MULTI, EXPLICIT)

LEGACY_KINDS = {"snac": NATIVE, "simple_ratio": MONOTONE, "max_min": SYMMETRIC}

# Monotone grading controls; the keys of the same name hold their values.
MONOTONE_CONTROLS = ("n", "ratio", "cellRatio", "widthStart", "widthEnd")
LEGACY_KEYS = {"cell_ratio": "cellRatio", "width_start": "widthStart", "width_end": "widthEnd"}


def canonical_kind(kind: Any) -> str:
    """Return the canonical grading kind for a current or legacy name."""

    text = str(kind if kind is not None else NATIVE)
    return LEGACY_KINDS.get(text, text)


def canonical_control(name: Any) -> str:
    """Return the canonical monotone control name for a current or legacy name."""

    text = str(name)
    return LEGACY_KEYS.get(text, text)


def canonical_axis(data: dict[str, Any] | None) -> dict[str, Any]:
    """Return a copy of an axis specification with canonical kind and keys.

    A canonical key wins over its legacy spelling when both are present.
    """

    if not data:
        return {}
    result: dict[str, Any] = {}
    for key, value in data.items():
        if key in LEGACY_KEYS:
            continue
        result[key] = value
    for legacy, key in LEGACY_KEYS.items():
        if legacy in data and key not in data:
            result[key] = data[legacy]
    if "kind" in data:
        result["kind"] = canonical_kind(data["kind"])
    controls = data.get("controls")
    if isinstance(controls, (list, tuple)):
        result["controls"] = [canonical_control(name) for name in controls]
    return result
