"""Project-file migrations for the grid generator."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from .axis_format import canonical_axis

PROJECT_SCHEMA_VERSION = 4


def migrate_project_dict(data: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of *data* migrated to the current project schema."""

    migrated = deepcopy(data)
    version = int(migrated.get("schemaVersion", 1))
    if version < 1 or version > PROJECT_SCHEMA_VERSION:
        raise ValueError(f"unsupported project schema version {version}")

    while version < PROJECT_SCHEMA_VERSION:
        if version == 1:
            decomposition = migrated.setdefault("decomposition", {})
            # The default when schema 2 introduced the setting; later default
            # changes must not alter how old files migrate.
            decomposition.setdefault("minLocalCells", 4)
            decomposition.setdefault("maxLocalAspect", 0.0)
        elif version == 2:
            migrated.setdefault("fftAxis", None)
        elif version == 3:
            # Canonical grading kinds and camelCase axis keys.
            for block in migrated.get("blocks") or []:
                axes = block.get("axes") or {}
                for name, axis in list(axes.items()):
                    if isinstance(axis, dict):
                        axes[name] = canonical_axis(axis)
        version += 1

    migrated["schemaVersion"] = PROJECT_SCHEMA_VERSION
    return migrated
