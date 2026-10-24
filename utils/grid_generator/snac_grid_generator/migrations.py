"""Project-file migrations for the grid generator."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

PROJECT_SCHEMA_VERSION = 2


def migrate_project_dict(data: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of *data* migrated to the current project schema."""

    migrated = deepcopy(data)
    version = int(migrated.get("schemaVersion", 1))
    if version < 1 or version > PROJECT_SCHEMA_VERSION:
        raise ValueError(f"unsupported project schema version {version}")

    while version < PROJECT_SCHEMA_VERSION:
        if version == 1:
            decomposition = migrated.setdefault("decomposition", {})
            decomposition.setdefault("minLocalCells", 4)
            decomposition.setdefault("maxLocalAspect", 0.0)
        version += 1

    migrated["schemaVersion"] = PROJECT_SCHEMA_VERSION
    return migrated
