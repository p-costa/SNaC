"""Interactive multi-block grid generation helpers for SNaC."""

from .export import export_project
from .grid import (
    AXIS_NAMES,
    GRID_FUNCTIONS,
    GridArrays,
    SpacingSolution,
    axis_grid_arrays,
    solve_spacing,
)
from .model import AxisSpec, Block, GradingSegment, Project
from .validation import CheckResult, check_project, infer_project_connectivity, update_project_structure

__all__ = [
    "AXIS_NAMES",
    "GRID_FUNCTIONS",
    "AxisSpec",
    "Block",
    "GradingSegment",
    "GridArrays",
    "Project",
    "SpacingSolution",
    "CheckResult",
    "axis_grid_arrays",
    "check_project",
    "export_project",
    "infer_project_connectivity",
    "solve_spacing",
    "update_project_structure",
]
