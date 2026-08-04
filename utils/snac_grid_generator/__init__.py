"""Interactive multi-block grid generation helpers for SNaC."""

from .export import export_project
from .grid import (
    AXIS_NAMES,
    GRID_FUNCTIONS,
    GridArrays,
    MinMaxFit,
    MonotoneFit,
    SpacingSolution,
    axis_grid_arrays,
    axis_grid_diagnostics,
    fit_min_max_cell_count,
    fit_monotone_spacing,
    solve_spacing,
)
from .model import PROJECT_SCHEMA_VERSION, AxisSpec, Block, GradingSegment, Project
from .validation import (
    CheckResult,
    apply_axis_to_aligned_blocks,
    check_project,
    infer_project_connectivity,
    update_project_structure,
)

__all__ = [
    "AXIS_NAMES",
    "GRID_FUNCTIONS",
    "AxisSpec",
    "Block",
    "GradingSegment",
    "GridArrays",
    "MinMaxFit",
    "MonotoneFit",
    "PROJECT_SCHEMA_VERSION",
    "Project",
    "SpacingSolution",
    "CheckResult",
    "axis_grid_arrays",
    "axis_grid_diagnostics",
    "apply_axis_to_aligned_blocks",
    "check_project",
    "export_project",
    "fit_min_max_cell_count",
    "fit_monotone_spacing",
    "infer_project_connectivity",
    "solve_spacing",
    "update_project_structure",
]
