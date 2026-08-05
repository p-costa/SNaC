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
from .model import PROJECT_SCHEMA_VERSION, AxisSpec, Block, DecompositionSpec, GradingSegment, Project
from .decomposition import (
    MIN_LOCAL_CELLS,
    BlockDecomposition,
    DecompositionResult,
    optimize_project_decomposition,
)
from .validation import (
    CheckResult,
    GridRepair,
    RepairResult,
    apply_axis_to_aligned_blocks,
    check_project,
    infer_project_connectivity,
    repair_project_grids,
    update_project_structure,
)

__all__ = [
    "AXIS_NAMES",
    "GRID_FUNCTIONS",
    "AxisSpec",
    "Block",
    "BlockDecomposition",
    "DecompositionResult",
    "DecompositionSpec",
    "GradingSegment",
    "GridArrays",
    "GridRepair",
    "MIN_LOCAL_CELLS",
    "MinMaxFit",
    "MonotoneFit",
    "PROJECT_SCHEMA_VERSION",
    "Project",
    "RepairResult",
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
    "optimize_project_decomposition",
    "repair_project_grids",
    "solve_spacing",
    "update_project_structure",
]
