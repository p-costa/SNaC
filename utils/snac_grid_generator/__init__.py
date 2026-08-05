"""Interactive multi-block grid generation helpers for SNaC."""

from .case_import import ImportResult, import_case
from .export import export_project
from .geometry import align_blocks, duplicate_blocks, mirror_blocks, snap_block_to_face
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
from .repair import GridRepair, RepairResult, repair_project_grids
from .snac_bc import apply_bc_preset
from .validation import (
    CheckResult,
    apply_axis_to_aligned_blocks,
    check_project,
    infer_project_connectivity,
    interface_spacing_metrics,
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
    "ImportResult",
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
    "apply_bc_preset",
    "align_blocks",
    "check_project",
    "export_project",
    "duplicate_blocks",
    "fit_min_max_cell_count",
    "fit_monotone_spacing",
    "infer_project_connectivity",
    "interface_spacing_metrics",
    "import_case",
    "optimize_project_decomposition",
    "mirror_blocks",
    "repair_project_grids",
    "solve_spacing",
    "snap_block_to_face",
    "update_project_structure",
]
