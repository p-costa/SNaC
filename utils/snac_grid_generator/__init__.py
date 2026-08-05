"""Interactive multi-block grid generation helpers for SNaC."""

from .case_import import ImportResult, import_case
from .diagnostics import Diagnostic, DiagnosticLocation
from .export import export_project
from .decomposition import (
    MIN_LOCAL_CELLS,
    BlockDecomposition,
    DecompositionResult,
    optimize_project_decomposition,
)
from .geometry import (
    align_blocks,
    duplicate_blocks,
    mirror_blocks,
    snap_block_to_face,
)
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
from .model import (
    PROJECT_SCHEMA_VERSION,
    AxisSpec,
    Block,
    DecompositionSpec,
    GradingSegment,
    Project,
)
from .repair import GridRepair, RepairResult, repair_project_grids
from .snac_bc import apply_bc_preset
from .snac_grid import BinaryGrid, read_grid_binary
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
    "BinaryGrid",
    "CheckResult",
    "DecompositionResult",
    "DecompositionSpec",
    "Diagnostic",
    "DiagnosticLocation",
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
    "align_blocks",
    "apply_axis_to_aligned_blocks",
    "apply_bc_preset",
    "axis_grid_arrays",
    "axis_grid_diagnostics",
    "check_project",
    "duplicate_blocks",
    "export_project",
    "fit_min_max_cell_count",
    "fit_monotone_spacing",
    "infer_project_connectivity",
    "interface_spacing_metrics",
    "import_case",
    "mirror_blocks",
    "optimize_project_decomposition",
    "repair_project_grids",
    "read_grid_binary",
    "solve_spacing",
    "snap_block_to_face",
    "update_project_structure",
]
