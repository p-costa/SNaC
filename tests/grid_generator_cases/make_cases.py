#!/usr/bin/env python
"""Write SNaC cases with the grid generator for the end-to-end solver test.

Usage: python make_cases.py <directory>
"""

import re
import shutil
import sys
from pathlib import Path

from snac_grid_generator import (
    Project,
    apply_balance_alternative,
    export_project,
    optimize_project_decomposition,
    update_project_structure,
)


def two_block_project(z_axis, x_axis=None, fft_axis=None):
    """Two blocks side by side in x, periodic in z, with generated grids."""

    blocks = []
    for index in range(2):
        axes = {"z": z_axis, "y": {"kind": "simple_ratio", "profile": "tanh", "ratio": 3.0}}
        if x_axis is not None:
            axes["x"] = x_axis
        blocks.append(
            {
                "id": index + 1,
                "ng": [8, 8, 8],
                "lmin": [float(index), 0.0, 0.0],
                "lmax": [float(index + 1), 1.0, 1.0],
                "axes": axes,
            }
        )
    return Project.from_dict({"periodicAxes": [False, False, True], "fftAxis": fft_axis, "blocks": blocks})


def balanced_step(ranks):
    """A backward-facing step whose cell counts were adjusted to balance *ranks* ranks."""

    project = Project.from_dict(
        {
            "name": "balanced-step",
            "fftAxis": "z",
            "periodicAxes": [False, False, True],
            "decomposition": {"targetRanks": ranks},
            "blocks": [
                {"id": 1, "name": "inlet", "ng": [32, 24, 16], "lmin": [0, 1, 0], "lmax": [2, 2, 2]},
                {"id": 2, "name": "upper", "ng": [96, 24, 16], "lmin": [2, 1, 0], "lmax": [10, 2, 2]},
                {"id": 3, "name": "lower", "ng": [96, 24, 16], "lmin": [2, 0, 0], "lmax": [10, 1, 2]},
            ],
        }
    )
    project, _ = update_project_structure(project)
    project, result = optimize_project_decomposition(project, alternatives=True)
    best = min(result.alternatives, key=lambda item: item["imbalance"])
    assert best["imbalance"] < result.imbalance - 0.1, (result.imbalance, best)
    return apply_balance_alternative(project, best)


def main(root):
    root = Path(root)
    uniform = {"kind": "simple_ratio", "ratio": 1.0}

    # A valid FFT case with generated, non-uniform x and y grids.
    export_project(two_block_project(uniform, {"kind": "simple_ratio", "ratio": 2.0}, fft_axis="z"), root / "valid")

    # A stretched z grid: SNaC's gr = 0 check cannot see it, so the solver must
    # reject the grid it loads. The generator would refuse with fftAxis = "z".
    export_project(two_block_project({"kind": "simple_ratio", "ratio": 4.0}), root / "stretched_fft")

    # Older generator versions wrote gt = 4 with gr = 0 next to a grid file for
    # that axis. The grid file defines the axis, so SNaC must accept the case.
    legacy = root / "legacy_geometric"
    export_project(two_block_project(uniform, {"kind": "snac", "gt": 0, "gr": 0.0}), legacy)
    namelist = legacy / "blocks.nml"
    namelist.write_text(re.sub(r"block_gt\(1:3,(\d)\) = 0, 0, 0", r"block_gt(1:3,\1) = 4, 0, 0", namelist.read_text()))

    # Without the grid file, initgrid would divide by gr = 0, so SNaC must stop.
    no_file = root / "geometric_without_file"
    shutil.copytree(legacy, no_file)
    for path in (no_file / "grid").glob("grid_x_*"):
        path.unlink()

    # Cell counts changed so that 9 ranks share the load almost evenly; SNaC
    # must run it on exactly 9 ranks with the generated grids.
    export_project(balanced_step(9), root / "balanced_step")


if __name__ == "__main__":
    main(sys.argv[1])
