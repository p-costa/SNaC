#!/usr/bin/env python
"""Write SNaC cases with the grid generator for the end-to-end solver test.

Usage: python make_cases.py <directory>
"""

import re
import shutil
import sys
from pathlib import Path

from snac_grid_generator import Project, export_project


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


if __name__ == "__main__":
    main(sys.argv[1])
