#!/usr/bin/env python
"""Check that the grids SNaC saved in data/ match the generator's files in grid/.

Usage: python check_grid.py <case directory>
"""

import sys
from pathlib import Path

import numpy as np


def main(case):
    case = Path(case)
    written = sorted((case / "grid").glob("grid_*_b_*.bin"))
    assert written, "no generated grid files found"
    for path in written:
        generated = np.fromfile(path, dtype="<f8")
        saved = np.fromfile(case / "data" / path.name, dtype=np.float64)
        np.testing.assert_allclose(saved, generated, rtol=1e-12, atol=1e-14, err_msg=path.name)
    print(f"{len(written)} saved grids match the generated grids")


if __name__ == "__main__":
    main(sys.argv[1])
