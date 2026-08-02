#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path


def main() -> None:
    log = Path(sys.argv[1]).read_text()
    values = [
        float(value)
        for value in re.findall(r"Maximum divergence\s*=\s*([+\-0-9.Ee]+)", log)
    ]
    assert len(values) == 3, f"expected three divergence checks, found {len(values)}"
    assert all(math.isfinite(value) for value in values), values
    assert max(values) < 1.0e-5, values
    assert "Aborting" not in log


if __name__ == "__main__":
    main()
