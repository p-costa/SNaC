#!/usr/bin/env python3
import math
import re
from pathlib import Path


log = Path("run.log").read_text()
values = [
    float(value)
    for value in re.findall(r"Maximum divergence\s*=\s*([+\-0-9.Ee]+)", log)
]
assert len(values) == 1, f"expected one divergence check, found {len(values)}"
assert all(math.isfinite(value) for value in values), values
assert max(values) < 1.0e-8, values
assert "*** Fim ***" in log
