"""Command-line entry point for SNaC grid project export."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .export import export_project


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export SNaC multi-block grid files from a project JSON file.")
    parser.add_argument("project", type=Path, help="Path to snac_grid_project.json")
    parser.add_argument("-o", "--output", type=Path, default=Path("."), help="Output case directory")
    args = parser.parse_args(argv)

    project = json.loads(args.project.read_text(encoding="utf-8"))
    result = export_project(project, args.output)
    print(f"Wrote {len(result.files)} files to {result.output_dir}")
    for warning in result.warnings:
        print(f"warning: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
