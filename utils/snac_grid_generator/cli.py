"""Command-line workflows for SNaC grid projects."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Callable

from .case_import import import_case
from .decomposition import optimize_project_decomposition
from .export import export_project
from .model import Project
from .repair import repair_project_grids
from .validation import check_project, update_project_structure

_WORKFLOW_COMMANDS = {"check", "decompose", "migrate", "repair", "update"}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        if argv and argv[0] == "import":
            return _import_main(argv[1:])
        if argv and argv[0] in _WORKFLOW_COMMANDS:
            return _workflow_main(argv[0], argv[1:])
        if argv and argv[0] == "export":
            argv = argv[1:]
        return _export_main(argv)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        if "--json" in argv:
            print(json.dumps({"ok": False, "errors": [str(exc)], "warnings": []}, indent=2))
        else:
            print(f"error: {exc}", file=sys.stderr)
        return 2


def _export_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Export SNaC multi-block files from project JSON.")
    parser.add_argument("project", type=Path, help="Path to snac_grid_project.json")
    parser.add_argument("-o", "--output", type=Path, default=Path("."), help="Output case directory")
    args = parser.parse_args(argv)

    result = export_project(_read_project(args.project), args.output)
    print(f"Wrote {len(result.files)} files to {result.output_dir}")
    for warning in result.warnings:
        print(f"warning: {warning}")
    return 0


def _import_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Import an existing SNaC case as project JSON.")
    parser.add_argument("case", type=Path, help="Case directory or blocks.nml path")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("snac_grid_project.json"),
        help="Output project JSON path",
    )
    args = parser.parse_args(argv)
    result = import_case(args.case)
    _write_project(args.output, result.project)
    print(f"Imported {len(result.project.blocks)} blocks from {result.case_dir} to {args.output}")
    for warning in result.warnings:
        print(f"warning: {warning}")
    return 0


def _workflow_main(command: str, argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=f"{command.capitalize()} a grid project.")
    parser.add_argument("project", type=Path, help="Path to snac_grid_project.json")
    parser.add_argument("--json", action="store_true", help="Print machine-readable diagnostics")
    if command != "check":
        parser.add_argument("-o", "--output", type=Path, help="Write the resulting project JSON")
    if command in {"repair", "update"}:
        parser.add_argument("--source-block", type=int, help="Preferred source block id")
    args = parser.parse_args(argv)

    project = _read_project(args.project)
    if command == "check":
        result = check_project(project)
        payload = result.to_dict()
    elif command == "migrate":
        result = None
        payload = {"ok": True, "project": project.to_dict()}
    else:
        operation: Callable[..., tuple[Project, Any]]
        kwargs: dict[str, Any] = {}
        if command == "update":
            operation = update_project_structure
            kwargs["source_block_id"] = args.source_block
        elif command == "repair":
            operation = repair_project_grids
            kwargs["source_block_id"] = args.source_block
        else:
            operation = optimize_project_decomposition
        project, result = operation(project, **kwargs)
        payload = {**result.to_dict(), "project": project.to_dict()}

    if command != "check" and args.output and payload.get("ok", False):
        _write_project(args.output, project)
        payload["output"] = str(args.output)
    _print_diagnostics(command, payload, args.json)
    return 0 if payload.get("ok", False) else 1


def _read_project(path: Path) -> Project:
    return Project.from_dict(json.loads(path.read_text(encoding="utf-8")))


def _write_project(path: Path, project: Project) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(project.to_dict(), indent=2) + "\n", encoding="utf-8")


def _print_diagnostics(command: str, payload: dict[str, Any], as_json: bool) -> None:
    if as_json:
        print(json.dumps(payload, indent=2))
        return
    state = "passed" if payload.get("ok", False) else "failed"
    print(f"{command.capitalize()} {state}.")
    for warning in payload.get("warnings", []):
        print(f"warning: {warning}")
    for error in payload.get("errors", []):
        print(f"error: {error}")
    if payload.get("output"):
        print(f"Wrote {payload['output']}")


if __name__ == "__main__":
    raise SystemExit(main())
