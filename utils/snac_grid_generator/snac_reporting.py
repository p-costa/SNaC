"""SNaC case reports built from format-neutral grid-quality metrics."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from .model import PROJECT_SCHEMA_VERSION, Project
from .quality import GridQualityReport, analyze_grid_quality
from .validation import check_project, copy_project

REPORT_SCHEMA_VERSION = 1


def build_case_report(project: Project | dict[str, Any]) -> dict[str, Any]:
    """Return deterministic validation, quality, and SNaC storage data."""

    project = copy_project(project)
    validation = check_project(project)
    quality: GridQualityReport | None = None
    quality_error: str | None = None
    try:
        quality = analyze_grid_quality(project)
    except (ValueError, ArithmeticError) as exc:
        quality_error = str(exc)

    binary_bytes_per_copy = 32 * sum(sum(block.ng) for block in project.blocks)
    copies = (
        0
        if not project.write_external_grid
        else 2 if project.external_grid_source == "both" else 1
    )
    report = {
        "schemaVersion": REPORT_SCHEMA_VERSION,
        "project": {
            "name": project.name,
            "projectSchemaVersion": PROJECT_SCHEMA_VERSION,
            "blocks": len(project.blocks),
            "scalars": project.nscal,
        },
        "validation": validation.to_dict(),
        "quality": quality.to_dict() if quality is not None else None,
        "storage": {
            "snacBinaryBytesPerCopy": binary_bytes_per_copy,
            "snacBinaryCopies": copies,
            "snacBinaryBytesTotal": binary_bytes_per_copy * copies,
        },
    }
    if quality_error is not None:
        report["qualityError"] = quality_error
    return report


def render_case_report_markdown(report: dict[str, Any]) -> str:
    """Render one deterministic, human-readable case report."""

    project = report["project"]
    validation = report["validation"]
    quality = report.get("quality")
    storage = report["storage"]
    status = "PASS" if validation.get("ok") else "FAIL"
    lines = [
        f"# SNaC grid report: {_markdown_text(project.get('name') or 'unnamed')}",
        "",
        f"**Validation:** {status}",
        "",
    ]

    if quality is not None:
        summary = quality["summary"]
        lines.extend(
            [
                "## Summary",
                "",
                "| Metric | Value |",
                "| --- | ---: |",
                f"| Blocks | {summary['blockCount']} |",
                f"| Cells | {summary['totalCells']} |",
                f"| MPI ranks | {summary['totalRanks']} |",
                f"| Cells per rank | {summary['cellsPerRankMin']} to {summary['cellsPerRankMax']} |",
                f"| Load imbalance | {_number(summary['loadImbalance'])} |",
                f"| Spacing range | {_number(summary['spacingMin'])} to {_number(summary['spacingMax'])} |",
                f"| Maximum adjacent-cell ratio | {_number(summary['maxAdjacentRatio'])} |",
                f"| Worst cell aspect ratio | {_number(summary['worstCellAspect'])} |",
                f"| Maximum interface ratio | {_number(summary['maxInterfaceRatio'])} |",
                f"| Rectilinear coordinate storage | {_bytes(summary['coordinateBytes'])} |",
                "",
                "## Blocks",
                "",
                "| Block | Cells | Ranks | Cells/rank | Spacing | Adjacent ratio | Cell aspect |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for block in quality["blocks"]:
            label = f"{block['blockId']}: {_markdown_text(block.get('name') or '')}".rstrip()
            lines.append(
                f"| {label} | {block['cells']} | {block['ranks']} | "
                f"{block['cellsPerRankMin']} to {block['cellsPerRankMax']} | "
                f"{_number(block['spacingMin'])} to {_number(block['spacingMax'])} | "
                f"{_number(block['maxAdjacentRatio'])} | {_number(block['worstCellAspect'])} |"
            )

        lines.extend(
            [
                "",
                "## Axis Quality",
                "",
                "| Block | Axis | Cells | Length | Spacing | Endpoint ratio | Adjacent ratio |",
                "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for block in quality["blocks"]:
            for axis in block["axes"]:
                lines.append(
                    f"| {block['blockId']} | {axis['axis'].upper()} | {axis['cells']} | "
                    f"{_number(axis['length'])} | {_number(axis['spacingMin'])} to "
                    f"{_number(axis['spacingMax'])} | {_number(axis['endpointRatio'])} | "
                    f"{_number(axis['maxAdjacentRatio'])} |"
                )

        lines.extend(
            [
                "",
                "## Interfaces",
                "",
                "| Face pair | Axis | Spacing | Neighbor spacing | Ratio |",
                "| --- | --- | ---: | ---: | ---: |",
            ]
        )
        if quality["interfaces"]:
            for interface in quality["interfaces"]:
                pair = (
                    f"B{interface['blockId']} {interface['face']} / "
                    f"B{interface['neighborId']} {interface['neighborFace']}"
                )
                lines.append(
                    f"| {pair} | {interface['axis'].upper()} | {_number(interface['spacing'])} | "
                    f"{_number(interface['neighborSpacing'])} | {_number(interface['ratio'])} |"
                )
        else:
            lines.append("| None | - | - | - | - |")
    else:
        lines.extend(["## Quality", "", f"Unavailable: {_markdown_text(report.get('qualityError', 'unknown error'))}"])

    lines.extend(["", "## Validation", ""])
    _append_messages(lines, "Errors", validation.get("errors", []))
    _append_messages(lines, "Warnings", validation.get("warnings", []))
    if not validation.get("errors") and not validation.get("warnings"):
        lines.append("No errors or warnings.")

    lines.extend(
        [
            "",
            "## SNaC Storage",
            "",
            "| Metric | Value |",
            "| --- | ---: |",
            f"| Binary grid bytes per copy | {storage['snacBinaryBytesPerCopy']} |",
            f"| Binary grid copies | {storage['snacBinaryCopies']} |",
            f"| Total binary grid bytes | {storage['snacBinaryBytesTotal']} |",
            "",
        ]
    )
    return "\n".join(lines)


def write_case_report(
    project: Project | dict[str, Any],
    path: str | Path,
    report_format: str | None = None,
) -> Path:
    """Write a JSON or Markdown report with single-file atomic replacement."""

    path = Path(path).expanduser().resolve()
    report_format = (report_format or path.suffix.lstrip(".") or "json").lower()
    if report_format in {"markdown", "md"}:
        content = render_case_report_markdown(build_case_report(project))
    elif report_format == "json":
        content = json.dumps(build_case_report(project), indent=2) + "\n"
    else:
        raise ValueError("report format must be json or markdown")

    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".snac_report_", dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        staged.write_text(content, encoding="utf-8")
        staged.replace(path)
    return path


def _append_messages(lines: list[str], label: str, messages: list[str]) -> None:
    if not messages:
        return
    lines.extend([f"### {label}", ""])
    lines.extend(f"- {_markdown_text(message)}" for message in messages)
    lines.append("")


def _number(value: Any) -> str:
    return f"{float(value):.6g}"


def _bytes(value: int) -> str:
    size = float(value)
    units = ("B", "KiB", "MiB", "GiB")
    unit = units[0]
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            break
        size /= 1024.0
    return f"{size:.4g} {unit}"


def _markdown_text(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")
