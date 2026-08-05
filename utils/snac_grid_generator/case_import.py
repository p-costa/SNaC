"""Import existing SNaC cases into grid-generator projects."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from .grid import AXIS_NAMES
from .model import AxisSpec, Project
from .numeric import coordinates_close, geometry_tolerance
from .snac_grid import read_grid_binary
from .topology import axis_extent


@dataclass
class ImportResult:
    """An imported project and any non-fatal format observations."""

    project: Project
    case_dir: Path
    source: str
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ok": True,
            "project": self.project.to_dict(),
            "caseDir": str(self.case_dir),
            "source": self.source,
            "warnings": self.warnings,
        }


_ASSIGNMENT_START_RE = re.compile(r"\b([A-Za-z]\w*)\s*(?:\(([^()]*)\))?\s*=")
_VALUE_RE = re.compile(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|[^,\s]+")


def import_case(case_dir: str | Path) -> ImportResult:
    """Load project JSON or reconstruct a project from ``blocks.nml``."""

    case_path = Path(case_dir).expanduser().resolve()
    if case_path.is_file():
        if case_path.name == "blocks.nml":
            case_path = case_path.parent
        else:
            raise ValueError("case import expects a case directory or blocks.nml")
    if not case_path.is_dir():
        raise ValueError(f"case directory does not exist: {case_path}")

    project_path = case_path / "snac_grid_project.json"
    blocks_path = case_path / "blocks.nml"
    if project_path.is_file():
        project = Project.from_dict(json.loads(project_path.read_text(encoding="utf-8")))
        source = project_path.name
    elif blocks_path.is_file():
        project = _project_from_blocks_nml(blocks_path)
        source = blocks_path.name
    else:
        raise ValueError(f"no snac_grid_project.json or blocks.nml found in {case_path}")

    warnings: list[str] = []
    imported_sources = _load_external_grids(project, case_path, warnings)
    if imported_sources:
        project.write_external_grid = True
        project.external_grid_source = (
            "both" if imported_sources == {"grid", "data"} else next(iter(imported_sources))
        )
    project.name = project.name or case_path.name
    project.validate()
    return ImportResult(project=project, case_dir=case_path, source=source, warnings=warnings)


def _project_from_blocks_nml(path: Path) -> Project:
    assignments = _parse_namelist(path.read_text(encoding="utf-8"))
    nblocks = int(_values(assignments, "nblocks", (), required=True)[0])
    scalar_indices = [
        indices[-2]
        for (name, indices) in assignments
        if name in {"block_cbcscal", "block_bcscal"} and len(indices) >= 2
    ]
    nscal = max(scalar_indices, default=0)
    blocks: list[dict[str, Any]] = []
    for block_id in range(1, nblocks + 1):
        dims = _vector(assignments, "block_dims", block_id, [1, 1, 1], int)
        ng = _vector(assignments, "block_ng", block_id, [32, 32, 2], int)
        lmin = _vector(assignments, "block_lmin", block_id, [0.0, 0.0, 0.0], float)
        lmax = _vector(assignments, "block_lmax", block_id, [1.0, 1.0, 0.05], float)
        gt = _vector(assignments, "block_gt", block_id, [0, 0, 0], int)
        gr = _vector(assignments, "block_gr", block_id, [0.0, 0.0, 0.0], float)
        blocks.append(
            {
                "id": block_id,
                "name": f"block-{block_id}",
                "dims": dims,
                "ng": ng,
                "lmin": lmin,
                "lmax": lmax,
                "axes": {
                    axis: AxisSpec(kind="snac", gt=gt[index], gr=gr[index]).to_dict()
                    for index, axis in enumerate(AXIS_NAMES)
                },
                "cbcvel": [
                    _indexed_vector(assignments, "block_cbcvel", (component, block_id), ["D"] * 6, str)
                    for component in range(1, 4)
                ],
                "cbcpre": _indexed_vector(assignments, "block_cbcpre", (block_id,), ["N"] * 6, str),
                "bcvel": [
                    _indexed_vector(assignments, "block_bcvel", (component, block_id), [0.0] * 6, float)
                    for component in range(1, 4)
                ],
                "bcpre": _indexed_vector(assignments, "block_bcpre", (block_id,), [0.0] * 6, float),
                "cbcscal": [
                    _indexed_vector(assignments, "block_cbcscal", (scalar, block_id), ["N"] * 6, str)
                    for scalar in range(1, nscal + 1)
                ],
                "bcscal": [
                    _indexed_vector(assignments, "block_bcscal", (scalar, block_id), [0.0] * 6, float)
                    for scalar in range(1, nscal + 1)
                ],
                "inflow": _indexed_vector(
                    assignments, "block_inflow_type", (block_id,), [0] * 6, int
                ),
                "inivel": str(
                    _values(assignments, "block_inivel", (block_id,), default=["zer"])[0]
                ),
            }
        )

    project = Project.from_dict(
        {
            "name": path.parent.name,
            "nscal": nscal,
            "blocks": blocks,
            "inferConnectivity": False,
            "writeExternalGrid": False,
        }
    )
    project.periodic_axes = _infer_periodic_axes(project)
    project.validate()
    return project


def _parse_namelist(text: str) -> dict[tuple[str, tuple[int, ...]], list[Any]]:
    clean = _strip_comments(text)
    group = re.search(r"(?is)&blocks\b(.*?)(?:^\s*/\s*$)", clean, re.MULTILINE)
    if group is None:
        raise ValueError("blocks.nml does not contain a complete &blocks namelist")
    body = group.group(1)
    masked = _mask_quoted_text(body)
    matches = list(_ASSIGNMENT_START_RE.finditer(masked))
    assignments: dict[tuple[str, tuple[int, ...]], list[Any]] = {}
    for position, match in enumerate(matches):
        name = match.group(1).lower()
        indices = tuple(
            int(item.strip())
            for item in (match.group(2) or "").split(",")
            if item.strip() and ":" not in item
        )
        end = matches[position + 1].start() if position + 1 < len(matches) else len(body)
        value_text = body[match.end() : end].strip().rstrip(",").rstrip()
        values = _parse_values(value_text)
        if not values:
            raise ValueError(f"blocks.nml assignment {name} has no value")
        assignments[(name, indices)] = values
    if ("nblocks", ()) not in assignments:
        raise ValueError("blocks.nml does not define nblocks")
    return assignments


def _mask_quoted_text(text: str) -> str:
    """Mask quoted values while retaining assignment offsets."""

    masked = list(text)
    quote: str | None = None
    index = 0
    while index < len(text):
        character = text[index]
        if quote is None:
            if character in {"'", '"'}:
                quote = character
                masked[index] = " "
        else:
            masked[index] = " "
            if character == quote:
                if index + 1 < len(text) and text[index + 1] == quote:
                    masked[index + 1] = " "
                    index += 1
                else:
                    quote = None
        index += 1
    if quote is not None:
        raise ValueError("blocks.nml contains an unterminated quoted value")
    return "".join(masked)


def _strip_comments(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        quote: str | None = None
        result: list[str] = []
        for character in line:
            if character in {"'", '"'}:
                quote = None if quote == character else character if quote is None else quote
            if character == "!" and quote is None:
                break
            result.append(character)
        lines.append("".join(result))
    return "\n".join(lines)


def _parse_values(text: str) -> list[Any]:
    values: list[Any] = []
    for token in _VALUE_RE.findall(text):
        repeated = re.fullmatch(r"(\d+)\*(.+)", token)
        if repeated:
            values.extend([_parse_scalar(repeated.group(2))] * int(repeated.group(1)))
        else:
            values.append(_parse_scalar(token))
    return values


def _parse_scalar(token: str) -> Any:
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] and token[0] in {"'", '"'}:
        return token[1:-1].replace(token[0] * 2, token[0])
    lowered = token.lower()
    if lowered in {".true.", "true"}:
        return True
    if lowered in {".false.", "false"}:
        return False
    number = token.replace("D", "e").replace("d", "e")
    try:
        return float(number) if any(character in number for character in ".eE") else int(number)
    except ValueError:
        return token


def _values(
    assignments: dict[tuple[str, tuple[int, ...]], list[Any]],
    name: str,
    indices: tuple[int, ...],
    *,
    default: list[Any] | None = None,
    required: bool = False,
) -> list[Any]:
    values = assignments.get((name, indices))
    if values is not None:
        return values
    if required:
        suffix = "" if not indices else f" at indices {indices}"
        raise ValueError(f"blocks.nml does not define {name}{suffix}")
    return list(default or [])


def _vector(
    assignments: dict[tuple[str, tuple[int, ...]], list[Any]],
    name: str,
    block_id: int,
    default: list[Any],
    converter,
) -> list[Any]:
    return _indexed_vector(assignments, name, (block_id,), default, converter)


def _indexed_vector(
    assignments: dict[tuple[str, tuple[int, ...]], list[Any]],
    name: str,
    indices: tuple[int, ...],
    default: list[Any],
    converter,
) -> list[Any]:
    values = _values(assignments, name, indices, default=default)
    result = [converter(value) for value in values[: len(default)]]
    result.extend(default[len(result) :])
    return result


def _load_external_grids(project: Project, case_dir: Path, warnings: list[str]) -> set[str]:
    sources: set[str] = set()
    for block in project.blocks:
        for axis_index, axis in enumerate(AXIS_NAMES):
            candidates = [
                (source, case_dir / source / f"grid_{axis}_b_{block.id:03d}.bin")
                for source in ("grid", "data")
            ]
            available = [(source, path) for source, path in candidates if path.is_file()]
            if not available:
                continue
            decoded = []
            failures = []
            for source, path in available:
                try:
                    grid = read_grid_binary(
                        path,
                        block.ng[axis_index],
                        block.lmin[axis_index],
                        block.lmax[axis_index],
                    )
                except (OSError, ValueError) as exc:
                    failures.append((source, path, str(exc)))
                    continue
                decoded.append((source, path, grid))
                sources.add(source)
                if grid.dtype != "<f8":
                    warnings.append(
                        f"block {block.id} {axis}: imported {source}/ {grid.dtype} grid; "
                        "export uses little-endian float64"
                    )

            if not decoded:
                details = "; ".join(f"{source}/: {message}" for source, _, message in failures)
                raise ValueError(f"block {block.id} {axis}: no valid external grid ({details})")
            for source, _, message in failures:
                warnings.append(
                    f"block {block.id} {axis}: ignored invalid {source}/ grid ({message})"
                )

            source, _, selected = decoded[0]
            block.axes[axis] = AxisSpec(kind="explicit", faces=selected.faces.tolist())
            for other_source, _, other in decoded[1:]:
                tolerance = geometry_tolerance(
                    block.lmin[axis_index],
                    block.lmax[axis_index],
                    scale=axis_extent(block, axis_index),
                )
                if not np.allclose(selected.faces, other.faces, rtol=0.0, atol=tolerance):
                    warnings.append(
                        f"block {block.id} {axis}: {source}/ and {other_source}/ grids differ; "
                        f"imported {source}/"
                    )
    return sources


def _infer_periodic_axes(project: Project) -> list[bool]:
    if not project.blocks:
        return [False, False, False]
    by_id = {block.id: block for block in project.blocks}
    global_min = [min(block.lmin[axis] for block in project.blocks) for axis in range(3)]
    global_max = [max(block.lmax[axis] for block in project.blocks) for axis in range(3)]
    result: list[bool] = []
    for axis in range(3):
        low_face = 2 * axis
        high_face = low_face + 1
        low_blocks = [
            block
            for block in project.blocks
            if coordinates_close(block.lmin[axis], global_min[axis], scale=axis_extent(block, axis))
        ]
        high_blocks = [
            block
            for block in project.blocks
            if coordinates_close(block.lmax[axis], global_max[axis], scale=axis_extent(block, axis))
        ]
        matched_high: set[int] = set()
        valid = bool(low_blocks and high_blocks)
        for low in low_blocks:
            if low.cbcpre[low_face] != "F":
                valid = False
                break
            friend_id = int(round(low.bcpre[low_face]))
            high = by_id.get(friend_id)
            if high is None or high not in high_blocks or not _same_transverse_extent(low, high, axis):
                valid = False
                break
            if high.cbcpre[high_face] != "F" or int(round(high.bcpre[high_face])) != low.id:
                valid = False
                break
            matched_high.add(high.id)
        result.append(valid and matched_high == {block.id for block in high_blocks})
    return result


def _same_transverse_extent(a, b, normal_axis: int) -> bool:
    for axis in range(3):
        if axis == normal_axis:
            continue
        scale = max(axis_extent(a, axis), axis_extent(b, axis))
        if not (
            coordinates_close(a.lmin[axis], b.lmin[axis], scale=scale)
            and coordinates_close(a.lmax[axis], b.lmax[axis], scale=scale)
        ):
            return False
    return True
