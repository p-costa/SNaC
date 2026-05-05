"""Writers for SNaC block namelists and external grid binaries."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .grid import AXIS_NAMES, axis_grid_arrays
from .model import Block, Project
from .validation import check_project, infer_project_connectivity

@dataclass
class ExportResult:
    output_dir: Path
    files: list[Path]
    warnings: list[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "outputDir": str(self.output_dir),
            "files": [str(path) for path in self.files],
            "warnings": self.warnings,
        }


def export_project(project: Project | dict[str, Any], output_dir: str | Path) -> ExportResult:
    """Write ``blocks.nml`` and optional external grid files."""

    if isinstance(project, dict):
        project = Project.from_dict(project)
    else:
        project = Project.from_dict(project.to_dict())

    warnings: list[str] = []
    if project.infer_connectivity:
        project, connectivity_check = infer_project_connectivity(project)
        if connectivity_check.errors:
            raise ValueError("; ".join(connectivity_check.errors))
        warnings.extend(connectivity_check.warnings)

    check = check_project(project)
    if check.errors:
        raise ValueError("; ".join(check.errors))
    warnings.extend(check.warnings)

    output_path = Path(output_dir).expanduser().resolve()
    grid_path = output_path / "grid"
    data_path = output_path / "data"
    output_path.mkdir(parents=True, exist_ok=True)
    if project.write_external_grid:
        grid_path.mkdir(parents=True, exist_ok=True)
        data_path.mkdir(parents=True, exist_ok=True)

    blocks = _renumber_blocks(project.blocks)
    extents, extent_warnings = _global_index_extents(blocks)
    warnings.extend(extent_warnings)

    files: list[Path] = []
    blocks_file = output_path / "blocks.nml"
    blocks_file.write_text(_format_blocks_nml(blocks, project.nscal), encoding="utf-8")
    files.append(blocks_file)

    for block in sorted(blocks, key=lambda item: item.id):
        if project.write_external_grid:
            for grid_dir in _grid_dirs(project.external_grid_source, grid_path, data_path):
                files.extend(_write_block_grids(block, grid_dir))
            files.append(_write_geometry(block, data_path, extents[block.id]))

    project_file = output_path / "snac_grid_project.json"
    project_file.write_text(json.dumps(project.to_dict(), indent=2) + "\n", encoding="utf-8")
    files.append(project_file)

    return ExportResult(output_dir=output_path, files=files, warnings=warnings)


def _renumber_blocks(blocks: list[Block]) -> list[Block]:
    id_map = {block.id: index for index, block in enumerate(sorted(blocks, key=lambda item: item.id), start=1)}
    renumbered = [Block.from_dict(block.to_dict()) for block in sorted(blocks, key=lambda item: item.id)]
    for index, block in enumerate(renumbered, start=1):
        block.id = index
        for face, cbc in enumerate(block.cbcpre):
            if cbc == "F":
                block.bcpre[face] = float(id_map.get(int(round(block.bcpre[face])), int(round(block.bcpre[face]))))
        for component in range(3):
            for face, cbc in enumerate(block.cbcvel[component]):
                if cbc == "F":
                    block.bcvel[component][face] = float(
                        id_map.get(int(round(block.bcvel[component][face])), int(round(block.bcvel[component][face])))
                    )
        for iscal in range(len(block.cbcscal)):
            for face, cbc in enumerate(block.cbcscal[iscal]):
                if cbc == "F":
                    block.bcscal[iscal][face] = float(
                        id_map.get(int(round(block.bcscal[iscal][face])), int(round(block.bcscal[iscal][face])))
                    )
    return renumbered


def _grid_dirs(source: str, grid_path: Path, data_path: Path) -> list[Path]:
    if source == "data":
        return [data_path]
    if source == "both":
        return [grid_path, data_path]
    return [grid_path]


def _write_block_grids(block: Block, grid_dir: Path) -> list[Path]:
    files: list[Path] = []
    suffix = f"_b_{block.id:03d}"
    for idx, axis in enumerate(AXIS_NAMES):
        arrays = axis_grid_arrays(block.axes[axis].to_dict(), block.lmin[idx], block.lmax[idx], block.ng[idx])
        base = grid_dir / f"grid_{axis}{suffix}"
        arrays.binary_payload.tofile(base.with_suffix(".bin"))
        _write_grid_out(base.with_suffix(".out"), arrays)
        files.append(base.with_suffix(".bin"))
        files.append(base.with_suffix(".out"))
    return files


def _write_grid_out(path: Path, arrays) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for values in zip(
            arrays.faces,
            arrays.centers,
            arrays.face_spacing,
            arrays.center_spacing,
        ):
            handle.write(f"{0.0:15.7E}{values[0]:15.7E}{values[1]:15.7E}{values[2]:15.7E}{values[3]:15.7E}\n")


def _write_geometry(block: Block, data_path: Path, extent: tuple[list[int], list[int]]) -> Path:
    path = data_path / f"geometry_b_{block.id:03d}.out"
    lo, hi = extent
    with path.open("w", encoding="utf-8") as handle:
        handle.write(_numbers(lo, integer=True) + "\n")
        handle.write(_numbers(hi, integer=True) + "\n")
        handle.write(_numbers(block.lmin) + "\n")
        handle.write(_numbers(block.lmax) + "\n")
    return path


def _format_blocks_nml(blocks: list[Block], nscal: int = 0) -> str:
    lines = ["&blocks", f"nblocks = {len(blocks)}"]
    for block in sorted(blocks, key=lambda item: item.id):
        gt = [block.axes[axis].gt if block.axes[axis].kind == "snac" else 0 for axis in AXIS_NAMES]
        gr = [block.axes[axis].gr if block.axes[axis].kind == "snac" else 0.0 for axis in AXIS_NAMES]
        ib = block.id
        lines.extend(
            [
                f"block_dims(1:3,{ib}) = {_nml_numbers(block.dims, integer=True)}",
                f"block_ng(1:3,{ib}) = {_nml_numbers(block.ng, integer=True)}",
                f"block_lmin(1:3,{ib}) = {_nml_numbers(block.lmin)}",
                f"block_lmax(1:3,{ib}) = {_nml_numbers(block.lmax)}",
                f"block_gt(1:3,{ib}) = {_nml_numbers(gt, integer=True)}",
                f"block_gr(1:3,{ib}) = {_nml_numbers(gr)}",
            ]
        )
        for component in range(3):
            lines.append(f"block_cbcvel(0:1,1:3,{component + 1},{ib}) = {_nml_strings(block.cbcvel[component])}")
        lines.append(f"block_cbcpre(0:1,1:3,{ib}) = {_nml_strings(block.cbcpre)}")
        for component in range(3):
            lines.append(f"block_bcvel(0:1,1:3,{component + 1},{ib}) = {_nml_numbers(block.bcvel[component])}")
        lines.append(f"block_bcpre(0:1,1:3,{ib}) = {_nml_numbers(block.bcpre)}")
        for iscal in range(nscal):
            lines.append(f"block_cbcscal(0:1,1:3,{iscal + 1},{ib}) = {_nml_strings(block.cbcscal[iscal])}")
        for iscal in range(nscal):
            lines.append(f"block_bcscal(0:1,1:3,{iscal + 1},{ib}) = {_nml_numbers(block.bcscal[iscal])}")
        lines.append(f"block_inflow_type(0:1,1:3,{ib}) = {_nml_numbers(block.inflow, integer=True)}")
        lines.append(f"block_inivel({ib}) = '{block.inivel}'")
        lines.append("")
    lines.append("/")
    return "\n".join(lines) + "\n"


def _global_index_extents(blocks: list[Block]) -> tuple[dict[int, tuple[list[int], list[int]]], list[str]]:
    warnings: list[str] = []
    tol = 1.0e-10
    by_id = {block.id: block for block in blocks}
    lo_by_id: dict[int, list[int]] = {}
    for seed in sorted(blocks, key=lambda item: item.id):
        if seed.id in lo_by_id:
            continue
        lo_by_id[seed.id] = [1, 1, 1]
        queue = [seed]
        while queue:
            current = queue.pop(0)
            current_lo = lo_by_id[current.id]
            for other in blocks:
                if other.id == current.id:
                    continue
                candidate = _neighbor_lo(current, other, current_lo, tol)
                if candidate is None:
                    continue
                if other.id in lo_by_id:
                    if lo_by_id[other.id] != candidate:
                        warnings.append(f"block {other.id} has conflicting inferred global indices")
                    continue
                lo_by_id[other.id] = candidate
                queue.append(other)
    extents = {}
    for block_id, block in by_id.items():
        lo = lo_by_id[block_id]
        hi = [lo[idx] + block.ng[idx] - 1 for idx in range(3)]
        extents[block_id] = (lo, hi)
    return extents, warnings


def _neighbor_lo(current: Block, other: Block, current_lo: list[int], tol: float) -> list[int] | None:
    for axis_index in range(3):
        if _touches(current.lmax[axis_index], other.lmin[axis_index], tol) and _same_cross_section(current, other, axis_index, tol):
            candidate = list(current_lo)
            candidate[axis_index] = current_lo[axis_index] + current.ng[axis_index]
            return candidate
        if _touches(other.lmax[axis_index], current.lmin[axis_index], tol) and _same_cross_section(current, other, axis_index, tol):
            candidate = list(current_lo)
            candidate[axis_index] = current_lo[axis_index] - other.ng[axis_index]
            return candidate
    return None


def _same_cross_section(a: Block, b: Block, axis_index: int, tol: float) -> bool:
    for idx in range(3):
        if idx == axis_index:
            continue
        if not _touches(a.lmin[idx], b.lmin[idx], tol):
            return False
        if not _touches(a.lmax[idx], b.lmax[idx], tol):
            return False
    return True


def _touches(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def _numbers(values: list[int] | list[float], *, integer: bool = False) -> str:
    if integer:
        return " ".join(str(int(value)) for value in values)
    return " ".join(_float_token(float(value)) for value in values)


def _nml_numbers(values: list[int] | list[float], *, integer: bool = False) -> str:
    if integer:
        return ", ".join(str(int(value)) for value in values)
    return ", ".join(_float_token(float(value)) for value in values)


def _nml_strings(values: list[str]) -> str:
    return ", ".join(f"'{str(value)[:1]}'" for value in values)


def _float_token(value: float) -> str:
    if value == 0.0:
        return "0."
    text = f"{value:.10g}"
    if "e" not in text.lower() and "." not in text:
        text += "."
    return text
