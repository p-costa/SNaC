"""Format-neutral quality metrics for rectilinear multi-block grids."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import prod
from typing import Any

import numpy as np

from .grid import AXIS_NAMES, axis_grid_arrays
from .model import Block, Project
from .topology import FACE_INFO, FACE_ORDER, axis_extent, build_topology


@dataclass(frozen=True)
class AxisQuality:
    """Spacing metrics for one block axis."""

    axis: str
    cells: int
    length: float
    spacing_min: float
    spacing_max: float
    spacing_mean: float
    endpoint_ratio: float
    max_adjacent_ratio: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "axis": self.axis,
            "cells": self.cells,
            "length": self.length,
            "spacingMin": self.spacing_min,
            "spacingMax": self.spacing_max,
            "spacingMean": self.spacing_mean,
            "endpointRatio": self.endpoint_ratio,
            "maxAdjacentRatio": self.max_adjacent_ratio,
        }


@dataclass(frozen=True)
class BlockQuality:
    """Quality and decomposition metrics for one block."""

    block_id: int
    name: str
    cells: int
    ranks: int
    cells_per_rank_min: int
    cells_per_rank_max: int
    spacing_min: float
    spacing_max: float
    max_adjacent_ratio: float
    worst_cell_aspect: float
    cell_volume_min: float
    cell_volume_max: float
    coordinate_bytes: int
    axes: tuple[AxisQuality, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "blockId": self.block_id,
            "name": self.name,
            "cells": self.cells,
            "ranks": self.ranks,
            "cellsPerRankMin": self.cells_per_rank_min,
            "cellsPerRankMax": self.cells_per_rank_max,
            "spacingMin": self.spacing_min,
            "spacingMax": self.spacing_max,
            "maxAdjacentRatio": self.max_adjacent_ratio,
            "worstCellAspect": self.worst_cell_aspect,
            "cellVolumeMin": self.cell_volume_min,
            "cellVolumeMax": self.cell_volume_max,
            "coordinateBytes": self.coordinate_bytes,
            "axes": [axis.to_dict() for axis in self.axes],
        }


@dataclass(frozen=True)
class InterfaceQuality:
    """Normal-spacing quality for a connected block face pair."""

    block_id: int
    face: str
    neighbor_id: int
    neighbor_face: str
    axis: str
    spacing: float
    neighbor_spacing: float
    ratio: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "blockId": self.block_id,
            "face": self.face,
            "neighborId": self.neighbor_id,
            "neighborFace": self.neighbor_face,
            "axis": self.axis,
            "spacing": self.spacing,
            "neighborSpacing": self.neighbor_spacing,
            "ratio": self.ratio,
        }


@dataclass(frozen=True)
class QualitySummary:
    """Project-wide grid and decomposition metrics."""

    block_count: int
    total_cells: int
    total_ranks: int
    cells_per_rank_min: int
    cells_per_rank_max: int
    mean_cells_per_rank: float
    load_imbalance: float
    spacing_min: float
    spacing_max: float
    max_adjacent_ratio: float
    worst_cell_aspect: float
    max_interface_ratio: float
    coordinate_bytes: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "blockCount": self.block_count,
            "totalCells": self.total_cells,
            "totalRanks": self.total_ranks,
            "cellsPerRankMin": self.cells_per_rank_min,
            "cellsPerRankMax": self.cells_per_rank_max,
            "meanCellsPerRank": self.mean_cells_per_rank,
            "loadImbalance": self.load_imbalance,
            "spacingMin": self.spacing_min,
            "spacingMax": self.spacing_max,
            "maxAdjacentRatio": self.max_adjacent_ratio,
            "worstCellAspect": self.worst_cell_aspect,
            "maxInterfaceRatio": self.max_interface_ratio,
            "coordinateBytes": self.coordinate_bytes,
        }


@dataclass(frozen=True)
class GridQualityReport:
    """Complete format-neutral grid-quality analysis."""

    summary: QualitySummary
    blocks: tuple[BlockQuality, ...] = field(default_factory=tuple)
    interfaces: tuple[InterfaceQuality, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "summary": self.summary.to_dict(),
            "blocks": [block.to_dict() for block in self.blocks],
            "interfaces": [interface.to_dict() for interface in self.interfaces],
        }


def analyze_grid_quality(project: Project | dict[str, Any]) -> GridQualityReport:
    """Compute deterministic grid and partition quality metrics."""

    project = (
        Project.from_dict(project.to_dict())
        if isinstance(project, Project)
        else Project.from_dict(project)
    )
    blocks = tuple(_block_quality(block) for block in sorted(project.blocks, key=lambda item: item.id))
    interfaces = _interface_quality(project)
    total_cells = sum(block.cells for block in blocks)
    total_ranks = sum(block.ranks for block in blocks)
    cells_per_rank_min = min((block.cells_per_rank_min for block in blocks), default=0)
    cells_per_rank_max = max((block.cells_per_rank_max for block in blocks), default=0)
    summary = QualitySummary(
        block_count=len(blocks),
        total_cells=total_cells,
        total_ranks=total_ranks,
        cells_per_rank_min=cells_per_rank_min,
        cells_per_rank_max=cells_per_rank_max,
        mean_cells_per_rank=total_cells / total_ranks if total_ranks else 0.0,
        load_imbalance=(
            cells_per_rank_max / cells_per_rank_min
            if cells_per_rank_min
            else 0.0
        ),
        spacing_min=min((block.spacing_min for block in blocks), default=0.0),
        spacing_max=max((block.spacing_max for block in blocks), default=0.0),
        max_adjacent_ratio=max((block.max_adjacent_ratio for block in blocks), default=1.0),
        worst_cell_aspect=max((block.worst_cell_aspect for block in blocks), default=1.0),
        max_interface_ratio=max((interface.ratio for interface in interfaces), default=1.0),
        coordinate_bytes=sum(block.coordinate_bytes for block in blocks),
    )
    return GridQualityReport(summary=summary, blocks=blocks, interfaces=interfaces)


def _block_quality(block: Block) -> BlockQuality:
    axes = tuple(_axis_quality(block, index, axis) for index, axis in enumerate(AXIS_NAMES))
    axis_minima = [axis.spacing_min for axis in axes]
    axis_maxima = [axis.spacing_max for axis in axes]
    partition_minima = [_partition_extent(n, parts)[0] for n, parts in zip(block.ng, block.dims)]
    partition_maxima = [_partition_extent(n, parts)[1] for n, parts in zip(block.ng, block.dims)]
    return BlockQuality(
        block_id=block.id,
        name=block.name,
        cells=prod(block.ng),
        ranks=prod(block.dims),
        cells_per_rank_min=prod(partition_minima),
        cells_per_rank_max=prod(partition_maxima),
        spacing_min=min(axis_minima),
        spacing_max=max(axis_maxima),
        max_adjacent_ratio=max(axis.max_adjacent_ratio for axis in axes),
        worst_cell_aspect=max(axis_maxima) / min(axis_minima),
        cell_volume_min=prod(axis_minima),
        cell_volume_max=prod(axis_maxima),
        coordinate_bytes=8 * sum(n + 1 for n in block.ng),
        axes=axes,
    )


def _axis_quality(block: Block, axis_index: int, axis: str) -> AxisQuality:
    arrays = axis_grid_arrays(
        block.axes[axis].to_dict(),
        block.lmin[axis_index],
        block.lmax[axis_index],
        block.ng[axis_index],
    )
    widths = arrays.face_spacing[1:-1]
    adjacent = widths[1:] / widths[:-1]
    max_adjacent = (
        float(np.max(np.maximum(adjacent, 1.0 / adjacent)))
        if adjacent.size
        else 1.0
    )
    return AxisQuality(
        axis=axis,
        cells=block.ng[axis_index],
        length=axis_extent(block, axis_index),
        spacing_min=float(np.min(widths)),
        spacing_max=float(np.max(widths)),
        spacing_mean=float(np.mean(widths)),
        endpoint_ratio=float(widths[-1] / widths[0]),
        max_adjacent_ratio=max_adjacent,
    )


def _interface_quality(project: Project) -> tuple[InterfaceQuality, ...]:
    topology = build_topology(project.blocks, project.periodic_axes)
    by_id = {block.id: block for block in project.blocks}
    result = []
    for connection in topology.connections:
        a = by_id[connection.a_id]
        b = by_id[connection.b_id]
        a_spacing = _face_spacing(a, connection.a_face)
        b_spacing = _face_spacing(b, connection.b_face)
        result.append(
            InterfaceQuality(
                block_id=a.id,
                face=FACE_ORDER[connection.a_face],
                neighbor_id=b.id,
                neighbor_face=FACE_ORDER[connection.b_face],
                axis=AXIS_NAMES[connection.axis_index],
                spacing=a_spacing,
                neighbor_spacing=b_spacing,
                ratio=max(a_spacing, b_spacing) / min(a_spacing, b_spacing),
            )
        )
    return tuple(result)


def _face_spacing(block: Block, face: int) -> float:
    axis_index, side = FACE_INFO[face]
    axis = AXIS_NAMES[axis_index]
    arrays = axis_grid_arrays(
        block.axes[axis].to_dict(),
        block.lmin[axis_index],
        block.lmax[axis_index],
        block.ng[axis_index],
    )
    widths = arrays.face_spacing[1:-1]
    return float(widths[0] if side == 0 else widths[-1])


def _partition_extent(cells: int, partitions: int) -> tuple[int, int]:
    base, remainder = divmod(cells, partitions)
    return base, base + (1 if remainder else 0)
