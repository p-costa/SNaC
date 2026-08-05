"""Serializable project model for the SNaC grid generator."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from typing import Any

from .grid import AXIS_NAMES, fit_min_max_cell_count, fit_monotone_spacing
from .migrations import PROJECT_SCHEMA_VERSION, migrate_project_dict
from .numeric import geometry_tolerance


@dataclass
class GradingSegment:
    """OpenFOAM-style multi-grading segment weights."""

    length: float = 1.0
    cells: float = 1.0
    ratio: float = 1.0
    continuous: bool = False

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GradingSegment":
        return cls(
            length=float(data.get("length", 1.0)),
            cells=float(data.get("cells", 1.0)),
            ratio=float(data.get("ratio", 1.0)),
            continuous=bool(data.get("continuous", False)),
        )

    def to_dict(self) -> dict[str, float | bool]:
        return {
            "length": self.length,
            "cells": self.cells,
            "ratio": self.ratio,
            "continuous": self.continuous,
        }


@dataclass
class AxisSpec:
    """Grid controls for one local block axis."""

    kind: str = "snac"
    gt: int = 0
    gr: float = 0.0
    ratio: float = 1.0
    cell_ratio: float = 1.0
    width_start: float = 0.01
    width_end: float = 0.02
    controls: list[str] = field(default_factory=lambda: ["n", "ratio"])
    profile: str = "geometric"
    min: float = 0.01
    max: float = 0.02
    side: str = "end"
    segments: list[GradingSegment] = field(default_factory=lambda: [GradingSegment()])
    faces: list[float] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "AxisSpec":
        if not data:
            return cls()
        kind = str(data.get("kind", "snac"))
        side = str(data.get("side", "end"))
        ratio = float(data.get("ratio", 1.0))
        width_start = float(data.get("width_start", data.get("widthStart", data.get("min", 0.01))))
        width_end = float(data.get("width_end", data.get("widthEnd", data.get("max", 0.02))))
        controls_data = data.get("controls")

        if kind == "simple_ratio" and controls_data is None:
            if side == "start":
                ratio = 1.0 / ratio
            controls_data = ["n", "ratio"]
            side = "end"
        elif kind == "max_min" and side in {"end", "start"}:
            minimum = float(data.get("min", 0.01))
            maximum = float(data.get("max", 0.02))
            width_start, width_end = (minimum, maximum) if side == "end" else (maximum, minimum)
            ratio = width_end / width_start
            controls_data = ["width_start", "width_end"]
            kind = "simple_ratio"
            side = "end"

        return cls(
            kind=kind,
            gt=int(data.get("gt", 0)),
            gr=float(data.get("gr", 0.0)),
            ratio=ratio,
            cell_ratio=float(data.get("cell_ratio", data.get("cellRatio", 1.0))),
            width_start=width_start,
            width_end=width_end,
            controls=[str(value) for value in (controls_data or ["n", "ratio"])],
            profile=str(data.get("profile", "geometric")),
            min=float(data.get("min", 0.01)),
            max=float(data.get("max", 0.02)),
            side=side,
            segments=[
                GradingSegment.from_dict(item)
                for item in data.get("segments", [{"length": 1, "cells": 1, "ratio": 1}])
            ],
            faces=[float(value) for value in data.get("faces", [])],
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "gt": self.gt,
            "gr": self.gr,
            "ratio": self.ratio,
            "cell_ratio": self.cell_ratio,
            "width_start": self.width_start,
            "width_end": self.width_end,
            "controls": self.controls,
            "profile": self.profile,
            "min": self.min,
            "max": self.max,
            "side": self.side,
            "segments": [segment.to_dict() for segment in self.segments],
            "faces": self.faces,
        }


@dataclass
class DecompositionSpec:
    """Project-level request for an automatic MPI decomposition."""

    target_ranks: int = 0
    mode: str = "auto"
    axes: list[bool] = field(default_factory=lambda: [True, True, True])
    min_local_cells: int = 4
    max_local_aspect: float = 0.0

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "DecompositionSpec":
        data = data or {}
        return cls(
            target_ranks=max(0, int(data.get("targetRanks", data.get("target_ranks", 0)))),
            mode=str(data.get("mode", "auto")).lower(),
            axes=_bool_list(data.get("axes"), 3, True),
            min_local_cells=int(data.get("minLocalCells", data.get("min_local_cells", 4))),
            max_local_aspect=float(data.get("maxLocalAspect", data.get("max_local_aspect", 0.0))),
        )

    def validate(self) -> None:
        if self.mode not in {"auto", "1d", "2d", "3d"}:
            raise ValueError("decomposition mode must be auto, 1d, 2d, or 3d")
        self.axes.extend([True] * (3 - len(self.axes)))
        self.axes = self.axes[:3]
        if self.target_ranks and not any(self.axes):
            raise ValueError("automatic decomposition requires at least one partition axis")
        if self.min_local_cells < 1:
            raise ValueError("minimum local cells must be at least 1")
        if self.max_local_aspect and self.max_local_aspect < 1.0:
            raise ValueError("maximum local aspect must be zero or at least 1")

    def to_dict(self) -> dict[str, Any]:
        return {
            "targetRanks": self.target_ranks,
            "mode": self.mode,
            "axes": self.axes,
            "minLocalCells": self.min_local_cells,
            "maxLocalAspect": self.max_local_aspect,
        }


@dataclass
class Block:
    """A rectilinear SNaC block."""

    id: int
    name: str = ""
    dims: list[int] = field(default_factory=lambda: [1, 1, 1])
    ng: list[int] = field(default_factory=lambda: [32, 32, 2])
    lmin: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    lmax: list[float] = field(default_factory=lambda: [1.0, 1.0, 0.05])
    axes: dict[str, AxisSpec] = field(default_factory=lambda: {axis: AxisSpec() for axis in AXIS_NAMES})
    axis_locks: list[bool] = field(default_factory=lambda: [False, False, False])
    cbcvel: list[list[str]] = field(default_factory=lambda: [["D", "D", "D", "D", "D", "D"] for _ in range(3)])
    cbcpre: list[str] = field(default_factory=lambda: ["N", "N", "N", "N", "N", "N"])
    bcvel: list[list[float]] = field(default_factory=lambda: [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in range(3)])
    bcpre: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    cbcscal: list[list[str]] = field(default_factory=list)
    bcscal: list[list[float]] = field(default_factory=list)
    inflow: list[int] = field(default_factory=lambda: [0, 0, 0, 0, 0, 0])
    inivel: str = "zer"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Block":
        axes_data = data.get("axes", {})
        scalar_rows = max(len(data.get("cbcscal") or []), len(data.get("bcscal") or []))
        block = cls(
            id=int(data.get("id", 1)),
            name=str(data.get("name", "")),
            dims=[int(v) for v in data.get("dims", [1, 1, 1])],
            ng=[int(v) for v in data.get("ng", [32, 32, 2])],
            lmin=[float(v) for v in data.get("lmin", [0.0, 0.0, 0.0])],
            lmax=[float(v) for v in data.get("lmax", [1.0, 1.0, 0.05])],
            axes={axis: AxisSpec.from_dict(axes_data.get(axis)) for axis in AXIS_NAMES},
            axis_locks=_bool_list(data.get("axisLocks", data.get("axis_locks")), 3, False),
            cbcvel=_string_rows(data.get("cbcvel"), 3, ["D"] * 6),
            cbcpre=_string_list(data.get("cbcpre"), ["N"] * 6),
            bcvel=_float_rows(data.get("bcvel"), 3, [0.0] * 6),
            bcpre=_float_list(data.get("bcpre"), [0.0] * 6),
            cbcscal=_string_rows(data.get("cbcscal"), scalar_rows, ["N"] * 6),
            bcscal=_float_rows(data.get("bcscal"), scalar_rows, [0.0] * 6),
            inflow=_int_list(data.get("inflow"), [0] * 6),
            inivel=str(data.get("inivel", "zer")),
        )
        block.validate()
        return block

    def validate(self) -> None:
        if self.id < 1:
            raise ValueError("block id must be positive")
        if len(self.dims) != 3 or any(value < 1 for value in self.dims):
            raise ValueError(f"block {self.id}: dims must contain three positive integers")
        if len(self.ng) != 3 or any(value < 1 for value in self.ng):
            raise ValueError(f"block {self.id}: ng must contain three positive integers")
        if len(self.lmin) != 3 or len(self.lmax) != 3:
            raise ValueError(f"block {self.id}: lmin/lmax must contain three coordinates")
        if len(self.axis_locks) != 3:
            raise ValueError(f"block {self.id}: axis locks must contain three flags")
        if not all(isfinite(value) for value in [*self.lmin, *self.lmax]):
            raise ValueError(f"block {self.id}: lmin/lmax coordinates must be finite")
        if any(high <= low for low, high in zip(self.lmin, self.lmax)):
            raise ValueError(f"block {self.id}: every lmax coordinate must exceed lmin")
        if not all(isfinite(value) for row in [*self.bcvel, *self.bcscal] for value in row):
            raise ValueError(f"block {self.id}: velocity/scalar boundary values must be finite")
        if not all(isfinite(value) for value in self.bcpre):
            raise ValueError(f"block {self.id}: pressure boundary values must be finite")
        for index, axis in enumerate(AXIS_NAMES):
            self.axes.setdefault(axis, AxisSpec())
            length = self.lmax[index] - self.lmin[index]
            if self.axes[axis].kind == "explicit":
                faces = self.axes[axis].faces
                if len(faces) < 2 or not all(isfinite(value) for value in faces):
                    raise ValueError(f"block {self.id}: explicit {axis} grid requires finite face coordinates")
                if any(high <= low for low, high in zip(faces, faces[1:])):
                    raise ValueError(f"block {self.id}: explicit {axis} grid faces must increase")
                tolerance = geometry_tolerance(
                    faces[0], faces[-1], self.lmin[index], self.lmax[index], scale=length
                )
                if (
                    abs(faces[0] - self.lmin[index]) > tolerance
                    or abs(faces[-1] - self.lmax[index]) > tolerance
                ):
                    raise ValueError(f"block {self.id}: explicit {axis} grid must span the block extent")
                self.ng[index] = len(faces) - 1
            elif self.axes[axis].kind == "max_min":
                self.ng[index] = fit_min_max_cell_count(self.axes[axis].to_dict(), length).n
            elif self.axes[axis].kind == "simple_ratio":
                self.ng[index] = fit_monotone_spacing(self.axes[axis].to_dict(), length, self.ng[index]).n

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "dims": self.dims,
            "ng": self.ng,
            "lmin": self.lmin,
            "lmax": self.lmax,
            "axes": {axis: self.axes[axis].to_dict() for axis in AXIS_NAMES},
            "axisLocks": self.axis_locks,
            "cbcvel": self.cbcvel,
            "cbcpre": self.cbcpre,
            "bcvel": self.bcvel,
            "bcpre": self.bcpre,
            "cbcscal": self.cbcscal,
            "bcscal": self.bcscal,
            "inflow": self.inflow,
            "inivel": self.inivel,
        }


@dataclass
class Project:
    """A complete multi-block grid project."""

    name: str = "snac-grid"
    nscal: int = 0
    blocks: list[Block] = field(default_factory=lambda: [Block(id=1)])
    periodic_axes: list[bool] = field(default_factory=lambda: [False, False, False])
    decomposition: DecompositionSpec = field(default_factory=DecompositionSpec)
    infer_connectivity: bool = True
    write_external_grid: bool = True
    external_grid_source: str = "grid"
    external_grid_precision: str = "double"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Project":
        data = migrate_project_dict(data)
        blocks = [Block.from_dict(item) for item in data.get("blocks", [Block(id=1).to_dict()])]
        inferred_nscal = max((len(block.cbcscal) for block in blocks), default=0)
        external_grid_source = str(
            data.get("externalGridSource", data.get("external_grid_source", "grid"))
        )
        if external_grid_source == "data":
            external_grid_source = "grid"
        project = cls(
            name=str(data.get("name", "snac-grid")),
            nscal=max(0, int(data.get("nscal", data.get("nScal", inferred_nscal)))),
            blocks=blocks,
            periodic_axes=_bool_list(data.get("periodicAxes", data.get("periodic_axes")), 3, False),
            decomposition=DecompositionSpec.from_dict(data.get("decomposition")),
            infer_connectivity=bool(data.get("inferConnectivity", data.get("infer_connectivity", True))),
            write_external_grid=bool(data.get("writeExternalGrid", data.get("write_external_grid", True))),
            external_grid_source=external_grid_source,
            external_grid_precision=str(
                data.get("externalGridPrecision", data.get("external_grid_precision", "double"))
            ).lower(),
        )
        project.validate()
        return project

    def validate(self) -> None:
        self.nscal = max(0, int(self.nscal))
        self.periodic_axes.extend([False] * (3 - len(self.periodic_axes)))
        self.periodic_axes = self.periodic_axes[:3]
        self.decomposition.validate()
        if self.external_grid_source not in {"grid", "both"}:
            raise ValueError("external grid source must be 'grid' or 'both'")
        if self.external_grid_precision not in {"single", "double"}:
            raise ValueError("external grid precision must be 'single' or 'double'")
        seen: set[int] = set()
        for block in self.blocks:
            block.cbcscal = _resize_string_rows(block.cbcscal, self.nscal, ["N"] * 6)
            block.bcscal = _resize_float_rows(block.bcscal, self.nscal, [0.0] * 6)
            block.validate()
            if block.id in seen:
                raise ValueError(f"duplicate block id {block.id}")
            seen.add(block.id)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schemaVersion": PROJECT_SCHEMA_VERSION,
            "name": self.name,
            "nscal": self.nscal,
            "periodicAxes": self.periodic_axes,
            "decomposition": self.decomposition.to_dict(),
            "inferConnectivity": self.infer_connectivity,
            "writeExternalGrid": self.write_external_grid,
            "externalGridSource": self.external_grid_source,
            "externalGridPrecision": self.external_grid_precision,
            "blocks": [block.to_dict() for block in self.blocks],
        }


def _string_rows(data: Any, rows: int, default: list[str]) -> list[list[str]]:
    if not data:
        return [list(default) for _ in range(rows)]
    result = []
    for row in data[:rows]:
        row_values = [str(value)[:1] for value in row[:6]]
        row_values.extend(default[len(row_values) :])
        result.append(row_values)
    while len(result) < rows:
        result.append(list(default))
    return result


def _float_rows(data: Any, rows: int, default: list[float]) -> list[list[float]]:
    if not data:
        return [list(default) for _ in range(rows)]
    result = []
    for row in data[:rows]:
        row_values = [float(value) for value in row[:6]]
        row_values.extend(default[len(row_values) :])
        result.append(row_values)
    while len(result) < rows:
        result.append(list(default))
    return result


def _string_list(data: Any, default: list[str]) -> list[str]:
    values = [str(value)[:1] for value in (data or [])[: len(default)]]
    values.extend(default[len(values) :])
    return values


def _float_list(data: Any, default: list[float]) -> list[float]:
    values = [float(value) for value in (data or [])[: len(default)]]
    values.extend(default[len(values) :])
    return values


def _int_list(data: Any, default: list[int]) -> list[int]:
    values = [int(value) for value in (data or [])[: len(default)]]
    values.extend(default[len(values) :])
    return values


def _resize_string_rows(data: list[list[str]], rows: int, default: list[str]) -> list[list[str]]:
    return _string_rows(data, rows, default)


def _resize_float_rows(data: list[list[float]], rows: int, default: list[float]) -> list[list[float]]:
    return _float_rows(data, rows, default)


def _bool_list(data: Any, length: int, default: bool) -> list[bool]:
    if not data:
        return [default] * length
    result = [bool(value) for value in data[:length]]
    result.extend([default] * (length - len(result)))
    return result
