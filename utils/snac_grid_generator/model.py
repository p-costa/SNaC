"""Serializable project model for the SNaC grid generator."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .grid import AXIS_NAMES


@dataclass
class GradingSegment:
    """OpenFOAM-style multi-grading segment weights."""

    length: float = 1.0
    cells: float = 1.0
    ratio: float = 1.0

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GradingSegment":
        return cls(
            length=float(data.get("length", 1.0)),
            cells=float(data.get("cells", 1.0)),
            ratio=float(data.get("ratio", 1.0)),
        )

    def to_dict(self) -> dict[str, float]:
        return {"length": self.length, "cells": self.cells, "ratio": self.ratio}


@dataclass
class AxisSpec:
    """Grid controls for one local block axis."""

    kind: str = "snac"
    gt: int = 0
    gr: float = 0.0
    ratio: float = 1.0
    profile: str = "geometric"
    min: float = 0.01
    max: float = 0.02
    side: str = "end"
    segments: list[GradingSegment] = field(default_factory=lambda: [GradingSegment()])

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "AxisSpec":
        if not data:
            return cls()
        return cls(
            kind=str(data.get("kind", "snac")),
            gt=int(data.get("gt", 0)),
            gr=float(data.get("gr", 0.0)),
            ratio=float(data.get("ratio", 1.0)),
            profile=str(data.get("profile", "geometric")),
            min=float(data.get("min", 0.01)),
            max=float(data.get("max", 0.02)),
            side=str(data.get("side", "end")),
            segments=[GradingSegment.from_dict(item) for item in data.get("segments", [{"length": 1, "cells": 1, "ratio": 1}])],
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "gt": self.gt,
            "gr": self.gr,
            "ratio": self.ratio,
            "profile": self.profile,
            "min": self.min,
            "max": self.max,
            "side": self.side,
            "segments": [segment.to_dict() for segment in self.segments],
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
            cbcvel=_string_rows(data.get("cbcvel"), 3, ["D"] * 6),
            cbcpre=[str(v)[:1] for v in data.get("cbcpre", ["N"] * 6)],
            bcvel=_float_rows(data.get("bcvel"), 3, [0.0] * 6),
            bcpre=[float(v) for v in data.get("bcpre", [0.0] * 6)],
            cbcscal=_string_rows(data.get("cbcscal"), scalar_rows, ["N"] * 6),
            bcscal=_float_rows(data.get("bcscal"), scalar_rows, [0.0] * 6),
            inflow=[int(v) for v in data.get("inflow", [0] * 6)],
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
        if any(high <= low for low, high in zip(self.lmin, self.lmax)):
            raise ValueError(f"block {self.id}: every lmax coordinate must exceed lmin")
        for axis in AXIS_NAMES:
            self.axes.setdefault(axis, AxisSpec())

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "dims": self.dims,
            "ng": self.ng,
            "lmin": self.lmin,
            "lmax": self.lmax,
            "axes": {axis: self.axes[axis].to_dict() for axis in AXIS_NAMES},
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
    infer_connectivity: bool = True
    write_external_grid: bool = True
    external_grid_source: str = "grid"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Project":
        blocks = [Block.from_dict(item) for item in data.get("blocks", [Block(id=1).to_dict()])]
        inferred_nscal = max((len(block.cbcscal) for block in blocks), default=0)
        project = cls(
            name=str(data.get("name", "snac-grid")),
            nscal=max(0, int(data.get("nscal", data.get("nScal", inferred_nscal)))),
            blocks=blocks,
            periodic_axes=_bool_list(data.get("periodicAxes", data.get("periodic_axes")), 3, False),
            infer_connectivity=bool(data.get("inferConnectivity", data.get("infer_connectivity", True))),
            write_external_grid=bool(data.get("writeExternalGrid", data.get("write_external_grid", True))),
            external_grid_source=str(data.get("externalGridSource", data.get("external_grid_source", "grid"))),
        )
        project.validate()
        return project

    def validate(self) -> None:
        self.nscal = max(0, int(self.nscal))
        self.periodic_axes.extend([False] * (3 - len(self.periodic_axes)))
        self.periodic_axes = self.periodic_axes[:3]
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
            "name": self.name,
            "nscal": self.nscal,
            "periodicAxes": self.periodic_axes,
            "inferConnectivity": self.infer_connectivity,
            "writeExternalGrid": self.write_external_grid,
            "externalGridSource": self.external_grid_source,
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
