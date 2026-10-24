"""Structured validation messages for the GUI, the CLI, and reports.

Every check creates its diagnostics here, with the affected blocks, faces, and
axis attached at the source. Nothing is recovered later by parsing message text,
so a message can be reworded without breaking GUI navigation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Union


@dataclass(frozen=True)
class DiagnosticLocation:
    """One block, and optionally one face, affected by a diagnostic."""

    block_id: int
    face: str | None = None

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"blockId": self.block_id}
        if self.face is not None:
            result["face"] = self.face
        return result


@dataclass(frozen=True)
class Diagnostic:
    """A machine-navigable validation message."""

    severity: str
    code: str
    message: str
    locations: tuple[DiagnosticLocation, ...] = ()
    axis: str | None = None

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
            "locations": [location.to_dict() for location in self.locations],
        }
        if self.axis is not None:
            result["axis"] = self.axis
        return result


# A place is a block id, a (block id, face label) pair, or a location.
Place = Union[int, "tuple[int, str]", DiagnosticLocation]


def error(code: str, message: str, *places: Place, axis: str | None = None) -> Diagnostic:
    """Create an error diagnostic located at *places*."""

    return Diagnostic("error", code, message, _locations(places), axis)


def warning(code: str, message: str, *places: Place, axis: str | None = None) -> Diagnostic:
    """Create a warning diagnostic located at *places*."""

    return Diagnostic("warning", code, message, _locations(places), axis)


def _locations(places: tuple[Place, ...]) -> tuple[DiagnosticLocation, ...]:
    """Normalize places, dropping a bare block that a face location already names."""

    locations: list[DiagnosticLocation] = []
    for place in places:
        if isinstance(place, DiagnosticLocation):
            location = place
        elif isinstance(place, tuple):
            location = DiagnosticLocation(int(place[0]), str(place[1]))
        else:
            location = DiagnosticLocation(int(place))
        if location not in locations:
            locations.append(location)
    with_faces = {location.block_id for location in locations if location.face is not None}
    return tuple(
        location
        for location in locations
        if location.face is not None or location.block_id not in with_faces
    )
