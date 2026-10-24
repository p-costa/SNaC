"""Structured locations for validation messages and GUI navigation."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Iterable


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


def diagnostics_for_messages(
    errors: Iterable[str],
    warnings: Iterable[str],
    interfaces: Iterable[dict[str, Any]] = (),
) -> list[Diagnostic]:
    """Attach stable categories and geometric locations to check messages."""

    interface_list = list(interfaces)
    return [
        *(_diagnostic(message, "error", interface_list) for message in errors),
        *(_diagnostic(message, "warning", interface_list) for message in warnings),
    ]


def _diagnostic(
    message: str,
    severity: str,
    interfaces: list[dict[str, Any]],
) -> Diagnostic:
    locations = _message_locations(message)
    block_ids = {location.block_id for location in locations}
    for interface in interfaces:
        interface_ids = {int(interface["blockId"]), int(interface["neighborId"])}
        if len(block_ids) >= 2 and interface_ids.issubset(block_ids):
            locations = (
                DiagnosticLocation(int(interface["blockId"]), str(interface["face"])),
                DiagnosticLocation(
                    int(interface["neighborId"]),
                    str(interface["neighborFace"]),
                ),
            )
            break
    return Diagnostic(
        severity=severity,
        code=_message_code(message),
        message=message,
        locations=locations,
        axis=_message_axis(message),
    )


def _message_locations(message: str) -> tuple[DiagnosticLocation, ...]:
    exact = re.findall(r"\bblock\s+(\d+)\s+face\s+([xyz][+-])", message, re.IGNORECASE)
    locations = [DiagnosticLocation(int(block_id), face.lower()) for block_id, face in exact]

    block_ids: list[int] = []
    for match in re.finditer(r"\bblocks?\s+\[([^]]+)\]", message, re.IGNORECASE):
        block_ids.extend(int(value) for value in re.findall(r"\d+", match.group(1)))
    for match in re.finditer(r"\bblocks?\s+(\d+)(?:\s+and\s+(\d+))?", message, re.IGNORECASE):
        block_ids.append(int(match.group(1)))
        if match.group(2):
            block_ids.append(int(match.group(2)))
    block_ids.extend(
        int(value)
        for value in re.findall(
            r"\b(?:touches|references|expected|from|to match)\s+block\s+(\d+)",
            message,
            re.IGNORECASE,
        )
    )

    located = {location.block_id for location in locations}
    for block_id in block_ids:
        if block_id not in located:
            locations.append(DiagnosticLocation(block_id))
            located.add(block_id)
    return tuple(locations)


def _message_axis(message: str) -> str | None:
    patterns = (
        r"\balong\s+([xyz])\b",
        r"\b(?:mpi|periodic|SNaC)\s+([xyz])\b",
        r"\binvalid\s+([xyz])\s+grid\b",
        r"\b([xyz])\s+interface\b",
    )
    for pattern in patterns:
        match = re.search(pattern, message, re.IGNORECASE)
        if match:
            return match.group(1).lower()
    face = re.search(r"\bface\s+([xyz])[+-]", message, re.IGNORECASE)
    return face.group(1).lower() if face else None


def _message_code(message: str) -> str:
    lowered = message.lower()
    categories = (
        ("project has no blocks", "no_blocks"),
        ("disconnected", "disconnected_blocks"),
        ("partial", "partial_face_contact"),
        ("overlap", "block_overlap"),
        ("different ng", "interface_cell_count"),
        ("touching grid lines", "interface_grid_mismatch"),
        ("touching mpi partitions", "interface_partition_mismatch"),
        ("interface spacing jumps", "interface_spacing_jump"),
        ("marked f", "friend_boundary"),
        ("references block", "friend_boundary"),
        ("not fully marked f", "friend_boundary"),
        ("bc pair", "boundary_condition"),
        ("inflow type", "boundary_condition"),
        ("mapping function", "grid_mapping"),
        ("grid growth parameter", "grid_mapping"),
        ("invalid", "invalid_grid"),
        ("decomposition", "decomposition"),
        ("mpi partitions", "decomposition"),
    )
    return next((code for phrase, code in categories if phrase in lowered), "validation")
