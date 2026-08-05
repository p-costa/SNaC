"""SNaC boundary-condition presets for grid projects."""

from __future__ import annotations

from typing import Any, Iterable, Sequence

from .model import Block, Project
from .topology import FACE_ORDER


def apply_bc_preset(
    project: Project | dict[str, Any],
    block_ids: Iterable[int],
    preset: str,
    *,
    face: str | None = None,
    velocity: Sequence[float] = (0.0, 0.0, 0.0),
    pressure: float = 0.0,
) -> tuple[Project, list[int]]:
    """Apply a common SNaC velocity/pressure BC pattern to selected blocks."""

    project = Project.from_dict(project.to_dict() if isinstance(project, Project) else project)
    selected_ids = {int(block_id) for block_id in block_ids}
    blocks = [block for block in project.blocks if block.id in selected_ids]
    if not blocks or len(blocks) != len(selected_ids):
        raise ValueError("select existing blocks before applying a boundary preset")
    if preset not in {"free_slip", "inlet", "moving_wall", "no_slip", "outlet"}:
        raise ValueError(f"unknown boundary preset {preset!r}")
    face_index = _face_index(face) if preset in {"inlet", "moving_wall", "outlet"} else None
    values = [float(value) for value in velocity]
    if len(values) != 3:
        raise ValueError("velocity must have three components")

    for block in blocks:
        if preset in {"no_slip", "moving_wall"}:
            _set_all_external(block, "D", (0.0, 0.0, 0.0), "N", 0.0)
        elif preset == "free_slip":
            _set_free_slip(block)
        if preset in {"inlet", "moving_wall"}:
            _set_face(block, face_index, "D", values, "N", 0.0)
        elif preset == "outlet":
            _set_face(block, face_index, "N", (0.0, 0.0, 0.0), "D", float(pressure))
    project.validate()
    return project, sorted(selected_ids)


def clear_friend_boundaries(blocks: Iterable[Block]) -> None:
    """Replace inferred SNaC friend boundaries with neutral physical BCs."""

    for block in blocks:
        for face in range(6):
            if block.cbcpre[face] == "F":
                block.cbcpre[face] = "N"
                block.bcpre[face] = 0.0
            for component in range(3):
                if block.cbcvel[component][face] == "F":
                    block.cbcvel[component][face] = "D"
                    block.bcvel[component][face] = 0.0
            for scalar in range(len(block.cbcscal)):
                if block.cbcscal[scalar][face] == "F":
                    block.cbcscal[scalar][face] = "N"
                    block.bcscal[scalar][face] = 0.0


def connect_friend_boundaries(a: Block, a_face: int, b: Block, b_face: int) -> None:
    """Encode one inferred block connection using SNaC friend boundaries."""

    for component in range(3):
        a.cbcvel[component][a_face] = "F"
        b.cbcvel[component][b_face] = "F"
        a.bcvel[component][a_face] = float(b.id)
        b.bcvel[component][b_face] = float(a.id)
    for scalar in range(len(a.cbcscal)):
        a.cbcscal[scalar][a_face] = "F"
        b.cbcscal[scalar][b_face] = "F"
        a.bcscal[scalar][a_face] = float(b.id)
        b.bcscal[scalar][b_face] = float(a.id)
    a.cbcpre[a_face] = "F"
    b.cbcpre[b_face] = "F"
    a.bcpre[a_face] = float(b.id)
    b.bcpre[b_face] = float(a.id)


def _face_index(face: str | None) -> int:
    if face not in FACE_ORDER:
        raise ValueError("choose a valid boundary face")
    return FACE_ORDER.index(face)


def _set_all_external(
    block: Block,
    velocity_code: str,
    velocity: Sequence[float],
    pressure_code: str,
    pressure: float,
) -> None:
    for face in range(6):
        if not _is_friend_face(block, face):
            _set_face(block, face, velocity_code, velocity, pressure_code, pressure)


def _set_free_slip(block: Block) -> None:
    for face in range(6):
        if _is_friend_face(block, face):
            continue
        normal_component = face // 2
        for component in range(3):
            block.cbcvel[component][face] = "D" if component == normal_component else "N"
            block.bcvel[component][face] = 0.0
        block.cbcpre[face] = "N"
        block.bcpre[face] = 0.0


def _set_face(
    block: Block,
    face: int | None,
    velocity_code: str,
    velocity: Sequence[float],
    pressure_code: str,
    pressure: float,
) -> None:
    if face is None or _is_friend_face(block, face):
        return
    for component in range(3):
        block.cbcvel[component][face] = velocity_code
        block.bcvel[component][face] = float(velocity[component])
    block.cbcpre[face] = pressure_code
    block.bcpre[face] = float(pressure)


def _is_friend_face(block: Block, face: int) -> bool:
    codes = [block.cbcpre[face], *[row[face] for row in block.cbcvel], *[row[face] for row in block.cbcscal]]
    return "F" in codes
