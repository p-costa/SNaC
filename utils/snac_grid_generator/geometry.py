"""Format-neutral geometry operations for rectilinear block projects."""

from __future__ import annotations

from typing import Any, Iterable

from .grid import AXIS_NAMES
from .model import Block, Project


def duplicate_blocks(
    project: Project | dict[str, Any],
    block_ids: Iterable[int],
    axis: str,
    count: int = 1,
    gap: float = 0.0,
) -> tuple[Project, list[int]]:
    """Duplicate a selected block group as a regularly spaced array."""

    project = _copy_project(project)
    selected = _selected_blocks(project, block_ids)
    axis_index = _axis_index(axis)
    count = int(count)
    if count < 1:
        raise ValueError("array count must be at least 1")
    group_min = min(block.lmin[axis_index] for block in selected)
    group_max = max(block.lmax[axis_index] for block in selected)
    pitch = group_max - group_min + float(gap)
    if pitch <= 0.0:
        raise ValueError("array spacing must keep copies ordered")

    next_id = max((block.id for block in project.blocks), default=0) + 1
    new_ids: list[int] = []
    for copy_index in range(1, count + 1):
        copies: list[Block] = []
        id_map: dict[int, int] = {}
        for source in selected:
            block = Block.from_dict(source.to_dict())
            block.id = next_id
            block.name = f"{source.name or f'block-{source.id}'}-{copy_index + 1}"
            id_map[source.id] = block.id
            next_id += 1
            _translate_block(block, axis_index, copy_index * pitch)
            copies.append(block)
            new_ids.append(block.id)
        for block in copies:
            _remap_friend_boundaries(block, id_map)
        project.blocks.extend(copies)
    project.validate()
    return project, new_ids


def mirror_blocks(
    project: Project | dict[str, Any],
    block_ids: Iterable[int],
    axis: str,
    plane: float,
) -> tuple[Project, list[int]]:
    """Create mirrored copies of selected blocks about a coordinate plane."""

    project = _copy_project(project)
    selected = _selected_blocks(project, block_ids)
    axis_index = _axis_index(axis)
    next_id = max((block.id for block in project.blocks), default=0) + 1
    copies: list[Block] = []
    id_map: dict[int, int] = {}
    for source in selected:
        block = Block.from_dict(source.to_dict())
        block.id = next_id
        block.name = f"{source.name or f'block-{source.id}'}-mirror"
        id_map[source.id] = block.id
        next_id += 1
        _mirror_block(block, axis_index, float(plane))
        copies.append(block)
    for block in copies:
        _remap_friend_boundaries(block, id_map)
    project.blocks.extend(copies)
    project.validate()
    return project, [block.id for block in copies]


def align_blocks(
    project: Project | dict[str, Any],
    block_ids: Iterable[int],
    source_block_id: int,
    axis: str,
    mode: str,
) -> tuple[Project, list[int]]:
    """Align lower, center, or upper coordinates to a source block."""

    project = _copy_project(project)
    selected = _selected_blocks(project, block_ids)
    source = _block_by_id(project, source_block_id)
    axis_index = _axis_index(axis)
    if mode not in {"lower", "center", "upper"}:
        raise ValueError("alignment mode must be lower, center, or upper")
    source_value = _alignment_value(source, axis_index, mode)
    changed: list[int] = []
    for block in selected:
        if block.id == source.id:
            continue
        shift = source_value - _alignment_value(block, axis_index, mode)
        _translate_block(block, axis_index, shift)
        changed.append(block.id)
    project.validate()
    return project, changed


def snap_block_to_face(
    project: Project | dict[str, Any],
    block_id: int,
    source_block_id: int,
    face: str,
) -> tuple[Project, list[int]]:
    """Snap one block to a source face and match its transverse extent."""

    project = _copy_project(project)
    block = _block_by_id(project, block_id)
    source = _block_by_id(project, source_block_id)
    if block.id == source.id:
        raise ValueError("choose a different block to snap")
    if face not in {f"{axis}{side}" for axis in AXIS_NAMES for side in ("-", "+")}:
        raise ValueError(f"unknown face {face!r}")
    axis_index = _axis_index(face[0])
    thickness = block.lmax[axis_index] - block.lmin[axis_index]
    if face[1] == "+":
        new_min = source.lmax[axis_index]
        new_max = new_min + thickness
    else:
        new_max = source.lmin[axis_index]
        new_min = new_max - thickness
    _rescale_block_axis(block, axis_index, new_min, new_max)
    for index in range(3):
        if index != axis_index:
            _rescale_block_axis(block, index, source.lmin[index], source.lmax[index])
    project.validate()
    return project, [block.id]


def _copy_project(project: Project | dict[str, Any]) -> Project:
    return Project.from_dict(project.to_dict() if isinstance(project, Project) else project)


def _selected_blocks(project: Project, block_ids: Iterable[int]) -> list[Block]:
    requested = {int(block_id) for block_id in block_ids}
    selected = [block for block in project.blocks if block.id in requested]
    if not selected:
        raise ValueError("select at least one block")
    if len(selected) != len(requested):
        raise ValueError("one or more selected blocks do not exist")
    return sorted(selected, key=lambda block: block.id)


def _block_by_id(project: Project, block_id: int) -> Block:
    block = next((item for item in project.blocks if item.id == int(block_id)), None)
    if block is None:
        raise ValueError(f"block {block_id} does not exist")
    return block


def _axis_index(axis: str) -> int:
    if axis not in AXIS_NAMES:
        raise ValueError(f"unknown axis {axis!r}")
    return AXIS_NAMES.index(axis)


def _alignment_value(block: Block, axis_index: int, mode: str) -> float:
    if mode == "lower":
        return block.lmin[axis_index]
    if mode == "upper":
        return block.lmax[axis_index]
    return 0.5 * (block.lmin[axis_index] + block.lmax[axis_index])


def _translate_block(block: Block, axis_index: int, shift: float) -> None:
    _rescale_block_axis(
        block,
        axis_index,
        block.lmin[axis_index] + shift,
        block.lmax[axis_index] + shift,
    )


def _rescale_block_axis(block: Block, axis_index: int, new_min: float, new_max: float) -> None:
    axis = block.axes[AXIS_NAMES[axis_index]]
    old_min = block.lmin[axis_index]
    old_max = block.lmax[axis_index]
    if axis.kind == "explicit":
        scale = (new_max - new_min) / (old_max - old_min)
        axis.faces = [new_min + (value - old_min) * scale for value in axis.faces]
    block.lmin[axis_index] = new_min
    block.lmax[axis_index] = new_max


def _mirror_block(block: Block, axis_index: int, plane: float) -> None:
    axis = block.axes[AXIS_NAMES[axis_index]]
    old_min = block.lmin[axis_index]
    old_max = block.lmax[axis_index]
    block.lmin[axis_index] = 2.0 * plane - old_max
    block.lmax[axis_index] = 2.0 * plane - old_min
    if axis.kind == "explicit":
        axis.faces = [2.0 * plane - value for value in reversed(axis.faces)]
    elif axis.kind == "snac":
        axis.gt = {1: 3, 3: 1, 4: 5, 5: 4}.get(axis.gt, axis.gt)
    elif axis.kind == "simple_ratio":
        if axis.ratio <= 0.0 or axis.cell_ratio <= 0.0:
            raise ValueError("cannot mirror a non-positive grading ratio")
        axis.ratio = 1.0 / axis.ratio
        axis.cell_ratio = 1.0 / axis.cell_ratio
        axis.width_start, axis.width_end = axis.width_end, axis.width_start
    elif axis.kind == "multi":
        old_segments = list(axis.segments)
        axis.segments = list(reversed(axis.segments))
        for index, segment in enumerate(axis.segments):
            if segment.ratio <= 0.0:
                raise ValueError("cannot mirror a non-positive segment ratio")
            segment.ratio = 1.0 / segment.ratio
            segment.continuous = index > 0 and old_segments[len(old_segments) - index].continuous
    lo_face = 2 * axis_index
    hi_face = lo_face + 1
    for values in [*block.cbcvel, block.cbcpre, *block.bcvel, block.bcpre, *block.cbcscal, *block.bcscal, block.inflow]:
        values[lo_face], values[hi_face] = values[hi_face], values[lo_face]


def _remap_friend_boundaries(copy: Block, id_map: dict[int, int]) -> None:
    rows = [
        *((codes, values, "D") for codes, values in zip(copy.cbcvel, copy.bcvel)),
        (copy.cbcpre, copy.bcpre, "N"),
        *((codes, values, "N") for codes, values in zip(copy.cbcscal, copy.bcscal)),
    ]
    for codes, values, external_code in rows:
        for face, code in enumerate(codes):
            if code != "F":
                continue
            friend = int(round(values[face]))
            if friend in id_map:
                values[face] = float(id_map[friend])
            else:
                codes[face] = external_code
                values[face] = 0.0
