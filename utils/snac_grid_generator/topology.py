"""Format-neutral topology helpers for rectilinear block projects."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, Sequence

from .numeric import coordinates_close, geometry_tolerance

AXIS_NAMES = ("x", "y", "z")
FACE_ORDER = ("x-", "x+", "y-", "y+", "z-", "z+")
FACE_INDEX = {
    (0, 0): 0,
    (0, 1): 1,
    (1, 0): 2,
    (1, 1): 3,
    (2, 0): 4,
    (2, 1): 5,
}
FACE_INFO = {
    0: (0, 0),
    1: (0, 1),
    2: (1, 0),
    3: (1, 1),
    4: (2, 0),
    5: (2, 1),
}


class RectilinearBlock(Protocol):
    """Geometry required by the format-neutral topology engine."""

    id: int
    lmin: Sequence[float]
    lmax: Sequence[float]


class StructuredBlock(RectilinearBlock, Protocol):
    """Rectilinear block with cell counts along each axis."""

    ng: Sequence[int]


@dataclass(frozen=True)
class FaceConnection:
    """Full-face connection between two blocks."""

    a_id: int
    a_face: int
    b_id: int
    b_face: int
    axis_index: int


@dataclass
class Topology:
    """Connections and geometry diagnostics for a block collection."""

    connections: list[FaceConnection] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def build_topology(
    blocks: Sequence[RectilinearBlock],
    periodic_axes: Sequence[bool] | None = None,
) -> Topology:
    """Find full-face and periodic connections for rectilinear blocks."""

    topology = Topology()
    face_owner: dict[tuple[int, int], FaceConnection] = {}
    periodic_axes = periodic_axes or [False, False, False]

    for left_index, right_index in _candidate_block_pairs(blocks):
        left = blocks[left_index]
        right = blocks[right_index]
        if volume_overlap(left, right):
            topology.errors.append(f"blocks {left.id} and {right.id} overlap")
            continue
        for axis_index in range(3):
            scale = max(axis_extent(left, axis_index), axis_extent(right, axis_index))
            if coordinates_close(left.lmax[axis_index], right.lmin[axis_index], scale=scale):
                _add_face_contact(topology, face_owner, left, 1, right, 0, axis_index)
            if coordinates_close(right.lmax[axis_index], left.lmin[axis_index], scale=scale):
                _add_face_contact(topology, face_owner, right, 1, left, 0, axis_index)

    for axis_index, is_periodic in enumerate(periodic_axes[:3]):
        if is_periodic:
            _add_periodic_contacts(topology, face_owner, blocks, axis_index)

    return topology


def infer_global_index_extents(
    blocks: Sequence[StructuredBlock],
    connections: Sequence[FaceConnection] | None = None,
) -> dict[int, tuple[list[int], list[int]]]:
    """Infer positive one-based index extents from full-face connections."""

    if not blocks:
        return {}
    by_id = {block.id: block for block in blocks}
    if len(by_id) != len(blocks):
        raise ValueError("global indices require unique block IDs")

    if connections is None:
        topology = build_topology(blocks)
        if topology.errors:
            raise ValueError("; ".join(topology.errors))
        connections = topology.connections

    adjacency: dict[int, list[tuple[int, list[int]]]] = {
        block.id: [] for block in blocks
    }
    for connection in connections:
        a = by_id.get(connection.a_id)
        b = by_id.get(connection.b_id)
        if a is None or b is None:
            raise ValueError("global-index connection references a missing block")
        axis = connection.axis_index
        a_axis, a_side = FACE_INFO[connection.a_face]
        b_axis, b_side = FACE_INFO[connection.b_face]
        if a_axis != axis or b_axis != axis:
            raise ValueError(
                f"blocks {a.id} and {b.id} have inconsistent connection axes"
            )
        delta = [0, 0, 0]
        if (a_side, b_side) == (1, 0):
            delta[axis] = int(a.ng[axis])
        elif (a_side, b_side) == (0, 1):
            delta[axis] = -int(b.ng[axis])
        else:
            raise ValueError(
                f"blocks {a.id} and {b.id} have incompatible connected faces"
            )
        adjacency[a.id].append((b.id, delta))
        adjacency[b.id].append((a.id, [-value for value in delta]))

    offsets: dict[int, list[int]] = {}
    seed = min(by_id)
    offsets[seed] = [0, 0, 0]
    queue = [seed]
    while queue:
        block_id = queue.pop(0)
        for neighbor_id, delta in adjacency[block_id]:
            candidate = [
                offsets[block_id][axis] + delta[axis]
                for axis in range(3)
            ]
            previous = offsets.get(neighbor_id)
            if previous is not None:
                if previous != candidate:
                    raise ValueError(
                        f"inconsistent global-index loop at block {neighbor_id}: "
                        f"{previous} conflicts with {candidate}"
                    )
                continue
            offsets[neighbor_id] = candidate
            queue.append(neighbor_id)

    missing = sorted(set(by_id) - set(offsets))
    if missing:
        raise ValueError(
            f"cannot infer global indices for disconnected blocks {missing}"
        )

    minima = [min(offset[axis] for offset in offsets.values()) for axis in range(3)]
    extents: dict[int, tuple[list[int], list[int]]] = {}
    for block_id, offset in offsets.items():
        lo = [offset[axis] - minima[axis] + 1 for axis in range(3)]
        hi = [lo[axis] + int(by_id[block_id].ng[axis]) - 1 for axis in range(3)]
        extents[block_id] = (lo, hi)

    _check_index_extent_overlaps(extents)
    return extents


def same_cross_section(
    a: RectilinearBlock,
    b: RectilinearBlock,
    axis_index: int,
) -> bool:
    """Return whether two blocks have equal transverse extents."""

    for index in range(3):
        if index == axis_index:
            continue
        scale = max(axis_extent(a, index), axis_extent(b, index))
        if not (
            coordinates_close(a.lmin[index], b.lmin[index], scale=scale)
            and coordinates_close(a.lmax[index], b.lmax[index], scale=scale)
        ):
            return False
    return True


def cross_section_overlap(
    a: RectilinearBlock,
    b: RectilinearBlock,
    axis_index: int,
) -> bool:
    """Return whether two block cross-sections overlap with positive area."""

    return all(
        interval_overlap(a.lmin[index], a.lmax[index], b.lmin[index], b.lmax[index])
        for index in range(3)
        if index != axis_index
    )


def volume_overlap(a: RectilinearBlock, b: RectilinearBlock) -> bool:
    """Return whether two blocks overlap with positive volume."""

    return all(
        interval_overlap(a.lmin[index], a.lmax[index], b.lmin[index], b.lmax[index])
        for index in range(3)
    )


def interval_overlap(a0: float, a1: float, b0: float, b1: float) -> bool:
    """Return whether two intervals overlap beyond geometric tolerance."""

    overlap = min(a1, b1) - max(a0, b0)
    scale = max(a1 - a0, b1 - b0)
    return overlap > geometry_tolerance(a0, a1, b0, b1, scale=scale)


def axis_extent(block: RectilinearBlock, axis_index: int) -> float:
    """Return one positive block extent."""

    return float(block.lmax[axis_index] - block.lmin[axis_index])


def _candidate_block_pairs(
    blocks: Sequence[RectilinearBlock],
) -> list[tuple[int, int]]:
    """Return pairs whose closed bounding boxes overlap within tolerance."""

    if len(blocks) < 2:
        return []
    sweep_tolerance = max(
        geometry_tolerance(
            block.lmin[0],
            block.lmax[0],
            scale=axis_extent(block, 0),
        )
        for block in blocks
    )
    order = sorted(
        range(len(blocks)),
        key=lambda index: (blocks[index].lmin[0], index),
    )
    active: list[int] = []
    candidates: list[tuple[int, int]] = []
    for current_index in order:
        current = blocks[current_index]
        active = [
            index
            for index in active
            if blocks[index].lmax[0]
            >= current.lmin[0] - sweep_tolerance
        ]
        for other_index in active:
            other = blocks[other_index]
            if all(
                _closed_intervals_overlap(other, current, axis_index)
                for axis_index in range(3)
            ):
                candidates.append(
                    (min(other_index, current_index), max(other_index, current_index))
                )
        active.append(current_index)
    return sorted(candidates)


def _closed_intervals_overlap(
    a: RectilinearBlock,
    b: RectilinearBlock,
    axis_index: int,
) -> bool:
    scale = max(axis_extent(a, axis_index), axis_extent(b, axis_index))
    tolerance = geometry_tolerance(
        a.lmin[axis_index],
        a.lmax[axis_index],
        b.lmin[axis_index],
        b.lmax[axis_index],
        scale=scale,
    )
    return min(a.lmax[axis_index], b.lmax[axis_index]) >= (
        max(a.lmin[axis_index], b.lmin[axis_index]) - tolerance
    )


def _check_index_extent_overlaps(
    extents: dict[int, tuple[list[int], list[int]]],
) -> None:
    block_ids = sorted(extents)
    for index, left_id in enumerate(block_ids):
        left_lo, left_hi = extents[left_id]
        for right_id in block_ids[index + 1 :]:
            right_lo, right_hi = extents[right_id]
            if all(
                max(left_lo[axis], right_lo[axis])
                <= min(left_hi[axis], right_hi[axis])
                for axis in range(3)
            ):
                raise ValueError(
                    f"blocks {left_id} and {right_id} have overlapping global indices"
                )


def _add_periodic_contacts(
    topology: Topology,
    face_owner: dict[tuple[int, int], FaceConnection],
    blocks: Sequence[RectilinearBlock],
    axis_index: int,
) -> None:
    if not blocks:
        return

    axis = AXIS_NAMES[axis_index]
    global_min = min(block.lmin[axis_index] for block in blocks)
    global_max = max(block.lmax[axis_index] for block in blocks)
    low_blocks = [
        block
        for block in blocks
        if coordinates_close(block.lmin[axis_index], global_min, scale=axis_extent(block, axis_index))
    ]
    high_blocks = [
        block
        for block in blocks
        if coordinates_close(block.lmax[axis_index], global_max, scale=axis_extent(block, axis_index))
    ]
    matched_high: set[int] = set()

    for low_block in low_blocks:
        matches = [block for block in high_blocks if same_cross_section(low_block, block, axis_index)]
        if len(matches) != 1:
            topology.errors.append(
                f"periodic {axis}- face of block {low_block.id} has {len(matches)} matching "
                "full faces on the opposite boundary"
            )
            continue
        high_block = matches[0]
        matched_high.add(high_block.id)
        _add_connection(
            topology,
            face_owner,
            low_block,
            FACE_INDEX[(axis_index, 0)],
            high_block,
            FACE_INDEX[(axis_index, 1)],
            axis_index,
        )

    for high_block in high_blocks:
        if high_block.id not in matched_high:
            topology.errors.append(
                f"periodic {axis}+ face of block {high_block.id} has no matching full face on the opposite boundary"
            )


def _add_face_contact(
    topology: Topology,
    face_owner: dict[tuple[int, int], FaceConnection],
    lower_block: RectilinearBlock,
    lower_side: int,
    upper_block: RectilinearBlock,
    upper_side: int,
    axis_index: int,
) -> None:
    if same_cross_section(lower_block, upper_block, axis_index):
        _add_connection(
            topology,
            face_owner,
            lower_block,
            FACE_INDEX[(axis_index, lower_side)],
            upper_block,
            FACE_INDEX[(axis_index, upper_side)],
            axis_index,
        )
    elif cross_section_overlap(lower_block, upper_block, axis_index):
        topology.errors.append(
            f"blocks {lower_block.id} and {upper_block.id} touch on a partial {AXIS_NAMES[axis_index]} face"
        )


def _add_connection(
    topology: Topology,
    face_owner: dict[tuple[int, int], FaceConnection],
    a: RectilinearBlock,
    a_face: int,
    b: RectilinearBlock,
    b_face: int,
    axis_index: int,
) -> None:
    connection = FaceConnection(a.id, a_face, b.id, b_face, axis_index)
    for block_id, face in ((connection.a_id, connection.a_face), (connection.b_id, connection.b_face)):
        old = face_owner.get((block_id, face))
        if old is not None:
            topology.errors.append(f"block {block_id} face {FACE_ORDER[face]} has multiple full-face neighbors")
        face_owner[(block_id, face)] = connection
    topology.connections.append(connection)
