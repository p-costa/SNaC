"""Cell-count alternatives that balance an exact MPI decomposition better."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import itertools
import json
from math import ceil, floor, prod
from pathlib import Path
import random
import tempfile
from time import perf_counter
import unittest

from snac_grid_generator import (
    Project,
    apply_balance_alternative,
    check_project,
    describe_project_decomposition,
    optimize_project_decomposition,
    update_project_structure,
)
from snac_grid_generator.cli import main as grid_cli_main
from snac_grid_generator.topology import FACE_INFO, axis_components, build_topology


def step_project(ranks: int, tolerance: float = 0.1) -> Project:
    """A backward-facing step: inlet 32x24x32, upper and lower 96x24x32, FFT along z."""

    project = Project.from_dict(
        {
            "name": "step",
            "fftAxis": "z",
            "periodicAxes": [False, False, True],
            "decomposition": {"targetRanks": ranks, "ngTolerance": tolerance},
            "blocks": [
                {"id": 1, "name": "inlet", "ng": [32, 24, 32], "lmin": [0, 1, 0], "lmax": [2, 2, 2]},
                {"id": 2, "name": "upper", "ng": [96, 24, 32], "lmin": [2, 1, 0], "lmax": [10, 2, 2]},
                {"id": 3, "name": "lower", "ng": [96, 24, 32], "lmin": [2, 0, 0], "lmax": [10, 1, 2]},
            ],
        }
    )
    project, result = update_project_structure(project)
    # Only the requested rank count is unmet until Balance runs.
    assert all(issue.code == "decomposition" for issue in result.issues if issue.severity == "error"), result.errors
    return project


def balance(project: Project):
    return optimize_project_decomposition(project, alternatives=True)


def imbalance_of(alternative: dict, ranks: int) -> float:
    blocks = alternative["blocks"]
    largest = max(prod(ceil(n / d) for n, d in zip(block["ng"], block["dims"])) for block in blocks)
    return ranks * largest / sum(prod(block["ng"]) for block in blocks)


class BalanceAlternativeTests(unittest.TestCase):
    def test_two_unequal_blocks_meet_in_the_middle(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 2},
                "blocks": [
                    {"id": 1, "ng": [30, 8, 4], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [34, 8, 4], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        _, result = balance(project)
        self.assertAlmostEqual(result.imbalance, 68 / 64)
        alternatives = result.alternatives
        self.assertEqual(
            [(round(item["maxChange"], 4), round(item["imbalance"], 4)) for item in alternatives],
            [(round(1 / 30, 4), round(66 / 64, 4)), (round(2 / 34, 4), round(64 / 63, 4)), (round(2 / 30, 4), 1.0)],
        )
        perfect = alternatives[-1]
        self.assertTrue(perfect["exact"])
        self.assertEqual(
            sorted((change["from"], change["to"]) for change in perfect["changes"]),
            [(30, 32), (34, 32)],
        )
        self.assertEqual([block["dims"] for block in perfect["blocks"]], [[1, 1, 1], [1, 1, 1]])
        self.assertTrue(result.alternatives_optimal)

    def test_step_case_trade_offs(self) -> None:
        # Confirmed by independent brute-force searches over every partition
        # and every cell count within 10%.
        expected = {
            6: (1.2857, [(0.0938, 1.2488)]),
            7: (1.0, []),
            9: (1.2857, [(0.0417, 1.1640), (0.0625, 1.1117), (0.0833, 1.0707), (0.0938, 1.0557)]),
            12: (1.0714, [(0.0104, 1.0270), (0.0625, 1.0093), (0.0938, 1.0047)]),
        }
        for ranks, (baseline, points) in expected.items():
            with self.subTest(ranks=ranks):
                _, result = balance(step_project(ranks))
                self.assertAlmostEqual(result.imbalance, baseline, places=4)
                found = [(round(item["maxChange"], 4), round(item["imbalance"], 4)) for item in result.alternatives]
                self.assertEqual(found, points)
                for alternative in result.alternatives:
                    self.assertAlmostEqual(alternative["imbalance"], imbalance_of(alternative, ranks))
                    self.assertTrue(all(block["ng"][2] == 32 for block in alternative["blocks"]))
                    self.assertTrue(all(block["dims"][2] == 1 for block in alternative["blocks"]))
                self.assertEqual(
                    result.fixed_components,
                    [] if baseline == 1.0 else [{"axis": "z", "blockIds": [1, 2, 3], "ng": 32, "reason": "FFT axis"}],
                )

    def test_counts_set_by_gradings_locks_and_fft_stay_fixed(self) -> None:
        project = Project.from_dict(
            {
                "fftAxis": "z",
                "decomposition": {"targetRanks": 4},
                "blocks": [
                    {
                        "id": 1,
                        "ng": [16, 16, 8],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "explicit", "faces": [index / 16 for index in range(17)]}},
                    },
                    {
                        "id": 2,
                        "ng": [20, 16, 8],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "axisLocks": [True, False, False],
                    },
                    {
                        "id": 3,
                        "ng": [20, 16, 8],
                        "lmin": [2, 0, 0],
                        "lmax": [3, 1, 1],
                        "axes": {
                            "x": {
                                "kind": "monotone",
                                "controls": ["widthStart", "widthEnd"],
                                "widthStart": 0.05,
                                "widthEnd": 0.05,
                            }
                        },
                    },
                    {"id": 4, "ng": [22, 16, 8], "lmin": [3, 0, 0], "lmax": [4, 1, 1]},
                ],
            }
        )
        project, structured = update_project_structure(project)
        self.assertTrue(structured.ok, structured.errors)
        _, result = balance(project)
        reasons = {(item["axis"], tuple(item["blockIds"])): item["reason"] for item in result.fixed_components}
        self.assertEqual(
            reasons,
            {
                ("x", (1,)): "the grading sets the cell count",
                ("x", (2,)): "locked",
                ("x", (3,)): "the grading sets the cell count",
                ("z", (1, 2, 3, 4)): "FFT axis",
            },
        )
        self.assertTrue(result.alternatives)
        for alternative in result.alternatives:
            self.assertEqual([change["blockIds"] for change in alternative["changes"]], [[4]])
            ng = {block["blockId"]: block["ng"] for block in alternative["blocks"]}
            self.assertEqual([ng[1][0], ng[2][0], ng[3][0]], [16, 20, 20])
        self.assertEqual(result.alternatives[-1]["changes"][0]["to"], 20)
        self.assertAlmostEqual(result.alternatives[-1]["imbalance"], 4 * 20 / 76)

    def test_tolerance_setting_migrates_validates_and_disables(self) -> None:
        self.assertEqual(Project.from_dict({"schemaVersion": 4, "blocks": []}).decomposition.ng_tolerance, 0.1)
        self.assertEqual(Project.from_dict({"blocks": []}).to_dict()["decomposition"]["ngTolerance"], 0.1)
        with self.assertRaisesRegex(ValueError, "cell-count tolerance"):
            Project.from_dict({"decomposition": {"ngTolerance": 0.6}, "blocks": []})
        _, result = balance(step_project(9, tolerance=0.0))
        self.assertEqual(result.alternatives, [])
        self.assertEqual(result.fixed_components, [])
        self.assertEqual(result.to_dict()["ngTolerance"], 0.0)
        # A smaller tolerance keeps every change within it.
        _, result = balance(step_project(9, tolerance=0.05))
        self.assertEqual([round(item["imbalance"], 4) for item in result.alternatives], [1.1640])

    def test_applied_alternative_is_a_valid_balanced_project(self) -> None:
        decomposed, result = balance(step_project(9))
        alternative = result.alternatives[2]
        applied = apply_balance_alternative(decomposed, alternative)
        self.assertEqual(sum(prod(block.dims) for block in applied.blocks), 9)
        for before, after in zip(decomposed.blocks, applied.blocks):
            for old, new in zip(before.ng, after.ng):
                self.assertLessEqual(abs(new - old), 0.1 * old)
        self.assertTrue(check_project(applied).ok)
        described = describe_project_decomposition(applied)
        self.assertTrue(described.ok)
        self.assertAlmostEqual(described.imbalance, alternative["imbalance"])
        self.assertEqual(described.total_ranks, 9)

        with self.assertRaisesRegex(ValueError, "missing block 7"):
            apply_balance_alternative(decomposed, {"blocks": [{"blockId": 7, "ng": [8, 8, 8], "dims": [1, 1, 1]}]})
        with self.assertRaisesRegex(ValueError, "must list the blocks"):
            apply_balance_alternative(decomposed, {})

    def test_alternatives_match_a_brute_force_search(self) -> None:
        rng = random.Random(20261024)
        checked = 0
        for case in range(60):
            kind = rng.choice(["line", "periodic", "corner", "step", "square"])
            project = _random_layout(rng, kind, ranks=rng.randint(2, 9))
            decomposed, result = balance(project)
            if not result.ok:
                continue
            checked += 1
            with self.subTest(case=case, kind=kind, ranks=project.decomposition.target_ranks):
                reference = min(result.imbalance, _brute_force(project, 0.0))
                best = _brute_force(project, project.decomposition.ng_tolerance)
                found = min((item["imbalance"] for item in result.alternatives), default=None)
                if best <= reference - 0.01:
                    self.assertIsNotNone(found)
                    self.assertAlmostEqual(found, best, places=12)
                else:
                    self.assertEqual(result.alternatives, [])
                _assert_valid_alternatives(self, decomposed, result)
        self.assertGreater(checked, 40)

    def test_large_lattice_stays_within_its_time_budget(self) -> None:
        rng = random.Random(3)
        counts = [[rng.randint(20, 60) for _ in range(4)] for _ in range(3)]
        blocks = [
            {
                "id": 1 + i * 16 + j * 4 + k,
                "ng": [counts[0][i], counts[1][j], counts[2][k]],
                "lmin": [i, j, k],
                "lmax": [i + 1, j + 1, k + 1],
            }
            for i in range(4)
            for j in range(4)
            for k in range(4)
        ]
        project = Project.from_dict({"blocks": blocks, "decomposition": {"targetRanks": 64}})
        start = perf_counter()
        decomposed, result = optimize_project_decomposition(
            project, time_limit=1.0, alternatives=True, alternatives_time_limit=1.0
        )
        self.assertLess(perf_counter() - start, 6.0)
        self.assertTrue(result.alternatives)
        self.assertLess(result.alternatives[-1]["imbalance"], result.imbalance - 0.1)
        _assert_valid_alternatives(self, decomposed, result)

    def test_cli_lists_and_applies_alternatives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "step.json"
            output = Path(tmp) / "balanced.json"
            source.write_text(json.dumps(step_project(9).to_dict()), encoding="utf-8")
            with redirect_stdout(StringIO()) as printed:
                self.assertEqual(grid_cli_main(["decompose", str(source)]), 0)
            text = printed.getvalue()
            self.assertIn("load imbalance 1.2857", text)
            self.assertIn("1. up to 4.17%: imbalance 1.1640, exact division", text)
            self.assertIn("x of upper, lower: 96 -> 100", text)

            with redirect_stdout(StringIO()) as printed:
                code = grid_cli_main(
                    ["decompose", str(source), "--ng-tolerance", "0.05", "--apply-alternative", "1", "-o", str(output)]
                )
            self.assertEqual(code, 0)
            self.assertIn("Applied alternative 1.", printed.getvalue())
            saved = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(saved["decomposition"]["ngTolerance"], 0.05)
            self.assertEqual([block["ng"][0] for block in saved["blocks"]], [31, 100, 100])
            self.assertEqual(sum(prod(block["dims"]) for block in saved["blocks"]), 9)

            with redirect_stdout(StringIO()):
                self.assertEqual(grid_cli_main(["decompose", str(source), "--apply-alternative", "9"]), 2)


def _random_layout(rng: random.Random, kind: str, ranks: int) -> Project:
    periodic = [False, False, False]
    if kind in {"line", "periodic"}:
        edges = [0.0]
        for _ in range(rng.choice([2, 3])):
            edges.append(edges[-1] + rng.choice([0.5, 1.0, 1.5, 2.0]))
        boxes = [([edges[i], 0, 0], [edges[i + 1], 1, 1]) for i in range(len(edges) - 1)]
        periodic[0] = kind == "periodic"
    elif kind == "corner":
        boxes = [([0, 0, 0], [1, 1, 1]), ([1, 0, 0], [2.5, 1, 1]), ([0, 1, 0], [1, 1.5, 1])]
    elif kind == "step":
        a, b = rng.choice([1.0, 2.0]), rng.choice([4.0, 8.0])
        boxes = [([0, 1, 0], [a, 2, 1]), ([a, 1, 0], [a + b, 2, 1]), ([a, 0, 0], [a + b, 1, 1])]
        periodic[2] = True
    else:
        boxes = [([i, j, 0], [i + 1, j + 1, 1]) for i in range(2) for j in range(2)]
    blocks = [
        {"id": index + 1, "ng": [8, 8, 4], "lmin": list(map(float, low)), "lmax": list(map(float, high))}
        for index, (low, high) in enumerate(boxes)
    ]
    shape = Project.from_dict({"blocks": blocks, "periodicAxes": periodic})
    connections = build_topology(shape.blocks, shape.periodic_axes).connections
    for axis in range(3):
        for group in axis_components(shape.blocks, connections, axis):
            # z counts are too small for a 10% change, so they stay fixed.
            count = rng.randint(8, 30) if axis < 2 else rng.choice([4, 6, 8])
            for block_id in group:
                blocks[block_id - 1]["ng"][axis] = count
    return Project.from_dict(
        {
            "blocks": blocks,
            "periodicAxes": periodic,
            "fftAxis": "z" if kind == "step" else None,
            "decomposition": {
                "targetRanks": ranks,
                "minLocalCells": rng.choice([2, 3, 4]),
                "maxLocalAspect": rng.choice([0.0, 0.0, 3.0, 5.0]),
            },
        }
    )


def _brute_force(project: Project, tolerance: float) -> float:
    """The lowest imbalance over every partition and every cell count within *tolerance*.

    Only for uniform grids, whose interface spacing is length / count.
    """

    connections = build_topology(project.blocks, project.periodic_axes).connections
    ranks = project.decomposition.target_ranks
    minimum = project.decomposition.min_local_cells
    aspect = project.decomposition.max_local_aspect
    groups = []
    group_of: dict[tuple[int, int], int] = {}
    for axis in range(3):
        for group in axis_components(project.blocks, connections, axis):
            n0 = next(block.ng[axis] for block in project.blocks if block.id in group)
            low = max(1, ceil(n0 * (1 - tolerance) - 1e-9))
            high = floor(n0 * (1 + tolerance) + 1e-9)
            counts = [n0] if axis == 2 else list(range(low, high + 1))
            for block_id in group:
                group_of[(block_id, axis)] = len(groups)
            groups.append((axis, counts))
    splittable = [
        project.decomposition.axes[axis] and axis != project.fft_axis_index for axis in range(3)
    ]
    lengths = {block.id: [block.lmax[a] - block.lmin[a] for a in range(3)] for block in project.blocks}
    pairs = []
    for connection in connections:
        axis = FACE_INFO[connection.a_face][0]
        pairs.append(
            (
                group_of[(connection.a_id, axis)],
                lengths[connection.a_id][axis],
                group_of[(connection.b_id, axis)],
                lengths[connection.b_id][axis],
            )
        )
    base = [
        next(block.ng[axis] for block in project.blocks if group_of[(block.id, axis)] == index)
        for index, (axis, _) in enumerate(groups)
    ]
    limits = [
        max(3.0, max(la / base[a], lb / base[b]) / min(la / base[a], lb / base[b])) * (1 + 1e-9)
        for a, la, b, lb in pairs
    ]
    ids = [block.id for block in project.blocks]
    best = float("inf")
    split_choices = [
        list(range(1, max(counts) // minimum + 1)) or [1] if splittable[axis] else [1] for axis, counts in groups
    ]
    for splits in itertools.product(*split_choices):
        if sum(prod(splits[group_of[(i, a)]] for a in range(3)) for i in ids) != ranks:
            continue
        options = [
            [n for n in counts if splits[g] == 1 or n // splits[g] >= minimum] for g, (_, counts) in enumerate(groups)
        ]
        for counts in itertools.product(*options):
            if any(
                max(la / counts[a], lb / counts[b]) / min(la / counts[a], lb / counts[b]) > limit
                for (a, la, b, lb), limit in zip(pairs, limits)
            ):
                continue
            if aspect and any(
                max(local) / min(local) > aspect * (1 + 1e-12)
                for local in (
                    [counts[group_of[(i, a)]] / splits[group_of[(i, a)]] for a in range(3)] for i in ids
                )
            ):
                continue
            largest = max(
                prod(-(-counts[group_of[(i, a)]] // splits[group_of[(i, a)]]) for a in range(3)) for i in ids
            )
            best = min(best, ranks * largest / sum(prod(counts[group_of[(i, a)]] for a in range(3)) for i in ids))
    return best


def _assert_valid_alternatives(test: unittest.TestCase, project: Project, result) -> None:
    ranks = project.decomposition.target_ranks
    tolerance = project.decomposition.ng_tolerance
    baseline_errors = set(check_project(project).errors)
    minimum = project.decomposition.min_local_cells
    for alternative in result.alternatives:
        blocks = {block["blockId"]: block for block in alternative["blocks"]}
        test.assertEqual(sum(prod(block["dims"]) for block in blocks.values()), ranks)
        test.assertLess(alternative["imbalance"], result.imbalance)
        test.assertAlmostEqual(alternative["imbalance"], imbalance_of(alternative, ranks))
        for block in project.blocks:
            entry = blocks[block.id]
            for axis in range(3):
                test.assertLessEqual(abs(entry["ng"][axis] - block.ng[axis]), tolerance * block.ng[axis] + 1e-9)
                if entry["dims"][axis] > 1:
                    test.assertGreaterEqual(entry["ng"][axis] // entry["dims"][axis], minimum)
        applied = apply_balance_alternative(project, alternative)
        test.assertLessEqual(set(check_project(applied).errors), baseline_errors)


if __name__ == "__main__":
    unittest.main()
