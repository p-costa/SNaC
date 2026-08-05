from __future__ import annotations

import json
from contextlib import redirect_stdout
from io import StringIO
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

from utils.snac_grid_generator import (
    MIN_LOCAL_CELLS,
    Project,
    align_blocks,
    apply_bc_preset,
    apply_axis_to_aligned_blocks,
    axis_grid_arrays,
    axis_grid_diagnostics,
    check_project,
    duplicate_blocks,
    export_project,
    fit_monotone_spacing,
    import_case,
    optimize_project_decomposition,
    mirror_blocks,
    repair_project_grids,
    solve_spacing,
    snap_block_to_face,
    update_project_structure,
)
from utils.snac_grid_generator.snac_grid import binary_payload
from utils.snac_grid_generator.cli import main as grid_cli_main
from utils.snac_grid_generator.topology import FaceConnection, infer_global_index_extents


class GridGeneratorTests(unittest.TestCase):
    def test_blockmesh_spacing_relations(self) -> None:
        solution = solve_spacing(length=7.0, n=3, ratio=4.0)
        self.assertAlmostEqual(solution.cell_ratio, 2.0)
        self.assertAlmostEqual(solution.width_start, 1.0)
        self.assertAlmostEqual(solution.width_end, 4.0)

    def test_spacing_solver_derives_cell_count_without_overflow(self) -> None:
        solution = solve_spacing(length=1.0, width_start=0.1, ratio=3.0)
        self.assertEqual(solution.n, 5)
        self.assertAlmostEqual(solution.width_end / solution.width_start, 3.0)
        geometric_sum = solution.width_start * sum(solution.cell_ratio**index for index in range(solution.n))
        self.assertAlmostEqual(geometric_sum, solution.length)

        single = solve_spacing(length=2.0, n=1)
        self.assertEqual(single.width_start, 2.0)
        self.assertEqual(single.width_end, 2.0)
        with self.assertRaisesRegex(ValueError, "n=1"):
            solve_spacing(length=2.0, n=1, ratio=2.0)

    def test_geometric_monotone_accepts_every_control_pair(self) -> None:
        for n, ratio in ((2, 0.1), (3, 0.7), (15, 1.3), (32, 10.0)):
            reference = solve_spacing(length=1.0, n=n, ratio=ratio)
            values = {
                "n": reference.n,
                "ratio": reference.ratio,
                "cell_ratio": reference.cell_ratio,
                "width_start": reference.width_start,
                "width_end": reference.width_end,
            }
            controls = tuple(values)
            for left_index, left in enumerate(controls):
                for right in controls[left_index + 1 :]:
                    axis = {
                        "profile": "geometric",
                        "controls": [left, right],
                        left: values[left],
                        right: values[right],
                    }
                    with self.subTest(n=n, ratio=ratio, controls=(left, right)):
                        fit = fit_monotone_spacing(axis, 1.0, reference.n)
                        self.assertEqual(fit.n, reference.n)
                        self.assertAlmostEqual(fit.ratio / reference.ratio, 1.0, places=8)
                        self.assertAlmostEqual(fit.cell_ratio / reference.cell_ratio, 1.0, places=8)
                        self.assertAlmostEqual(fit.width_start / reference.width_start, 1.0, places=8)
                        self.assertAlmostEqual(fit.width_end / reference.width_end, 1.0, places=8)

    def test_mapped_monotone_inverts_endpoint_controls(self) -> None:
        for profile in ("tanh", "erf"):
            for n, ratio in ((4, 0.25), (8, 2.0), (32, 8.0)):
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    reference = axis_grid_diagnostics(
                        {
                            "kind": "simple_ratio",
                            "profile": profile,
                            "controls": ["n", "ratio"],
                            "ratio": ratio,
                        },
                        0.0,
                        1.0,
                        n,
                    )
                    values = {
                        "n": n,
                        "ratio": reference["expansion"],
                        "width_start": reference["lower"],
                        "width_end": reference["upper"],
                    }
                    pairs = (
                        ("n", "width_start"),
                        ("n", "width_end"),
                        ("ratio", "width_start"),
                        ("ratio", "width_end"),
                        ("width_start", "width_end"),
                    )
                    for controls in pairs:
                        axis = {
                            "profile": profile,
                            "controls": list(controls),
                            **{control: values[control] for control in controls},
                        }
                        with self.subTest(profile=profile, n=n, ratio=ratio, controls=controls):
                            fit = fit_monotone_spacing(axis, 1.0, n)
                            self.assertEqual(fit.n, n)
                            self.assertAlmostEqual(fit.ratio / ratio, 1.0, places=7)
                            self.assertAlmostEqual(fit.width_start / reference["lower"], 1.0, places=7)
                            self.assertAlmostEqual(fit.width_end / reference["upper"], 1.0, places=7)

    def test_monotone_rejects_underdetermined_or_profile_specific_controls(self) -> None:
        with self.assertRaisesRegex(ValueError, "do not determine n"):
            fit_monotone_spacing(
                {"profile": "geometric", "controls": ["ratio", "cell_ratio"], "ratio": 1.0, "cell_ratio": 1.0},
                1.0,
                32,
            )
        with self.assertRaisesRegex(ValueError, "only for geometric"):
            fit_monotone_spacing(
                {"profile": "tanh", "controls": ["n", "cell_ratio"], "cell_ratio": 1.01},
                1.0,
                32,
            )

    def test_monotone_diagnostics_report_achieved_fit(self) -> None:
        metrics = axis_grid_diagnostics(
            {
                "kind": "simple_ratio",
                "profile": "geometric",
                "controls": ["width_start", "width_end"],
                "width_start": 0.03,
                "width_end": 0.12,
            },
            0.0,
            1.0,
            32,
        )
        self.assertEqual(metrics["n"], 15)
        self.assertAlmostEqual(metrics["expansion"], 4.0)
        self.assertAlmostEqual(metrics["idealN"], 15.236778, places=5)
        self.assertAlmostEqual(metrics["residuals"]["width_start"], 0.015599, places=5)
        self.assertAlmostEqual(metrics["residuals"]["width_end"], 0.015599, places=5)
        arrays = axis_grid_arrays(
            {
                "kind": "simple_ratio",
                "profile": "geometric",
                "controls": ["width_start", "width_end"],
                "width_start": 0.03,
                "width_end": 0.12,
            },
            0.0,
            1.0,
            32,
        )
        self.assertEqual(arrays.faces[-2], 1.0)

    def test_multi_grading_normalizes_weights(self) -> None:
        arrays = axis_grid_arrays(
            {
                "kind": "multi",
                "segments": [
                    {"length": 20, "cells": 30, "ratio": 4},
                    {"length": 60, "cells": 40, "ratio": 1},
                    {"length": 20, "cells": 30, "ratio": 0.25},
                ],
            },
            0.0,
            10.0,
            30,
        )
        self.assertEqual(binary_payload(arrays).shape, (120,))
        self.assertAlmostEqual(arrays.faces[0], 0.0)
        self.assertAlmostEqual(arrays.faces[-2], 10.0)
        self.assertTrue(np.all(np.diff(arrays.faces[:-1]) > 0.0))

    def test_multi_grading_distributes_remainder_across_segments(self) -> None:
        arrays = axis_grid_arrays(
            {
                "kind": "multi",
                "segments": [
                    {"length": 1, "cells": 1, "ratio": 1},
                    {"length": 1, "cells": 1, "ratio": 1},
                    {"length": 1, "cells": 1, "ratio": 1},
                ],
            },
            0.0,
            1.0,
            5,
        )
        widths = arrays.face_spacing[1:-1]
        np.testing.assert_allclose(widths[:4], np.full(4, 1.0 / 6.0))
        self.assertAlmostEqual(widths[4], 1.0 / 3.0)

    def test_multi_grading_can_link_segment_spacing(self) -> None:
        base = {
            "kind": "multi",
            "profile": "geometric",
            "segments": [
                {"length": 1, "cells": 1, "ratio": 4},
                {"length": 1, "cells": 1, "ratio": 4},
            ],
        }
        unlinked = axis_grid_diagnostics(base, 0.0, 1.0, 32)
        self.assertGreater(unlinked["segments"][1]["jump"], 1.1)

        base["segments"][1]["continuous"] = True
        linked = axis_grid_diagnostics(base, 0.0, 1.0, 32)
        self.assertAlmostEqual(linked["segments"][1]["jump"], 1.0, places=10)
        self.assertTrue(linked["segments"][1]["continuous"])

        project = Project.from_dict(
            {"blocks": [{"id": 1, "ng": [32, 2, 2], "axes": {"x": base}}]}
        )
        saved = project.to_dict()
        self.assertTrue(saved["blocks"][0]["axes"]["x"]["segments"][1]["continuous"])
        self.assertEqual(Project.from_dict(saved).to_dict(), saved)

    def test_ratio_profiles_preserve_requested_direction(self) -> None:
        for profile in ("geometric", "tanh", "erf"):
            arrays = axis_grid_arrays({"kind": "simple_ratio", "ratio": 4.0, "profile": profile}, 0.0, 1.0, 32)
            widths = arrays.face_spacing[1:-1]
            self.assertAlmostEqual(widths[-1] / widths[0], 4.0, places=7)

    def test_single_ratio_can_be_biased_to_upper_end(self) -> None:
        for profile in ("geometric", "tanh", "erf"):
            arrays = axis_grid_arrays({"kind": "simple_ratio", "ratio": 4.0, "side": "start", "profile": profile}, 0.0, 1.0, 32)
            widths = arrays.face_spacing[1:-1]
            self.assertAlmostEqual(widths[-1] / widths[0], 0.25, places=7)

    def test_min_max_clustering_options(self) -> None:
        arrays = axis_grid_arrays({"kind": "max_min", "min": 0.01, "max": 0.04, "side": "end"}, 0.0, 1.0, 32)
        widths = arrays.face_spacing[1:-1]
        self.assertLess(widths[0], widths[-1])

        arrays = axis_grid_arrays({"kind": "max_min", "min": 0.01, "max": 0.04, "side": "both"}, 0.0, 1.0, 32)
        widths = arrays.face_spacing[1:-1]
        self.assertLess(widths[0], widths[15])
        self.assertLess(widths[-1], widths[16])

        arrays = axis_grid_arrays({"kind": "max_min", "min": 0.01, "max": 0.04, "side": "middle"}, 0.0, 1.0, 32)
        widths = arrays.face_spacing[1:-1]
        self.assertGreater(widths[0], widths[15])
        self.assertGreater(widths[-1], widths[16])

    def test_min_max_spacing_derives_cell_count_and_preserves_odd_symmetry(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "ng": [31, 2, 2],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "max_min", "min": 0.01, "max": 0.04, "side": "both"}},
                    }
                ]
            }
        )
        block = project.blocks[0]
        arrays = axis_grid_arrays(block.axes["x"].to_dict(), 0.0, 1.0, block.ng[0])
        widths = arrays.face_spacing[1:-1]
        self.assertEqual(block.ng[0], 46)
        np.testing.assert_allclose(widths, widths[::-1], rtol=1.0e-12, atol=1.0e-14)
        self.assertAlmostEqual(float(widths.min()), 0.01, delta=3.0e-5)
        self.assertAlmostEqual(float(widths.max()), 0.04, delta=1.1e-4)

        odd = axis_grid_arrays(
            {"kind": "max_min", "min": 0.01, "max": 0.04, "side": "middle"},
            0.0,
            1.0,
            31,
        ).face_spacing[1:-1]
        np.testing.assert_allclose(odd, odd[::-1], rtol=1.0e-12, atol=1.0e-14)

        with self.assertRaisesRegex(ValueError, "max width"):
            axis_grid_arrays({"kind": "max_min", "min": 0.04, "max": 0.01}, 0.0, 1.0, 32)

    def test_legacy_one_sided_min_max_migrates_to_monotone_controls(self) -> None:
        lower = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "ng": [32, 2, 2],
                        "axes": {"x": {"kind": "max_min", "min": 0.01, "max": 0.04, "side": "end"}},
                    }
                ]
            }
        ).blocks[0].axes["x"]
        self.assertEqual(lower.kind, "simple_ratio")
        self.assertEqual(lower.controls, ["width_start", "width_end"])
        self.assertEqual((lower.width_start, lower.width_end), (0.01, 0.04))

        upper = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "ng": [32, 2, 2],
                        "axes": {"x": {"kind": "max_min", "min": 0.01, "max": 0.04, "side": "start"}},
                    }
                ]
            }
        ).blocks[0].axes["x"]
        self.assertEqual((upper.width_start, upper.width_end), (0.04, 0.01))

        symmetric = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "axes": {"x": {"kind": "max_min", "min": 0.01, "max": 0.04, "side": "both"}},
                    }
                ]
            }
        ).blocks[0].axes["x"]
        self.assertEqual(symmetric.kind, "max_min")
        self.assertEqual(symmetric.side, "both")

    def test_export_writes_block_and_external_grids(self) -> None:
        project = Project.from_dict(
            {
                "name": "unit",
                "inferConnectivity": True,
                "writeExternalGrid": True,
                "externalGridSource": "grid",
                "blocks": [
                    {
                        "id": 1,
                        "ng": [4, 5, 2],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 2, 0.1],
                        "axes": {
                            "x": {"kind": "simple_ratio", "ratio": 2},
                            "y": {"kind": "snac", "gt": 0, "gr": 0},
                            "z": {"kind": "snac", "gt": 0, "gr": 0},
                        },
                    }
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            result = export_project(project, tmp)
            root = Path(tmp)
            self.assertTrue((root / "blocks.nml").exists())
            self.assertTrue((root / "grid" / "grid_x_b_001.bin").exists())
            payload = np.fromfile(root / "grid" / "grid_x_b_001.bin", dtype="<f8")
            self.assertEqual(payload.shape, (16,))
            saved = json.loads((root / "snac_grid_project.json").read_text())
            self.assertEqual(saved["name"], "unit")
            self.assertGreaterEqual(len(result.files), 8)

    def test_imports_an_existing_blocks_namelist(self) -> None:
        result = import_case(Path("examples/__MULTI_BLOCK_GEOMETRY/L_channel"))
        project = result.project
        self.assertEqual(result.source, "blocks.nml")
        self.assertEqual(len(project.blocks), 3)
        self.assertEqual(project.blocks[2].dims, [2, 1, 1])
        self.assertEqual(project.blocks[0].cbcvel[1], ["D", "D", "D", "F", "D", "D"])
        self.assertEqual(project.blocks[1].bcpre[1], 3.0)
        self.assertFalse(project.write_external_grid)
        self.assertTrue(check_project(project).ok)

        periodic = import_case(Path("tests/differentially_heated_cavity")).project
        self.assertEqual(periodic.nscal, 1)
        self.assertEqual(periodic.periodic_axes, [False, True, False])
        self.assertTrue(check_project(periodic).ok)

    def test_imports_same_line_namelist_assignments(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp)
            (case_dir / "blocks.nml").write_text(
                "&blocks\n"
                "nblocks=1, block_dims(1)=2,3,4, block_ng(1)=8,9,10, "
                "block_inivel(1)='zer'\n"
                "/\n",
                encoding="utf-8",
            )
            project = import_case(case_dir).project

        self.assertEqual(project.blocks[0].dims, [2, 3, 4])
        self.assertEqual(project.blocks[0].ng, [8, 9, 10])
        self.assertEqual(project.blocks[0].inivel, "zer")

    def test_imports_external_binary_grids_exactly(self) -> None:
        project = Project.from_dict(
            {
                "name": "external-round-trip",
                "writeExternalGrid": True,
                "blocks": [
                    {
                        "id": 1,
                        "ng": [12, 8, 4],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 7, "profile": "erf"}},
                    }
                ],
            }
        )
        expected = axis_grid_arrays(project.blocks[0].axes["x"].to_dict(), 0.0, 1.0, 12).faces[:-1]
        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            (Path(tmp) / "snac_grid_project.json").unlink()
            imported = import_case(tmp).project
        axis = imported.blocks[0].axes["x"]
        self.assertEqual(axis.kind, "explicit")
        np.testing.assert_allclose(axis.faces, expected, rtol=0.0, atol=1.0e-14)
        self.assertTrue(imported.write_external_grid)

    def test_explicit_grids_require_external_grid_writing(self) -> None:
        project = Project.from_dict(
            {
                "writeExternalGrid": False,
                "blocks": [
                    {
                        "id": 1,
                        "ng": [2, 2, 2],
                        "axes": {"x": {"kind": "explicit", "faces": [0.0, 0.25, 1.0]}},
                    }
                ],
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertIn("non-native grids require external grid writing", result.errors)

    def test_export_preserves_per_component_velocity_bcs(self) -> None:
        project = Project.from_dict(
            {
                "name": "component-bcs",
                "blocks": [
                    {
                        "id": 1,
                        "ng": [4, 4, 4],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "cbcvel": [
                            ["D", "N", "D", "D", "D", "D"],
                            ["N", "D", "D", "D", "D", "D"],
                            ["D", "D", "N", "D", "D", "D"],
                        ],
                        "bcvel": [
                            [1.25, 2.5, 0, 0, 0, 0],
                            [3.75, 4.5, 0, 0, 0, 0],
                            [0, 0, 5.25, 6.5, 0, 0],
                        ],
                    }
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            blocks = (Path(tmp) / "blocks.nml").read_text()
            self.assertIn("block_cbcvel(0:1,1:3,1,1) = 'D', 'N', 'D', 'D', 'D', 'D'", blocks)
            self.assertIn("block_cbcvel(0:1,1:3,2,1) = 'N', 'D', 'D', 'D', 'D', 'D'", blocks)
            self.assertIn("block_cbcvel(0:1,1:3,3,1) = 'D', 'D', 'N', 'D', 'D', 'D'", blocks)
            self.assertIn("block_bcvel(0:1,1:3,1,1) = 1.25, 2.5, 0., 0., 0., 0.", blocks)
            self.assertIn("block_bcvel(0:1,1:3,2,1) = 3.75, 4.5, 0., 0., 0., 0.", blocks)
            self.assertIn("block_bcvel(0:1,1:3,3,1) = 0., 0., 5.25, 6.5, 0., 0.", blocks)

    def test_export_writes_arbitrary_scalar_bcs_and_renumbers_friends(self) -> None:
        project = Project.from_dict(
            {
                "name": "scalar-bcs",
                "nscal": 2,
                "blocks": [
                    {
                        "id": 10,
                        "ng": [4, 4, 4],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "cbcscal": [
                            ["D", "F", "N", "N", "N", "N"],
                            ["N", "F", "D", "D", "N", "N"],
                        ],
                        "bcscal": [
                            [1.0, 20.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 20.0, -0.5, 0.5, 0.0, 0.0],
                        ],
                    },
                    {
                        "id": 20,
                        "ng": [4, 4, 4],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "cbcscal": [
                            ["F", "N", "N", "N", "N", "N"],
                            ["F", "N", "D", "D", "N", "N"],
                        ],
                        "bcscal": [
                            [10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [10.0, 0.0, -0.5, 0.5, 0.0, 0.0],
                        ],
                    },
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            blocks = (Path(tmp) / "blocks.nml").read_text()
            self.assertIn("block_cbcscal(0:1,1:3,1,1) = 'D', 'F', 'N', 'N', 'N', 'N'", blocks)
            self.assertIn("block_bcscal(0:1,1:3,1,1) = 1., 2., 0., 0., 0., 0.", blocks)
            self.assertIn("block_cbcscal(0:1,1:3,2,2) = 'F', 'N', 'D', 'D', 'N', 'N'", blocks)
            self.assertIn("block_bcscal(0:1,1:3,2,2) = 1., 0., -0.5, 0.5, 0., 0.", blocks)

    def test_export_geometry_uses_inferred_global_indices(self) -> None:
        project = Project.from_dict(
            {
                "name": "indices",
                "blocks": [
                    {"id": 1, "ng": [4, 5, 2], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [6, 5, 2], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            lines = (Path(tmp) / "data" / "geometry_b_002.out").read_text().splitlines()
            self.assertEqual(lines[0].split(), ["5", "1", "1"])
            self.assertEqual(lines[1].split(), ["10", "5", "2"])

    def test_export_geometry_indices_do_not_depend_on_block_ids(self) -> None:
        project = Project.from_dict(
            {
                "name": "reverse-indices",
                "blocks": [
                    {"id": 1, "ng": [6, 5, 2], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                    {"id": 2, "ng": [4, 5, 2], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                ],
            }
        )

        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            right = (Path(tmp) / "data" / "geometry_b_001.out").read_text().splitlines()
            left = (Path(tmp) / "data" / "geometry_b_002.out").read_text().splitlines()

        self.assertEqual(right[0].split(), ["5", "1", "1"])
        self.assertEqual(right[1].split(), ["10", "5", "2"])
        self.assertEqual(left[0].split(), ["1", "1", "1"])
        self.assertEqual(left[1].split(), ["4", "5", "2"])

    def test_global_index_inference_rejects_disconnected_and_conflicting_graphs(self) -> None:
        blocks = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "ng": [4, 5, 2]},
                    {"id": 2, "ng": [6, 5, 2], "lmin": [1, 0, 0], "lmax": [2, 1, 0.05]},
                ]
            }
        ).blocks
        with self.assertRaisesRegex(ValueError, "disconnected blocks"):
            infer_global_index_extents(blocks, [])

        conflicting = [
            FaceConnection(1, 1, 2, 0, 0),
            FaceConnection(1, 0, 2, 1, 0),
        ]
        with self.assertRaisesRegex(ValueError, "inconsistent global-index loop"):
            infer_global_index_extents(blocks, conflicting)

    def test_empty_project_can_exist_but_cannot_be_exported(self) -> None:
        project = Project.from_dict({"name": "empty", "blocks": []})
        check = check_project(project)
        self.assertFalse(check.ok)
        self.assertTrue(any("no blocks" in error for error in check.errors))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                export_project(project, tmp)
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_model_rejects_non_finite_extents_and_normalizes_face_arrays(self) -> None:
        with self.assertRaisesRegex(ValueError, "finite"):
            Project.from_dict({"blocks": [{"id": 1, "lmin": [float("nan"), 0, 0]}]})

        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "cbcpre": ["D"],
                        "bcpre": [2.0],
                        "inflow": [3],
                    }
                ]
            }
        )
        block = project.blocks[0]
        self.assertEqual(block.cbcpre, ["D", "N", "N", "N", "N", "N"])
        self.assertEqual(block.bcpre, [2.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.assertEqual(block.inflow, [3, 0, 0, 0, 0, 0])

    def test_project_json_has_an_explicit_schema_version(self) -> None:
        project = Project.from_dict({"blocks": []})
        self.assertEqual(project.to_dict()["schemaVersion"], 2)
        self.assertEqual(project.decomposition.min_local_cells, 4)
        with self.assertRaisesRegex(ValueError, "schema version"):
            Project.from_dict({"schemaVersion": 3, "blocks": []})

    def test_project_json_preserves_decomposition_and_axis_locks(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {
                    "targetRanks": 64,
                    "mode": "2d",
                    "axes": [True, True, False],
                    "minLocalCells": 6,
                    "maxLocalAspect": 8.0,
                },
                "blocks": [{"id": 1, "axisLocks": [False, True, False]}],
            }
        )
        saved = project.to_dict()
        self.assertEqual(saved["decomposition"]["targetRanks"], 64)
        self.assertEqual(saved["decomposition"]["mode"], "2d")
        self.assertEqual(saved["decomposition"]["axes"], [True, True, False])
        self.assertEqual(saved["decomposition"]["minLocalCells"], 6)
        self.assertEqual(saved["decomposition"]["maxLocalAspect"], 8.0)
        self.assertEqual(saved["blocks"][0]["axisLocks"], [False, True, False])
        self.assertEqual(Project.from_dict(saved).to_dict(), saved)

    def test_check_rejects_disconnected_blocks(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [2, 0, 0], "lmax": [3, 1, 1]},
                ]
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("disconnected" in error for error in result.errors))

    def test_export_removes_stale_owned_grids_and_preserves_other_files(self) -> None:
        project = Project.from_dict({"name": "managed", "blocks": [{"id": 1, "ng": [4, 4, 4]}]})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_project(project, root)
            stale_grid = root / "grid" / "grid_x_b_001.bin"
            self.assertTrue(stale_grid.exists())
            unrelated = root / "data" / "simulation.out"
            unrelated.write_text("keep\n", encoding="utf-8")

            updated = project.to_dict()
            updated["writeExternalGrid"] = False
            export_project(updated, root)

            self.assertFalse(stale_grid.exists())
            self.assertFalse((root / "data" / "geometry_b_001.out").exists())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep\n")
            manifest = json.loads((root / ".snac_grid_generator_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["version"], 1)
            self.assertNotIn("grid/grid_x_b_001.bin", manifest["files"])

    def test_failed_export_does_not_replace_an_existing_case(self) -> None:
        valid = Project.from_dict({"name": "valid", "blocks": [{"id": 1, "ng": [4, 4, 4]}]})
        invalid = Project.from_dict(
            {
                "name": "invalid",
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [2, 0, 0], "lmax": [3, 1, 1]},
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_project(valid, root)
            old_blocks = (root / "blocks.nml").read_bytes()
            old_project = (root / "snac_grid_project.json").read_bytes()

            with self.assertRaises(ValueError):
                export_project(invalid, root)

            self.assertEqual((root / "blocks.nml").read_bytes(), old_blocks)
            self.assertEqual((root / "snac_grid_project.json").read_bytes(), old_project)

    def test_check_rejects_partial_face_contact(self) -> None:
        project = Project.from_dict(
            {
                "name": "partial",
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [1, 0.5, 0], "lmax": [2, 1.5, 1]},
                ],
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("partial" in error for error in result.errors))

    def test_tiny_blocks_have_only_their_true_face_connection(self) -> None:
        size = 1.0e-12
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "lmin": [0.0, 0.0, 0.0], "lmax": [size, size, size]},
                    {"id": 2, "lmin": [size, 0.0, 0.0], "lmax": [2.0 * size, size, size]},
                ]
            }
        )

        updated, result = update_project_structure(project, source_block_id=1)

        self.assertTrue(result.ok, result.errors)
        self.assertEqual(len(result.interfaces), 1)
        self.assertEqual(result.interfaces[0]["face"], "x+")
        self.assertEqual(result.interfaces[0]["neighborFace"], "x-")
        self.assertEqual(updated.blocks[0].cbcpre, ["N", "F", "N", "N", "N", "N"])
        self.assertEqual(updated.blocks[1].cbcpre, ["F", "N", "N", "N", "N", "N"])

    def test_tiny_explicit_grid_must_span_its_actual_extent(self) -> None:
        size = 1.0e-12
        with self.assertRaisesRegex(ValueError, "must span the block extent"):
            Project.from_dict(
                {
                    "blocks": [
                        {
                            "id": 1,
                            "lmin": [0.0, 0.0, 0.0],
                            "lmax": [size, 1.0, 1.0],
                            "axes": {
                                "x": {
                                    "kind": "explicit",
                                    "faces": [0.0, 0.5 * size, 1.5 * size],
                                }
                            },
                        }
                    ]
                }
            )

    def test_check_mirrors_snac_grid_and_initialization_contracts(self) -> None:
        bad_growth = Project.from_dict(
            {"blocks": [{"id": 1, "axes": {"x": {"kind": "snac", "gt": 1, "gr": -2.0}}}]}
        )
        self.assertTrue(any("non-negative" in error for error in check_project(bad_growth).errors))

        bad_mapping = Project.from_dict(
            {"blocks": [{"id": 1, "axes": {"x": {"kind": "snac", "gt": 99, "gr": 1.0}}}]}
        )
        self.assertTrue(any("mapping function 99" in error for error in check_project(bad_mapping).errors))

        bad_initialization = Project.from_dict({"blocks": [{"id": 1, "inivel": "unknown"}]})
        self.assertTrue(
            any("initial velocity field" in error for error in check_project(bad_initialization).errors)
        )
        for field in ("kov", "pdc"):
            with self.subTest(field=field):
                self.assertTrue(check_project(Project.from_dict({"blocks": [{"id": 1, "inivel": field}]})).ok)

    def test_check_rejects_inflow_on_a_periodic_axis(self) -> None:
        project = Project.from_dict(
            {
                "periodicAxes": [True, False, False],
                "blocks": [{"id": 1, "inflow": [1, 0, 0, 0, 0, 0]}],
            }
        )

        _, result = update_project_structure(project)

        self.assertFalse(result.ok)
        self.assertTrue(any("inflow type on periodic face x-" in error for error in result.errors))

    def test_update_structure_infers_friends_and_propagates_mpi(self) -> None:
        project = Project.from_dict(
            {
                "name": "repair",
                "blocks": [
                    {"id": 1, "dims": [1, 3, 2], "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "dims": [4, 1, 1], "ng": [8, 8, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        updated, result = update_project_structure(project, source_block_id=1)
        self.assertTrue(result.ok)
        self.assertEqual(updated.blocks[1].dims, [4, 3, 2])
        self.assertEqual(updated.blocks[0].cbcpre[1], "F")
        self.assertEqual(updated.blocks[1].cbcpre[0], "F")
        self.assertEqual(updated.blocks[0].bcpre[1], 2.0)
        self.assertEqual(updated.blocks[1].bcpre[0], 1.0)

    def test_apply_axis_propagates_only_across_aligned_block_faces(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "ng": [12, 17, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [20, 9, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                    {"id": 3, "ng": [12, 5, 8], "lmin": [0, 1, 0], "lmax": [1, 2, 1]},
                ],
            }
        )
        project.blocks[0].axes["y"].kind = "simple_ratio"
        project.blocks[0].axes["y"].ratio = 3.0

        updated, changed = apply_axis_to_aligned_blocks(project, 1, "y")
        self.assertEqual(changed, [2])
        self.assertEqual(updated.blocks[1].ng[1], 17)
        self.assertEqual(updated.blocks[1].axes["y"].ratio, 3.0)
        self.assertEqual(updated.blocks[2].ng[1], 5)

        transverse, changed = apply_axis_to_aligned_blocks(project, 1, "x")
        self.assertEqual(changed, [3])
        self.assertEqual(transverse.blocks[1].ng[0], 20)
        self.assertEqual(transverse.blocks[2].ng[0], 12)

    def test_repair_propagates_selected_tangential_grid(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "ng": [12, 17, 8],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"y": {"kind": "simple_ratio", "ratio": 3}},
                    },
                    {"id": 2, "ng": [20, 9, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        repaired, result = repair_project_grids(project, source_block_id=1)
        self.assertTrue(result.ok, result.errors)
        self.assertEqual([change.to_dict() for change in result.changes], [
            {"blockId": 2, "axis": "y", "sourceBlockId": 1, "oldN": 9, "newN": 17}
        ])
        self.assertEqual(repaired.blocks[1].ng, [20, 17, 8])
        self.assertEqual(repaired.blocks[1].axes["y"].ratio, 3.0)

    def test_repair_respects_locks_and_rejects_locked_conflicts(self) -> None:
        raw = {
            "blocks": [
                {
                    "id": 1,
                    "ng": [8, 12, 8],
                    "lmin": [0, 0, 0],
                    "lmax": [1, 1, 1],
                    "axes": {"y": {"kind": "simple_ratio", "ratio": 2}},
                },
                {
                    "id": 2,
                    "ng": [8, 20, 8],
                    "lmin": [1, 0, 0],
                    "lmax": [2, 1, 1],
                    "axes": {"y": {"kind": "simple_ratio", "ratio": 4}},
                    "axisLocks": [False, True, False],
                },
            ],
        }
        repaired, result = repair_project_grids(raw, source_block_id=1)
        self.assertTrue(result.ok, result.errors)
        self.assertEqual(repaired.blocks[0].ng[1], 20)
        self.assertEqual(result.changes[0].source_block_id, 2)

        raw["blocks"][0]["axisLocks"] = [False, True, False]
        raw["blocks"][1]["ng"][2] = 10
        failed, result = repair_project_grids(raw, source_block_id=1)
        self.assertFalse(result.ok)
        self.assertTrue(any("locked y grids conflict" in error for error in result.errors))
        self.assertEqual(failed.blocks[1].ng[2], 10)
        self.assertEqual(result.changes, [])

    def test_repair_matches_normal_interface_spacing(self) -> None:
        project = Project.from_dict(
            {
                "writeExternalGrid": False,
                "blocks": [
                    {"id": 1, "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {
                        "id": 2,
                        "ng": [16, 8, 8],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 100}},
                    },
                ],
            }
        )
        repaired, result = repair_project_grids(project, source_block_id=1)
        self.assertTrue(result.ok, result.errors)
        changes = [change for change in result.changes if change.kind == "spacing"]
        self.assertEqual(len(changes), 1)
        self.assertEqual((changes[0].block_id, changes[0].face, changes[0].source_block_id), (2, "x-", 1))
        self.assertAlmostEqual(changes[0].new_spacing, 0.125)
        self.assertEqual(repaired.blocks[1].axes["x"].kind, "simple_ratio")
        self.assertEqual(repaired.blocks[1].axes["x"].controls, ["n", "width_start"])
        self.assertTrue(repaired.write_external_grid)
        self.assertFalse(any("interface spacing jumps" in warning for warning in check_project(repaired).warnings))

    def test_repair_uses_a_locked_normal_grid_as_authority(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {
                        "id": 2,
                        "ng": [16, 8, 8],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 100}},
                        "axisLocks": [True, False, False],
                    },
                ],
            }
        )
        repaired, result = repair_project_grids(project, source_block_id=1)
        self.assertTrue(result.ok, result.errors)
        change = next(change for change in result.changes if change.kind == "spacing")
        self.assertEqual((change.block_id, change.face, change.source_block_id), (1, "x+", 2))
        self.assertEqual(repaired.blocks[0].axes["x"].kind, "explicit")
        self.assertEqual(repaired.blocks[1].axes["x"].kind, "simple_ratio")

    def test_repair_rejects_a_normal_jump_between_locked_grids(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "ng": [8, 8, 8],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axisLocks": [True, False, False],
                    },
                    {
                        "id": 2,
                        "ng": [16, 8, 8],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 100}},
                        "axisLocks": [True, False, False],
                    },
                ],
            }
        )
        unchanged, result = repair_project_grids(project, source_block_id=1)
        self.assertFalse(result.ok)
        self.assertTrue(any("locked x grids have an interface spacing jump" in error for error in result.errors))
        self.assertEqual(unchanged.to_dict(), project.to_dict())

    def test_decomposition_balances_unequal_blocks_exactly(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 20, "mode": "1d", "axes": [True, False, False]},
                "blocks": [
                    {"id": 1, "ng": [128, 64, 64], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [32, 64, 64], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        decomposed, result = optimize_project_decomposition(project)
        self.assertTrue(result.ok, result.errors)
        self.assertEqual([block.dims for block in decomposed.blocks], [[16, 1, 1], [4, 1, 1]])
        self.assertEqual(result.total_ranks, 20)
        self.assertEqual(result.active_axes, ("x",))
        self.assertAlmostEqual(result.imbalance, 1.0)
        _, structured = update_project_structure(decomposed)
        self.assertTrue(structured.ok, structured.errors)

    def test_decomposition_modes_obey_topology(self) -> None:
        for mode, expected_axes in (("1d", 1), ("2d", 2), ("3d", 3)):
            project = Project.from_dict(
                {
                    "decomposition": {"targetRanks": 16, "mode": mode, "axes": [True, True, True]},
                    "blocks": [
                        {"id": 1, "ng": [64, 32, 16], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                        {"id": 2, "ng": [64, 32, 16], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                    ],
                }
            )
            with self.subTest(mode=mode):
                decomposed, result = optimize_project_decomposition(project)
                self.assertTrue(result.ok, result.errors)
                self.assertEqual(len(result.active_axes), expected_axes)
                self.assertEqual(sum(np.prod(block.dims) for block in decomposed.blocks), 16)
                self.assertEqual(decomposed.blocks[0].dims[1:], decomposed.blocks[1].dims[1:])

    def test_decomposition_allows_one_serial_rank_per_block_in_any_mode(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 2, "mode": "3d", "axes": [True, True, True]},
                "blocks": [
                    {"id": 1, "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [8, 8, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        decomposed, result = optimize_project_decomposition(project)
        self.assertTrue(result.ok, result.errors)
        self.assertEqual(result.active_axes, ())
        self.assertEqual([block.dims for block in decomposed.blocks], [[1, 1, 1], [1, 1, 1]])

    def test_decomposition_rejects_disconnected_blocks(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 4, "mode": "auto"},
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [2, 0, 0], "lmax": [3, 1, 1]},
                ],
            }
        )
        _, result = optimize_project_decomposition(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("disconnected" in error for error in result.errors))

    def test_decomposition_reports_nearby_counts_when_exact_is_impossible(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 5, "mode": "1d", "axes": [False, True, False]},
                "blocks": [
                    {"id": 1, "ng": [32, 32, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [32, 32, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        _, result = optimize_project_decomposition(project)
        self.assertFalse(result.ok)
        self.assertIn(4, result.nearby_ranks)
        self.assertIn(6, result.nearby_ranks)

    def test_large_decompositions_remain_balanced_and_well_shaped(self) -> None:
        for target in (512, 1024, 2048):
            project = Project.from_dict(
                {
                    "decomposition": {"targetRanks": target, "mode": "3d"},
                    "blocks": [
                        {
                            "id": block_id,
                            "ng": [nx, 128, 128],
                            "lmin": [block_id - 1, 0, 0],
                            "lmax": [block_id, 1, 1],
                        }
                        for block_id, nx in enumerate((128, 64, 128, 64), 1)
                    ],
                }
            )
            with self.subTest(target=target):
                decomposed, result = optimize_project_decomposition(project, time_limit=0.5)
                self.assertTrue(result.ok, result.errors)
                self.assertEqual(sum(np.prod(block.dims) for block in decomposed.blocks), target)
                self.assertLessEqual(result.imbalance, 1.1)
                self.assertLessEqual(result.max_local_aspect, 2.0)
                for block in decomposed.blocks:
                    for cells, partitions in zip(block.ng, block.dims):
                        if partitions > 1:
                            self.assertGreaterEqual(cells // partitions, MIN_LOCAL_CELLS)

    def test_check_enforces_requested_rank_count_and_disabled_axes(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 8, "mode": "1d", "axes": [True, False, False]},
                "blocks": [{"id": 1, "dims": [2, 2, 1], "ng": [16, 16, 8]}],
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("uses 4 ranks, expected 8" in error for error in result.errors))
        self.assertTrue(any("disabled y" in error for error in result.errors))

    def test_export_rejects_stale_requested_rank_count(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {"targetRanks": 8, "mode": "auto"},
                "blocks": [{"id": 1, "dims": [2, 2, 1], "ng": [16, 16, 8]}],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "uses 4 ranks, expected 8"):
                export_project(project, tmp)

    def test_check_rejects_broken_manual_friend_value(self) -> None:
        project = Project.from_dict(
            {
                "name": "broken-friend",
                "blocks": [
                    {"id": 1, "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [8, 8, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        updated, result = update_project_structure(project)
        self.assertTrue(result.ok)
        updated.blocks[0].bcpre[1] = 99.0
        result = check_project(updated)
        self.assertFalse(result.ok)
        self.assertTrue(any("references missing block 99" in error for error in result.errors))

    def test_periodic_axis_sets_self_friend_boundaries(self) -> None:
        project = Project.from_dict(
            {
                "name": "periodic",
                "periodicAxes": [False, False, True],
                "blocks": [
                    {"id": 1, "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [8, 8, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        updated, result = update_project_structure(project)
        self.assertTrue(result.ok)
        for block in updated.blocks:
            self.assertEqual(block.cbcpre[4:6], ["F", "F"])
            self.assertEqual(block.bcpre[4:6], [float(block.id), float(block.id)])
            for component in range(3):
                self.assertEqual(block.cbcvel[component][4:6], ["F", "F"])
                self.assertEqual(block.bcvel[component][4:6], [float(block.id), float(block.id)])

    def test_periodic_axis_wraps_a_multi_block_chain(self) -> None:
        project = Project.from_dict(
            {
                "name": "periodic-chain",
                "periodicAxes": [True, False, False],
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        updated, result = update_project_structure(project)
        self.assertTrue(result.ok)
        self.assertEqual(updated.blocks[0].bcpre[0:2], [2.0, 2.0])
        self.assertEqual(updated.blocks[1].bcpre[0:2], [1.0, 1.0])

    def test_periodic_axis_rejects_unmatched_outer_faces(self) -> None:
        project = Project.from_dict(
            {
                "name": "bad-periodic",
                "periodicAxes": [True, False, False],
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [2, 1, 0], "lmax": [3, 2, 1]},
                ],
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("periodic x" in error for error in result.errors))

    def test_periodic_axis_can_be_interrupted_by_a_hole(self) -> None:
        extents = [
            (1, [0, 0, 0], [0.5, 1, 1], [4, 4, 4]),
            (2, [0.5, 0, 0], [1.5, 1, 1], [8, 4, 4]),
            (3, [1.5, 0, 0], [2, 1, 1], [4, 4, 4]),
            (4, [0, 1, 0], [0.5, 2, 1], [4, 4, 4]),
            (5, [1.5, 1, 0], [2, 2, 1], [4, 4, 4]),
            (6, [0, 2, 0], [0.5, 3, 1], [4, 4, 4]),
            (7, [0.5, 2, 0], [1.5, 3, 1], [8, 4, 4]),
            (8, [1.5, 2, 0], [2, 3, 1], [4, 4, 4]),
        ]
        project = Project.from_dict(
            {
                "name": "periodic-hole",
                "periodicAxes": [True, False, False],
                "blocks": [
                    {"id": block_id, "lmin": lmin, "lmax": lmax, "ng": ng}
                    for block_id, lmin, lmax, ng in extents
                ],
            }
        )
        updated, result = update_project_structure(project)
        self.assertTrue(result.ok, result.errors)
        by_id = {block.id: block for block in updated.blocks}
        self.assertEqual(by_id[1].bcpre[0], 3.0)
        self.assertEqual(by_id[4].bcpre[0], 5.0)
        self.assertEqual(by_id[6].bcpre[0], 8.0)
        self.assertEqual(by_id[4].cbcpre[1], "N")
        self.assertEqual(by_id[5].cbcpre[0], "N")
        self.assertEqual(by_id[4].cbcvel[0][1], "D")
        self.assertEqual(by_id[5].cbcvel[0][0], "D")

    def test_export_rejects_unstructured_partitioning(self) -> None:
        project = Project.from_dict(
            {
                "name": "bad-partition",
                "inferConnectivity": True,
                "blocks": [
                    {"id": 1, "dims": [1, 2, 1], "ng": [8, 8, 8], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "dims": [1, 3, 1], "ng": [8, 8, 8], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                export_project(project, tmp)
            self.assertFalse((Path(tmp) / "blocks.nml").exists())

    def test_spacing_jump_becomes_warning(self) -> None:
        project = Project.from_dict(
            {
                "name": "spacing-warning",
                "inferConnectivity": True,
                "blocks": [
                    {
                        "id": 1,
                        "ng": [8, 8, 8],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 8}},
                    },
                    {
                        "id": 2,
                        "ng": [8, 8, 8],
                        "lmin": [1, 0, 0],
                        "lmax": [2, 1, 1],
                        "axes": {"x": {"kind": "simple_ratio", "ratio": 8}},
                    },
                ],
            }
        )
        updated, result = update_project_structure(project)
        self.assertTrue(result.ok)
        self.assertTrue(any("interface spacing jumps" in warning for warning in result.warnings))
        self.assertEqual(len(result.interfaces), 1)
        self.assertGreater(result.interfaces[0]["ratio"], 3.0)

    def test_geometry_operations_preserve_explicit_grids(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "explicit", "faces": [0.0, 0.2, 1.0]}},
                    }
                ]
            }
        )
        duplicated, new_ids = duplicate_blocks(project, [1], "x", count=2, gap=0.5)
        self.assertEqual(new_ids, [2, 3])
        self.assertEqual(duplicated.blocks[1].axes["x"].faces, [1.5, 1.7, 2.5])
        mirrored, mirror_ids = mirror_blocks(project, [1], "x", 0.0)
        self.assertEqual(mirror_ids, [2])
        self.assertEqual(mirrored.blocks[1].axes["x"].faces, [-1.0, -0.2, 0.0])

        aligned, changed = align_blocks(duplicated, [1, 2], 1, "y", "upper")
        self.assertEqual(changed, [2])
        snapped, changed = snap_block_to_face(aligned, 2, 1, "x+")
        self.assertEqual(changed, [2])
        self.assertEqual(snapped.blocks[1].lmin, [1.0, 0.0, 0.0])

    def test_mirror_reverses_native_and_monotone_grading(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                        "axes": {"x": {"kind": "snac", "gt": 1, "gr": 2.0}},
                    },
                    {
                        "id": 2,
                        "lmin": [2, 0, 0],
                        "lmax": [3, 1, 1],
                        "axes": {
                            "x": {
                                "kind": "simple_ratio",
                                "controls": ["n", "ratio"],
                                "ratio": 4.0,
                            }
                        },
                    },
                ]
            }
        )
        mirrored, _ = mirror_blocks(project, [1, 2], "x", 0.0)
        self.assertEqual(mirrored.blocks[2].axes["x"].gt, 3)
        self.assertAlmostEqual(mirrored.blocks[3].axes["x"].ratio, 0.25)

    def test_snac_boundary_presets_keep_friend_faces(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "cbcvel": [["D", "F", "D", "D", "D", "D"]] * 3,
                        "cbcpre": ["N", "F", "N", "N", "N", "N"],
                        "bcvel": [[0, 2, 0, 0, 0, 0]] * 3,
                        "bcpre": [0, 2, 0, 0, 0, 0],
                    },
                    {"id": 2},
                ]
            }
        )
        updated, changed = apply_bc_preset(project, [1], "moving_wall", face="z+", velocity=[1, 0, 0])
        block = updated.blocks[0]
        self.assertEqual(changed, [1])
        self.assertEqual([block.cbcvel[index][5] for index in range(3)], ["D", "D", "D"])
        self.assertEqual([block.bcvel[index][5] for index in range(3)], [1.0, 0.0, 0.0])
        self.assertEqual(block.cbcpre[1], "F")
        self.assertEqual(block.bcpre[1], 2.0)

    def test_decomposition_respects_minimum_local_cells(self) -> None:
        project = Project.from_dict(
            {
                "decomposition": {
                    "targetRanks": 8,
                    "mode": "1d",
                    "axes": [True, False, False],
                    "minLocalCells": 10,
                },
                "blocks": [{"id": 1, "ng": [64, 8, 8]}],
            }
        )
        _, result = optimize_project_decomposition(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("no exact" in error for error in result.errors))

        aspect_limited = Project.from_dict(
            {
                "decomposition": {
                    "targetRanks": 8,
                    "mode": "1d",
                    "axes": [True, False, False],
                    "maxLocalAspect": 4.0,
                },
                "blocks": [{"id": 1, "ng": [64, 64, 4]}],
            }
        )
        _, result = optimize_project_decomposition(aspect_limited)
        self.assertFalse(result.ok)
        self.assertTrue(any("no exact" in error for error in result.errors))

    def test_headless_cli_checks_and_migrates_projects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            project_path = Path(tmp) / "old.json"
            migrated_path = Path(tmp) / "migrated.json"
            project_path.write_text(json.dumps({"schemaVersion": 1, "blocks": [{"id": 1}]}), encoding="utf-8")
            with redirect_stdout(StringIO()) as output:
                self.assertEqual(grid_cli_main(["check", str(project_path), "--json"]), 0)
                self.assertTrue(json.loads(output.getvalue())["ok"])
            with redirect_stdout(StringIO()):
                self.assertEqual(
                    grid_cli_main(["migrate", str(project_path), "-o", str(migrated_path), "--json"]),
                    0,
                )
            migrated = json.loads(migrated_path.read_text(encoding="utf-8"))
            self.assertEqual(migrated["schemaVersion"], 2)


if __name__ == "__main__":
    unittest.main()
