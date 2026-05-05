from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from utils.snac_grid_generator import (
    Project,
    axis_grid_arrays,
    check_project,
    export_project,
    solve_spacing,
    update_project_structure,
)


class GridGeneratorTests(unittest.TestCase):
    def test_blockmesh_spacing_relations(self) -> None:
        solution = solve_spacing(length=7.0, n=3, ratio=4.0)
        self.assertAlmostEqual(solution.cell_ratio, 2.0)
        self.assertAlmostEqual(solution.width_start, 1.0)
        self.assertAlmostEqual(solution.width_end, 4.0)

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
        self.assertEqual(arrays.binary_payload.shape, (120,))
        self.assertAlmostEqual(arrays.faces[0], 0.0)
        self.assertAlmostEqual(arrays.faces[-2], 10.0)
        self.assertTrue(np.all(np.diff(arrays.faces[:-1]) > 0.0))

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

    def test_empty_project_exports_project_file(self) -> None:
        project = Project.from_dict({"name": "empty", "blocks": []})
        with tempfile.TemporaryDirectory() as tmp:
            result = export_project(project, tmp)
            root = Path(tmp)
            self.assertTrue((root / "blocks.nml").exists())
            self.assertTrue((root / "snac_grid_project.json").exists())
            self.assertEqual(len(result.files), 2)

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

    def test_periodic_axis_requires_common_extent(self) -> None:
        project = Project.from_dict(
            {
                "name": "bad-periodic",
                "periodicAxes": [True, False, False],
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ],
            }
        )
        result = check_project(project)
        self.assertFalse(result.ok)
        self.assertTrue(any("periodic x" in error for error in result.errors))

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


if __name__ == "__main__":
    unittest.main()
