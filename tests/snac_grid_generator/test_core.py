from __future__ import annotations

import json
from contextlib import redirect_stdout
from io import StringIO
import tempfile
import unittest
from unittest.mock import patch
import warnings
from pathlib import Path

import numpy as np

from utils.snac_grid_generator import (
    MIN_LOCAL_CELLS,
    Project,
    align_blocks,
    analyze_grid_quality,
    apply_bc_preset,
    apply_axis_to_aligned_blocks,
    axis_grid_arrays,
    axis_grid_diagnostics,
    check_project,
    build_case_report,
    duplicate_blocks,
    export_project,
    fit_monotone_spacing,
    import_case,
    optimize_project_decomposition,
    mirror_blocks,
    repair_project_grids,
    render_case_report_markdown,
    solve_spacing,
    snap_block_to_face,
    update_project_structure,
    write_case_report,
)
from utils.snac_grid_generator.snac_grid import binary_payload, read_grid_binary
from utils.snac_grid_generator.cli import main as grid_cli_main
from utils.snac_grid_generator.topology import (
    FaceConnection,
    build_topology,
    infer_global_index_extents,
)


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
            arrays = axis_grid_arrays(
                {"kind": "simple_ratio", "ratio": 4.0, "side": "start", "profile": profile},
                0.0,
                1.0,
                32,
            )
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

    def test_export_matches_selected_snac_precision(self) -> None:
        for precision, dtype in (("single", "<f4"), ("double", "<f8")):
            with self.subTest(precision=precision), tempfile.TemporaryDirectory() as tmp:
                project = Project.from_dict(
                    {
                        "externalGridPrecision": precision,
                        "blocks": [{"id": 1, "ng": [4, 3, 2]}],
                    }
                )
                export_project(project, tmp)
                path = Path(tmp) / "grid" / "grid_x_b_001.bin"
                payload = np.fromfile(path, dtype=dtype)
                decoded = read_grid_binary(path, 4, 0.0, 1.0)

                self.assertEqual(payload.shape, (16,))
                self.assertEqual(path.stat().st_size, 16 * np.dtype(dtype).itemsize)
                self.assertEqual(decoded.dtype, dtype)
                self.assertEqual(
                    json.loads((Path(tmp) / "snac_grid_project.json").read_text())["externalGridPrecision"],
                    precision,
                )

    def test_single_precision_export_rejects_collapsed_coordinates(self) -> None:
        project = Project.from_dict(
            {
                "externalGridPrecision": "single",
                "blocks": [
                    {
                        "id": 1,
                        "ng": [4, 3, 2],
                        "lmin": [1.0e9, 0.0, 0.0],
                        "lmax": [1.0e9 + 1.0e-3, 1.0, 1.0],
                    }
                ],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "case"
            with self.assertRaisesRegex(ValueError, "cannot represent all face coordinates"):
                export_project(project, output)
            self.assertTrue(output.is_dir())
            self.assertFalse(any(output.iterdir()))

    def test_legacy_data_only_destination_exports_startup_grids(self) -> None:
        project = Project.from_dict(
            {
                "externalGridSource": "data",
                "blocks": [{"id": 1, "ng": [4, 3, 2]}],
            }
        )
        self.assertEqual(project.external_grid_source, "grid")

        with tempfile.TemporaryDirectory() as tmp:
            export_project(project, tmp)
            self.assertTrue((Path(tmp) / "grid" / "grid_x_b_001.bin").is_file())
            self.assertFalse((Path(tmp) / "data" / "grid_x_b_001.bin").exists())

    def test_data_only_case_import_exports_startup_grids(self) -> None:
        project = Project.from_dict(
            {
                "externalGridSource": "both",
                "externalGridPrecision": "single",
                "blocks": [{"id": 1, "ng": [4, 3, 2]}],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            destination = root / "destination"
            export_project(project, source)
            (source / "snac_grid_project.json").unlink()
            for path in (source / "grid").iterdir():
                path.unlink()
            (source / "grid").rmdir()

            imported = import_case(source).project
            export_project(imported, destination)

            self.assertEqual(imported.external_grid_source, "grid")
            self.assertEqual(imported.external_grid_precision, "single")
            startup_grid = destination / "grid" / "grid_x_b_001.bin"
            self.assertTrue(startup_grid.is_file())
            self.assertEqual(startup_grid.stat().st_size, 16 * np.dtype("<f4").itemsize)

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

    def test_import_reads_scalar_count_from_dns_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp)
            (case_dir / "dns.nml").write_text("&dns\nnscal = 2\n/\n", encoding="utf-8")
            (case_dir / "blocks.nml").write_text(
                "&blocks\n"
                "nblocks = 1\n"
                "block_ng(1:3,1) = 4, 4, 4\n"
                "block_lmin(1:3,1) = 0., 0., 0.\n"
                "block_lmax(1:3,1) = 1., 1., 1.\n"
                "/\n",
                encoding="utf-8",
            )
            project = import_case(case_dir).project

        self.assertEqual(project.nscal, 2)
        self.assertEqual(project.blocks[0].cbcscal, [["N"] * 6, ["N"] * 6])
        self.assertEqual(project.blocks[0].bcscal, [[0.0] * 6, [0.0] * 6])

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

    def test_binary_grid_reader_validates_all_stored_arrays(self) -> None:
        project = Project.from_dict({"blocks": [{"id": 1, "ng": [4, 4, 4]}]})
        arrays = axis_grid_arrays(project.blocks[0].axes["x"].to_dict(), 0.0, 1.0, 4)
        payload = binary_payload(arrays).copy()
        payload[4] += 0.1

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "grid_x_b_001.bin"
            payload.tofile(path)
            with self.assertRaisesRegex(ValueError, "cell centers"):
                read_grid_binary(path, 4, 0.0, 1.0)

    def test_binary_grid_reader_accepts_big_endian_single_precision(self) -> None:
        project = Project.from_dict({"blocks": [{"id": 1, "ng": [4, 4, 4]}]})
        arrays = axis_grid_arrays(project.blocks[0].axes["x"].to_dict(), 0.0, 1.0, 4)

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "grid_x_b_001.bin"
            binary_payload(arrays).astype(">f4").tofile(path)
            decoded = read_grid_binary(path, 4, 0.0, 1.0)

        self.assertEqual(decoded.dtype, ">f4")
        np.testing.assert_allclose(decoded.faces, arrays.faces[:-1], rtol=0.0, atol=1.0e-6)

    def test_case_import_falls_back_to_valid_duplicate_grid(self) -> None:
        project = Project.from_dict(
            {
                "writeExternalGrid": True,
                "externalGridSource": "both",
                "blocks": [{"id": 1, "ng": [4, 4, 4]}],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_project(project, root)
            (root / "snac_grid_project.json").unlink()
            grid_path = root / "grid" / "grid_x_b_001.bin"
            payload = np.fromfile(grid_path, dtype="<f8")
            payload[4] += 0.1
            payload.tofile(grid_path)

            result = import_case(root)

        self.assertEqual(result.project.blocks[0].axes["x"].kind, "explicit")
        self.assertTrue(any("ignored invalid grid/ grid" in warning for warning in result.warnings))

    def test_case_round_trip_preserves_namelist_and_binary_grids(self) -> None:
        project = Project.from_dict(
            {
                "name": "round-trip",
                "writeExternalGrid": True,
                "blocks": [
                    {
                        "id": 1,
                        "ng": [5, 4, 3],
                        "axes": {
                            "x": {
                                "kind": "explicit",
                                "faces": [0.0, 0.08, 0.23, 0.47, 0.72, 1.0],
                            }
                        },
                    }
                ],
            }
        )
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_path = Path(first)
            second_path = Path(second)
            export_project(project, first_path)
            (first_path / "snac_grid_project.json").unlink()
            imported = import_case(first_path).project
            export_project(imported, second_path)

            self.assertEqual(
                (first_path / "blocks.nml").read_bytes(),
                (second_path / "blocks.nml").read_bytes(),
            )
            for axis in ("x", "y", "z"):
                relative = Path("grid") / f"grid_{axis}_b_001.bin"
                self.assertEqual(
                    (first_path / relative).read_bytes(),
                    (second_path / relative).read_bytes(),
                )

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

    def test_generated_block_lattice_indices_ignore_ids_and_input_order(self) -> None:
        x_cells = [2, 5, 3, 4]
        y_cells = [3, 6, 2]
        ids = [9, 2, 12, 5, 7, 1, 11, 4, 10, 3, 8, 6]
        records = []
        expected = {}
        for j, ny in enumerate(y_cells):
            for i, nx in enumerate(x_cells):
                block_id = ids[j * len(x_cells) + i]
                records.append(
                    {
                        "id": block_id,
                        "ng": [nx, ny, 2],
                        "lmin": [float(i), float(j), 0.0],
                        "lmax": [float(i + 1), float(j + 1), 1.0],
                    }
                )
                lo = [sum(x_cells[:i]) + 1, sum(y_cells[:j]) + 1, 1]
                expected[block_id] = (lo, [lo[0] + nx - 1, lo[1] + ny - 1, 2])

        blocks = Project.from_dict({"blocks": list(reversed(records))}).blocks
        self.assertEqual(infer_global_index_extents(blocks), expected)

    def test_generated_grading_profiles_remain_monotone(self) -> None:
        for profile in ("geometric", "tanh", "erf"):
            for n in (2, 3, 8, 31):
                for ratio in (0.1, 1.0, 10.0):
                    with self.subTest(profile=profile, n=n, ratio=ratio):
                        arrays = axis_grid_arrays(
                            {
                                "kind": "simple_ratio",
                                "profile": profile,
                                "controls": ["n", "ratio"],
                                "ratio": ratio,
                            },
                            -2.0,
                            3.0,
                            n,
                        )
                        faces = arrays.faces[:-1]
                        widths = np.diff(faces)
                        self.assertTrue(np.all(np.isfinite(widths)))
                        self.assertTrue(np.all(widths > 0.0))
                        self.assertAlmostEqual(float(faces[0]), -2.0)
                        self.assertAlmostEqual(float(faces[-1]), 3.0)
                        self.assertAlmostEqual(float(widths[-1] / widths[0]), ratio, places=8)

    def test_large_coordinate_offset_with_small_cells_remains_structured(self) -> None:
        origin = 1.0e9
        length = 1.0e-3
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 8,
                        "ng": [32, 4, 4],
                        "lmin": [origin, 0.0, 0.0],
                        "lmax": [origin + length, 1.0, 1.0],
                    },
                    {
                        "id": 3,
                        "ng": [48, 4, 4],
                        "lmin": [origin + length, 0.0, 0.0],
                        "lmax": [origin + 2.0 * length, 1.0, 1.0],
                    },
                ]
            }
        )
        extents = infer_global_index_extents(project.blocks)
        self.assertEqual(extents[8][0], [1, 1, 1])
        self.assertEqual(extents[3][0], [33, 1, 1])

        arrays = axis_grid_arrays(
            {"kind": "simple_ratio", "ratio": 4.0},
            origin,
            origin + length,
            32,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "offset-grid.bin"
            binary_payload(arrays).tofile(path)
            decoded = read_grid_binary(path, 32, origin, origin + length)
            export_project(project, root / "case")
            (root / "case" / "snac_grid_project.json").unlink()
            imported = import_case(root / "case").project
        np.testing.assert_array_equal(decoded.faces, arrays.faces[:-1])
        self.assertEqual(
            sorted((block.lmin, block.lmax) for block in imported.blocks),
            sorted((block.lmin, block.lmax) for block in project.blocks),
        )

    def test_malformed_namelists_fail_with_clear_errors(self) -> None:
        cases = {
            "missing group": "nblocks=1\n/\n",
            "incomplete group": "&blocks\nnblocks=1\n",
            "missing nblocks": "&blocks\nblock_ng(1)=4,4,4\n/\n",
            "empty assignment": "&blocks\nnblocks=\n/\n",
            "unterminated quote": "&blocks\nnblocks=1\nblock_inivel(1)='zer\n/\n",
        }
        for label, content in cases.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                (Path(tmp) / "blocks.nml").write_text(content, encoding="utf-8")
                with self.assertRaises(ValueError):
                    import_case(tmp)

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

    def test_check_returns_navigable_interface_diagnostics(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "ng": [4, 4, 4], "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "ng": [4, 5, 4], "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                ]
            }
        )
        diagnostics = check_project(project).to_dict()["diagnostics"]
        diagnostic = next(item for item in diagnostics if item["code"] == "interface_cell_count")

        self.assertEqual(diagnostic["axis"], "y")
        self.assertEqual(
            diagnostic["locations"],
            [{"blockId": 1, "face": "x+"}, {"blockId": 2, "face": "x-"}],
        )

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

    def test_export_refuses_unmanaged_collisions_and_bad_manifests(self) -> None:
        project = Project.from_dict({"blocks": [{"id": 1, "ng": [4, 4, 4]}]})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            blocks = root / "blocks.nml"
            blocks.write_text("manual case\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unmanaged file"):
                export_project(project, root)
            self.assertEqual(blocks.read_text(encoding="utf-8"), "manual case\n")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".snac_grid_generator_manifest.json").write_text(
                "not json\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "cannot read generator manifest"):
                export_project(project, root)

    def test_export_rolls_back_after_publication_failure(self) -> None:
        original = Project.from_dict(
            {
                "name": "before",
                "writeExternalGrid": True,
                "blocks": [{"id": 1, "ng": [4, 4, 4]}],
            }
        )
        updated = Project.from_dict(
            {
                "name": "after",
                "writeExternalGrid": True,
                "blocks": [{"id": 1, "ng": [5, 4, 4]}],
            }
        )
        real_replace = Path.replace

        def fail_during_publish(source: Path, target: Path) -> Path:
            if ".snac_grid_stage_" in source.as_posix() and source.name == "snac_grid_project.json":
                raise OSError("injected publication failure")
            return real_replace(source, target)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_project(original, root)
            before = {
                path.relative_to(root): path.read_bytes()
                for path in root.rglob("*")
                if path.is_file()
            }
            with patch.object(Path, "replace", fail_during_publish):
                with self.assertRaisesRegex(OSError, "injected publication failure"):
                    export_project(updated, root)
            after = {
                path.relative_to(root): path.read_bytes()
                for path in root.rglob("*")
                if path.is_file()
            }

        self.assertEqual(after, before)

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

    def test_topology_candidate_sweep_preserves_contacts_and_overlaps(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {"id": 1, "lmin": [0, 0, 0], "lmax": [1, 1, 1]},
                    {"id": 2, "lmin": [1, 0, 0], "lmax": [2, 1, 1]},
                    {"id": 3, "lmin": [0, 1, 0], "lmax": [1, 2, 1]},
                    {"id": 4, "lmin": [0, 0, 1], "lmax": [1, 1, 2]},
                    {"id": 5, "lmin": [1, 0.5, 0], "lmax": [2, 1.5, 1]},
                    {"id": 6, "lmin": [0.25, 0.25, 0.25], "lmax": [0.75, 0.75, 0.75]},
                    {"id": 7, "lmin": [3, 3, 3], "lmax": [4, 4, 4]},
                ]
            }
        )

        topology = build_topology(project.blocks)

        connected = {
            (connection.a_id, connection.b_id, connection.axis_index)
            for connection in topology.connections
        }
        self.assertIn((1, 2, 0), connected)
        self.assertIn((1, 3, 1), connected)
        self.assertIn((1, 4, 2), connected)
        self.assertTrue(any("blocks 1 and 5 touch on a partial x face" == error for error in topology.errors))
        self.assertTrue(any("blocks 1 and 6 overlap" == error for error in topology.errors))
        self.assertFalse(any("block 7" in error for error in topology.errors))

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
                    "decomposition": {
                        "targetRanks": target,
                        "mode": "3d",
                        "maxLocalAspect": 2.0,
                    },
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

    def test_grid_quality_reports_exact_block_and_interface_metrics(self) -> None:
        project = Project.from_dict(
            {
                "blocks": [
                    {
                        "id": 1,
                        "name": "fine",
                        "ng": [4, 2, 2],
                        "dims": [3, 1, 1],
                        "lmin": [0, 0, 0],
                        "lmax": [1, 1, 1],
                    },
                    {
                        "id": 2,
                        "name": "coarse",
                        "ng": [2, 2, 2],
                        "lmin": [1, 0, 0],
                        "lmax": [3, 1, 1],
                    },
                ]
            }
        )
        project, check = update_project_structure(project, source_block_id=1)
        self.assertTrue(check.ok)

        quality = analyze_grid_quality(project).to_dict()
        summary = quality["summary"]
        self.assertEqual(summary["totalCells"], 24)
        self.assertEqual(summary["totalRanks"], 4)
        self.assertEqual((summary["cellsPerRankMin"], summary["cellsPerRankMax"]), (4, 8))
        self.assertEqual(summary["meanCellsPerRank"], 6.0)
        self.assertEqual(summary["loadImbalance"], 2.0)
        self.assertEqual((summary["spacingMin"], summary["spacingMax"]), (0.25, 1.0))
        self.assertEqual(summary["worstCellAspect"], 2.0)
        self.assertEqual(summary["maxInterfaceRatio"], 4.0)
        self.assertEqual(summary["coordinateBytes"], 160)
        self.assertEqual(quality["interfaces"][0]["ratio"], 4.0)

    def test_case_reports_are_deterministic_and_exportable(self) -> None:
        project = Project.from_dict({"name": "report-case", "blocks": [{"id": 1, "ng": [4, 3, 2]}]})
        report = build_case_report(project)
        markdown = render_case_report_markdown(report)

        self.assertTrue(report["validation"]["ok"])
        self.assertEqual(report["storage"]["snacBinaryBytesPerCopy"], 288)
        self.assertEqual(report["storage"]["snacBinaryPrecision"], "double")
        self.assertEqual(markdown, render_case_report_markdown(build_case_report(project)))
        self.assertIn("# SNaC grid report: report-case", markdown)
        self.assertIn("## Axis Quality", markdown)

        single = Project.from_dict(
            {"externalGridPrecision": "single", "blocks": [{"id": 1, "ng": [4, 3, 2]}]}
        )
        self.assertEqual(build_case_report(single)["storage"]["snacBinaryBytesPerCopy"], 144)

        with tempfile.TemporaryDirectory() as tmp:
            json_path = write_case_report(project, Path(tmp) / "quality.json")
            markdown_path = write_case_report(project, Path(tmp) / "quality.md")
            self.assertEqual(json.loads(json_path.read_text(encoding="utf-8")), report)
            self.assertEqual(markdown_path.read_text(encoding="utf-8"), markdown)

    def test_headless_cli_writes_quality_reports(self) -> None:
        project = Project.from_dict({"name": "cli-report", "blocks": [{"id": 1, "ng": [4, 3, 2]}]})
        with tempfile.TemporaryDirectory() as tmp:
            project_path = Path(tmp) / "project.json"
            report_path = Path(tmp) / "report.md"
            project_path.write_text(json.dumps(project.to_dict()), encoding="utf-8")

            with redirect_stdout(StringIO()) as output:
                self.assertEqual(
                    grid_cli_main(["report", str(project_path), "--format", "json"]),
                    0,
                )
                payload = json.loads(output.getvalue())
            self.assertEqual(payload["quality"]["summary"]["totalCells"], 24)

            with redirect_stdout(StringIO()):
                self.assertEqual(
                    grid_cli_main(
                        [
                            "report",
                            str(project_path),
                            "--format",
                            "markdown",
                            "-o",
                            str(report_path),
                        ]
                    ),
                    0,
                )
            self.assertIn("## SNaC Storage", report_path.read_text(encoding="utf-8"))

            invalid_path = Path(tmp) / "invalid.json"
            invalid_report_path = Path(tmp) / "invalid-report.json"
            invalid_path.write_text(json.dumps({"blocks": []}), encoding="utf-8")
            with redirect_stdout(StringIO()):
                self.assertEqual(
                    grid_cli_main(
                        [
                            "report",
                            str(invalid_path),
                            "--format",
                            "json",
                            "-o",
                            str(invalid_report_path),
                        ]
                    ),
                    1,
                )
            invalid_report = json.loads(invalid_report_path.read_text(encoding="utf-8"))
            self.assertFalse(invalid_report["validation"]["ok"])
            self.assertEqual(invalid_report["quality"]["summary"]["totalCells"], 0)


if __name__ == "__main__":
    unittest.main()
