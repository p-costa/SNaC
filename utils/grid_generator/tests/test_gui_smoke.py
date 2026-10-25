from __future__ import annotations

from io import BytesIO
from http.client import HTTPConnection
import json
import os
from pathlib import Path
import shutil
import tempfile
from threading import Thread
import unittest


from snac_grid_generator import Project, export_project
from snac_grid_generator.server import GridServer, _Handler


def np_prod(values) -> int:
    result = 1
    for value in values:
        result *= value
    return result


class _QuietHandler(_Handler):
    def log_message(self, fmt: str, *args) -> None:
        pass


def _browser_executable() -> str | None:
    candidates = [
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
        shutil.which("chromium"),
        shutil.which("chromium-browser"),
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
    ]
    return next((str(path) for path in candidates if path and Path(path).is_file()), None)


class GridGeneratorGuiSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        try:
            from PIL import Image, ImageStat
            from playwright.sync_api import expect, sync_playwright
        except ImportError as exc:
            if os.environ.get("SNAC_REQUIRE_BROWSER_TEST"):
                raise RuntimeError("browser smoke-test dependencies are required") from exc
            raise unittest.SkipTest("browser smoke-test dependencies are not installed")

        executable = _browser_executable()
        if executable is None:
            if os.environ.get("SNAC_REQUIRE_BROWSER_TEST"):
                raise RuntimeError("Chrome or Chromium is required for the browser smoke test")
            raise unittest.SkipTest("Chrome or Chromium is not installed")

        cls.Image = Image
        cls.ImageStat = ImageStat
        cls.expect = expect
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch(
            executable_path=executable,
            headless=True,
            args=["--no-sandbox"],
        )
        cls.server = GridServer(("127.0.0.1", 0), _QuietHandler)
        cls.thread = Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f"http://127.0.0.1:{cls.server.server_port}/"

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "server"):
            cls.server.shutdown()
            cls.server.server_close()
        if hasattr(cls, "thread"):
            cls.thread.join(timeout=2.0)
        if hasattr(cls, "browser"):
            cls.browser.close()
        if hasattr(cls, "playwright"):
            cls.playwright.stop()

    def test_editor_workflow_desktop_and_mobile(self) -> None:
        context = self.browser.new_context(viewport={"width": 1440, "height": 900})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))

        with tempfile.TemporaryDirectory() as output_dir:
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            self.assertEqual(page.locator(".block-item").count(), 2)
            self._assert_canvas_is_nonblank(page)
            self.assertEqual(page.locator("#grid-precision").input_value(), "double")
            page.locator("#grid-precision").select_option("single")
            page.locator("#grid-precision").select_option("double")

            page.evaluate(
                """() => {
                  const nativeFetch = window.fetch.bind(window);
                  let delayUpdate = true;
                  window.fetch = (input, options) => {
                    const url = typeof input === 'string' ? input : input.url;
                    if (delayUpdate && url.endsWith('/api/update')) {
                      delayUpdate = false;
                      return new Promise((resolve, reject) => {
                        setTimeout(() => nativeFetch(input, options).then(resolve, reject), 300);
                      });
                    }
                    return nativeFetch(input, options);
                  };
                }"""
            )
            page.locator("#update-structure").click()
            page.locator("#project-name").fill("edited-during-update")
            self.expect(page.locator("#status")).to_contain_text(
                "discarded because the project changed"
            )
            self.assertEqual(page.locator("#project-name").input_value(), "edited-during-update")

            self.assertEqual(page.locator("[data-axis-view]").count(), 6)
            for view in ("+x", "-x", "+y", "-y", "+z", "-z"):
                button = page.locator(f'[data-axis-view="{view}"]')
                button.click()
                self.expect(page.locator("#scene")).to_have_attribute("data-projection", "orthographic")
                self.expect(page.locator("#scene")).to_have_attribute("data-view", view)
                self.expect(button).to_have_attribute("aria-pressed", "true")
            axis_button = page.locator('[data-axis-view="+x"]')
            axis_button.click()
            canvas = page.locator("#scene").bounding_box()
            self.assertIsNotNone(canvas)
            page.mouse.move(canvas["x"] + canvas["width"] * 0.5, canvas["y"] + canvas["height"] * 0.5)
            page.mouse.down()
            page.mouse.move(
                canvas["x"] + canvas["width"] * 0.65,
                canvas["y"] + canvas["height"] * 0.6,
                steps=6,
            )
            page.mouse.up()
            self.expect(page.locator("#scene")).to_have_attribute("data-projection", "orthographic")
            self.expect(page.locator("#scene")).to_have_attribute("data-view", "free-orthographic")
            self.expect(page.locator("#scene")).to_have_attribute("data-selected-blocks", "1")
            self.expect(axis_button).to_have_attribute("aria-pressed", "false")
            page.locator("#focus-selection").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-view", "free-orthographic")
            page.locator("#fit-view").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-projection", "perspective")
            self.expect(page.locator("#scene")).to_have_attribute("data-view", "3d")

            page.locator("#focus-selection").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-framed-blocks", "1")
            page.locator("#hide-selected").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-hidden-blocks", "1")
            self.assertIn("hidden-block", page.locator(".block-item").first.get_attribute("class"))
            self.assertTrue(page.locator("#show-all-blocks").is_enabled())
            page.locator("#show-all-blocks").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-hidden-blocks", "")
            page.locator("#isolate-selected").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-hidden-blocks", "2")
            page.locator("#show-all-blocks").click()

            page.locator("#toggle-clipping").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-clipping", "true")
            self.assertTrue(page.locator("#clip-tools").is_visible())
            page.locator('[data-clip-axis="y"]').click()
            self.expect(page.locator("#scene")).to_have_attribute("data-clip-axis", "y")
            page.locator("#clip-position").fill("0.25")
            page.locator("#flip-clipping").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-clip-direction", "-1")
            page.locator("#toggle-clipping").click()
            self.expect(page.locator("#scene")).to_have_attribute("data-clipping", "false")
            self.assertFalse(page.locator("#clip-tools").is_visible())

            page.get_by_label("Select 2: right", exact=True).check()
            page.get_by_role("button", name="Duplicate array", exact=True).click()
            self.expect(page.locator(".block-item")).to_have_count(4)

            page.locator("#check-project").click()
            self.expect(page.locator("#status")).to_contain_text("passed")
            page.locator("#bc-preset").select_option("moving_wall")
            page.locator("#apply-bc-preset").click()
            self.expect(page.locator("#status")).to_contain_text("Applied preset")
            page.locator("#check-project").click()
            self.expect(page.locator("#status")).to_contain_text("passed")

            page.locator("#target-ranks").fill("8")
            self.assertTrue(page.locator("#balance-decomposition").is_enabled())
            page.locator("#balance-decomposition").click()
            self.expect(page.locator("#status")).to_contain_text("Balanced 8")

            page.reload(wait_until="load")
            page.locator("#recovery-dialog").wait_for(state="visible")
            page.locator("#restore-recovery").click()
            self.expect(page.locator(".block-item")).to_have_count(4)

            page.locator("#output-dir").fill(output_dir)
            page.locator("#export-case").click()
            self.expect(page.locator("#status")).to_contain_text("Wrote")
            self.assertTrue((Path(output_dir) / "blocks.nml").is_file())

            page.set_viewport_size({"width": 390, "height": 844})
            self.assertLessEqual(
                page.evaluate("document.documentElement.scrollWidth"),
                page.evaluate("window.innerWidth") + 1,
            )
            self.assertTrue(page.locator("#duplicate-array").is_visible())
            self._assert_canvas_is_nonblank(page)

            page.reload(wait_until="load")
            self.assertFalse(page.locator("#recovery-dialog").is_visible())

        context.close()
        self.assertEqual(errors, [])

    def test_block_creation_preserves_explicit_grid_coordinates(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))

        with tempfile.TemporaryDirectory() as import_dir, tempfile.TemporaryDirectory() as output_dir:
            project = Project.from_dict(
                {
                    "name": "explicit-clone",
                    "writeExternalGrid": True,
                    "blocks": [
                        {
                            "id": 1,
                            "ng": [2, 4, 2],
                            "axes": {"x": {"kind": "explicit", "faces": [0.0, 0.25, 1.0]}},
                        }
                    ],
                }
            )
            export_project(project, import_dir)

            page.goto(self.url, wait_until="load")
            page.locator("#output-dir").fill(import_dir)
            page.locator("#import-case").click()
            self.expect(page.locator("#status")).to_contain_text("Imported")
            self.assertEqual(page.locator("#inivel option").all_text_contents()[-2:], ["kov", "pdc"])

            page.locator("#add-block").click()
            self.expect(page.locator(".block-item")).to_have_count(2)
            page.locator("#adjacent-size").fill("2")
            page.locator("#add-adjacent").click()
            self.expect(page.locator(".block-item")).to_have_count(3)

            page.locator("#output-dir").fill(output_dir)
            page.locator("#export-case").click()
            self.expect(page.locator("#status")).to_contain_text("Wrote")
            saved = json.loads((Path(output_dir) / "snac_grid_project.json").read_text(encoding="utf-8"))

        np_faces = [block["axes"]["x"]["faces"] for block in saved["blocks"]]
        self.assertEqual(np_faces[0], [0.0, 0.25, 1.0])
        self.assertEqual(np_faces[1], [1.0, 1.25, 2.0])
        self.assertEqual(np_faces[2], [2.0, 2.5, 4.0])
        context.close()
        self.assertEqual(errors, [])

    def test_multiblock_selection_grid_scope_tools_and_ordering(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(self.url, wait_until="load")
        page.wait_for_selector(".block-item")

        friend_remap = page.evaluate(
            """async () => {
              const {defaultBlock, renumberBlocksByOrder} = await import('/block-model.js');
              const lower = defaultBlock(4);
              const upper = defaultBlock(9);
              for (const [block, face, friend] of [[lower, 1, 9], [upper, 0, 4]]) {
                block.cbcvel[0][face] = 'F';
                block.bcvel[0][face] = friend;
                block.cbcpre[face] = 'F';
                block.bcpre[face] = friend;
                block.cbcscal = [['N', 'N', 'N', 'N', 'N', 'N']];
                block.bcscal = [[0, 0, 0, 0, 0, 0]];
                block.cbcscal[0][face] = 'F';
                block.bcscal[0][face] = friend;
              }
              const blocks = [upper, lower];
              renumberBlocksByOrder(blocks);
              return blocks.map((block, index) => ({
                id: block.id,
                name: block.name,
                refs: [block.bcvel[0][index], block.bcpre[index], block.bcscal[0][index]],
              }));
            }"""
        )
        self.assertEqual(friend_remap[0], {"id": 1, "name": "block-1", "refs": [2, 2, 2]})
        self.assertEqual(friend_remap[1], {"id": 2, "name": "block-2", "refs": [1, 1, 1]})

        scene = page.locator("#scene")
        first = page.locator('.block-item[data-block-id="1"]')
        second = page.locator('.block-item[data-block-id="2"]')
        self.expect(scene).to_have_attribute("data-selected-blocks", "1")
        self.expect(scene).to_have_attribute("data-primary-block", "1")

        second.locator(".block-main").click(modifiers=["Shift"])
        self.expect(scene).to_have_attribute("data-selected-blocks", "1,2")
        self.expect(scene).to_have_attribute("data-primary-block", "2")
        self.expect(scene).to_have_attribute("data-grid-blocks", "1,2")
        self.expect(scene).to_have_attribute("data-boundary-blocks", "1,2")
        self.expect(scene).to_have_attribute("data-partition-blocks", "1,2")
        self.assertTrue(first.locator('input[type="checkbox"]').is_checked())
        self.assertTrue(second.locator('input[type="checkbox"]').is_checked())

        second.locator(".block-main").click(modifiers=["Shift"])
        self.expect(scene).to_have_attribute("data-selected-blocks", "1")
        self.expect(scene).to_have_attribute("data-primary-block", "1")
        second.locator(".block-main").click()
        self.expect(scene).to_have_attribute("data-selected-blocks", "2")
        self.expect(scene).to_have_attribute("data-primary-block", "2")

        page.locator("#grid-lines-mode").select_option("all")
        self.expect(scene).to_have_attribute("data-grid-mode", "all")
        self.expect(scene).to_have_attribute("data-grid-blocks", "1,2")
        page.locator("#grid-lines-mode").select_option("off")
        self.expect(scene).to_have_attribute("data-grid-blocks", "")
        page.locator("#grid-lines-mode").select_option("selected")
        self.expect(scene).to_have_attribute("data-grid-blocks", "2")

        move = page.locator("#tool-move")
        move.click()
        self.expect(move).to_have_attribute("aria-pressed", "true")
        move.click()
        self.expect(move).to_have_attribute("aria-pressed", "false")
        self.expect(page.locator("#tool-select")).to_have_attribute("aria-pressed", "true")
        scale = page.locator("#tool-scale")
        scale.click()
        self.expect(scale).to_have_attribute("aria-pressed", "true")
        scale.click()
        self.expect(scale).to_have_attribute("aria-pressed", "false")

        second_box = second.bounding_box()
        self.assertIsNotNone(second_box)
        first.locator(".block-drag-handle").drag_to(
            second,
            target_position={"x": second_box["width"] / 2, "y": second_box["height"] - 2},
        )
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "2")
        self.expect(page.locator(".block-item").nth(1)).to_have_attribute("data-block-id", "1")
        page.wait_for_timeout(300)
        self.assertTrue(page.locator("#undo-project").is_enabled())
        page.locator(".block-drag-handle").nth(1).press("ArrowUp")
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "1")
        page.wait_for_timeout(300)
        page.locator("#undo-project").click()
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "2")
        page.locator("#redo-project").click()
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "1")
        page.locator('.block-item[data-block-id="1"] .block-drag-handle').press("ArrowDown")
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "2")

        page.locator('.block-item[data-block-id="2"] .block-main').click()
        page.locator("#block-name").fill("custom-right")
        page.locator("#renumber-blocks").click()
        self.expect(page.locator(".block-item").first).to_have_attribute("data-block-id", "1")
        self.expect(page.locator(".block-item").first.locator("strong")).to_have_text("1: block-1")
        self.expect(page.locator(".block-item").nth(1)).to_have_attribute("data-block-id", "2")
        self.expect(page.locator(".block-item").nth(1).locator("strong")).to_have_text("2: block-2")
        self.expect(scene).to_have_attribute("data-primary-block", "1")
        self.assertEqual(page.locator('[data-field="lmin"][data-index="0"]').input_value(), "1")
        with page.expect_download() as download_info:
            page.locator("#download-json").click()
        renumbered = json.loads(Path(download_info.value.path()).read_text(encoding="utf-8"))
        self.assertEqual([block["id"] for block in renumbered["blocks"]], [1, 2])
        self.assertEqual([block["name"] for block in renumbered["blocks"]], ["block-1", "block-2"])
        self.assertEqual([block["lmin"][0] for block in renumbered["blocks"]], [1.0, 0.0])
        self.assertEqual(renumbered["blocks"][0]["bcpre"][0], 2)
        self.assertEqual(renumbered["blocks"][1]["bcpre"][1], 1)
        page.wait_for_timeout(300)
        page.locator("#add-block").click()
        self.expect(page.locator(".block-item")).to_have_count(3)
        page.wait_for_timeout(300)
        page.locator("#undo-project").click()
        self.expect(page.locator(".block-item")).to_have_count(2)
        page.reload(wait_until="load")
        self.expect(page.locator("#recovery-dialog")).to_be_visible()
        page.locator("#restore-recovery").click()
        self.expect(page.locator(".block-item")).to_have_count(2)
        self._assert_canvas_is_nonblank(page)

        context.close()
        self.assertEqual(errors, [])

    def test_local_api_rejects_cross_origin_and_unmarked_requests(self) -> None:
        payload = json.dumps({"project": {"blocks": []}})
        headers = {
            "Content-Type": "application/json",
            "X-SNaC-Grid-Request": "1",
            "Origin": self.url.rstrip("/"),
        }
        status, body = self._post_api(payload, headers)
        self.assertEqual(status, 200)
        self.assertFalse(json.loads(body)["ok"])

        status, _ = self._post_api(payload, {**headers, "Origin": "https://example.com"})
        self.assertEqual(status, 403)

        missing_marker = {key: value for key, value in headers.items() if key != "X-SNaC-Grid-Request"}
        status, _ = self._post_api(payload, missing_marker)
        self.assertEqual(status, 403)

        status, body = self._post_api(payload, {**headers, "Content-Type": "text/plain"})
        self.assertEqual(status, 415)
        self.assertIn("application/json", json.loads(body)["error"])

        evil_host = {**headers, "Host": f"example.com:{self.server.server_port}"}
        status, _ = self._post_api(payload, evil_host)
        self.assertEqual(status, 403)

        oversized = {**headers, "Content-Length": str(16 * 1024 * 1024 + 1)}
        status, body = self._post_api("{}", oversized)
        self.assertEqual(status, 413)
        self.assertIn("too large", json.loads(body)["error"])

    def test_validation_diagnostic_selects_the_affected_interface(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        page = context.new_page()
        page.goto(self.url, wait_until="load")
        page.locator(".block-item .block-main").nth(1).click()
        page.locator('[data-field="ng"][data-index="1"]').fill("33")
        page.locator("#check-project").click()

        diagnostic = page.locator(".check-line.actionable").filter(
            has_text="different ng along y"
        )
        self.expect(diagnostic).to_have_count(1)
        diagnostic.click()

        self.assertTrue(page.get_by_label("Select 1: left", exact=True).is_checked())
        self.assertTrue(page.get_by_label("Select 2: right", exact=True).is_checked())
        self.expect(page.locator("#status")).to_contain_text("interface cell count")
        context.close()

    def test_quality_summary_and_report_downloads(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(self.url, wait_until="load")
        page.locator("#check-project").click()

        self.expect(page.locator("#quality-summary")).to_contain_text("Cells")
        self.expect(page.locator("#quality-summary")).to_contain_text("Cell aspect")
        self.assertTrue(page.locator("#download-report-json").is_enabled())

        with page.expect_download() as json_download_info:
            page.locator("#download-report-json").click()
        json_download = json_download_info.value
        report = json.loads(Path(json_download.path()).read_text(encoding="utf-8"))
        self.assertTrue(report["validation"]["ok"])
        self.assertGreater(report["quality"]["summary"]["totalCells"], 0)
        self.assertTrue(json_download.suggested_filename.endswith("-quality.json"))

        with page.expect_download() as markdown_download_info:
            page.locator("#download-report-markdown").click()
        markdown_download = markdown_download_info.value
        markdown = Path(markdown_download.path()).read_text(encoding="utf-8")
        self.assertIn("## Axis Quality", markdown)
        self.assertIn("## SNaC Storage", markdown)
        self.assertTrue(markdown_download.suggested_filename.endswith("-quality.md"))

        context.close()
        self.assertEqual(errors, [])

    def test_grid_library_presets_templates_and_import_export(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(self.url, wait_until="load")

        explicit = page.evaluate(
            """async () => {
              const library = await import('/library.js');
              const preset = library.captureAxisPreset(
                'Imported faces',
                {kind: 'explicit', faces: [2, 2.5, 4]},
                2,
                2,
                4,
                'axis-test'
              );
              return {
                stored: preset.axis.faces,
                applied: library.applyAxisPreset(preset, 10, 14).axis.faces,
              };
            }"""
        )
        self.assertEqual(explicit["stored"], [0, 0.25, 1])
        self.assertEqual(explicit["applied"], [10, 11, 14])
        schema_error = page.evaluate(
            """async () => {
              const library = await import('/library.js');
              try {
                library.normalizeLibrary({schemaVersion: 99});
                return '';
              } catch (error) {
                return error.message;
              }
            }"""
        )
        self.assertIn("Unsupported library schema version", schema_error)

        page.locator("#open-library").click()
        self.expect(page.locator("#library-dialog")).to_be_visible()
        page.locator("#library-items").select_option("builtin:builtin-uniform-32")
        page.locator("#apply-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Applied Uniform 32")
        self.assertEqual(page.locator("#grid-kind").input_value(), "monotone")
        page.locator("#library-name").fill("My uniform")
        page.locator("#save-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Saved My uniform")
        page.locator("#library-name").fill("My grading")
        page.locator("#rename-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Renamed")
        page.locator("#apply-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Applied My grading")
        self.assertEqual(page.locator("#grid-kind").input_value(), "monotone")
        self.assertEqual(page.locator('[data-field="ng"][data-index="0"]').input_value(), "32")

        page.get_by_role("button", name="Case templates", exact=True).click()
        page.locator("#library-items").select_option("builtin:builtin-periodic-channel")
        page.locator("#apply-library-item").click()
        self.expect(page.locator("#library-dialog")).not_to_be_visible()
        self.expect(page.locator("#status")).to_contain_text("Applied Periodic channel")
        self.expect(page.locator(".block-item")).to_have_count(1)
        self.assertTrue(page.locator('#periodic-grid input[data-axis-index="0"]').is_checked())
        self.assertEqual(page.locator('#boundary-table select[data-face="0"]').first.input_value(), "F")
        self.assertEqual(page.locator('#boundary-table input[data-face="0"]').first.input_value(), "1")
        self.assertEqual(page.locator('#boundary-table select[data-face="1"]').first.input_value(), "F")

        page.locator("#open-library").click()
        page.locator("#library-name").fill("Periodic working case")
        page.locator("#save-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Saved Periodic working case")
        with page.expect_download() as download_info:
            page.locator("#export-library").click()
        download = download_info.value
        exported = json.loads(Path(download.path()).read_text(encoding="utf-8"))
        self.assertEqual(exported["schemaVersion"], 1)
        self.assertEqual([item["name"] for item in exported["axisPresets"]], ["My grading"])
        self.assertEqual([item["name"] for item in exported["projectTemplates"]], ["Periodic working case"])

        page.locator("#delete-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Deleted Periodic working case")
        page.locator("#import-library").set_input_files(download.path())
        self.expect(page.locator("#status")).to_contain_text("Imported 2 library items")
        self.expect(page.locator("#library-items")).to_contain_text("Periodic working case")

        page.reload(wait_until="load")
        page.locator("#recovery-dialog").wait_for(state="visible")
        page.locator("#restore-recovery").click()
        page.locator("#open-library").click()
        page.get_by_role("button", name="Axis presets", exact=True).click()
        self.expect(page.locator("#library-items")).to_contain_text("My grading")
        page.get_by_role("button", name="Case templates", exact=True).click()
        self.expect(page.locator("#library-items")).to_contain_text("Periodic working case")

        context.close()
        self.assertEqual(errors, [])

    def test_check_is_read_only_and_editor_inputs_are_guarded(self) -> None:
        context = self.browser.new_context(viewport={"width": 1440, "height": 900})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))

        with tempfile.TemporaryDirectory() as tmp:
            # Two touching blocks whose shared faces are not yet marked as friends.
            project_path = Path(tmp) / "unconnected.json"
            project_path.write_text(
                json.dumps(
                    {
                        "name": "unconnected",
                        "blocks": [
                            {"id": 1, "ng": [8, 8, 2], "lmin": [0, 0, 0], "lmax": [1, 1, 0.1]},
                            {"id": 2, "ng": [8, 8, 2], "lmin": [1, 0, 0], "lmax": [2, 1, 0.1]},
                        ],
                    }
                ),
                encoding="utf-8",
            )
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            page.locator("#open-json").set_input_files(str(project_path))
            self.expect(page.locator("#status")).to_contain_text("Loaded")

            page.locator("#check-project").click()
            self.expect(page.locator("#status")).to_contain_text("Update applies the inferred friend boundaries")
            self.assertEqual(self._download_project(page)["blocks"][0]["cbcpre"][1], "N")

            page.locator("#update-structure").click()
            self.expect(page.locator("#status")).to_contain_text("Structure updated")
            self.assertEqual(self._download_project(page)["blocks"][0]["cbcpre"][1], "F")

            # An empty or partial number keeps the last valid value and is restored on blur.
            lmax_x = page.locator('#geometry-grid input[data-field="lmax"][data-index="0"]')
            lmax_x.fill("")
            self.assertEqual(self._download_project(page)["blocks"][0]["lmax"][0], 1.0)
            lmax_x.blur()
            self.assertEqual(lmax_x.input_value(), "1")

            # The FFT axis is never partitioned.
            page.locator("#fft-axis").select_option("z")
            z_partition = page.locator('#decomposition-axes input[data-axis="2"]')
            self.expect(z_partition).to_be_disabled()
            self.expect(z_partition).not_to_be_checked()
            saved = self._download_project(page)
            self.assertEqual(saved["fftAxis"], "z")
            self.assertEqual(saved["decomposition"]["axes"], [True, True, False])

            # Geometric native functions start from their uniform gr instead of gr = 0.
            page.locator("#grid-kind").select_option("native")
            page.locator("#snac-gt").select_option("6")
            self.assertEqual(page.locator("#snac-gr").input_value(), "2")
            # ...and returning to a tanh function drops that automatic value again.
            page.locator("#snac-gt").select_option("0")
            self.assertEqual(page.locator("#snac-gr").input_value(), "0")

            # The FFT axis never counts as an enabled partition axis.
            page.locator('#decomposition-axes input[data-axis="0"]').uncheck()
            page.locator('#decomposition-axes input[data-axis="1"]').click()
            self.expect(page.locator("#status")).to_contain_text("Keep at least one partition axis enabled")
            self.assertEqual(self._download_project(page)["decomposition"]["axes"], [False, True, False])

            # Clearing other number fields keeps their last valid values.
            page.locator("#grid-kind").select_option("monotone")
            cells = page.locator('#axis-editor input[data-control="n"]')
            self.expect(cells).to_have_value("8")
            cells.fill("")
            pressure = page.locator('#boundary-table input[data-kind="pre"][data-face="2"]')
            pressure.fill("2.5")
            pressure.fill("")
            aspect = page.locator("#max-local-aspect")
            aspect.fill("0.5")
            saved = self._download_project(page)
            self.assertEqual(saved["blocks"][0]["ng"][0], 8)
            self.assertEqual(saved["blocks"][0]["bcpre"][2], 2.5)
            self.assertEqual(saved["decomposition"]["maxLocalAspect"], 0)
            cells.blur()
            self.assertEqual(cells.input_value(), "8")

        context.close()
        self.assertEqual(errors, [])

    def test_legacy_library_axis_preset_applies_with_canonical_names(self) -> None:
        context = self.browser.new_context(viewport={"width": 1440, "height": 900})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))

        with tempfile.TemporaryDirectory() as tmp:
            # A library saved before schema 4, with legacy axis kind and keys.
            library_path = Path(tmp) / "legacy-library.json"
            library_path.write_text(
                json.dumps(
                    {
                        "schemaVersion": 1,
                        "axisPresets": [
                            {
                                "id": "legacy-axis",
                                "name": "Legacy monotone",
                                "ng": 12,
                                "coordinateSpace": "native",
                                "axis": {
                                    "kind": "simple_ratio",
                                    "profile": "geometric",
                                    "controls": ["n", "width_start"],
                                    "width_start": 0.02,
                                },
                            }
                        ],
                        "projectTemplates": [],
                    }
                ),
                encoding="utf-8",
            )
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            page.locator("#open-library").click()
            page.locator("#import-library").set_input_files(str(library_path))
            self.expect(page.locator("#status")).to_contain_text("Imported 1 library item")
            page.locator("#library-items").select_option("user:legacy-axis")
            page.locator("#apply-library-item").click()
            self.expect(page.locator("#status")).to_contain_text("Applied Legacy monotone")
            page.locator("#close-library").click()
            self.expect(page.locator("#grid-kind")).to_have_value("monotone")
            axis = self._download_project(page)["blocks"][0]["axes"]["x"]

        self.assertEqual(axis["kind"], "monotone")
        self.assertEqual(axis["controls"], ["n", "widthStart"])
        self.assertEqual(axis["widthStart"], 0.02)
        self.assertNotIn("width_start", axis)
        context.close()
        self.assertEqual(errors, [])

    def test_each_tab_keeps_its_own_autosave(self) -> None:
        context = self.browser.new_context(viewport={"width": 1280, "height": 800})
        errors: list[str] = []

        def open_tab():
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            return page

        alpha = open_tab()
        alpha.locator("#project-name").fill("alpha-tab")
        alpha.wait_for_timeout(700)  # autosave writes are batched

        # Work in a tab that is still open is not offered to another tab.
        beta = open_tab()
        self.assertFalse(beta.locator("#recovery-dialog").is_visible())
        beta.locator("#project-name").fill("beta-tab")
        beta.wait_for_timeout(700)

        # Reloading a tab offers that tab's own work, not the other tab's.
        alpha.reload(wait_until="load")
        self.expect(alpha.locator("#recovery-dialog")).to_be_visible()
        self.expect(alpha.locator("#recovery-summary")).to_contain_text("alpha-tab")
        alpha.locator("#restore-recovery").click()
        self.expect(alpha.locator("#project-name")).to_have_value("alpha-tab")

        # Work left by a closed tab is offered to the next new tab.
        beta.close()
        fresh = open_tab()
        self.expect(fresh.locator("#recovery-dialog")).to_be_visible()
        self.expect(fresh.locator("#recovery-summary")).to_contain_text("beta-tab")

        context.close()
        self.assertEqual(errors, [])

    def test_issues_repair_choice_balancing_and_new_project(self) -> None:
        context = self.browser.new_context(viewport={"width": 1440, "height": 900})
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))

        with tempfile.TemporaryDirectory() as tmp:
            # A backward-facing step whose upper block was stretched along x, so
            # its x grid no longer matches the lower block below it.
            path = Path(tmp) / "step.json"
            path.write_text(
                json.dumps(
                    {
                        "name": "step",
                        "fftAxis": "z",
                        "periodicAxes": [False, False, True],
                        "blocks": [
                            {"id": 1, "name": "inlet", "ng": [32, 24, 32], "lmin": [0, 1, 0], "lmax": [2, 2, 2]},
                            {
                                "id": 2,
                                "name": "upper",
                                "ng": [96, 24, 32],
                                "lmin": [2, 1, 0],
                                "lmax": [10, 2, 2],
                                "axes": {"x": {"kind": "monotone", "controls": ["n", "ratio"], "ratio": 4}},
                            },
                            {"id": 3, "name": "lower", "ng": [96, 24, 32], "lmin": [2, 0, 0], "lmax": [10, 1, 2]},
                        ],
                    }
                ),
                encoding="utf-8",
            )
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            page.locator("#open-json").set_input_files(str(path))
            self.expect(page.locator("#status")).to_contain_text("Loaded")

        # A failing check opens the issues panel in the viewport.
        page.locator("#check-project").click()
        panel = page.locator("#issues-panel")
        self.expect(panel).to_be_visible()
        self.expect(page.locator("#issues-title")).to_have_text("1 error")
        self.expect(page.locator("#export-case")).to_be_disabled()
        self.assertIn("Fix the check errors", page.locator("#export-case").get_attribute("title"))
        page.locator("#toggle-issues").click()
        self.expect(page.locator("#check-results")).to_be_hidden()
        page.locator("#toggle-issues").click()

        # Clicking the issue frames the two affected blocks.
        page.locator("#check-results .check-line.actionable").first.click()
        self.expect(page.locator("#scene")).to_have_attribute("data-framed-blocks", "2,3")

        # Repair offers both grids and keeps the one picked, here the stretched one.
        page.locator("#repair-grids").click()
        self.expect(page.locator("#repair-dialog")).to_be_visible()
        options = page.locator('#repair-summary input[type="radio"]')
        self.expect(options).to_have_count(2)
        self.expect(page.locator("#repair-summary")).to_contain_text("monotone geometric, upper/lower ratio 4")
        page.locator('#repair-summary label', has_text="2: upper").locator("input").check()
        self.expect(page.locator("#repair-summary .repair-change")).to_contain_text("3: lower")
        self.expect(page.locator("#repair-summary .repair-change")).to_contain_text("→ monotone geometric")
        page.locator("#confirm-repair").click()
        self.expect(page.locator("#status")).to_contain_text("Applied 1 grid repair")
        self.assertEqual(self._download_project(page)["blocks"][2]["axes"]["x"]["kind"], "monotone")

        # Six ranks balance poorly here, so better totals are offered.
        page.locator("#target-ranks").fill("6")
        page.locator("#balance-decomposition").click()
        self.expect(page.locator("#decomposition-summary")).to_contain_text("Better balanced totals")
        page.locator("#decomposition-summary .nearby-rank", has_text="7 ·").click()
        self.expect(page.locator("#status")).to_contain_text("Balanced 7 MPI ranks")
        self.expect(page.locator("#decomposition-summary")).to_contain_text("imbalance 1.0000")

        # A new empty project resets the name; Undo brings the old project back.
        page.locator("#reset-project").click()
        self.expect(page.locator(".block-item")).to_have_count(0)
        self.expect(page.locator("#project-name")).to_have_value("snac-grid")
        page.locator("#undo-project").click()
        self.expect(page.locator(".block-item")).to_have_count(3)
        self.expect(page.locator("#project-name")).to_have_value("step")

        # Adding a block that ends up out of view frames the model again.
        page.locator("#adjacent-face").select_option("x+")
        page.locator("#adjacent-size").fill("40")
        page.locator("#add-adjacent").click()
        self.expect(page.locator(".block-item")).to_have_count(4)
        self.assertEqual(len(page.locator("#scene").get_attribute("data-framed-blocks").split(",")), 4)

        context.close()
        self.assertEqual(errors, [])

    def test_balance_offers_and_applies_cell_count_changes(self) -> None:
        context = self.browser.new_context(viewport={"width": 1440, "height": 900}, accept_downloads=True)
        page = context.new_page()
        errors: list[str] = []
        page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: errors.append(str(error)))
        blocks = [
            {"id": 1, "name": "inlet", "ng": [32, 24, 32], "lmin": [0, 1, 0], "lmax": [2, 2, 2]},
            {"id": 2, "name": "upper", "ng": [96, 24, 32], "lmin": [2, 1, 0], "lmax": [10, 2, 2]},
            {"id": 3, "name": "lower", "ng": [96, 24, 32], "lmin": [2, 0, 0], "lmax": [10, 1, 2]},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "step.json"
            path.write_text(
                json.dumps({"name": "step", "fftAxis": "z", "periodicAxes": [False, False, True], "blocks": blocks}),
                encoding="utf-8",
            )
            page.goto(self.url, wait_until="load")
            page.wait_for_selector(".block-item")
            page.locator("#open-json").set_input_files(str(path))
            self.expect(page.locator("#status")).to_contain_text("Loaded")

        self.expect(page.locator("#ng-tolerance")).to_have_value("10")
        page.locator("#target-ranks").fill("9")
        page.locator("#balance-decomposition").click()
        summary = page.locator("#decomposition-summary")
        self.expect(summary).to_contain_text("Adjust cell counts (within 10%)")
        self.expect(page.locator(".balance-alternative")).to_have_count(4)
        self.expect(page.locator(".balance-alternative").nth(2)).to_contain_text("≤8.33% · imbalance 1.0707 · exact")
        self.expect(page.locator(".balance-alternative").nth(2)).to_contain_text("upper, lower: 96 → 104")
        page.locator(".balance-fixed summary").click()
        self.expect(page.locator(".balance-fixed")).to_contain_text("z · all blocks · 32 · FFT axis")

        page.locator(".apply-alternative").nth(2).click()
        self.expect(page.locator("#status")).to_contain_text("Changed 4 cell counts; imbalance 1.0707")
        self.expect(summary).to_contain_text("Cell counts changed")
        saved = self._download_project(page)
        self.assertEqual(sum(np_prod(block["dims"]) for block in saved["blocks"]), 9)
        for before, after in zip(blocks, saved["blocks"]):
            for old, new in zip(before["ng"], after["ng"]):
                self.assertLessEqual(abs(new - old), 0.1 * old)
        self.assertEqual([block["ng"][0] for block in saved["blocks"]], [30, 104, 104])

        # Later proposals measure the tolerance from the original counts.
        page.locator("#target-ranks").fill("6")
        page.locator("#balance-decomposition").click()
        self.expect(summary).to_contain_text("Adjust cell counts (within 10% of the original counts)")
        self.expect(page.locator(".balance-alternative").first).to_contain_text("Original counts · imbalance 1.2857")
        self.expect(summary).to_contain_text("upper, lower: 104 → 87 (original 96)")
        page.locator("#target-ranks").fill("9")

        # A zero tolerance turns the proposals off.
        page.locator("#ng-tolerance").fill("0")
        page.locator("#balance-decomposition").click()
        self.expect(page.locator("#status")).to_contain_text("Balanced 9 MPI ranks")
        self.expect(page.locator(".balance-alternative")).to_have_count(0)
        self.assertEqual(self._download_project(page)["decomposition"]["ngTolerance"], 0)
        self.assertEqual(errors, [])
        context.close()

    def _download_project(self, page) -> dict:
        with page.expect_download() as download_info:
            page.locator("#download-json").click()
        return json.loads(Path(download_info.value.path()).read_text(encoding="utf-8"))

    def _post_api(self, payload: str, headers: dict[str, str]) -> tuple[int, str]:
        connection = HTTPConnection("127.0.0.1", self.server.server_port, timeout=5)
        try:
            connection.request("POST", "/api/check", body=payload, headers=headers)
            response = connection.getresponse()
            return response.status, response.read().decode("utf-8")
        finally:
            connection.close()

    def _assert_canvas_is_nonblank(self, page) -> None:
        image = self.Image.open(BytesIO(page.locator("#scene").screenshot())).convert("RGB")
        variation = max(self.ImageStat.Stat(image).var)
        self.assertGreater(variation, 40.0)


if __name__ == "__main__":
    unittest.main()
