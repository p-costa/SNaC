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

from http.server import ThreadingHTTPServer

from utils.snac_grid_generator import Project, export_project
from utils.snac_grid_generator.server import _Handler


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
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), _QuietHandler)
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
        page.locator("#reset-block-names").click()
        self.expect(page.locator('.block-item[data-block-id="2"] strong')).to_have_text("2: block-1")
        self.expect(page.locator('.block-item[data-block-id="1"] strong')).to_have_text("1: block-2")
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
        self.assertEqual(page.locator("#grid-kind").input_value(), "simple_ratio")
        page.locator("#library-name").fill("My uniform")
        page.locator("#save-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Saved My uniform")
        page.locator("#library-name").fill("My grading")
        page.locator("#rename-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Renamed")
        page.locator("#apply-library-item").click()
        self.expect(page.locator("#status")).to_contain_text("Applied My grading")
        self.assertEqual(page.locator("#grid-kind").input_value(), "simple_ratio")
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
