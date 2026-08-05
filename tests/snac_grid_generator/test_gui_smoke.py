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
