"""Small local web server for the SNaC grid GUI."""

from __future__ import annotations

import argparse
import json
import mimetypes
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

from .export import export_project
from .model import Project
from .validation import check_project, update_project_structure

ROOT = Path(__file__).resolve().parent
STATIC = ROOT / "static"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the SNaC multi-block grid GUI.")
    parser.add_argument("--host", default="127.0.0.1", help="Bind address")
    parser.add_argument("--port", type=int, default=8765, help="Bind port")
    parser.add_argument("--no-open", action="store_true", help="Do not open a browser tab")
    args = parser.parse_args(argv)

    server = ThreadingHTTPServer((args.host, args.port), _Handler)
    url = f"http://{args.host}:{server.server_port}/"
    print(f"SNaC grid generator listening on {url}")
    if not args.no_open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping SNaC grid generator.")
    return 0


class _Handler(BaseHTTPRequestHandler):
    server_version = "SNaCGridGUI/0.1"

    def do_GET(self) -> None:
        path = urlparse(self.path).path
        if path == "/api/sample":
            self._send_json(_sample_project().to_dict())
            return
        self._serve_static(path)

    def do_POST(self) -> None:
        path = urlparse(self.path).path
        if path not in {"/api/check", "/api/export", "/api/update"}:
            self.send_error(HTTPStatus.NOT_FOUND, "unknown API route")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            project = Project.from_dict(payload.get("project", payload))
            if path == "/api/check":
                result = check_project(project)
                self._send_json(result.to_dict())
                return
            if path == "/api/update":
                source_block_id = payload.get("sourceBlockId")
                project, result = update_project_structure(project, int(source_block_id) if source_block_id else None)
                self._send_json({**result.to_dict(), "project": project.to_dict()})
                return
            output_dir = payload.get("outputDir") or "generated/snac_grid_case"
            result = export_project(project, output_dir)
        except Exception as exc:
            self._send_json({"ok": False, "error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
            return
        self._send_json({"ok": True, **result.to_dict()})

    def log_message(self, fmt: str, *args) -> None:
        print(f"{self.address_string()} - {fmt % args}")

    def _serve_static(self, raw_path: str) -> None:
        path = "/index.html" if raw_path in {"", "/"} else raw_path
        relative = Path(unquote(path).lstrip("/"))
        candidate = (STATIC / relative).resolve()
        if STATIC not in candidate.parents and candidate != STATIC:
            self.send_error(HTTPStatus.FORBIDDEN, "invalid path")
            return
        if not candidate.exists() or not candidate.is_file():
            self.send_error(HTTPStatus.NOT_FOUND, "file not found")
            return
        content_type = mimetypes.guess_type(candidate.name)[0] or "application/octet-stream"
        body = candidate.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, data: dict, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = (json.dumps(data, indent=2) + "\n").encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _sample_project() -> Project:
    project = Project.from_dict(
        {
            "name": "two-block-demo",
            "inferConnectivity": True,
            "writeExternalGrid": True,
            "externalGridSource": "grid",
            "blocks": [
                {
                    "id": 1,
                    "name": "left",
                    "dims": [1, 1, 1],
                    "ng": [32, 32, 4],
                    "lmin": [0.0, 0.0, 0.0],
                    "lmax": [1.0, 1.0, 0.1],
                    "axes": {
                        "x": {"kind": "snac", "gt": 0, "gr": 0.0},
                        "y": {
                            "kind": "multi",
                            "segments": [
                                {"length": 20, "cells": 30, "ratio": 4},
                                {"length": 60, "cells": 40, "ratio": 1},
                                {"length": 20, "cells": 30, "ratio": 0.25},
                            ],
                        },
                        "z": {"kind": "snac", "gt": 0, "gr": 0.0},
                    },
                    "inivel": "zer",
                },
                {
                    "id": 2,
                    "name": "right",
                    "dims": [1, 1, 1],
                    "ng": [32, 32, 4],
                    "lmin": [1.0, 0.0, 0.0],
                    "lmax": [2.0, 1.0, 0.1],
                    "axes": {
                        "x": {"kind": "simple_ratio", "ratio": 1.5},
                        "y": {
                            "kind": "multi",
                            "segments": [
                                {"length": 20, "cells": 30, "ratio": 4},
                                {"length": 60, "cells": 40, "ratio": 1},
                                {"length": 20, "cells": 30, "ratio": 0.25},
                            ],
                        },
                        "z": {"kind": "snac", "gt": 0, "gr": 0.0},
                    },
                    "inivel": "zer",
                },
            ],
        }
    )
    project, _ = update_project_structure(project)
    return project


if __name__ == "__main__":
    raise SystemExit(main())
