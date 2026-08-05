"""Small local web server for the SNaC grid GUI."""

from __future__ import annotations

import argparse
from ipaddress import ip_address
import json
import mimetypes
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

from .case_import import import_case
from .decomposition import optimize_project_decomposition
from .export import export_project
from .geometry import align_blocks, duplicate_blocks, mirror_blocks, snap_block_to_face
from .grid import AXIS_NAMES, axis_grid_arrays, axis_grid_diagnostics
from .model import Project
from .repair import repair_project_grids
from .snac_bc import apply_bc_preset
from .validation import (
    apply_axis_to_aligned_blocks,
    check_project,
    update_project_structure,
)

ROOT = Path(__file__).resolve().parent
STATIC = ROOT / "static"
_MAX_REQUEST_BYTES = 16 * 1024 * 1024
_REQUEST_HEADER = "X-SNaC-Grid-Request"
_SECURITY_HEADERS = {
    "Content-Security-Policy": (
        "default-src 'self'; script-src 'self' 'unsafe-inline'; "
        "style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; "
        "object-src 'none'; base-uri 'none'; frame-ancestors 'none'"
    ),
    "Cross-Origin-Resource-Policy": "same-origin",
    "Referrer-Policy": "no-referrer",
    "X-Content-Type-Options": "nosniff",
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the SNaC multi-block grid GUI.")
    parser.add_argument("--host", default="127.0.0.1", help="Bind address")
    parser.add_argument("--port", type=int, default=8765, help="Bind port")
    parser.add_argument("--no-open", action="store_true", help="Do not open a browser tab")
    parser.add_argument(
        "--allow-remote",
        action="store_true",
        help="Allow a non-loopback bind address",
    )
    args = parser.parse_args(argv)
    if not args.allow_remote and not _is_loopback_host(args.host):
        parser.error("non-loopback --host requires --allow-remote")

    server = ThreadingHTTPServer((args.host, args.port), _Handler)
    server.allow_remote = args.allow_remote
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
        if not self._request_origin_allowed():
            return
        path = urlparse(self.path).path
        if path == "/api/sample":
            self._send_json(_sample_project().to_dict())
            return
        self._serve_static(path)

    def do_POST(self) -> None:
        if not self._request_origin_allowed(require_marker=True):
            return
        path = urlparse(self.path).path
        if path not in {
            "/api/apply-axis",
            "/api/check",
            "/api/bc-preset",
            "/api/decompose",
            "/api/export",
            "/api/geometry",
            "/api/import",
            "/api/migrate",
            "/api/preview",
            "/api/repair",
            "/api/update",
        }:
            self.send_error(HTTPStatus.NOT_FOUND, "unknown API route")
            return
        if self.headers.get_content_type() != "application/json":
            self._send_json(
                {"ok": False, "error": "API requests require application/json"},
                status=HTTPStatus.UNSUPPORTED_MEDIA_TYPE,
            )
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length < 1:
                raise ValueError("API request body is empty")
            if length > _MAX_REQUEST_BYTES:
                self._send_json(
                    {"ok": False, "error": "API request body is too large"},
                    status=HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
                )
                return
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if path == "/api/import":
                result = import_case(str(payload.get("caseDir", "")))
                self._send_json(result.to_dict())
                return
            project = Project.from_dict(payload.get("project", payload))
            if path == "/api/migrate":
                self._send_json({"ok": True, "project": project.to_dict()})
                return
            if path == "/api/geometry":
                project, changed = _apply_geometry_operation(project, payload)
                self._send_json({"ok": True, "project": project.to_dict(), "changedBlockIds": changed})
                return
            if path == "/api/bc-preset":
                project, changed = apply_bc_preset(
                    project,
                    payload.get("blockIds", []),
                    str(payload.get("preset", "")),
                    face=payload.get("face"),
                    velocity=payload.get("velocity", [0.0, 0.0, 0.0]),
                    pressure=float(payload.get("pressure", 0.0)),
                )
                self._send_json({"ok": True, "project": project.to_dict(), "changedBlockIds": changed})
                return
            if path == "/api/check":
                result = check_project(project)
                self._send_json(result.to_dict())
                return
            if path == "/api/update":
                source_block_id = payload.get("sourceBlockId")
                project, result = update_project_structure(project, int(source_block_id) if source_block_id else None)
                self._send_json({**result.to_dict(), "project": project.to_dict()})
                return
            if path == "/api/apply-axis":
                source_block_id = int(payload.get("sourceBlockId", 0))
                axis = str(payload.get("axis", ""))
                project, changed = apply_axis_to_aligned_blocks(project, source_block_id, axis)
                self._send_json({"ok": True, "project": project.to_dict(), "changedBlockIds": changed})
                return
            if path == "/api/repair":
                source_block_id = payload.get("sourceBlockId")
                project, result = repair_project_grids(
                    project,
                    int(source_block_id) if source_block_id else None,
                )
                self._send_json({**result.to_dict(), "project": project.to_dict()})
                return
            if path == "/api/decompose":
                project, result = optimize_project_decomposition(project)
                self._send_json({**result.to_dict(), "project": project.to_dict()})
                return
            if path == "/api/preview":
                block_id = int(payload.get("blockId", 0))
                block = next((item for item in project.blocks if item.id == block_id), None)
                if block is None:
                    raise ValueError(f"block {block_id} does not exist")
                axes = {}
                for index, name in enumerate(AXIS_NAMES):
                    arrays = axis_grid_arrays(
                        block.axes[name].to_dict(),
                        block.lmin[index],
                        block.lmax[index],
                        block.ng[index],
                    )
                    axes[name] = {
                        "faces": arrays.faces[:-1].tolist(),
                        **axis_grid_diagnostics(
                            block.axes[name].to_dict(),
                            block.lmin[index],
                            block.lmax[index],
                            block.ng[index],
                        ),
                    }
                self._send_json({"ok": True, "blockId": block.id, "ng": block.ng, "axes": axes})
                return
            output_dir = payload.get("outputDir") or "generated/snac_grid_case"
            result = export_project(project, output_dir)
        except Exception as exc:
            self._send_json({"ok": False, "error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
            return
        self._send_json({"ok": True, **result.to_dict()})

    def do_OPTIONS(self) -> None:
        self.send_error(HTTPStatus.METHOD_NOT_ALLOWED, "cross-origin API access is disabled")

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
        self._send_security_headers()
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, data: dict, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = (json.dumps(data, indent=2) + "\n").encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self._send_security_headers()
        self.end_headers()
        self.wfile.write(body)

    def _request_origin_allowed(self, *, require_marker: bool = False) -> bool:
        host_header = self.headers.get("Host", "")
        host, port = _parse_authority(host_header, self.server.server_port)
        allow_remote = bool(getattr(self.server, "allow_remote", False))
        if host is None or (not allow_remote and not _is_loopback_host(host)):
            self.send_error(HTTPStatus.FORBIDDEN, "request host is not allowed")
            return False

        origin = self.headers.get("Origin")
        if origin:
            try:
                parsed = urlparse(origin)
                origin_port = parsed.port or (80 if parsed.scheme == "http" else 443)
            except ValueError:
                self.send_error(HTTPStatus.FORBIDDEN, "invalid request origin")
                return False
            if (
                parsed.scheme != "http"
                or _normalize_host(parsed.hostname) != host
                or origin_port != port
            ):
                self.send_error(HTTPStatus.FORBIDDEN, "cross-origin requests are disabled")
                return False
        if self.headers.get("Sec-Fetch-Site") == "cross-site":
            self.send_error(HTTPStatus.FORBIDDEN, "cross-origin requests are disabled")
            return False
        if require_marker and self.headers.get(_REQUEST_HEADER) != "1":
            self.send_error(HTTPStatus.FORBIDDEN, "missing local API request marker")
            return False
        return True

    def _send_security_headers(self) -> None:
        for name, value in _SECURITY_HEADERS.items():
            self.send_header(name, value)


def _normalize_host(host: str | None) -> str | None:
    if host is None:
        return None
    return host.rstrip(".").lower()


def _parse_authority(authority: str, default_port: int) -> tuple[str | None, int]:
    try:
        parsed = urlparse(f"//{authority}")
        return _normalize_host(parsed.hostname), parsed.port or default_port
    except ValueError:
        return None, default_port


def _is_loopback_host(host: str) -> bool:
    normalized = _normalize_host(host)
    if normalized == "localhost":
        return True
    try:
        return bool(normalized and ip_address(normalized).is_loopback)
    except ValueError:
        return False


def _apply_geometry_operation(project: Project, payload: dict) -> tuple[Project, list[int]]:
    operation = str(payload.get("operation", ""))
    block_ids = payload.get("blockIds", [])
    axis = str(payload.get("axis", "x"))
    if operation == "duplicate":
        return duplicate_blocks(
            project,
            block_ids,
            axis,
            count=int(payload.get("count", 1)),
            gap=float(payload.get("gap", 0.0)),
        )
    if operation == "mirror":
        return mirror_blocks(project, block_ids, axis, float(payload.get("plane", 0.0)))
    if operation == "align":
        return align_blocks(
            project,
            block_ids,
            int(payload.get("sourceBlockId", 0)),
            axis,
            str(payload.get("mode", "lower")),
        )
    if operation == "snap":
        if len(block_ids) != 1:
            raise ValueError("select exactly one block to snap")
        return snap_block_to_face(
            project,
            int(block_ids[0]),
            int(payload.get("sourceBlockId", 0)),
            str(payload.get("face", "")),
        )
    raise ValueError(f"unknown geometry operation {operation!r}")


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
