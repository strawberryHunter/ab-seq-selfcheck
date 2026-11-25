"""The browser app (web/abcheck.html) must be built from the current template and engine.

The engine itself is checked against Python in tests/test_js_parity.py.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "web" / "abcheck.html"


def test_committed_app_is_current():
    result = subprocess.run(
        [sys.executable, "web/build.py", "--check"],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr + " — run: python3 web/build.py"


def test_app_loads_nothing_but_fonts_from_outside():
    """Sequences must never leave the browser: no fetch/XHR, no external scripts."""
    html = APP.read_text(encoding="utf-8")
    hosts = set(re.findall(r'(?:src|href)="(https?://[^/"]+)', html))
    assert hosts <= {"https://fonts.googleapis.com", "https://fonts.gstatic.com", "https://github.com"}, hosts
    for api in ["fetch(", "XMLHttpRequest", "WebSocket", "sendBeacon"]:
        assert api not in html, api


def test_app_carries_the_limitations_and_no_verdict_notice():
    html = APP.read_text(encoding="utf-8")
    assert "<title>AbCheck 항체 서열 점검</title>" in html
    assert "numbering 엔진이 아니다" in html
    assert "종합 판정을 내리지 않는다" in html


def test_example_uses_only_synthetic_sequences():
    html = APP.read_text(encoding="utf-8")
    example = re.search(r"const EXAMPLE = (\".*?\");", html).group(1)
    assert ">SYN-" in example and "OAS" not in example
