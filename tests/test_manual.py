"""Check that the PDF user manual has not drifted from the current code.

The manual contains real terminal output, error messages and JSON. If the code changes
and the manual is not rebuilt, it silently goes stale. So there are three layers of defense:
  1. Is the committed manual.html identical to one generated from the current code?
  2. Was the committed PDF built from that manual.html? (source hash)
  3. Do README / TUTORIAL link to the manual?

To rebuild:  bash docs/manual/build.sh
"""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MANUAL_DIR = ROOT / "docs" / "manual"
PDF = MANUAL_DIR / "abcheck_manual.pdf"
HTML = MANUAL_DIR / "manual.html"
STAMP = MANUAL_DIR / "abcheck_manual.pdf.source-sha256"
REBUILD = "bash docs/manual/build.sh 로 다시 만들어 커밋하라."


@pytest.fixture(scope="module")
def builder():
    spec = importlib.util.spec_from_file_location("build_manual", MANUAL_DIR / "build_manual.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pdf_exists_and_is_a_pdf():
    assert PDF.exists(), f"설명서 PDF 가 없다. {REBUILD}"
    data = PDF.read_bytes()
    assert data.startswith(b"%PDF"), "PDF 형식이 아니다"
    assert len(data) > 100_000, "PDF 가 너무 작다 — 스크린샷이 빠졌을 수 있다"


def test_screenshots_exist():
    for name in ["report_summary.png", "report_detail.png", "report_failed.png"]:
        assert (MANUAL_DIR / "img" / name).exists(), f"{name} 가 없다. {REBUILD}"


def test_committed_manual_matches_current_code(builder):
    """Do the manual's terminal output, error messages and JSON match the current code's actual output?"""
    assert HTML.exists(), f"manual.html 이 없다. {REBUILD}"
    committed = HTML.read_text(encoding="utf-8")
    current = builder.build_html()
    assert committed == current, f"설명서가 현재 코드의 출력과 다르다. {REBUILD}"


def test_pdf_was_built_from_current_manual():
    """Does the source HTML hash of the PDF match the current manual.html -- i.e. is the PDF not stale?"""
    assert STAMP.exists(), f"원본 해시 파일이 없다. {REBUILD}"
    recorded = STAMP.read_text(encoding="utf-8").strip()
    actual = hashlib.sha256(HTML.read_bytes()).hexdigest()
    assert recorded == actual, f"PDF 가 지금 manual.html 로 만든 것이 아니다. {REBUILD}"


def test_manual_is_linked_from_readme_and_tutorial():
    link = "docs/manual/abcheck_manual.pdf"
    assert link in (ROOT / "README.md").read_text(encoding="utf-8")
    assert link in (ROOT / "TUTORIAL.md").read_text(encoding="utf-8")


def test_manual_covers_input_output_and_usage():
    """Are all three things the user asked about -- input, results, usage -- covered?"""
    text = HTML.read_text(encoding="utf-8")
    for section in ["입력 — 무엇을 어떤 형식으로 넣는가", "실행 — 명령 네 가지",
                    "터미널 출력 읽는 법", "HTML 리포트", "JSON",
                    "'신뢰 불가'가 나왔을 때", "오류 메시지와 대처", "실전"]:
        assert section in text, f"설명서에 '{section}' 절이 없다"


def test_manual_uses_only_synthetic_examples():
    """Every example sequence in the manual must come from the synthetic corpus."""
    import json
    import re

    corpus = json.loads((ROOT / "data/synthetic/sequences.json").read_text(encoding="utf-8"))
    known = {r["sequence"] for r in corpus["sequences"]}
    text = HTML.read_text(encoding="utf-8")
    for seq in set(re.findall(r"[ACDEFGHIKLMNPQRSTVWY]{100,}", text)):
        assert any(seq in k or k in seq for k in known), "설명서에 합성 코퍼스 밖의 긴 서열이 있다"
