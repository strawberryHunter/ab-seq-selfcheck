"""PDF 사용 설명서가 현재 코드와 어긋나지 않는지 확인한다.

설명서에는 실제 터미널 출력·오류 메시지·JSON 이 들어 있다. 코드를 바꾸고 설명서를
다시 만들지 않으면 조용히 낡는다. 그래서 세 겹으로 막는다.
  1. 커밋된 manual.html 이 지금 코드로 생성한 것과 같은가
  2. 커밋된 PDF 가 그 manual.html 로 만든 것인가 (원본 해시)
  3. README / TUTORIAL 에서 설명서로 가는 링크가 있는가

다시 만들려면:  bash docs/manual/build.sh
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
    """설명서의 터미널 출력·오류 메시지·JSON 이 지금 코드의 실제 출력과 같은지."""
    assert HTML.exists(), f"manual.html 이 없다. {REBUILD}"
    committed = HTML.read_text(encoding="utf-8")
    current = builder.build_html()
    assert committed == current, f"설명서가 현재 코드의 출력과 다르다. {REBUILD}"


def test_pdf_was_built_from_current_manual():
    """PDF 를 만든 원본 HTML 해시가 지금 manual.html 과 같은지 — PDF 가 낡지 않았는지."""
    assert STAMP.exists(), f"원본 해시 파일이 없다. {REBUILD}"
    recorded = STAMP.read_text(encoding="utf-8").strip()
    actual = hashlib.sha256(HTML.read_bytes()).hexdigest()
    assert recorded == actual, f"PDF 가 지금 manual.html 로 만든 것이 아니다. {REBUILD}"


def test_manual_is_linked_from_readme_and_tutorial():
    link = "docs/manual/abcheck_manual.pdf"
    assert link in (ROOT / "README.md").read_text(encoding="utf-8")
    assert link in (ROOT / "TUTORIAL.md").read_text(encoding="utf-8")


def test_manual_covers_input_output_and_usage():
    """사용자가 물어본 세 가지 — 입력·결과·사용법 — 가 다 들어 있는지."""
    text = HTML.read_text(encoding="utf-8")
    for section in ["입력 — 무엇을 어떤 형식으로 넣는가", "실행 — 명령 네 가지",
                    "터미널 출력 읽는 법", "HTML 리포트", "JSON",
                    "'신뢰 불가'가 나왔을 때", "오류 메시지와 대처", "실전"]:
        assert section in text, f"설명서에 '{section}' 절이 없다"


def test_manual_uses_only_synthetic_examples():
    """설명서의 예시 서열은 전부 합성 코퍼스에 있는 것이어야 한다."""
    import json
    import re

    corpus = json.loads((ROOT / "data/synthetic/sequences.json").read_text(encoding="utf-8"))
    known = {r["sequence"] for r in corpus["sequences"]}
    text = HTML.read_text(encoding="utf-8")
    for seq in set(re.findall(r"[ACDEFGHIKLMNPQRSTVWY]{100,}", text)):
        assert any(seq in k or k in seq for k in known), "설명서에 합성 코퍼스 밖의 긴 서열이 있다"
