"""커밋된 `sample_output/report.html` 이 현재 코드의 출력과 같은지 확인한다.

리포트에는 생성 날짜가 들어가므로 그대로 비교하면 매일 깨진다. 날짜만 치환해
놓고 나머지를 비교한다. 이렇게 하면 "출력 양식을 바꿨는데 예시 리포트를 다시
만들지 않은 채 커밋" 하는 것을 pytest 가 잡고, pre-push 훅과 CI 가 그걸 막는다.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from abcheck.render import as_html
from abcheck.report import build_many

REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_HTML = REPO_ROOT / "sample_output" / "report.html"

#: 리포트 안의 ISO 날짜. 이것만 변동값이다.
_ISO_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")

REGENERATE_HINT = (
    "sample_output/report.html 이 현재 코드의 출력과 다르다. "
    "`python3 -m abcheck.cli data/synthetic/sequences.json --html` 로 다시 만들어 커밋하라."
)


def _normalize(html: str) -> str:
    """날짜를 고정 문자열로 바꿔 비교 가능하게 만든다."""
    return _ISO_DATE.sub("DATE", html)


def _current_html(records: list[dict]) -> str:
    """CLI 가 코퍼스로 만드는 것과 같은 리포트를 메모리에서 생성한다."""
    return as_html(build_many([(r["id"], r["sequence"]) for r in records]))


def test_sample_output_exists():
    assert SAMPLE_HTML.exists(), (
        f"{SAMPLE_HTML.relative_to(REPO_ROOT)} 가 없다. {REGENERATE_HINT}"
    )


def test_committed_sample_output_matches_current_render(records):
    """날짜를 뺀 내용이 현재 코드의 출력과 완전히 같아야 한다."""
    committed = _normalize(SAMPLE_HTML.read_text(encoding="utf-8"))
    current = _normalize(_current_html(records))

    if committed != current:
        # 어디가 달라졌는지 한 줄이라도 보여준다
        import difflib

        diff = list(difflib.unified_diff(
            committed.splitlines(), current.splitlines(),
            fromfile="커밋된 sample_output/report.html",
            tofile="현재 코드의 출력",
            lineterm="", n=1,
        ))
        pytest.fail(REGENERATE_HINT + "\n\n" + "\n".join(diff[:40]))


def test_normalization_only_touches_dates():
    """날짜 치환이 다른 내용을 건드리지 않는지 — 비교가 헐거워지는 것을 막는다."""
    sample = '<p>생성 2026-11-24 · 서열 16건</p><td>31-35</td><td>95-104</td>'
    assert _normalize(sample) == '<p>생성 DATE · 서열 16건</p><td>31-35</td><td>95-104</td>'


def test_stale_sample_output_is_detected(records):
    """양식이 바뀌면 실제로 잡히는지 — 이 테스트 자체가 헛돌지 않는지 확인한다."""
    committed = _normalize(SAMPLE_HTML.read_text(encoding="utf-8"))
    tampered = committed.replace("CDR별 liability 집계", "집계", 1)

    assert tampered != committed, "치환이 적용되지 않았다"
    assert tampered != _normalize(_current_html(records))


def test_sample_output_covers_whole_corpus(records):
    html = SAMPLE_HTML.read_text(encoding="utf-8")
    for rec in records:
        assert rec["id"] in html, f"{rec['id']} 가 예시 리포트에 없다. {REGENERATE_HINT}"
