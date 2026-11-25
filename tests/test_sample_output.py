"""Check that the committed `sample_output/report.html` matches the current code's output.

The report includes the generation date, so comparing it as-is would break daily. Only the
date is substituted and the rest is compared. This way pytest catches "changed the output
format but committed without regenerating the sample report", and the pre-push hook and CI block it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from abcheck.render import as_html
from abcheck.report import build_many

REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_HTML = REPO_ROOT / "sample_output" / "report.html"

#: ISO date inside the report. This is the only varying value.
_ISO_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")

REGENERATE_HINT = (
    "sample_output/report.html 이 현재 코드의 출력과 다르다. "
    "`python3 -m abcheck.cli data/synthetic/sequences.json --html` 로 다시 만들어 커밋하라."
)


def _normalize(html: str) -> str:
    """Replace dates with a fixed string so the reports can be compared."""
    return _ISO_DATE.sub("DATE", html)


def _current_html(records: list[dict]) -> str:
    """Generate in memory the same report the CLI builds from the corpus."""
    return as_html(build_many([(r["id"], r["sequence"]) for r in records]))


def test_sample_output_exists():
    assert SAMPLE_HTML.exists(), (
        f"{SAMPLE_HTML.relative_to(REPO_ROOT)} 가 없다. {REGENERATE_HINT}"
    )


def test_committed_sample_output_matches_current_render(records):
    """Apart from the date, the content must be identical to the current code's output."""
    committed = _normalize(SAMPLE_HTML.read_text(encoding="utf-8"))
    current = _normalize(_current_html(records))

    if committed != current:
        # Show at least one line of what changed
        import difflib

        diff = list(difflib.unified_diff(
            committed.splitlines(), current.splitlines(),
            fromfile="커밋된 sample_output/report.html",
            tofile="현재 코드의 출력",
            lineterm="", n=1,
        ))
        pytest.fail(REGENERATE_HINT + "\n\n" + "\n".join(diff[:40]))


def test_normalization_only_touches_dates():
    """Date substitution must not touch other content -- keeps the comparison from getting loose."""
    sample = '<p>생성 2026-11-24 · 서열 16건</p><td>31-35</td><td>95-104</td>'
    assert _normalize(sample) == '<p>생성 DATE · 서열 16건</p><td>31-35</td><td>95-104</td>'


def test_stale_sample_output_is_detected(records):
    """Is a format change actually caught? Checks that this test itself is not a no-op."""
    committed = _normalize(SAMPLE_HTML.read_text(encoding="utf-8"))
    tampered = committed.replace("CDR별 liability 집계", "집계", 1)

    assert tampered != committed, "치환이 적용되지 않았다"
    assert tampered != _normalize(_current_html(records))


def test_sample_output_covers_whole_corpus(records):
    html = SAMPLE_HTML.read_text(encoding="utf-8")
    for rec in records:
        assert rec["id"] in html, f"{rec['id']} 가 예시 리포트에 없다. {REGENERATE_HINT}"
