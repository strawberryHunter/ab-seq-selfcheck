"""Output tests -- does the limitations notice appear in all three outputs, and does HTML respect the style constraints?"""

from __future__ import annotations

import html as html_module
import json
import re

from abcheck import LIMITATION_NOTICE, NO_VERDICT_NOTICE
from abcheck.render import as_html, as_json, as_text
from abcheck.report import SUMMARY_NOTE, build, build_many

from .conftest import by_id


def _reports(records, *ids):
    return build_many([(r["id"], r["sequence"]) for r in records if r["id"] in ids])


# --------------------------------------------------------------------------
# Limitations notice -- must appear in all three outputs
# --------------------------------------------------------------------------

def test_text_output_carries_limitation_notice(records):
    text = as_text(_reports(records, "SYN-VH-001"))
    #  Line wrapping prevents comparing it whole, so check distinctive fragments
    assert "numbering 엔진이" in text
    assert "ANARCI" in text
    assert "종합 판정" in text


def test_html_output_carries_limitation_notice(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    # HTML escapes quotes too, so compare against the escaped form
    for note in LIMITATION_NOTICE:
        assert html_module.escape(note) in out, note[:30]
    assert html_module.escape(NO_VERDICT_NOTICE) in out


def test_json_output_carries_limitation_notice(records):
    payload = json.loads(as_json(_reports(records, "SYN-VH-001")))
    assert payload["limitations"] == list(LIMITATION_NOTICE)
    assert payload["no_verdict"] == NO_VERDICT_NOTICE
    # Must also be in each per-sequence record -- the notice travels with any single record
    assert payload["sequences"][0]["limitations"] == list(LIMITATION_NOTICE)
    assert payload["sequences"][0]["no_verdict"] == NO_VERDICT_NOTICE


def test_no_output_contains_overall_verdict_wording(records):
    """Never emit overall verdict phrases like 'high/low developability'.

    The notice text that declares "we do not make verdicts" is the exception. Whitespace in the
    output is normalized and the notice removed before checking -- the terminal output wraps the notice.
    """
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    # The two notices declaring that no verdict is made are removed as exceptions
    notices = [re.sub(r"\s+", " ", n) for n in (NO_VERDICT_NOTICE, SUMMARY_NOTE)]

    for text in [as_text(reports), as_html(reports), as_json(reports)]:
        flat = re.sub(r"\s+", " ", text)
        for notice in notices:
            flat = flat.replace(notice, "").replace(html_module.escape(notice), "")
        assert "개발 가능성" not in flat
        assert "developable" not in flat.lower()
        assert "종합 판정" not in flat
        assert "점수" not in flat
        assert "순위" not in flat


# --------------------------------------------------------------------------
# Is the confident=False reason visible in the output?
# --------------------------------------------------------------------------

def test_failure_reason_appears_in_every_output(records):
    reports = _reports(records, "SYN-VH-901")
    reason = reports[0].annotation.reasons[0]

    html_out = as_html(reports)
    assert reason in html_out
    assert "신뢰 불가" in html_out
    assert "위치를 추측하지 않는다" in html_out

    payload = json.loads(as_json(reports))
    assert payload["sequences"][0]["annotation"]["reasons"] == [reason]
    assert payload["sequences"][0]["annotation"]["confident"] is False

    text_out = as_text(reports)
    assert "신뢰 불가" in text_out
    assert "말단 모티프를 찾지 못했다" in text_out


def test_length_outlier_reason_appears(records):
    text_out = as_text(_reports(records, "SYN-VH-903"))
    assert "통상 범위" in text_out
    assert "ANARCI/IMGT 로 경계를 재확인" in text_out


def test_invalid_input_is_reported_not_crashed():
    reports = [build("bad", "ACDEF")]
    for out in [as_text(reports), as_html(reports), as_json(reports)]:
        assert "너무 짧다" in out
    assert "입력을 쓸 수 없다" in as_text(reports)


# --------------------------------------------------------------------------
# HTML style constraints -- Arial, black header, gray zebra, no color highlighting
# --------------------------------------------------------------------------

def test_html_uses_arial(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    assert "Arial" in out


def test_html_header_is_black_and_rows_are_gray_zebra(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    assert "background: #000000" in out          # Black header
    assert "nth-child(even) { background: #f0f0f0" in out   # Gray zebra


def test_html_has_no_chromatic_colors(records):
    """No colors other than gray/black-and-white -- a constraint so risk is never color-coded."""
    out = as_html(build_many([(r["id"], r["sequence"]) for r in records]))

    hex_colors = set(re.findall(r"#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})\b", out))
    for value in hex_colors:
        if len(value) == 3:
            value = "".join(ch * 2 for ch in value)
        r, g, b = (int(value[i:i + 2], 16) for i in (0, 2, 4))
        assert r == g == b, f"무채색이 아닌 색이 있다: #{value}"

    # No color names or other color notations either
    for banned in ["rgb(", "hsl(", "red", "green", "blue", "orange", "yellow", "crimson"]:
        assert banned not in out.lower(), banned


def test_html_marks_cdr_with_bold_only(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    ann = _reports(records, "SYN-VH-001")[0].annotation

    # CDR sequences appear wrapped in <b>
    for cdr in ann.cdrs.values():
        assert f"<b>{cdr.seq}</b>" in out, cdr.name

    # Emphasis is bold only -- no underline, italics or highlighter
    for banned in ["<u>", "<em>", "<i>", "<mark>", "text-decoration", "font-style: italic"]:
        assert banned not in out, banned


def test_html_sequence_block_splits_into_lines(records):
    from abcheck.render import SEQUENCE_LINE_WIDTH

    rec = by_id(records, "SYN-VH-001")
    out = as_html(_reports(records, "SYN-VH-001"))
    expected_rows = -(-len(rec["sequence"]) // SEQUENCE_LINE_WIDTH)
    assert out.count('<td class="num">') == expected_rows


def test_html_is_single_self_contained_file(records):
    """Single HTML file -- loads no external CSS/JS/images."""
    out = as_html(build_many([(r["id"], r["sequence"]) for r in records]))
    assert "<style>" in out
    for banned in ["<link", "<script", "<img", "http://", "https://", "@import"]:
        assert banned not in out, banned


def test_html_escapes_content():
    reports = [build("<script>alert(1)</script>", "A" * 120)]
    out = as_html(reports)
    assert "<script>" not in out
    assert "&lt;script&gt;" in out


# --------------------------------------------------------------------------
# Terminal table -- no color, aligned
# --------------------------------------------------------------------------

def test_text_output_has_no_ansi_escapes(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    assert "\x1b[" not in out


def test_text_table_columns_are_aligned(records):
    out = as_text(_reports(records, "SYN-VH-001"))
    # The CDR table has 4 columns. The property table also has rows starting with CDR-H, so filter by column count.
    rows = [
        ln for ln in out.splitlines()
        if ln.strip().startswith("CDR-H") and ln.count("|") == 3
    ]
    assert len(rows) == 3
    # The '|' separator must be at the same position on every row
    positions = {tuple(i for i, ch in enumerate(row) if ch == "|") for row in rows}
    assert len(positions) == 1, f"열 정렬이 틀어졌다: {positions}"


def test_text_marks_cdr_positions_under_sequence(records):
    out = as_text(_reports(records, "SYN-VH-001"))
    assert "^ 표시가 CDR 구간" in out
    assert "^^^^^" in out
    assert "CDR-H1" in out


def test_text_wraps_without_repeating_bullet(records):
    """Bullets are not repeated in wrapped description text."""
    out = as_text(_reports(records, "SYN-VH-903"))
    bullet_lines = [ln for ln in out.splitlines() if ln.strip().startswith("* ")]
    for line in bullet_lines:
        assert line.count("*") == 1


# --------------------------------------------------------------------------
# JSON structure
# --------------------------------------------------------------------------

def test_json_round_trips_cdr_coordinates(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    payload = json.loads(as_json(reports))

    assert payload["count"] == len(records)
    for rec, entry in zip(records, payload["sequences"]):
        assert entry["id"] == rec["id"]
        expected = rec["expected"]["cdrs"]
        if expected is None:
            assert entry["annotation"]["cdrs"] == {}
            continue
        for name, want in expected.items():
            got = entry["annotation"]["cdrs"][name]
            assert (got["start"], got["end"], got["seq"]) == (
                want["start"], want["end"], want["seq"]
            )


def test_json_includes_liabilities_and_properties(records):
    payload = json.loads(as_json(_reports(records, "SYN-VH-007")))
    entry = payload["sequences"][0]

    dg = [h for h in entry["liabilities"]["hits"] if h["match"] == "DG"]
    assert dg and dg[0]["risk"] == "높음" and dg[0]["region"] == "CDR-H3"
    assert "rationale" in dg[0] and "risk_note" in dg[0]

    assert "도메인 전체" in entry["properties"]
    assert set(entry["properties"]) >= {"CDR-H1", "CDR-H2", "CDR-H3"}
    assert entry["properties"]["도메인 전체"]["pH"] == 7.4


def test_json_includes_domain_flags(records):
    payload = json.loads(as_json(_reports(records, "SYN-VH-005")))
    flags = payload["sequences"][0]["liabilities"]["domain_flags"]
    assert len(flags) == 1
    assert flags[0]["kind"] == "unpaired_cysteine"


def test_json_states_coordinate_convention(records):
    payload = json.loads(as_json(_reports(records, "SYN-VH-001")))
    assert "1-based" in payload["coordinates"]
    assert "1-based" in payload["sequences"][0]["annotation"]["coordinates"]


# --------------------------------------------------------------------------
# Usability pass (phone layout, "what to open first")
# --------------------------------------------------------------------------

def test_summary_names_the_sequences_to_open_first(corpus_reports):
    """The triage box used to give only a count; the user had to scan every row to find them."""
    from abcheck.render import as_html, as_text

    text = as_text(corpus_reports)
    assert "먼저 볼 서열" in text
    assert "검토 필요 6건: SYN-VH-003" in text
    assert "ANARCI/IMGT 재확인 필요 (주석 신뢰 불가) 5건" in text

    page = as_html(corpus_reports)
    assert "먼저 볼 서열" in page
    assert 'href="#seq-SYN-VH-003"' in page


def test_every_html_table_scrolls_inside_its_own_box(corpus_reports):
    """On a phone, wide tables must scroll in place instead of widening the whole page."""
    import re

    from abcheck.render import as_html

    page = as_html(corpus_reports)
    tables = re.findall(r"(.{0,40})<table", page)
    assert tables
    for before in tables:
        assert before.endswith('<div class="scroll">'), before
    assert "word-break: keep-all" in page        # Korean breaks at spaces, not per syllable


def test_triage_reason_does_not_repeat_the_motif(corpus_reports):
    from abcheck.render import as_text

    text = as_text(corpus_reports)
    for motif in ["NG", "DG", "DP"]:
        assert f"({motif}) ({motif})" not in text


import pytest  # noqa: E402


@pytest.fixture(scope="module")
def corpus_reports(records):
    return build_many([(r["id"], r["sequence"]) for r in records])
