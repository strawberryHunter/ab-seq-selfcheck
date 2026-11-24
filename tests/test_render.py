"""출력 테스트 — 한계 고지가 세 출력 모두에 실리는지, HTML 양식 제약을 지키는지."""

from __future__ import annotations

import html as html_module
import json
import re

from abcheck import LIMITATION_NOTICE, NO_VERDICT_NOTICE
from abcheck.render import as_html, as_json, as_text
from abcheck.report import build, build_many

from .conftest import by_id


def _reports(records, *ids):
    return build_many([(r["id"], r["sequence"]) for r in records if r["id"] in ids])


# --------------------------------------------------------------------------
# 한계 고지 — 세 출력 모두에 반드시 들어간다
# --------------------------------------------------------------------------

def test_text_output_carries_limitation_notice(records):
    text = as_text(_reports(records, "SYN-VH-001"))
    #  줄바꿈 때문에 통째로 비교할 수 없으니 특징적인 조각으로 확인한다
    assert "numbering 엔진이" in text
    assert "ANARCI" in text
    assert "종합 판정" in text


def test_html_output_carries_limitation_notice(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    # HTML 은 따옴표까지 이스케이프하므로 이스케이프한 형태로 비교한다
    for note in LIMITATION_NOTICE:
        assert html_module.escape(note) in out, note[:30]
    assert html_module.escape(NO_VERDICT_NOTICE) in out


def test_json_output_carries_limitation_notice(records):
    payload = json.loads(as_json(_reports(records, "SYN-VH-001")))
    assert payload["limitations"] == list(LIMITATION_NOTICE)
    assert payload["no_verdict"] == NO_VERDICT_NOTICE
    # 서열별 레코드에도 들어 있어야 한다 — 한 건만 떼어 써도 고지가 따라간다
    assert payload["sequences"][0]["limitations"] == list(LIMITATION_NOTICE)
    assert payload["sequences"][0]["no_verdict"] == NO_VERDICT_NOTICE


def test_no_output_contains_overall_verdict_wording(records):
    """'개발 가능성 높음/낮음' 같은 종합 판정 문구를 내지 않는다.

    "내리지 않는다"고 선언하는 고지문 자체는 예외다. 출력에서 공백을 정규화한 뒤
    고지문을 떼어내고 나서 검사한다 — 터미널 출력은 고지문을 줄바꿈하기 때문이다.
    """
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    notice = re.sub(r"\s+", " ", NO_VERDICT_NOTICE)

    for text in [as_text(reports), as_html(reports), as_json(reports)]:
        flat = re.sub(r"\s+", " ", text)
        flat = flat.replace(notice, "").replace(html_module.escape(notice), "")
        assert "개발 가능성" not in flat
        assert "developable" not in flat.lower()
        assert "종합 판정" not in flat


# --------------------------------------------------------------------------
# confident=False 사유가 출력에 보이는지
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
# HTML 양식 제약 — Arial, 검정 헤더, 회색 지브라, 컬러 강조 없음
# --------------------------------------------------------------------------

def test_html_uses_arial(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    assert "Arial" in out


def test_html_header_is_black_and_rows_are_gray_zebra(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    assert "background: #000000" in out          # 검정 헤더
    assert "nth-child(even) { background: #f0f0f0" in out   # 회색 지브라


def test_html_has_no_chromatic_colors(records):
    """회색/흑백 외의 색을 쓰지 않는다 — 위험도를 색으로 칠하지 않기 위한 제약."""
    out = as_html(build_many([(r["id"], r["sequence"]) for r in records]))

    hex_colors = set(re.findall(r"#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})\b", out))
    for value in hex_colors:
        if len(value) == 3:
            value = "".join(ch * 2 for ch in value)
        r, g, b = (int(value[i:i + 2], 16) for i in (0, 2, 4))
        assert r == g == b, f"무채색이 아닌 색이 있다: #{value}"

    # 색 이름이나 다른 색 표기도 쓰지 않는다
    for banned in ["rgb(", "hsl(", "red", "green", "blue", "orange", "yellow", "crimson"]:
        assert banned not in out.lower(), banned


def test_html_marks_cdr_with_bold_only(records):
    out = as_html(_reports(records, "SYN-VH-001"))
    ann = _reports(records, "SYN-VH-001")[0].annotation

    # CDR 서열이 <b> 로 감싸여 나타난다
    for cdr in ann.cdrs.values():
        assert f"<b>{cdr.seq}</b>" in out, cdr.name

    # 강조는 굵게만 — 밑줄·기울임·형광펜 없음
    for banned in ["<u>", "<em>", "<i>", "<mark>", "text-decoration", "font-style: italic"]:
        assert banned not in out, banned


def test_html_sequence_block_splits_into_lines(records):
    from abcheck.render import SEQUENCE_LINE_WIDTH

    rec = by_id(records, "SYN-VH-001")
    out = as_html(_reports(records, "SYN-VH-001"))
    expected_rows = -(-len(rec["sequence"]) // SEQUENCE_LINE_WIDTH)
    assert out.count('<td class="num">') == expected_rows


def test_html_is_single_self_contained_file(records):
    """단일 HTML — 외부 CSS/JS/이미지를 불러오지 않는다."""
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
# 터미널 표 — 컬러 없음, 정렬
# --------------------------------------------------------------------------

def test_text_output_has_no_ansi_escapes(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    assert "\x1b[" not in out


def test_text_table_columns_are_aligned(records):
    out = as_text(_reports(records, "SYN-VH-001"))
    # CDR 표는 4열이다. 물성 표에도 CDR-H 로 시작하는 행이 있으므로 열 수로 가린다.
    rows = [
        ln for ln in out.splitlines()
        if ln.strip().startswith("CDR-H") and ln.count("|") == 3
    ]
    assert len(rows) == 3
    # 구분자 '|' 위치가 모든 행에서 같아야 한다
    positions = {tuple(i for i, ch in enumerate(row) if ch == "|") for row in rows}
    assert len(positions) == 1, f"열 정렬이 틀어졌다: {positions}"


def test_text_marks_cdr_positions_under_sequence(records):
    out = as_text(_reports(records, "SYN-VH-001"))
    assert "^ 표시가 CDR 구간" in out
    assert "^^^^^" in out
    assert "CDR-H1" in out


def test_text_wraps_without_repeating_bullet(records):
    """줄바꿈된 설명문에 글머리표가 반복되지 않는다."""
    out = as_text(_reports(records, "SYN-VH-903"))
    bullet_lines = [ln for ln in out.splitlines() if ln.strip().startswith("* ")]
    for line in bullet_lines:
        assert line.count("*") == 1


# --------------------------------------------------------------------------
# JSON 구조
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
