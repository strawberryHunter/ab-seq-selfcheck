"""CDR별 liability 집계 표 테스트.

집계는 후보 여러 건을 나란히 비교하기 위한 것이다. 점수를 매기거나 순위를
내지 않고 건수만 센다 — 그 선을 지키는지도 함께 확인한다.
"""

from __future__ import annotations

import json
import re

from abcheck.render import _display_width, as_html, as_json, as_text
from abcheck.report import SUMMARY_NOTE, SUMMARY_REGIONS, build, build_many, summarize

from .conftest import by_id


def _reports(records, *ids):
    return build_many([(r["id"], r["sequence"]) for r in records if r["id"] in ids])


def _row(summary: dict, seq_id: str) -> dict:
    return next(r for r in summary["rows"] if r["id"] == seq_id)


# --------------------------------------------------------------------------
# 집계 내용
# --------------------------------------------------------------------------

def test_summary_has_one_row_per_sequence(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    summary = summarize(reports)

    assert len(summary["rows"]) == len(records)
    assert [r["id"] for r in summary["rows"]] == [r["id"] for r in records]
    assert summary["regions"] == list(SUMMARY_REGIONS)


def test_cdr_hit_lands_in_the_right_cdr_column(records):
    """CDR-H2 의 NG 는 CDR2 칸에, CDR-H3 의 DG 는 CDR3 칸에 들어간다."""
    summary = summarize(_reports(records, "SYN-VH-003", "SYN-VH-007"))

    h003 = _row(summary, "SYN-VH-003")["regions"]
    assert h003["CDR2"]["높음"] == 1
    assert h003["CDR1"] == {"높음": 0, "중간": 0, "낮음": 0}
    assert h003["CDR3"] == {"높음": 0, "중간": 0, "낮음": 0}

    h007 = _row(summary, "SYN-VH-007")["regions"]
    assert h007["CDR3"]["높음"] == 1
    assert h007["CDR2"] == {"높음": 0, "중간": 0, "낮음": 0}


def test_light_chain_cdrs_share_the_same_columns(records):
    """CDR-L1 은 CDR-H1 과 같은 CDR1 칸을 쓴다 — 중쇄·경쇄를 한 표에서 비교하기 위함."""
    summary = summarize(_reports(records, "SYN-VL-004"))
    regions = _row(summary, "SYN-VL-004")["regions"]
    assert regions["CDR1"]["높음"] == 1          # CDR-L1 의 DP


def test_framework_hits_are_separated_from_cdr(records):
    summary = summarize(_reports(records, "SYN-VH-006"))
    regions = _row(summary, "SYN-VH-006")["regions"]

    assert regions["framework"]["중간"] == 1     # framework DG
    assert regions["framework"]["낮음"] == 3     # Met + 보존 Trp 2개
    for cdr in ["CDR1", "CDR2", "CDR3"]:
        assert sum(regions[cdr].values()) == 0


def test_normal_sequences_have_no_cdr_counts(records):
    summary = summarize(_reports(records, "SYN-VH-001", "SYN-VL-001"))
    for seq_id in ["SYN-VH-001", "SYN-VL-001"]:
        regions = _row(summary, seq_id)["regions"]
        for cdr in ["CDR1", "CDR2", "CDR3"]:
            assert sum(regions[cdr].values()) == 0, seq_id


def test_domain_flags_are_listed(records):
    summary = summarize(_reports(records, "SYN-VH-005", "SYN-VH-001"))

    flags = _row(summary, "SYN-VH-005")["domain_flags"]
    assert len(flags) == 1
    assert flags[0]["risk"] == "높음"
    assert "Cys" in flags[0]["label"]

    assert _row(summary, "SYN-VH-001")["domain_flags"] == []
    assert summary["domain_flag_total"] == 1


def test_totals_match_the_sum_of_rows(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    summary = summarize(reports)

    for region in SUMMARY_REGIONS:
        for level in summary["risk_levels"]:
            expected = sum(r["regions"][region][level] for r in summary["rows"])
            assert summary["totals"][region][level] == expected, (region, level)


def test_totals_account_for_every_hit(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    summary = summarize(reports)

    counted = sum(
        n for region in SUMMARY_REGIONS for n in summary["totals"][region].values()
    )
    actual = sum(len(r.liabilities.hits) for r in reports if r.ok)
    assert counted == actual


# --------------------------------------------------------------------------
# 주석이 실패한 서열 / 입력 오류
# --------------------------------------------------------------------------

def test_unannotated_sequence_is_marked_and_not_attributed_to_cdr(records):
    """구간을 가를 수 없는 서열의 모티프를 CDR 칸에 넣지 않는다."""
    summary = summarize(_reports(records, "SYN-VH-901"))
    row = _row(summary, "SYN-VH-901")

    assert row["annotated"] is False
    assert row["regions_unassigned"] is True
    assert row["confident"] is False
    for cdr in ["CDR1", "CDR2", "CDR3"]:
        assert sum(row["regions"][cdr].values()) == 0
    assert sum(row["regions"]["framework"].values()) == 1


def test_length_outlier_is_counted_normally(records):
    """앵커가 멀쩡한 서열은 confident=False 여도 CDR별로 집계된다."""
    summary = summarize(_reports(records, "SYN-VH-903"))
    row = _row(summary, "SYN-VH-903")
    assert row["annotated"] is True
    assert row["confident"] is False
    assert row["regions_unassigned"] is False


def test_invalid_input_row_is_present_but_empty():
    summary = summarize([build("bad", "ACDEF")])
    row = _row(summary, "bad")

    assert row["input_ok"] is False
    assert row["domain_flags"] == []
    for region in SUMMARY_REGIONS:
        assert sum(row["regions"][region].values()) == 0


# --------------------------------------------------------------------------
# 출력에 실리는지
# --------------------------------------------------------------------------

def test_summary_appears_before_sequence_sections_in_html(records):
    out = as_html(_reports(records, "SYN-VH-001", "SYN-VH-007"))

    summary_pos = out.index("CDR별 liability 집계")
    first_section = out.index('class="seq-section"')
    assert summary_pos < first_section, "집계 표가 서열 섹션보다 앞에 와야 한다"


def test_summary_appears_before_sequence_sections_in_text(records):
    out = as_text(_reports(records, "SYN-VH-001", "SYN-VH-007"))
    assert out.index("CDR별 liability 집계") < out.index("[SYN-VH-001]")


def test_html_summary_links_to_each_sequence_section(records):
    out = as_html(_reports(records, "SYN-VH-001", "SYN-VH-007"))
    for seq_id in ["SYN-VH-001", "SYN-VH-007"]:
        anchor = "seq-" + seq_id.replace("-", "-")
        assert f'href="#{anchor}"' in out
        assert f'id="{anchor}"' in out


def test_text_summary_table_is_aligned(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    section = out.split("CDR별 liability 집계")[1].split("[SYN-VH-001]")[0]
    rows = [ln for ln in section.splitlines() if ln.strip().startswith("SYN-")]

    assert len(rows) == len(records)
    # 한글은 표시 폭이 2칸이므로 문자열 인덱스가 아니라 표시 폭으로 비교해야 한다
    positions = {
        tuple(_display_width(row[:i]) for i, ch in enumerate(row) if ch == "|")
        for row in rows
    }
    assert len(positions) == 1, f"집계 표 열 정렬이 틀어졌다: {positions}"


def test_text_summary_has_totals_row(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    assert "합계" in out


def test_json_carries_summary(records):
    payload = json.loads(as_json(build_many([(r["id"], r["sequence"]) for r in records])))

    assert "summary" in payload
    assert payload["summary"]["note"] == SUMMARY_NOTE
    assert len(payload["summary"]["rows"]) == len(records)
    assert payload["summary"]["regions"] == list(SUMMARY_REGIONS)
    assert payload["summary"]["risk_levels"] == ["높음", "중간", "낮음"]


def test_summary_note_warns_it_is_not_a_verdict(records):
    for out in [as_text(_reports(records, "SYN-VH-001")),
                as_html(_reports(records, "SYN-VH-001"))]:
        flat = re.sub(r"\s+", " ", out)
        assert "종합 판정이 아니다" in flat
        assert "탈락시킬 근거도 되지 않는다" in flat


def test_summary_does_not_rank_or_score(records):
    """집계는 건수만 낸다 — 점수·순위·정렬 재배치를 하지 않는다."""
    records_order = [r["id"] for r in records]
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    assert [r["id"] for r in summary["rows"]] == records_order
    assert not any(k in summary for k in ["score", "rank", "ranking", "점수", "순위"])
