"""Tests for the per-CDR liability summary table.

The summary is for comparing several candidates side by side. It only counts; it does not
score or rank -- these tests also check that it stays within that line.
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
# Summary contents
# --------------------------------------------------------------------------

def test_summary_has_one_row_per_sequence(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])
    summary = summarize(reports)

    assert len(summary["rows"]) == len(records)
    assert [r["id"] for r in summary["rows"]] == [r["id"] for r in records]
    assert summary["regions"] == list(SUMMARY_REGIONS)


def test_cdr_hit_lands_in_the_right_cdr_column(records):
    """NG in CDR-H2 goes in the CDR2 cell, DG in CDR-H3 in the CDR3 cell."""
    summary = summarize(_reports(records, "SYN-VH-003", "SYN-VH-007"))

    h003 = _row(summary, "SYN-VH-003")["regions"]
    assert h003["CDR2"]["높음"] == 1
    assert h003["CDR1"] == {"높음": 0, "중간": 0, "낮음": 0}
    assert h003["CDR3"] == {"높음": 0, "중간": 0, "낮음": 0}

    h007 = _row(summary, "SYN-VH-007")["regions"]
    assert h007["CDR3"]["높음"] == 1
    assert h007["CDR2"] == {"높음": 0, "중간": 0, "낮음": 0}


def test_light_chain_cdrs_share_the_same_columns(records):
    """CDR-L1 uses the same CDR1 cell as CDR-H1 -- so heavy and light chains can be compared in one table."""
    summary = summarize(_reports(records, "SYN-VL-004"))
    regions = _row(summary, "SYN-VL-004")["regions"]
    assert regions["CDR1"]["높음"] == 1          # DP in CDR-L1


def test_framework_hits_are_separated_from_cdr(records):
    summary = summarize(_reports(records, "SYN-VH-006"))
    regions = _row(summary, "SYN-VH-006")["regions"]

    assert regions["framework"]["중간"] == 1     # framework DG
    assert regions["framework"]["낮음"] == 3     # Met + 2 conserved Trp
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
# Sequences with failed annotation / input errors
# --------------------------------------------------------------------------

def test_unannotated_sequence_is_marked_and_not_attributed_to_cdr(records):
    """Motifs from sequences that cannot be split into regions are not placed in CDR cells."""
    summary = summarize(_reports(records, "SYN-VH-901"))
    row = _row(summary, "SYN-VH-901")

    assert row["annotated"] is False
    assert row["regions_unassigned"] is True
    assert row["confident"] is False
    for cdr in ["CDR1", "CDR2", "CDR3"]:
        assert sum(row["regions"][cdr].values()) == 0
    assert sum(row["regions"]["framework"].values()) == 1


def test_length_outlier_is_counted_normally(records):
    """Sequences with intact anchors are summarized per CDR even when confident=False."""
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
# Does it appear in the output?
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


def test_text_summary_tables_are_aligned(records):
    """Tables are split by chain, so check alignment for each table separately."""
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    section = out.split("CDR별 liability 집계")[1].split("[SYN-VH-001]")[0]

    tables = section.split("■")[1:]
    assert len(tables) == 3, "중쇄 / 경쇄 / 판별 실패 세 표가 나와야 한다"

    total_rows = 0
    for table in tables:
        rows = [ln for ln in table.splitlines() if ln.strip().startswith("SYN-")]
        assert rows
        total_rows += len(rows)
        # Korean characters are 2 columns wide, so compare by display width, not string index
        positions = {
            tuple(_display_width(row[:i]) for i, ch in enumerate(row) if ch == "|")
            for row in rows
        }
        assert len(positions) == 1, f"집계 표 열 정렬이 틀어졌다: {positions}"

    assert total_rows == len(records), "모든 서열이 어느 한 표에는 들어가야 한다"


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
    """The summary only reports counts -- no scores, ranks or reordering."""
    records_order = [r["id"] for r in records]
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    assert [r["id"] for r in summary["rows"]] == records_order
    assert not any(k in summary for k in ["score", "rank", "ranking", "점수", "순위"])


# --------------------------------------------------------------------------
# Separate heavy / light chain tables
# --------------------------------------------------------------------------

def test_groups_split_heavy_light_and_unknown(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    groups = {g["chain"]: g for g in summary["groups"]}

    assert set(groups) == {"H", "L", None}
    assert [g["chain"] for g in summary["groups"]] == ["H", "L", None], "중쇄 → 경쇄 → 실패 순서"

    assert len(groups["H"]["rows"]) == 16
    assert len(groups["L"]["rows"]) == 7
    assert [r["id"] for r in groups[None]["rows"]] == ["SYN-VH-901"]


def test_every_sequence_lands_in_exactly_one_group(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    seen = [r["id"] for g in summary["groups"] for r in g["rows"]]
    assert sorted(seen) == sorted(r["id"] for r in records)
    assert len(seen) == len(set(seen)), "같은 서열이 두 표에 들어갔다"


def test_group_columns_use_chain_specific_cdr_names(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    groups = {g["chain"]: g for g in summary["groups"]}

    assert groups["H"]["cdr_labels"] == ["CDR-H1", "CDR-H2", "CDR-H3"]
    assert groups["L"]["cdr_labels"] == ["CDR-L1", "CDR-L2", "CDR-L3"]


def test_unannotated_but_chain_known_stays_in_its_chain_group(records):
    """SYN-VH-902 is still detected as chain H, so it stays in the heavy chain table."""
    summary = summarize(_reports(records, "SYN-VH-902", "SYN-VH-901"))
    groups = {g["chain"]: g for g in summary["groups"]}

    assert [r["id"] for r in groups["H"]["rows"]] == ["SYN-VH-902"]
    assert [r["id"] for r in groups[None]["rows"]] == ["SYN-VH-901"]


def test_invalid_input_goes_to_unknown_group():
    summary = summarize([build("bad", "ACDEF")])
    assert [g["chain"] for g in summary["groups"]] == [None]
    assert summary["groups"][0]["rows"][0]["input_ok"] is False


def test_empty_groups_are_omitted(records):
    """With only light chains as input, no heavy chain table is built."""
    summary = summarize(_reports(records, "SYN-VL-001", "SYN-VL-003"))
    assert [g["chain"] for g in summary["groups"]] == ["L"]


def test_group_totals_match_their_own_rows(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    for group in summary["groups"]:
        for region in SUMMARY_REGIONS:
            for level in summary["risk_levels"]:
                expected = sum(r["regions"][region][level] for r in group["rows"])
                assert group["totals"][region][level] == expected, (
                    group["chain"], region, level
                )


def test_group_totals_add_up_to_overall_totals(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    for region in SUMMARY_REGIONS:
        for level in summary["risk_levels"]:
            per_group = sum(g["totals"][region][level] for g in summary["groups"])
            assert per_group == summary["totals"][region][level], (region, level)

    assert sum(g["domain_flag_total"] for g in summary["groups"]) == summary["domain_flag_total"]


def test_html_renders_one_table_per_chain(records):
    out = as_html(build_many([(r["id"], r["sequence"]) for r in records]))

    assert "중쇄 (VH)" in out
    assert "경쇄 (VL)" in out
    assert "체인 판별 실패" in out
    # The heavy table uses CDR-H columns, the light table CDR-L columns
    assert "<th>CDR-H2</th>" in out
    assert "<th>CDR-L2</th>" in out
    assert "<th>CDR2</th>" not in out


def test_text_renders_one_table_per_chain(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    section = out.split("CDR별 liability 집계")[1].split("[SYN-VH-001]")[0]

    assert "중쇄 (VH) — 16건" in section
    assert "경쇄 (VL) — 7건" in section
    assert "CDR-H2" in section and "CDR-L2" in section


def test_unknown_group_has_no_cdr_columns(records):
    """If the chain is unknown, no CDR cells exist at all -- writing 0 would be misread as 'clean'."""
    out = as_text(_reports(records, "SYN-VH-901"))
    section = out.split("CDR별 liability 집계")[1].split("[SYN-VH-901]")[0]

    assert "검출 모티프" in section
    assert "CDR-H1" not in section
    assert "CDR1" not in section
