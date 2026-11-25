"""Triage rule tests.

Three things matter:
  - a 'high' inside a CDR is flagged 'needs review'
  - a 'high' in the framework or a 'medium' in a CDR is **not** flagged (so the rule does not get loose)
  - a sequence whose annotation failed is not passed as 'criteria not met' (what wasn't seen must not be reported as fine)
"""

from __future__ import annotations

import json
import re

from abcheck.render import as_html, as_json, as_text
from abcheck.report import build, build_many, summarize
from abcheck.triage import (
    DEFAULT_RULES,
    STATUS_CLEAR,
    STATUS_REVIEW,
    STATUS_UNKNOWN,
    TRIAGE_NOTICE,
    evaluate,
    rule_descriptions,
)

from .conftest import by_id


def _report(records, seq_id):
    rec = by_id(records, seq_id)
    return build(rec["id"], rec["sequence"])


# --------------------------------------------------------------------------
# What the rule must catch
# --------------------------------------------------------------------------

def test_high_risk_in_cdr_is_flagged_for_review(records):
    """A 'high' inside a CDR needs review, wherever it is in CDR1/2/3."""
    for seq_id, region in [
        ("SYN-VH-004", "CDR-H1"),     # Glycosylation
        ("SYN-VH-003", "CDR-H2"),     # NG
        ("SYN-VH-007", "CDR-H3"),     # DG
        ("SYN-VL-004", "CDR-L1"),     # DP
        ("SYN-VL-003", "CDR-L3"),     # NG
    ]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_REVIEW, seq_id
        assert report.triage.reasons, seq_id
        assert region in report.triage.reasons[0], (seq_id, report.triage.reasons)
        assert "높음" in report.triage.reasons[0]


def test_reason_names_position_and_motif(records):
    """The reason alone must tell you where to look."""
    report = _report(records, "SYN-VH-004")
    reason = report.triage.reasons[0]

    assert "CDR-H1" in reason
    assert "32번" in reason
    assert "NLS" in reason
    assert "N-결합 당화" in reason


def test_unpaired_cysteine_at_high_risk_is_flagged(records):
    """An unpaired Cys cannot be pinned to one position, but inside a CDR it is 'high', so it is caught too."""
    report = _report(records, "SYN-VH-005")
    assert report.triage.status == STATUS_REVIEW
    assert any("비결합 Cys" in r for r in report.triage.reasons)


# --------------------------------------------------------------------------
# What the rule must not catch -- check the rule does not get loose
# --------------------------------------------------------------------------

def test_clean_sequences_are_not_flagged(records):
    for seq_id in ["SYN-VH-001", "SYN-VH-002", "SYN-VL-001", "SYN-VL-002"]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_CLEAR, (seq_id, report.triage.reasons)
        assert report.triage.reasons == []


def test_framework_high_risk_is_not_flagged(records):
    """DG in the framework is adjusted down to 'medium', so it is not flagged."""
    report = _report(records, "SYN-VH-006")
    assert report.triage.status == STATUS_CLEAR
    # DG is clearly present in the framework, yet it must not be flagged
    assert any(h.match == "DG" and not h.in_cdr for h in report.liabilities.hits)


def test_weak_class_high_risk_in_cdr_is_not_flagged(records):
    """Even if **displayed** as 'high' inside a CDR, a weak-tier motif is not flagged by triage.

    The CDR upgrade makes NS/NT/NN/NH and DS/DT/DD 'high' too, but in real repertoires these
    are the majority, so flagging them all would give no discrimination. Triage therefore only looks
    at the top-risk motif of each category (N-linked glycosylation, NG, DG, DP). The display is unchanged, so no information is lost.
    """
    report = _report(records, "SYN-VH-010")          # NT in CDR-H2

    cdr_high = [h for h in report.liabilities.cdr_hits if h.risk == "높음"]
    assert cdr_high, "CDR 안에 '높음'으로 표시된 항목이 있어야 한다"
    assert all(not h.top_within_kind for h in cdr_high)

    assert report.triage.status == STATUS_CLEAR, "약한 등급은 선별에 걸리면 안 된다"
    assert report.triage.reasons == []


def test_top_class_motifs_are_the_ones_that_trigger(records):
    """Every item that triggered triage must be a top-risk motif."""
    for seq_id in ["SYN-VH-003", "SYN-VH-004", "SYN-VH-007", "SYN-VL-003", "SYN-VL-004"]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_REVIEW, seq_id
        triggering = [
            h for h in report.liabilities.cdr_hits
            if h.risk == "높음" and h.top_within_kind
        ]
        assert triggering, seq_id
        assert len(report.triage.reasons) >= len(triggering)


def test_medium_risk_in_cdr_is_not_flagged(records):
    """Even inside a CDR, a 'medium' does not trip this rule."""
    report = _report(records, "SYN-VH-008")      # Met/Trp in CDR-H3 -> medium
    assert report.triage.status == STATUS_CLEAR
    assert any(h.in_cdr and h.risk == "중간" for h in report.liabilities.hits)


def test_medium_domain_flag_is_not_flagged(records):
    """A 'medium' domain flag is not flagged."""
    report = _report(records, "SYN-VH-902")
    assert report.liabilities.domain_flags[0]["risk"] == "중간"
    # Though this sequence failed annotation, so it cannot be assessed in the first place
    assert report.triage.status == STATUS_UNKNOWN


# --------------------------------------------------------------------------
# Cases that cannot be assessed -- the most important part
# --------------------------------------------------------------------------

def test_failed_annotation_is_unknown_not_clear(records):
    """If CDRs cannot be split, the result is 'cannot assess', not 'criteria not met'."""
    for seq_id in ["SYN-VH-901", "SYN-VH-902"]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_UNKNOWN, seq_id
        assert report.triage.status != STATUS_CLEAR, "못 본 것을 괜찮다고 적으면 안 된다"
        assert report.triage.evaluable is False
        assert "사람이 직접 확인" in report.triage.note


def test_invalid_input_is_unknown():
    report = build("bad", "ACDEF")
    assert report.triage.status == STATUS_UNKNOWN
    assert report.triage.evaluable is False
    assert "검증을 통과하지 못해" in report.triage.note


def test_unconfident_sequence_is_never_reported_as_no_rule_matched(records):
    """Rules still run on an unconfident annotation, but "no match" cannot be claimed.

    SYN-VH-903 (CDR-H1 too long) and SYN-VH-013/014 (FR3 back-calculated) used to show
    "기준 해당 없음" in the summary table, with the warning only in a note that the table
    does not show.
    """
    for seq_id in ["SYN-VH-903", "SYN-VH-013", "SYN-VH-014"]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_UNKNOWN, seq_id
        assert report.triage.evaluable is True          # the rules did run
        assert "ANARCI" in report.triage.note


def test_unconfident_sequence_with_a_cdr_hit_still_needs_review(records):
    """Uncertain boundaries must not hide a hit the rules did find."""
    import dataclasses

    report = _report(records, "SYN-VH-003")                 # NG in CDR-H2
    report.annotation = dataclasses.replace(report.annotation, confident=False)
    from abcheck.triage import evaluate

    assert evaluate(report).status == STATUS_REVIEW


# --------------------------------------------------------------------------
# The rule itself
# --------------------------------------------------------------------------

def test_rules_are_described_in_korean():
    descriptions = rule_descriptions()
    assert len(descriptions) == len(DEFAULT_RULES)
    for desc in descriptions:
        assert any("가" <= ch <= "힣" for ch in desc)
        assert "검토 필요" in desc


def test_evaluate_with_empty_rule_set_never_flags(records):
    """An empty rule flags nothing -- confirms the rule is what produces the result."""
    report = _report(records, "SYN-VH-004")
    assert report.triage.status == STATUS_REVIEW

    result = evaluate(report, rules=())
    assert result.status == STATUS_CLEAR
    assert result.reasons == []


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def test_triage_appears_in_all_three_outputs(records):
    reports = build_many([(r["id"], r["sequence"]) for r in records])

    text = as_text(reports)
    assert "선별 분류: 검토 필요" in text
    assert "적용한 선별 규칙" in text

    html = as_html(reports)
    assert "선별 분류" in html
    assert "적용한 선별 규칙" in html

    payload = json.loads(as_json(reports))
    entry = next(e for e in payload["sequences"] if e["id"] == "SYN-VH-004")
    assert entry["triage"]["status"] == STATUS_REVIEW
    assert entry["triage"]["reasons"]


def test_summary_table_has_triage_column(records):
    out = as_text(build_many([(r["id"], r["sequence"]) for r in records]))
    section = out.split("CDR별 liability 집계")[1].split("[SYN-VH-001]")[0]

    assert "| 분류 " in section or "분류" in section
    assert "검토 필요" in section
    assert "기준 해당 없음" in section
    assert "판정 불가" in section


def test_summary_counts_reviews(records):
    summary = summarize(build_many([(r["id"], r["sequence"]) for r in records]))
    expected = sum(
        1 for r in summary["rows"] if r["triage"]["status"] == STATUS_REVIEW
    )
    assert summary["review_count"] == expected
    assert summary["review_count"] == 6


def test_every_output_carries_the_triage_notice(records):
    """The notice that 'needs review' is not a verdict must not be dropped."""
    reports = build_many([(r["id"], r["sequence"]) for r in records])

    flat_text = re.sub(r"\s+", " ", as_text(reports))
    assert re.sub(r"\s+", " ", TRIAGE_NOTICE) in flat_text

    import html as html_module
    assert html_module.escape(TRIAGE_NOTICE) in as_html(reports)

    payload = json.loads(as_json(reports))
    assert payload["summary"]["triage_notice"] == TRIAGE_NOTICE
    assert payload["sequences"][0]["triage"]["notice"] == TRIAGE_NOTICE


def test_notice_states_it_is_not_a_verdict():
    assert "나쁘다는 뜻이 아니다" in TRIAGE_NOTICE
    assert "안전하다는 뜻이 아니라" in TRIAGE_NOTICE
    assert "연구자가 정한 것" in TRIAGE_NOTICE
