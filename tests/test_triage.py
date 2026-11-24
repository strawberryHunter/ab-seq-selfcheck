"""선별 규칙 테스트.

핵심은 세 가지다.
  - CDR 안의 '높음'이 '검토 필요'로 잡히는지
  - framework 의 '높음'이나 CDR 의 '중간'은 잡히지 **않는지** (규칙이 헐거워지지 않게)
  - 주석이 실패한 서열을 '기준 해당 없음'으로 넘기지 않는지 (못 본 것을 괜찮다고 적으면 안 된다)
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
# 규칙이 잡아야 하는 것
# --------------------------------------------------------------------------

def test_high_risk_in_cdr_is_flagged_for_review(records):
    """CDR 안의 '높음'은 CDR1/2/3 어디에 있든 검토 필요."""
    for seq_id, region in [
        ("SYN-VH-004", "CDR-H1"),     # 당화
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
    """사유만 보고도 어디를 봐야 할지 알 수 있어야 한다."""
    report = _report(records, "SYN-VH-004")
    reason = report.triage.reasons[0]

    assert "CDR-H1" in reason
    assert "32번" in reason
    assert "NLS" in reason
    assert "N-결합 당화" in reason


def test_unpaired_cysteine_at_high_risk_is_flagged(records):
    """비결합 Cys 는 위치 하나로 특정할 수 없지만 CDR 안이면 '높음'이라 같이 잡는다."""
    report = _report(records, "SYN-VH-005")
    assert report.triage.status == STATUS_REVIEW
    assert any("비결합 Cys" in r for r in report.triage.reasons)


# --------------------------------------------------------------------------
# 규칙이 잡으면 안 되는 것 — 규칙이 헐거워지는지 확인
# --------------------------------------------------------------------------

def test_clean_sequences_are_not_flagged(records):
    for seq_id in ["SYN-VH-001", "SYN-VH-002", "SYN-VL-001", "SYN-VL-002"]:
        report = _report(records, seq_id)
        assert report.triage.status == STATUS_CLEAR, (seq_id, report.triage.reasons)
        assert report.triage.reasons == []


def test_framework_high_risk_is_not_flagged(records):
    """framework 의 DG 는 등급 조정으로 '중간'이 되므로 걸리지 않는다."""
    report = _report(records, "SYN-VH-006")
    assert report.triage.status == STATUS_CLEAR
    # framework 에 DG 가 분명히 있는데도 걸리지 않아야 한다
    assert any(h.match == "DG" and not h.in_cdr for h in report.liabilities.hits)


def test_medium_risk_in_cdr_is_not_flagged(records):
    """CDR 안이어도 '중간'이면 이 규칙에는 걸리지 않는다."""
    report = _report(records, "SYN-VH-008")      # CDR-H3 의 Met/Trp → 중간
    assert report.triage.status == STATUS_CLEAR
    assert any(h.in_cdr and h.risk == "중간" for h in report.liabilities.hits)


def test_medium_domain_flag_is_not_flagged(records):
    """도메인 플래그가 '중간'이면 걸리지 않는다."""
    report = _report(records, "SYN-VH-902")
    assert report.liabilities.domain_flags[0]["risk"] == "중간"
    # 다만 이 서열은 주석 실패라 애초에 판정 불가다
    assert report.triage.status == STATUS_UNKNOWN


# --------------------------------------------------------------------------
# 판정할 수 없는 경우 — 가장 중요한 부분
# --------------------------------------------------------------------------

def test_failed_annotation_is_unknown_not_clear(records):
    """CDR 을 못 가르면 '기준 해당 없음'이 아니라 '판정 불가'다."""
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


def test_length_outlier_is_evaluated_but_carries_a_warning(records):
    """앵커가 멀쩡하면 규칙은 적용하되, 신뢰도가 낮다는 경고를 남긴다."""
    report = _report(records, "SYN-VH-903")
    assert report.triage.status == STATUS_CLEAR
    assert report.triage.evaluable is True
    assert "ANARCI" in report.triage.note


# --------------------------------------------------------------------------
# 규칙 자체
# --------------------------------------------------------------------------

def test_rules_are_described_in_korean():
    descriptions = rule_descriptions()
    assert len(descriptions) == len(DEFAULT_RULES)
    for desc in descriptions:
        assert any("가" <= ch <= "힣" for ch in desc)
        assert "검토 필요" in desc


def test_evaluate_with_empty_rule_set_never_flags(records):
    """규칙을 비우면 아무것도 걸리지 않는다 — 규칙이 결과를 만든다는 확인."""
    report = _report(records, "SYN-VH-004")
    assert report.triage.status == STATUS_REVIEW

    result = evaluate(report, rules=())
    assert result.status == STATUS_CLEAR
    assert result.reasons == []


# --------------------------------------------------------------------------
# 출력
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
    """'검토 필요'가 판정이 아니라는 고지가 빠지지 않는지."""
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
