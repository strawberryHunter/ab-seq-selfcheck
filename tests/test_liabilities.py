"""liability 검출 테스트 — CDR 과 framework 의 위험도 구분이 핵심이다."""

from __future__ import annotations

from abcheck.annotate import annotate
from abcheck.liabilities import MOTIF_SPECS, RISK_LEVELS, find_liabilities

from .conftest import by_id


def _scan(records, seq_id):
    rec = by_id(records, seq_id)
    ann = annotate(rec["sequence"])
    return rec, ann, find_liabilities(ann)


def _kinds(hits) -> set[str]:
    return {h.kind for h in hits}


# --------------------------------------------------------------------------
# 구간별 위험도 구분 — 이 모듈의 존재 이유
# --------------------------------------------------------------------------

def test_same_motif_is_riskier_inside_cdr(records):
    """같은 DG 모티프가 CDR 안이면 '높음', framework 면 '중간' 이어야 한다."""
    _, _, fr_report = _scan(records, "SYN-VH-006")       # framework 에 DG
    _, _, cdr_report = _scan(records, "SYN-VH-007")      # CDR-H3 에 DG

    fr_dg = [h for h in fr_report.hits if h.match == "DG"]
    cdr_dg = [h for h in cdr_report.hits if h.match == "DG"]

    assert len(fr_dg) == 1 and len(cdr_dg) == 1
    assert fr_dg[0].in_cdr is False
    assert cdr_dg[0].in_cdr is True

    assert fr_dg[0].base_risk == cdr_dg[0].base_risk == "높음"
    assert fr_dg[0].risk == "중간", "framework DG 는 한 단계 내려야 한다"
    assert cdr_dg[0].risk == "높음", "CDR DG 는 최고 등급이어야 한다"

    # 조정 근거가 한글로 남는다
    assert "framework" in fr_dg[0].note
    assert "CDR-H3" in cdr_dg[0].note


def test_framework_met_is_lowest_but_cdr_met_is_raised(records):
    _, _, fr_report = _scan(records, "SYN-VH-006")       # framework 에 Met
    _, _, cdr_report = _scan(records, "SYN-VH-008")      # CDR-H3 에 Met

    fr_met = [h for h in fr_report.hits if h.match == "M"]
    cdr_met = [h for h in cdr_report.hits if h.match == "M" and h.in_cdr]

    assert fr_met and cdr_met
    assert fr_met[0].risk == "낮음"
    assert cdr_met[0].risk == "중간"


def test_risk_never_leaves_the_three_levels(records):
    for rec in records:
        ann = annotate(rec["sequence"])
        report = find_liabilities(ann)
        for hit in report.hits:
            assert hit.risk in RISK_LEVELS
            assert hit.base_risk in RISK_LEVELS


# --------------------------------------------------------------------------
# 정상 케이스 — CDR liability 가 없어야 한다
# --------------------------------------------------------------------------

def test_normal_cases_have_no_cdr_liabilities(records):
    """주입한 liability 가 없는 케이스는 CDR 안이 깨끗하다."""
    for seq_id in ["SYN-VH-001", "SYN-VH-002", "SYN-VL-001", "SYN-VL-002"]:
        _, _, report = _scan(records, seq_id)
        assert report.cdr_hits == [], f"{seq_id}: CDR liability 가 있다 — {report.cdr_hits}"
        assert report.domain_flags == [], f"{seq_id}: 도메인 플래그가 있다"


def test_normal_cases_only_flag_conserved_trp_in_framework(records):
    """정상 케이스의 framework 플래그는 보존 Trp 앵커 산화뿐이다.

    framework filler 알파벳에서 C/W/M/N/D 를 뺐으므로, 남는 것은 앵커 Trp 밖에 없다.
    중쇄는 Trp36 + WGxG 의 Trp 로 2건, 경쇄는 Trp35 로 1건.
    """
    for seq_id, expected in [("SYN-VH-001", 2), ("SYN-VH-002", 2),
                             ("SYN-VL-001", 1), ("SYN-VL-002", 1)]:
        _, ann, report = _scan(records, seq_id)
        assert _kinds(report.framework_hits) == {"oxidation"}, seq_id
        assert len(report.framework_hits) == expected, (
            f"{seq_id}: framework 플래그 {len(report.framework_hits)}건 "
            f"(기대 {expected}) — {[(h.match, h.start) for h in report.framework_hits]}"
        )
        for hit in report.framework_hits:
            assert hit.risk == "낮음"
        # 앵커 위치와 일치하는지
        anchor_positions = {ann.anchors["trp"]}
        if ann.chain == "H":
            anchor_positions.add(ann.anchors["terminal_motif"])
        assert {h.start for h in report.framework_hits} == anchor_positions


# --------------------------------------------------------------------------
# 주입한 liability 를 정확히 검출하는지
# --------------------------------------------------------------------------

def test_cdr_deamidation_ng_detected(records):
    for seq_id, cdr in [("SYN-VH-003", "CDR-H2"), ("SYN-VL-003", "CDR-L3")]:
        _, _, report = _scan(records, seq_id)
        ng = [h for h in report.cdr_hits if h.match == "NG"]
        assert len(ng) == 1, f"{seq_id}: NG {len(ng)}건"
        assert ng[0].kind == "deamidation"
        assert ng[0].region == cdr
        assert ng[0].risk == "높음"
        assert "deamidation" in ng[0].rationale


def test_cdr_glycosylation_detected(records):
    _, _, report = _scan(records, "SYN-VH-004")
    glyc = [h for h in report.cdr_hits if h.kind == "glycosylation"]
    assert len(glyc) == 1
    assert glyc[0].match == "NLS"
    assert glyc[0].region == "CDR-H1"
    assert glyc[0].risk == "높음"
    assert "불균일성" in glyc[0].rationale


def test_cdr_fragmentation_dp_detected(records):
    _, _, report = _scan(records, "SYN-VL-004")
    dp = [h for h in report.cdr_hits if h.match == "DP"]
    assert len(dp) == 1
    assert dp[0].kind == "fragmentation"
    assert dp[0].region == "CDR-L1"
    assert dp[0].risk == "높음", "framework 기본 중간 → CDR 에서 한 단계 상향"


def test_cdr_oxidation_met_and_trp_detected(records):
    _, _, report = _scan(records, "SYN-VH-008")
    ox = [h for h in report.cdr_hits if h.kind == "oxidation"]
    assert {h.match for h in ox} == {"M", "W"}
    assert all(h.region == "CDR-H3" for h in ox)
    assert all(h.risk == "중간" for h in ox)


# --------------------------------------------------------------------------
# 비결합 Cys — 도메인 단위 플래그
# --------------------------------------------------------------------------

def test_unpaired_cysteine_flagged_when_count_is_odd(records):
    rec, _, report = _scan(records, "SYN-VH-005")

    assert rec["sequence"].count("C") == 3
    assert len(report.domain_flags) == 1

    flag = report.domain_flags[0]
    assert flag["kind"] == "unpaired_cysteine"
    assert flag["cysteine_count"] == 3
    assert flag["risk"] == "높음", "CDR 안에 여분 Cys 가 있으면 위험도가 더 높다"
    assert "CDR-H3" in flag["regions"]
    assert "홀수" in flag["rationale"]
    assert "응집" in flag["rationale"]


def test_even_cysteine_count_is_not_flagged(records):
    for rec in records:
        if rec["sequence"].count("C") % 2 == 0:
            ann = annotate(rec["sequence"])
            report = find_liabilities(ann)
            assert report.domain_flags == [], rec["id"]


def test_cdr_cysteine_also_appears_as_regular_hit(records):
    """CDR 안의 여분 Cys 는 도메인 플래그로 보고하고, 위치 기반 모티프로는 잡지 않는다."""
    _, _, report = _scan(records, "SYN-VH-005")
    assert "cysteine" not in _kinds(report.hits)
    assert report.domain_flags[0]["positions"] == [22, 95, 103]


# --------------------------------------------------------------------------
# 주석이 실패한 서열에서도 동작해야 한다
# --------------------------------------------------------------------------

def test_works_without_annotation_and_warns(records):
    _, ann, report = _scan(records, "SYN-VH-901")

    assert ann.cdrs == {}
    assert report.hits, "주석 실패와 무관하게 모티프는 검출해야 한다"
    assert all(h.region == "구간 미확정" for h in report.hits)
    assert all(h.in_cdr is False for h in report.hits)
    # 등급 조정을 하지 않았으므로 기본 등급 그대로다
    assert all(h.risk == h.base_risk for h in report.hits)

    assert report.warnings
    assert "구간을 구분할 수 없다" in " ".join(report.warnings)
    assert "과소평가하지 말 것" in report.hits[0].note


def test_odd_cysteine_without_annotation_still_flagged(records):
    rec, ann, report = _scan(records, "SYN-VH-902")

    assert rec["sequence"].count("C") == 1
    assert ann.cdrs == {}
    assert len(report.domain_flags) == 1
    flag = report.domain_flags[0]
    assert flag["risk"] == "중간", "구간을 모르면 CDR 상향을 적용하지 않는다"
    assert flag["regions"] == ["구간 미확정"]


# --------------------------------------------------------------------------
# 겹치는 모티프 처리
# --------------------------------------------------------------------------

def test_higher_risk_motif_wins_within_same_kind():
    """NG 를 잡은 자리를 같은 분류의 약한 패턴으로 다시 잡지 않는다."""
    from abcheck.annotate import Annotation

    ann = Annotation(sequence="AANGSAA")
    report = find_liabilities(ann)

    deamid = [h for h in report.hits if h.kind == "deamidation"]
    assert len(deamid) == 1
    assert deamid[0].match == "NG"


def test_glycosylation_and_deamidation_both_reported_for_ngs():
    """NGS 는 당화 수용 서열이면서 NG deamidation 이기도 하다 — 분류가 다르므로 둘 다 본다."""
    from abcheck.annotate import Annotation

    ann = Annotation(sequence="AANGSAA")
    report = find_liabilities(ann)

    assert _kinds(report.hits) >= {"glycosylation", "deamidation"}


def test_motif_on_cdr_boundary_counts_as_cdr(records):
    """경계에 걸친 모티프는 위험도를 과소평가하지 않도록 CDR 쪽으로 본다."""
    from abcheck.annotate import CDR, Annotation

    seq = "A" * 10 + "NG" + "A" * 10
    ann = Annotation(sequence=seq, chain="H",
                     cdrs={"CDR-H3": CDR("CDR-H3", 12, 20, seq[11:20])})
    report = find_liabilities(ann)

    ng = [h for h in report.hits if h.match == "NG"]
    assert ng[0].start == 11 and ng[0].end == 12     # 11번은 framework, 12번은 CDR
    assert ng[0].in_cdr is True
    assert ng[0].region == "CDR-H3"


# --------------------------------------------------------------------------
# 집계
# --------------------------------------------------------------------------

def test_counts_by_risk_is_ordered_high_first(records):
    _, _, report = _scan(records, "SYN-VH-007")
    counts = report.counts_by_risk()
    assert list(counts) == ["높음", "중간", "낮음"]
    assert sum(counts.values()) == len(report.hits)


def test_every_spec_has_korean_rationale():
    """근거는 반드시 한글이어야 한다.

    라벨은 'deamidation (NG)' 처럼 영문 용어를 그대로 쓴다 — 연구실에서 통용되는
    표현이라 억지로 번역하면 오히려 알아보기 어렵다.
    """
    for spec in MOTIF_SPECS:
        assert any("가" <= ch <= "힣" for ch in spec.rationale), spec.kind
        assert spec.label.strip(), spec.kind
        assert spec.base_risk in RISK_LEVELS, spec.kind
