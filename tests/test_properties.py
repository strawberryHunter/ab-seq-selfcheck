"""물성 계산 테스트 — 손으로 검증할 수 있는 작은 케이스로 확인한다."""

from __future__ import annotations

import math

from abcheck.annotate import annotate
from abcheck.properties import (
    DEFAULT_PH,
    KYTE_DOOLITTLE,
    PKA_C_TERMINUS,
    PKA_N_TERMINUS,
    SIDE_CHAIN_PKA,
    compute,
    compute_all,
    gravy,
    isoelectric_point,
    net_charge,
)

from .conftest import by_id


# --------------------------------------------------------------------------
# net charge — Henderson-Hasselbalch 를 직접 검산한다
# --------------------------------------------------------------------------

def test_net_charge_of_single_alanine_is_termini_only():
    """Ala 는 이온화 측쇄가 없다 — 말단 두 그룹만 남는다."""
    expected = (
        1.0 / (1.0 + 10 ** (DEFAULT_PH - PKA_N_TERMINUS))
        - 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - DEFAULT_PH))
    )
    assert math.isclose(net_charge("A"), expected, rel_tol=1e-12)


def test_net_charge_of_lysine_adds_one_positive_group():
    base = net_charge("A")
    lys_contribution = 1.0 / (1.0 + 10 ** (DEFAULT_PH - SIDE_CHAIN_PKA["K"]))
    assert math.isclose(net_charge("AK"), base + lys_contribution, rel_tol=1e-12)


def test_net_charge_of_aspartate_adds_one_negative_group():
    base = net_charge("A")
    asp_contribution = -1.0 / (1.0 + 10 ** (SIDE_CHAIN_PKA["D"] - DEFAULT_PH))
    assert math.isclose(net_charge("AD"), base + asp_contribution, rel_tol=1e-12)


def test_net_charge_decreases_monotonically_with_ph():
    seq = "KRHDECYAAA"
    charges = [net_charge(seq, ph) for ph in [2, 4, 6, 8, 10, 12, 14]]
    assert charges == sorted(charges, reverse=True)


def test_poly_lysine_is_positive_and_poly_aspartate_negative():
    assert net_charge("K" * 10) > 5
    assert net_charge("D" * 10) < -5


def test_termini_can_be_excluded():
    assert net_charge("AAAA", include_termini=False) == 0.0
    assert net_charge("AAAA", include_termini=True) != 0.0


# --------------------------------------------------------------------------
# pI — net charge 가 0 이 되는 지점인지 되짚어 확인한다
# --------------------------------------------------------------------------

def test_pi_is_the_ph_where_net_charge_is_zero():
    for seq in ["AKDEAA", "KRKRKR", "DDEEDD", "AAHAAYAA", "CYHKRDE"]:
        pi = isoelectric_point(seq)
        assert pi is not None
        assert abs(net_charge(seq, pi)) < 1e-4, f"{seq}: pI={pi} 에서 전하가 0 이 아니다"


def test_poly_lysine_has_basic_pi_and_poly_aspartate_acidic():
    assert isoelectric_point("K" * 10) > 10
    assert isoelectric_point("D" * 10) < 4


def test_neutral_peptide_pi_sits_between_terminal_pkas():
    """이온화 측쇄가 없으면 pI 는 두 말단 pKa 사이에 있다."""
    pi = isoelectric_point("A" * 10)
    assert PKA_C_TERMINUS < pi < PKA_N_TERMINUS


def test_pi_without_any_ionizable_group_is_none():
    assert isoelectric_point("AAAA", include_termini=False) is None
    assert isoelectric_point("") is None


def test_pi_is_within_search_range():
    for seq in ["K" * 50, "D" * 50, "R" * 50]:
        pi = isoelectric_point(seq)
        assert 0.0 <= pi <= 14.0


# --------------------------------------------------------------------------
# GRAVY
# --------------------------------------------------------------------------

def test_gravy_is_arithmetic_mean_of_kyte_doolittle():
    seq = "AILVF"
    expected = sum(KYTE_DOOLITTLE[r] for r in seq) / len(seq)
    assert math.isclose(gravy(seq), expected, rel_tol=1e-12)


def test_gravy_sign_matches_hydropathy():
    assert gravy("IIIVVV") > 0, "분지 소수성 잔기는 양수여야 한다"
    assert gravy("RRRKKK") < 0, "하전 잔기는 음수여야 한다"


def test_gravy_of_single_residues_matches_table():
    for residue, value in KYTE_DOOLITTLE.items():
        assert math.isclose(gravy(residue), value, rel_tol=1e-12)


def test_gravy_of_empty_is_none():
    assert gravy("") is None


# --------------------------------------------------------------------------
# 구간별 계산
# --------------------------------------------------------------------------

def test_compute_all_covers_domain_and_every_cdr(records):
    rec = by_id(records, "SYN-VH-001")
    ann = annotate(rec["sequence"])
    props = compute_all(ann.sequence, ann.cdrs)

    assert "도메인 전체" in props
    assert set(props) == {"도메인 전체"} | set(ann.cdrs)

    domain = props["도메인 전체"]
    assert domain.length == len(rec["sequence"])
    for name, cdr in ann.cdrs.items():
        assert props[name].sequence == cdr.seq
        assert props[name].length == cdr.length


def test_compute_all_is_ordered_by_position(records):
    ann = annotate(by_id(records, "SYN-VH-001")["sequence"])
    props = compute_all(ann.sequence, ann.cdrs)
    names = [n for n in props if n != "도메인 전체"]
    assert names == ["CDR-H1", "CDR-H2", "CDR-H3"]


def test_compute_all_without_cdrs_returns_domain_only(records):
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    props = compute_all(ann.sequence, ann.cdrs)
    assert list(props) == ["도메인 전체"]


def test_basic_cdr_has_higher_pi_than_other_cdrs(records):
    """K/R 을 몰아넣은 CDR 은 같은 서열의 다른 CDR 보다 pI 가 높다."""
    rec = by_id(records, "SYN-VL-005")
    ann = annotate(rec["sequence"])
    props = compute_all(ann.sequence, ann.cdrs)

    loaded = props["CDR-L1"].pi          # DYKRKRKGRGG
    assert loaded is not None
    assert loaded > 10, f"염기성 CDR-L1 의 pI 가 {loaded}"
    assert props["CDR-L1"].net_charge_at_ph > 2
    assert loaded > props["CDR-L2"].pi


def test_basic_sequence_has_higher_pi_than_neutral_one(records):
    """K/R 집중 서열은 도메인 전체 pI 도 정상 VL 보다 높다."""
    loaded = annotate(by_id(records, "SYN-VL-005")["sequence"])
    normal = annotate(by_id(records, "SYN-VL-001")["sequence"])

    pi_loaded = compute("loaded", loaded.sequence).pi
    pi_normal = compute("normal", normal.sequence).pi
    assert pi_loaded > pi_normal


# --------------------------------------------------------------------------
# 전체 코퍼스에 대한 건전성 검사
# --------------------------------------------------------------------------

def test_every_corpus_sequence_gets_finite_properties(records):
    for rec in records:
        p = compute(rec["id"], rec["sequence"])
        assert p.pi is not None and 0.0 <= p.pi <= 14.0, rec["id"]
        assert p.gravy is not None and -5.0 < p.gravy < 5.0, rec["id"]
        assert math.isfinite(p.net_charge_at_ph), rec["id"]


def test_net_charge_sign_is_consistent_with_pi(records):
    """pH 7.4 가 pI 보다 낮으면 양전하, 높으면 음전하여야 한다."""
    for rec in records:
        p = compute(rec["id"], rec["sequence"])
        if abs(p.pi - DEFAULT_PH) < 0.2:
            continue                      # pI 근처는 부호가 흔들린다
        if p.pi > DEFAULT_PH:
            assert p.net_charge_at_ph > 0, f"{rec['id']}: pI={p.pi}"
        else:
            assert p.net_charge_at_ph < 0, f"{rec['id']}: pI={p.pi}"


def test_as_dict_rounds_for_display(records):
    p = compute("t", by_id(records, "SYN-VH-001")["sequence"])
    d = p.as_dict()
    assert d["pH"] == DEFAULT_PH
    assert isinstance(d["pI"], float)
    assert round(d["pI"], 2) == d["pI"]
    assert round(d["GRAVY"], 3) == d["GRAVY"]
