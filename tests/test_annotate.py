"""CDR 주석 테스트.

핵심은 `test_ground_truth_cdr_positions` 다 — 합성 서열의 정답 CDR 위치를 전부
복원해야 한다. 그 다음이 `test_broken_anchors_are_not_confident` — 앵커가 깨진
서열에서 confident=False 가 나와야 하고, 위치를 추측해서는 안 된다.
"""

from __future__ import annotations

import pytest

from abcheck.annotate import TYPICAL_LENGTHS, annotate, detect_chain

from .conftest import by_id


# --------------------------------------------------------------------------
# 핵심: 정답 CDR 위치 복원
# --------------------------------------------------------------------------

def test_ground_truth_cdr_positions(records):
    """정답 CDR 좌표가 정의된 모든 합성 서열에서 위치·서열을 정확히 복원한다."""
    checked = 0
    failures: list[str] = []

    for rec in records:
        expected = rec["expected"]["cdrs"]
        if expected is None:          # 앵커가 일부러 깨진 케이스는 여기서 제외
            continue

        ann = annotate(rec["sequence"])

        if set(ann.cdrs) != set(expected):
            failures.append(
                f"{rec['id']}: CDR 집합 불일치 — 기대 {sorted(expected)}, 결과 {sorted(ann.cdrs)}"
            )
            continue

        for name, want in expected.items():
            got = ann.cdrs[name]
            if (got.start, got.end, got.seq) != (want["start"], want["end"], want["seq"]):
                failures.append(
                    f"{rec['id']} {name}: 기대 {want['start']}-{want['end']} {want['seq']!r}, "
                    f"결과 {got.start}-{got.end} {got.seq!r}"
                )
        checked += 1

    assert not failures, "정답 CDR 위치 복원 실패:\n" + "\n".join(failures)
    assert checked >= 12, f"검증한 서열이 {checked}개뿐이다 — 코퍼스가 줄어든 것 아닌지 확인하라"


def test_ground_truth_anchor_positions(records):
    """의도한 앵커 위치(Cys·Trp·FR3 모티프·말단 모티프)를 그대로 집어낸다."""
    for rec in records:
        intended = rec["anchors_intended"]
        if intended is None:
            continue
        ann = annotate(rec["sequence"])
        for key, want in intended.items():
            assert ann.anchors.get(key) == want, (
                f"{rec['id']} 앵커 {key}: 기대 {want}, 결과 {ann.anchors.get(key)}"
            )


def test_expected_chain_is_detected(records):
    for rec in records:
        want = rec["expected"]["chain"]
        if want is None:
            continue
        ann = annotate(rec["sequence"])
        assert ann.chain == want, f"{rec['id']}: 체인 기대 {want}, 결과 {ann.chain}"


def test_expected_confidence_matches(records):
    """코퍼스가 기록한 expected.confident 와 주석기 판단이 일치한다."""
    for rec in records:
        want = rec["expected"]["confident"]
        ann = annotate(rec["sequence"])
        assert ann.confident is want, (
            f"{rec['id']}: confident 기대 {want}, 결과 {ann.confident} / 사유={ann.reasons}"
        )


# --------------------------------------------------------------------------
# 앵커가 깨진 경우 — 추측 금지
# --------------------------------------------------------------------------

def test_broken_anchors_are_not_confident(records):
    """앵커가 깨진 서열은 confident=False 이고, 한글 사유가 남아야 한다."""
    broken = [r for r in records if r["expected"]["anchors_broken"]]
    assert broken, "앵커가 깨진 케이스가 코퍼스에 없다"

    for rec in broken:
        ann = annotate(rec["sequence"])
        assert ann.confident is False, f"{rec['id']}: confident=False 여야 한다"
        assert ann.reasons, f"{rec['id']}: 사유가 비어 있다"
        # 사유는 한글로 남긴다
        assert any(any("가" <= ch <= "힣" for ch in r) for r in ann.reasons), (
            f"{rec['id']}: 사유에 한글이 없다 — {ann.reasons}"
        )


def test_broken_anchors_do_not_guess_positions(records):
    """앵커를 못 찾은 지점 이후의 CDR 을 추측해서 채우지 않는다."""
    # 말단 WGxG 가 깨짐 → 체인 판별부터 실패, CDR 하나도 없어야 한다
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    assert ann.chain is None
    assert ann.cdrs == {}
    assert ann.annotated is False

    # Cys22 가 깨짐 → 체인은 판별되지만 CDR 은 하나도 끊지 못해야 한다
    ann = annotate(by_id(records, "SYN-VH-902")["sequence"])
    assert ann.chain == "H"
    assert ann.cdrs == {}
    assert "Cys22" in " ".join(ann.reasons)


def test_broken_terminal_motif_reason_mentions_both_motifs(records):
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    joined = " ".join(ann.reasons)
    assert "WGxG" in joined and "FGxG" in joined


# --------------------------------------------------------------------------
# 길이 이탈 — 앵커는 찾았지만 confident=False
# --------------------------------------------------------------------------

def test_length_outlier_keeps_positions_but_drops_confidence(records):
    """CDR-H1 이 9 잔기인 케이스: 위치는 정확히 복원하면서 confident=False."""
    rec = by_id(records, "SYN-VH-903")
    ann = annotate(rec["sequence"])

    assert ann.confident is False
    assert ann.chain == "H"
    # 앵커가 멀쩡하므로 위치는 정답과 같아야 한다
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end) == (want["start"], want["end"])

    assert ann.cdrs["CDR-H1"].length == 9
    joined = " ".join(ann.reasons)
    assert "CDR-H1" in joined and "통상 범위" in joined
    assert "ANARCI" in joined


def test_normal_cases_have_typical_lengths(records):
    """confident=True 인 케이스는 모든 CDR 길이가 통상 범위 안이다."""
    for rec in records:
        if not rec["expected"]["confident"]:
            continue
        ann = annotate(rec["sequence"])
        for name, cdr in ann.cdrs.items():
            lo, hi = TYPICAL_LENGTHS[name]
            assert lo <= cdr.length <= hi, (
                f"{rec['id']} {name}: 길이 {cdr.length} 가 통상 범위({lo}-{hi}) 밖"
            )


# --------------------------------------------------------------------------
# 앵커 상대 거리로 찾는지 (고정 번호로 찍지 않는지)
# --------------------------------------------------------------------------

def test_shifted_anchors_are_followed(records):
    """CDR-L1 이 짧아 뒤 앵커가 전부 밀린 서열도 올바르게 주석한다."""
    rec = by_id(records, "SYN-VL-002")
    ann = annotate(rec["sequence"])

    # 통상 kappa 와 달리 Cys2 가 88번이 아니다 — 고정 번호 구현이면 여기서 깨진다
    assert ann.anchors["cys2"] != 88
    assert ann.anchors["cys2"] == rec["anchors_intended"]["cys2"]
    assert ann.confident is True
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end) == (want["start"], want["end"])


def test_unpaired_cys_in_cdr_h3_does_not_break_anchor(records):
    """CDR-H3 안의 여분 Cys 를 Cys92 앵커로 잘못 집지 않는다.

    'FR3 모티프 이후 첫 Cys' 규칙의 근거가 되는 테스트다. '마지막 Cys' 를 쓰면
    CDR-H3 안의 비결합 Cys 를 앵커로 집어 주석이 망가진다.
    """
    rec = by_id(records, "SYN-VH-005")
    ann = annotate(rec["sequence"])

    assert rec["sequence"].count("C") == 3
    assert ann.anchors["cys2"] == rec["anchors_intended"]["cys2"]
    assert "C" in ann.cdrs["CDR-H3"].seq
    assert ann.confident is True


# --------------------------------------------------------------------------
# 구간 분할 보조 기능
# --------------------------------------------------------------------------

def test_framework_spans_cover_sequence_without_overlap(records):
    for rec in records:
        ann = annotate(rec["sequence"])
        covered = sorted(
            [(c.start, c.end) for c in ann.cdrs.values()] + ann.framework_spans()
        )
        # 1 부터 끝까지 빈틈 없이, 겹침 없이 덮어야 한다
        cursor = 1
        for start, end in covered:
            assert start == cursor, f"{rec['id']}: 구간이 {cursor} 에서 끊겼다 ({start}-{end})"
            assert start <= end
            cursor = end + 1
        assert cursor == len(rec["sequence"]) + 1, f"{rec['id']}: 끝까지 덮지 못했다"


def test_region_of_matches_cdr_boundaries(records):
    rec = by_id(records, "SYN-VH-001")
    ann = annotate(rec["sequence"])
    h1 = ann.cdrs["CDR-H1"]

    assert ann.region_of(h1.start) == "CDR-H1"
    assert ann.region_of(h1.end) == "CDR-H1"
    assert ann.region_of(h1.start - 1) == "framework"
    assert ann.region_of(h1.end + 1) == "framework"


def test_framework_spans_cover_everything_when_annotation_failed(records):
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    assert ann.framework_spans() == [(1, len(ann.sequence))]


# --------------------------------------------------------------------------
# 체인 판별
# --------------------------------------------------------------------------

def test_detect_chain_prefers_motif_closer_to_c_terminus():
    # FGxG 가 앞, WGxG 가 뒤 → 중쇄
    seq = "A" * 40 + "FGQG" + "A" * 40 + "WGQG" + "A" * 10
    chain, pos = detect_chain(seq)
    assert chain == "H"
    assert pos == seq.rfind("WGQG") + 1

    # 반대 순서 → 경쇄
    seq = "A" * 40 + "WGQG" + "A" * 40 + "FGQG" + "A" * 10
    chain, pos = detect_chain(seq)
    assert chain == "L"


def test_detect_chain_returns_none_without_motifs():
    chain, pos = detect_chain("A" * 100)
    assert chain is None
    assert pos is None


def test_early_terminal_motif_raises_warning():
    """말단 모티프가 전반부에 있으면 경고를 남긴다."""
    rec_seq = "WGQG" + "A" * 100
    ann = annotate(rec_seq, validate_input=False)
    assert any("전반부" in w for w in ann.warnings)


# --------------------------------------------------------------------------
# 입력 검증 연동
# --------------------------------------------------------------------------

def test_annotate_rejects_invalid_input():
    from abcheck.validate import SequenceError

    with pytest.raises(SequenceError):
        annotate("ACDEFG")                       # 너무 짧다

    with pytest.raises(SequenceError):
        annotate(">seq1\n" + "A" * 120)          # FASTA 헤더


# --------------------------------------------------------------------------
# FR3 시작 모티프의 germline 변형
# --------------------------------------------------------------------------

def test_fr3_motif_variants_are_found(records):
    """실제 레퍼토리에서 흔한 FR3 변형(QVTIS, RITIN)도 찾아야 한다.

    좁은 패턴 [RK].[TSAV][IMLFV][STA] 로는 둘 다 놓친다. OAS VH 1,000건에서
    26건(2.6%)이 이 때문에 CDR-H2/H3 주석에 실패했다.
    """
    for seq_id, motif in [("SYN-VH-011", "QVTIS"), ("SYN-VH-012", "RITIN")]:
        rec = by_id(records, seq_id)
        ann = annotate(rec["sequence"])

        assert ann.confident is True, f"{seq_id}: {ann.reasons}"
        assert set(ann.cdrs) == {"CDR-H1", "CDR-H2", "CDR-H3"}, seq_id
        fr3 = ann.anchors["fr3_motif"]
        assert ann.sequence[fr3 - 1:fr3 + 4] == motif, seq_id
        for name, want in rec["expected"]["cdrs"].items():
            got = ann.cdrs[name]
            assert (got.start, got.end) == (want["start"], want["end"]), (seq_id, name)


def test_strict_fr3_pattern_would_miss_those_variants(records):
    """넓힌 패턴이 실제로 필요한지 — 좁은 패턴으로는 못 찾는 것을 확인한다.

    이 테스트가 없으면 누가 패턴을 좁혀도 아무도 모른다.
    """
    import re

    from abcheck.annotate import FR3_MOTIF_PATTERN_STRICT

    strict = re.compile(FR3_MOTIF_PATTERN_STRICT)
    for motif in ["QVTIS", "RITIN"]:
        assert not strict.fullmatch(motif), f"{motif} 가 좁은 패턴에 걸린다 — 테스트 전제가 틀렸다"


def test_broken_fr3_motif_falls_back_to_cys_distance(records):
    """FR3 모티프가 과변이로 깨져도 보존 Cys 에서 역산해 복원한다.

    OAS 마우스 plasmablast(과변이가 심한 세포)에서 8건이 이 경우였다. 모티프가
    EATLT / RLNIS / KAILT 처럼 바뀌어 정규식에 안 걸렸다. 실제 VH 992건 전부에서
    FR3 모티프→보존 Cys 거리가 29 로 일정했으므로 그 거리를 쓴다.
    """
    rec = by_id(records, "SYN-VH-013")
    ann = annotate(rec["sequence"])

    assert ann.confident is True, ann.reasons
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end, got.seq) == (want["start"], want["end"], want["seq"]), name

    # 추측이 아니라 역산이라는 것을 사용자에게 알려야 한다
    assert any("역산" in w for w in ann.warnings), ann.warnings
    assert any("ANARCI" in w for w in ann.warnings)


def test_fr3_to_cys_gap_holds_across_corpus(records):
    """FR3 모티프와 보존 Cys 의 거리가 코퍼스 전체에서 일정한지.

    이 거리가 흔들리면 역산 경로가 조용히 틀어진다.
    """
    from abcheck.annotate import FR3_TO_CYS_GAP

    checked = 0
    for rec in records:
        ann = annotate(rec["sequence"])
        if "fr3_motif" in ann.anchors and "cys2" in ann.anchors:
            assert ann.anchors["cys2"] - ann.anchors["fr3_motif"] == FR3_TO_CYS_GAP, rec["id"]
            checked += 1
    assert checked >= 10


def test_fallback_is_not_used_when_motif_is_present(records):
    """모티프가 멀쩡하면 역산 경고가 붙지 않아야 한다 — 대체 경로가 남발되면 곤란하다."""
    for seq_id in ["SYN-VH-001", "SYN-VH-002", "SYN-VH-011", "SYN-VH-012"]:
        ann = annotate(by_id(records, seq_id)["sequence"])
        assert not any("역산" in w for w in ann.warnings), (seq_id, ann.warnings)
