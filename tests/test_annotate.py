"""CDR annotation tests.

The core is `test_ground_truth_cdr_positions`: every ground-truth CDR position in
the synthetic sequences must be recovered. Next is `test_broken_anchors_are_not_confident`:
sequences with broken anchors must yield confident=False, and positions must not be guessed.
"""

from __future__ import annotations

import pytest

from abcheck.annotate import TYPICAL_LENGTHS, annotate, detect_chain

from .conftest import by_id


# --------------------------------------------------------------------------
# Core: recover ground-truth CDR positions
# --------------------------------------------------------------------------

def test_ground_truth_cdr_positions(records):
    """Recover positions and sequences exactly for every synthetic sequence with ground-truth CDR coordinates."""
    checked = 0
    failures: list[str] = []

    for rec in records:
        expected = rec["expected"]["cdrs"]
        if expected is None:          # Cases with deliberately broken anchors are excluded here
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
    """Pick out the intended anchor positions (Cys, Trp, FR3 motif, terminal motif) exactly."""
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
    """The annotator's verdict matches expected.confident recorded in the corpus."""
    for rec in records:
        want = rec["expected"]["confident"]
        ann = annotate(rec["sequence"])
        assert ann.confident is want, (
            f"{rec['id']}: confident 기대 {want}, 결과 {ann.confident} / 사유={ann.reasons}"
        )


# --------------------------------------------------------------------------
# Broken anchors: no guessing
# --------------------------------------------------------------------------

def test_broken_anchors_are_not_confident(records):
    """Sequences with broken anchors are confident=False and must leave a reason in Korean."""
    broken = [r for r in records if r["expected"]["anchors_broken"]]
    assert broken, "앵커가 깨진 케이스가 코퍼스에 없다"

    for rec in broken:
        ann = annotate(rec["sequence"])
        assert ann.confident is False, f"{rec['id']}: confident=False 여야 한다"
        assert ann.reasons, f"{rec['id']}: 사유가 비어 있다"
        # Reasons are written in Korean
        assert any(any("가" <= ch <= "힣" for ch in r) for r in ann.reasons), (
            f"{rec['id']}: 사유에 한글이 없다 — {ann.reasons}"
        )


def test_broken_anchors_do_not_guess_positions(records):
    """Do not guess and fill in CDRs past the point where an anchor was not found."""
    # Terminal WGxG broken -> chain detection itself fails; there must be no CDRs at all
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    assert ann.chain is None
    assert ann.cdrs == {}
    assert ann.annotated is False

    # Cys22 broken -> chain is detected, but no CDR can be cut
    ann = annotate(by_id(records, "SYN-VH-902")["sequence"])
    assert ann.chain == "H"
    assert ann.cdrs == {}
    assert "Cys22" in " ".join(ann.reasons)


def test_broken_terminal_motif_reason_mentions_both_motifs(records):
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    joined = " ".join(ann.reasons)
    assert "WGxG" in joined and "FGxG" in joined


# --------------------------------------------------------------------------
# Length outliers: anchors found, but confident=False
# --------------------------------------------------------------------------

def test_length_outlier_keeps_positions_but_drops_confidence(records):
    """Case where CDR-H1 is 9 residues: positions are recovered exactly, yet confident=False."""
    rec = by_id(records, "SYN-VH-903")
    ann = annotate(rec["sequence"])

    assert ann.confident is False
    assert ann.chain == "H"
    # Anchors are intact, so positions must match the ground truth
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end) == (want["start"], want["end"])

    assert ann.cdrs["CDR-H1"].length == 9
    joined = " ".join(ann.reasons)
    assert "CDR-H1" in joined and "통상 범위" in joined
    assert "ANARCI" in joined


def test_normal_cases_have_typical_lengths(records):
    """In confident=True cases every CDR length is within the typical range."""
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
# Check that CDRs are found by distance from anchors (not by fixed numbering)
# --------------------------------------------------------------------------

def test_shifted_anchors_are_followed(records):
    """Annotate correctly even when a short CDR-L1 shifts all downstream anchors."""
    rec = by_id(records, "SYN-VL-002")
    ann = annotate(rec["sequence"])

    # Unlike typical kappa, Cys2 is not at 88 -- a fixed-numbering implementation breaks here
    assert ann.anchors["cys2"] != 88
    assert ann.anchors["cys2"] == rec["anchors_intended"]["cys2"]
    assert ann.confident is True
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end) == (want["start"], want["end"])


def test_unpaired_cys_in_cdr_h3_does_not_break_anchor(records):
    """Do not mistake an extra Cys inside CDR-H3 for the Cys92 anchor.

    This test backs the 'first Cys after the FR3 motif' rule. Using the 'last Cys'
    would pick an unpaired Cys inside CDR-H3 as the anchor and break the annotation.
    """
    rec = by_id(records, "SYN-VH-005")
    ann = annotate(rec["sequence"])

    assert rec["sequence"].count("C") == 3
    assert ann.anchors["cys2"] == rec["anchors_intended"]["cys2"]
    assert "C" in ann.cdrs["CDR-H3"].seq
    assert ann.confident is True


# --------------------------------------------------------------------------
# Region-splitting helpers
# --------------------------------------------------------------------------

def test_framework_spans_cover_sequence_without_overlap(records):
    for rec in records:
        ann = annotate(rec["sequence"])
        covered = sorted(
            [(c.start, c.end) for c in ann.cdrs.values()] + ann.framework_spans()
        )
        # Must cover 1 through the end with no gaps and no overlaps
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
# Chain detection
# --------------------------------------------------------------------------

def test_detect_chain_prefers_motif_closer_to_c_terminus():
    # FGxG first, WGxG later -> heavy chain
    seq = "A" * 40 + "FGQG" + "A" * 40 + "WGQG" + "A" * 10
    chain, pos = detect_chain(seq)
    assert chain == "H"
    assert pos == seq.rfind("WGQG") + 1

    # Reverse order -> light chain
    seq = "A" * 40 + "WGQG" + "A" * 40 + "FGQG" + "A" * 10
    chain, pos = detect_chain(seq)
    assert chain == "L"


def test_detect_chain_returns_none_without_motifs():
    chain, pos = detect_chain("A" * 100)
    assert chain is None
    assert pos is None


def test_early_terminal_motif_raises_warning():
    """Leave a warning when the terminal motif is in the first half."""
    rec_seq = "WGQG" + "A" * 100
    ann = annotate(rec_seq, validate_input=False)
    assert any("전반부" in w for w in ann.warnings)


# --------------------------------------------------------------------------
# Integration with input validation
# --------------------------------------------------------------------------

def test_annotate_rejects_invalid_input():
    from abcheck.validate import SequenceError

    with pytest.raises(SequenceError):
        annotate("ACDEFG")                       # Too short

    with pytest.raises(SequenceError):
        annotate(">seq1\n" + "A" * 120)          # FASTA header


# --------------------------------------------------------------------------
# Germline variants of the FR3 start motif
# --------------------------------------------------------------------------

def test_fr3_motif_variants_are_found(records):
    """Must also find FR3 variants common in real repertoires (QVTIS, RITIN).

    The narrow pattern [RK].[TSAV][IMLFV][STA] misses both. In 1,000 OAS VH,
    26 (2.6%) failed CDR-H2/H3 annotation because of this.
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
    """Is the widened pattern actually needed? Confirm the narrow pattern cannot find these.

    Without this test, nobody would notice if someone narrowed the pattern again.
    """
    import re

    from abcheck.annotate import FR3_MOTIF_PATTERN_STRICT

    strict = re.compile(FR3_MOTIF_PATTERN_STRICT)
    for motif in ["QVTIS", "RITIN"]:
        assert not strict.fullmatch(motif), f"{motif} 가 좁은 패턴에 걸린다 — 테스트 전제가 틀렸다"


def test_broken_fr3_motif_falls_back_to_cys_distance(records):
    """Even if hypermutation breaks the FR3 motif, recover it by back-calculating from the conserved Cys.

    This happened in 8 OAS mouse plasmablasts (cells with heavy hypermutation). The motif
    had changed to EATLT / RLNIS / KAILT and no longer matched the regex. Across all 992 real VH,
    the FR3 motif -> conserved Cys distance was a constant 29, so that distance is used.
    """
    rec = by_id(records, "SYN-VH-013")
    ann = annotate(rec["sequence"])

    # Positions are still reported -- the result is not discarded
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end, got.seq) == (want["start"], want["end"], want["seq"]), name

    # But the evidence is weaker than seeing the motif directly, so confidence is lowered
    assert ann.confident is False
    assert any("역산" in r for r in ann.reasons), ann.reasons
    assert any("ANARCI" in r for r in ann.reasons)


def test_fr3_to_cys_gap_holds_across_corpus(records):
    """Is the FR3 motif to conserved Cys distance constant across the whole corpus?

    If this distance drifts, the back-calculation path silently goes wrong.
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
    """No back-calculation warning when the motif is intact -- the fallback path must not fire freely."""
    for seq_id in ["SYN-VH-001", "SYN-VH-002", "SYN-VH-011", "SYN-VH-012"]:
        ann = annotate(by_id(records, seq_id)["sequence"])
        assert not any("역산" in r for r in ann.reasons + ann.warnings), (seq_id, ann.reasons)
        assert ann.confident is True, (seq_id, ann.reasons)


def test_decoy_fr3_motif_is_corrected_by_cys_distance(records):
    """Correct a fake FR3 motif inside CDR-H2 by cross-checking against the conserved Cys distance.

    In 7,000 real sequences, the FR3 motif <-> conserved Cys distance was exactly 29 in all 3,429,
    and the 14 outliers were all cases where the motif matched in the wrong place (CDR-H2 came out as 15 or 29).
    """
    rec = by_id(records, "SYN-VH-014")
    ann = annotate(rec["sequence"])

    # Even if a fake motif matched, the final CDR-H2 must be correct
    for name, want in rec["expected"]["cdrs"].items():
        got = ann.cdrs[name]
        assert (got.start, got.end, got.seq) == (want["start"], want["end"], want["seq"]), name

    # Do not hide the fact that a correction was made
    assert ann.confident is False
    assert any("통상값" in r and "고쳐 잡았다" in r for r in ann.reasons), ann.reasons


# --------------------------------------------------------------------------
# Trp anchor based on FR2 context
# --------------------------------------------------------------------------

def test_trp_anchor_uses_fr2_context_for_light_chain():
    """In the IGKV1 form with two Trp at the end of CDR-L1 and in FR2, pick the later one (the real FR2).

    `...RASQGISS W LA W YQQ...` -- picking the earlier Trp shifts every later boundary by 3 residues.
    Compared against OAS's ANARCI annotations as ground truth, this fix brought CDR-L1/L2 mismatches to 0.
    """
    seq = ("DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSG"
           "TDFTLTISSLQPEDFATYYCQQYNSYPPTFGQGTKVEIK")
    ann = annotate(seq)

    assert ann.anchors["trp"] == 35, "FR2 의 Trp35 를 집어야 한다"
    assert (ann.cdrs["CDR-L1"].start, ann.cdrs["CDR-L1"].end) == (24, 34)
    assert (ann.cdrs["CDR-L2"].start, ann.cdrs["CDR-L2"].end) == (50, 56)
    assert (ann.cdrs["CDR-L3"].start, ann.cdrs["CDR-L3"].end) == (89, 97)
    assert ann.cdrs["CDR-L3"].seq == "QQYNSYPPT"


def test_trp_anchor_uses_fr2_context_for_heavy_chain():
    """Even with a Trp inside CDR-H1, FR2 context picks the real Trp36."""
    seq = ("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMQWVRQAPGQGLEWMGAIYPGDGDTRYTQKFKG"
           "RVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGGLRRFDYWGQGTLVTVSS")
    ann = annotate(seq)

    assert ann.anchors["trp"] == 36
    assert ann.cdrs["CDR-H1"].seq == "SYWMQ"
    assert not any("문맥" in w for w in ann.warnings), ann.warnings


def test_fallback_warns_when_fr2_context_is_absent(records):
    """If the context was not found and it fell back to the first Trp, say so."""
    from abcheck.annotate import FR2_CONTEXT_H, _find_conserved_trp

    seq = "A" * 30 + "W" + "A" * 30
    pos, by_context = _find_conserved_trp(seq, 20, FR2_CONTEXT_H)
    assert pos == 31
    assert by_context is False


# --------------------------------------------------------------------------
# Terminal motif mutated, J-gene tail intact (found on OAS data)
# --------------------------------------------------------------------------

def _with_end(seq: str, motif_pos: int, new_end: str) -> str:
    """Replace everything from the terminal motif (1-based) with new_end."""
    return seq[: motif_pos - 1] + new_end


@pytest.mark.parametrize("end", ["WAQGTLVTVSS", "WSQGTLVTVSS", "WGRASLVTVSS", "SGQGTMVTVSS"])
def test_heavy_with_mutated_wgxg_is_rescued_by_j_tail(records, end):
    """Complete VH whose WGxG is substituted used to fail chain detection outright."""
    rec = by_id(records, "SYN-VH-001")
    ref = annotate(rec["sequence"])
    seq = _with_end(rec["sequence"], ref.anchors["terminal_motif"], end)

    ann = annotate(seq)
    assert ann.chain == "H"
    assert ann.anchors["terminal_motif"] == ref.anchors["terminal_motif"]
    assert ann.cdrs["CDR-H3"].seq == ref.cdrs["CDR-H3"].seq
    assert any("J 유전자 꼬리" in w for w in ann.warnings)


@pytest.mark.parametrize("end", ["FAQGTRLEIK", "LGQGTKLEIK", "FAGGTKVTVL", "LGGGTKLTVL", "YGGGTKLTVL"])
def test_light_with_mutated_fgxg_is_rescued_by_j_tail(records, end):
    rec = by_id(records, "SYN-VL-001")
    ref = annotate(rec["sequence"])
    seq = _with_end(rec["sequence"], ref.anchors["terminal_motif"], end)

    ann = annotate(seq)
    assert ann.chain == "L"
    assert ann.anchors["terminal_motif"] == ref.anchors["terminal_motif"]
    assert ann.cdrs["CDR-L3"].seq == ref.cdrs["CDR-L3"].seq


def test_j_tails_do_not_cross_chains():
    """A heavy J tail must never read as light, and vice versa."""
    from abcheck.annotate import detect_chain_by_j_tail

    pad = "S" * 100
    assert detect_chain_by_j_tail(pad + "WAQGTLVTVSS")[0] == "H"
    assert detect_chain_by_j_tail(pad + "FAQGTRLEIK")[0] == "L"
    assert detect_chain_by_j_tail(pad + "LGGGTKLTVL")[0] == "L"
    assert detect_chain_by_j_tail(pad + "SSSSSSSSSS") == (None, None)


def test_no_motif_and_no_j_tail_still_fails(records):
    """SYN-VH-901 has a broken WGxG and a random FR4: still no chain, no guessed CDRs."""
    ann = annotate(by_id(records, "SYN-VH-901")["sequence"])
    assert ann.chain is None
    assert not ann.cdrs
