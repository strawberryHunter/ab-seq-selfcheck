#!/usr/bin/env python3
"""Synthetic antibody sequence generator (standard library only).

Purpose
-------
Build a test corpus for validating abcheck's CDR annotation / liability detection logic.
Sequences are assembled as "germline-like framework scaffold + inserted CDRs", so the
ground-truth CDR positions (1-based, inclusive) are known exactly at insertion time and are stored alongside.

Important -- statement on data provenance
-----------------------------------------
Every sequence this file produces is synthetic. Nothing is copied from real therapeutic/approved
antibody sequences or from specific germline database records. The framework parts are merely
"germline-like" scaffolds: weighted random residues mimicking the composition common in V domains
(hydrophobic/small residues dominate), with only conserved residues (Cys, Trp, WGxG/FGxG, etc.)
fixed at set positions. They do not correspond to any real sequence.
This data is therefore for regression-testing the annotation logic, not for biological interpretation.

Scaffold coordinate design (VH, 1-based linear positions; mimicking Kabat boundaries)
-------------------------------------------------------------------------------------
    1..21      FR1 filler
    22         Cys  (anchor 1)
    23..30     filler (8)
    31..       CDR-H1  (= starts at Cys22 + 9)
    after H1   Trp     (anchor 2, typically Trp36)
    Trp+1..+13 filler (13)
    Trp+14..   CDR-H2
    after H2   FR3 motif, 5 residues  [RK].[TSAV][IMLFV][STA]  (anchor 3)
    motif+5..  filler (21) -> then Cys (anchor 4, typically Cys92)
    Cys+1..+2  filler (2)
    Cys+3..    CDR-H3
    after H3   WGxG (anchor 5)
    then       FR4 filler (7)

VL
--
    1..22      FR1 filler
    23         Cys  (anchor 1)
    24..       CDR-L1 (= Cys23 + 1)
    after L1   Trp    (anchor 2, typically Trp35)
    Trp+1..+14 filler (14)
    Trp+15..   CDR-L2 (fixed at 7 residues)
    after L2   filler (31) -> then Cys (anchor 3, typically Cys88)
    Cys+1..    CDR-L3
    after L3   FGxG (anchor 4)
    then       FR4 filler (7)

After assembly, the generator checks itself that "the search rules the annotator will use" pick
exactly the intended anchors (see ASSERTIONS). If not, it rebuilds with a different seed. So this
corpus contains only controlled input where the search rules are not tripped by accidental duplicate
motifs. Cases with deliberately broken anchors are marked separately with expected.confident = False.
"""

from __future__ import annotations

import argparse
import json
import random
import re
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Alphabets / motifs
# --------------------------------------------------------------------------

#: Residues and weights used for framework filler.
#: Cys/Trp/Met/Asn/Asp are deliberately excluded -- this keeps the conserved Cys/Trp anchors unique
#: and keeps accidental liability motifs (NG, DG, DP, N[^P][ST], oxidizable M/W) out of the framework.
#: Cases that need a framework liability inject it explicitly below.
FRAMEWORK_WEIGHTS: dict[str, int] = {
    "A": 7, "G": 7, "V": 7, "L": 8, "I": 5, "P": 5, "F": 4,
    "S": 9, "T": 8, "Y": 4, "Q": 6, "E": 6, "K": 6, "R": 5, "H": 2,
}

#: CDR filler alphabet. Cys/Trp excluded (added only when injected on purpose).
CDR_WEIGHTS: dict[str, int] = {
    "A": 5, "G": 7, "V": 4, "L": 4, "I": 4, "P": 3, "F": 3,
    "S": 7, "T": 6, "Y": 7, "Q": 3, "E": 3, "K": 3, "R": 4,
    "H": 2, "N": 4, "D": 5,
}

#: Liability motifs used for the CDR cleanliness check.
#: They let the generator guarantee "no liabilities inside CDRs other than the intended ones",
#: and are kept separate from the real risk-assessment logic (abcheck.liabilities).
CDR_LIABILITY_PATTERNS: list[tuple[str, str]] = [
    ("glycosylation", r"N[^P][ST]"),
    ("deamidation", r"N[GSTNH]"),
    ("isomerization", r"D[GSTD]"),
    ("fragmentation", r"DP"),
    ("oxidation", r"[MW]"),
    ("cysteine", r"C"),
]

#: Must match abcheck.annotate.FR3_MOTIF_PATTERN.
FR3_MOTIF_RE = re.compile(r"[RKQ].[TSAV][IMLFV][STAN]")

#: Distance from the FR3 start motif to the conserved Cys. Must match abcheck.annotate.FR3_TO_CYS_GAP.
FR3_TO_CYS_GAP = 29

#: Offset at which to start searching for the conserved Trp anchor. Same value as abcheck.annotate.
MIN_CDR1_LENGTH_H = 5
MIN_CDR1_LENGTH_L = 7
WGXG_RE = re.compile(r"WG.G")
FGXG_RE = re.compile(r"FG.G")

# CDR lengths the annotator treats as "typical range" (a screening heuristic, not a strict numbering standard)
TYPICAL_LENGTHS = {
    "H1": (5, 7),
    "H2": (16, 20),
    "H3": (3, 30),
    "L1": (7, 17),
    "L2": (7, 7),
    "L3": (5, 15),
}


# --------------------------------------------------------------------------
# Scaffolds
# --------------------------------------------------------------------------

@dataclass
class HeavyScaffold:
    """VH framework pieces. All VH cases share the same scaffold."""

    fr1: str          # 21
    post_cys1: str    # 8
    fr2: str          # 13
    fr3_motif: str    # 5
    fr3_rest: str     # 21
    post_cys2: str    # 2
    fr4_rest: str     # 7


@dataclass
class LightScaffold:
    """VL framework pieces."""

    fr1: str          # 22
    fr2: str          # 14
    fr3_rest: str     # 31
    fr4_rest: str     # 7


def _draw(rng: random.Random, weights: dict[str, int], n: int) -> str:
    pool = "".join(res * w for res, w in weights.items())
    return "".join(rng.choice(pool) for _ in range(n))


def build_heavy_scaffold(rng: random.Random) -> HeavyScaffold:
    motif = (
        rng.choice("RK")
        + rng.choice("AGVLSTQE")
        + rng.choice("TSAV")
        + rng.choice("IMLFV")
        + rng.choice("STA")
    )
    return HeavyScaffold(
        fr1=_draw(rng, FRAMEWORK_WEIGHTS, 21),
        post_cys1=_draw(rng, FRAMEWORK_WEIGHTS, 8),
        # In real sequences VRQ/IRQ overwhelmingly follow the conserved Trp. The annotator uses this
        # context for the anchor search, so match it. Length (13) is unchanged.
        fr2="VRQ" + _draw(rng, FRAMEWORK_WEIGHTS, 10),
        fr3_motif=motif,
        # In 83% of real sequences the conserved Cys is preceded by YY. The annotator uses [YF]C context, so match it.
        # Length 22 + "YY" = 24. Adding the 5-residue motif makes the FR3->conserved Cys distance 29.
        # All 992 real VH were exactly 29, so match that (the annotator back-calculates with this distance).
        fr3_rest=_draw(rng, FRAMEWORK_WEIGHTS, 22) + "YY",
        post_cys2=_draw(rng, FRAMEWORK_WEIGHTS, 2),
        fr4_rest=_draw(rng, FRAMEWORK_WEIGHTS, 7),
    )


def build_light_scaffold(rng: random.Random) -> LightScaffold:
    return LightScaffold(
        fr1=_draw(rng, FRAMEWORK_WEIGHTS, 22),
        fr2="YQQ" + _draw(rng, FRAMEWORK_WEIGHTS, 11),
        fr3_rest=_draw(rng, FRAMEWORK_WEIGHTS, 29) + "YY",
        fr4_rest=_draw(rng, FRAMEWORK_WEIGHTS, 7),
    )


# --------------------------------------------------------------------------
# Case definitions
# --------------------------------------------------------------------------

@dataclass
class Case:
    """Blueprint for one synthetic sequence."""

    seq_id: str
    chain: str                       # "H" | "L"
    description: str                 # Korean description
    cdr_lengths: tuple[int, int, int]
    #: Fragments to force into CDRs. {"H1": [(offset, "NLS")], ...} (offset is 0-based)
    cdr_inserts: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    #: Fragments to force into framework pieces. {"fr3_rest": [(offset, "DG")], ...}
    fr_inserts: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    #: Substitutions applied to the assembled sequence (1-based position, new residue). Used to break anchors.
    breakages: list[tuple[int, str]] = field(default_factory=list)
    #: Case whose anchors are broken so annotation must fail
    anchors_broken: bool = False
    #: Case with intact anchors but a CDR length outside the typical range, so it must be confident=False
    length_outlier: bool = False
    #: Pin the FR3 start motif to this value (to reproduce variants of real germline families).
    fr3_motif: str | None = None
    #: Case with a deliberately inserted FR3-motif lookalike inside CDR-H2.
    #: Checks that the annotator corrects it by cross-checking the Cys distance.
    decoy_fr3_motif: bool = False
    #: Case where the motif does not match the regex, so the annotator must take the Cys back-calculation path.
    #: Back-calculation is weaker evidence than seeing the motif directly, so it must be confident=False.
    fr3_motif_broken: bool = False
    notes: list[str] = field(default_factory=list)


def _apply_inserts(seq: str, inserts: list[tuple[int, str]]) -> str:
    """Overwrite a fragment at offset while keeping the length (preserving length is what preserves coordinates)."""
    buf = list(seq)
    for offset, frag in inserts:
        if offset < 0 or offset + len(frag) > len(buf):
            raise ValueError(f"삽입 조각이 구간을 벗어난다: offset={offset} frag={frag} len={len(buf)}")
        buf[offset:offset + len(frag)] = list(frag)
    return "".join(buf)


def _make_cdr(rng: random.Random, name: str, length: int,
              inserts: dict[str, list[tuple[int, str]]]) -> str:
    seq = _draw(rng, CDR_WEIGHTS, length)
    return _apply_inserts(seq, inserts.get(name, []))


def _cdr_liability_hits(cdr: str) -> list[tuple[int, int, str]]:
    """Return positions of liability motifs inside a CDR as (0-based start, end+1, label)."""
    hits: list[tuple[int, int, str]] = []
    for label, pattern in CDR_LIABILITY_PATTERNS:
        for m in re.finditer(pattern, cdr):
            hits.append((m.start(), m.end(), label))
    return hits


#: Maximum liability motif length (N[^P][ST] = 3) minus 1.
#: To see motifs straddling a CDR boundary, scan this much wider on each side.
_MOTIF_PAD = 2


def unintended_cdr_liabilities_in_context(
    seq: str, cdr_start: int, cdr_end: int, inserts: list[tuple[int, str]]
) -> list[tuple[int, int, str]]:
    """Liabilities that **overlap** a CDR but are not contained in an injected fragment.

    Looking only at the CDR string misses motifs that straddle a boundary. For example, if a CDR
    ends in ...SNA and the following framework starts with S, NAS becomes a glycosylation sequon,
    but scanning only inside the CDR misses it. So scan 2 residues wider on each side and keep only hits that overlap the CDR.

    All coordinates are 1-based and inclusive.
    """
    lo = max(1, cdr_start - _MOTIF_PAD)
    hi = min(len(seq), cdr_end + _MOTIF_PAD)
    window = seq[lo - 1:hi]

    #: Global coordinate ranges of injected fragments
    spans = [(cdr_start + off, cdr_start + off + len(frag) - 1) for off, frag in inserts]

    out: list[tuple[int, int, str]] = []
    for w_start, w_end, label in _cdr_liability_hits(window):
        start = lo + w_start           # 1-based start
        end = lo + w_end - 1           # 1-based end (inclusive)
        if end < cdr_start or start > cdr_end:
            continue                   # If it doesn't overlap the CDR, it's the framework's business
        if any(s_lo <= start and end <= s_hi for s_lo, s_hi in spans):
            continue                   # Injected on purpose
        out.append((start, end, label))
    return out


def _assemble_heavy(case: Case, scaf: HeavyScaffold, rng: random.Random):
    n1, n2, n3 = case.cdr_lengths
    h1 = _make_cdr(rng, "H1", n1, case.cdr_inserts)
    h2 = _make_cdr(rng, "H2", n2, case.cdr_inserts)
    h3 = _make_cdr(rng, "H3", n3, case.cdr_inserts)

    fr1 = _apply_inserts(scaf.fr1, case.fr_inserts.get("fr1", []))
    post_cys1 = _apply_inserts(scaf.post_cys1, case.fr_inserts.get("post_cys1", []))
    fr2 = _apply_inserts(scaf.fr2, case.fr_inserts.get("fr2", []))
    fr3_rest = _apply_inserts(scaf.fr3_rest, case.fr_inserts.get("fr3_rest", []))
    fr4_rest = _apply_inserts(scaf.fr4_rest, case.fr_inserts.get("fr4_rest", []))

    motif = case.fr3_motif or scaf.fr3_motif
    parts = [fr1, "C", post_cys1, h1, "W", fr2, h2,
             motif, fr3_rest, "C", scaf.post_cys2, h3, "WGQG", fr4_rest]
    seq = "".join(parts)

    # Compute 1-based coordinates
    cys1 = len(fr1) + 1                                   # 22
    h1_start = cys1 + 9                                   # 31
    h1_end = h1_start + n1 - 1
    trp = h1_end + 1
    h2_start = trp + 14
    h2_end = h2_start + n2 - 1
    fr3_start = h2_end + 1
    cys2 = fr3_start + len(motif) + len(fr3_rest)
    h3_start = cys2 + 3
    h3_end = h3_start + n3 - 1
    wgxg = h3_end + 1

    anchors = {"cys1": cys1, "trp": trp, "fr3_motif": fr3_start, "cys2": cys2, "terminal_motif": wgxg}
    cdrs = {
        "CDR-H1": (h1_start, h1_end),
        "CDR-H2": (h2_start, h2_end),
        "CDR-H3": (h3_start, h3_end),
    }
    return seq, cdrs, anchors


def _assemble_light(case: Case, scaf: LightScaffold, rng: random.Random):
    n1, n2, n3 = case.cdr_lengths
    l1 = _make_cdr(rng, "L1", n1, case.cdr_inserts)
    l2 = _make_cdr(rng, "L2", n2, case.cdr_inserts)
    l3 = _make_cdr(rng, "L3", n3, case.cdr_inserts)

    fr1 = _apply_inserts(scaf.fr1, case.fr_inserts.get("fr1", []))
    fr2 = _apply_inserts(scaf.fr2, case.fr_inserts.get("fr2", []))
    fr3_rest = _apply_inserts(scaf.fr3_rest, case.fr_inserts.get("fr3_rest", []))
    fr4_rest = _apply_inserts(scaf.fr4_rest, case.fr_inserts.get("fr4_rest", []))

    parts = [fr1, "C", l1, "W", fr2, l2, fr3_rest, "C", l3, "FGQG", fr4_rest]
    seq = "".join(parts)

    cys1 = len(fr1) + 1                                   # 23
    l1_start = cys1 + 1                                   # 24
    l1_end = l1_start + n1 - 1
    trp = l1_end + 1
    l2_start = trp + 15
    l2_end = l2_start + n2 - 1
    cys2 = l2_end + len(fr3_rest) + 1
    l3_start = cys2 + 1
    l3_end = l3_start + n3 - 1
    fgxg = l3_end + 1

    anchors = {"cys1": cys1, "trp": trp, "cys2": cys2, "terminal_motif": fgxg}
    cdrs = {
        "CDR-L1": (l1_start, l1_end),
        "CDR-L2": (l2_start, l2_end),
        "CDR-L3": (l3_start, l3_end),
    }
    return seq, cdrs, anchors


# --------------------------------------------------------------------------
# Self-check -- confirm the annotator's search rules pick the intended anchors
# --------------------------------------------------------------------------

def _first_index(pattern: re.Pattern[str], seq: str, min_pos1: int) -> int | None:
    """1-based start of the first match starting at or after min_pos1 (1-based)."""
    m = pattern.search(seq, min_pos1 - 1)
    return m.start() + 1 if m else None


CONSERVED_CYS_RE = re.compile(r"[YF]C")

#: Must match the FR2 context in abcheck.annotate.
FR2_CONTEXT_H = re.compile(r"W[VIA][RKQ]Q")
FR2_CONTEXT_L = re.compile(r"W[YF][QL]Q")


def _first_fr2_trp(seq: str, min_pos: int, context: re.Pattern[str]) -> int | None:
    """Position of the Trp in the first FR2 context at or after min_pos (1-based). Falls back to the first Trp."""
    m = context.search(seq, max(min_pos - 1, 0))
    if m is not None and m.start() + 1 >= min_pos:
        return m.start() + 1
    idx = seq.find("W", max(min_pos - 1, 0))
    return idx + 1 if idx >= 0 else None


def _first_conserved_cys(seq: str, min_pos: int) -> int | None:
    """Position of the Cys in the first [YF]C at or after min_pos (1-based). Same rule as abcheck.annotate."""
    for m in CONSERVED_CYS_RE.finditer(seq, max(min_pos - 2, 0)):
        pos = m.start() + 2
        if pos >= min_pos:
            return pos
    return None


def _last_index(pattern: re.Pattern[str], seq: str) -> int | None:
    last = None
    for m in pattern.finditer(seq):
        last = m.start() + 1
    return last


def verify_heavy(seq: str, cdrs: dict, anchors: dict,
                 allow_motif_mismatch: bool = False) -> list[str]:
    """Check that the search rules uniquely recover the intended anchors. Returns problems as a list of strings."""
    problems: list[str] = []

    # The terminal WGxG must be the last WGxG match in the sequence; an FGxG after it would confuse chain detection
    if _last_index(WGXG_RE, seq) != anchors["terminal_motif"]:
        problems.append("말단 WGxG 가 마지막 매치가 아니다")
    last_fgxg = _last_index(FGXG_RE, seq)
    if last_fgxg is not None and last_fgxg > anchors["terminal_motif"]:
        problems.append("말단 WGxG 뒤에 FGxG 가 있어 체인 판별이 모호하다")

    # Cys22: must be unique in the 15..30 window
    window = [i + 1 for i, r in enumerate(seq[:30]) if r == "C" and i + 1 >= 15]
    if window != [anchors["cys1"]]:
        problems.append(f"15..30 창의 Cys 가 유일하지 않다: {window}")

    # Trp36: first FR2 context (W[VIA][RKQ]Q) after CDR-H1 start + minimum length
    h1_start = cdrs["CDR-H1"][0]
    first_w = _first_fr2_trp(seq, h1_start + MIN_CDR1_LENGTH_H, FR2_CONTEXT_H)
    if first_w != anchors["trp"]:
        problems.append(f"CDR-H1 이후 첫 FR2 문맥 Trp 위치가 다르다: {first_w} != {anchors['trp']}")

    # FR3 motif: first match after CDR-H2 start + 9.
    # Cases with a deliberately broken motif (reproducing hypermutation) must take the annotator's fallback path,
    # so there we check that "back-calculating 29 residues from the conserved Cys" lands on the intended position.
    h2_start = cdrs["CDR-H2"][0]
    search_from = h2_start + 9
    found = _first_index(FR3_MOTIF_RE, seq, search_from)

    if found is None:
        cys_fb = _first_conserved_cys(seq, search_from + FR3_TO_CYS_GAP)
        if cys_fb is None or cys_fb - FR3_TO_CYS_GAP != anchors["fr3_motif"]:
            problems.append(
                f"모티프도 없고 Cys 역산도 어긋난다: {cys_fb} - {FR3_TO_CYS_GAP} "
                f"!= {anchors['fr3_motif']}"
            )
    elif found != anchors["fr3_motif"] and not allow_motif_mismatch:
        # Normally this is an accidental match inside CDR-H2, so redraw.
        # But "fake motif" cases are built that way on purpose, so allow it --
        # those cases check that the annotator corrects it via the Cys distance.
        problems.append(f"FR3 모티프 첫 매치가 CDR-H2 안에 걸렸다: {found} != {anchors['fr3_motif']}")

    # Cys92: first [YF]C after the FR3 motif, and its distance matches reality
    first_c = _first_conserved_cys(seq, anchors["fr3_motif"])
    if first_c != anchors["cys2"]:
        problems.append(f"FR3 모티프 이후 첫 [YF]C 위치가 다르다: {first_c} != {anchors['cys2']}")
    gap = anchors["cys2"] - anchors["fr3_motif"]
    if gap != FR3_TO_CYS_GAP:
        problems.append(f"FR3→보존Cys 거리가 실제({FR3_TO_CYS_GAP})와 다르다: {gap}")

    return problems


def verify_light(seq: str, cdrs: dict, anchors: dict) -> list[str]:
    problems: list[str] = []

    if _last_index(FGXG_RE, seq) != anchors["terminal_motif"]:
        problems.append("말단 FGxG 가 마지막 매치가 아니다")
    last_wgxg = _last_index(WGXG_RE, seq)
    if last_wgxg is not None and last_wgxg > anchors["terminal_motif"]:
        problems.append("말단 FGxG 뒤에 WGxG 가 있어 체인 판별이 모호하다")

    window = [i + 1 for i, r in enumerate(seq[:30]) if r == "C" and i + 1 >= 16]
    if window != [anchors["cys1"]]:
        problems.append(f"16..30 창의 Cys 가 유일하지 않다: {window}")

    l1_start = cdrs["CDR-L1"][0]
    first_w = _first_fr2_trp(seq, l1_start + MIN_CDR1_LENGTH_L, FR2_CONTEXT_L)
    if first_w != anchors["trp"]:
        problems.append(f"CDR-L1 이후 첫 FR2 문맥 Trp 위치가 다르다: {first_w} != {anchors['trp']}")

    l2_end = cdrs["CDR-L2"][1]
    first_c = _first_conserved_cys(seq, l2_end + 1)
    if first_c != anchors["cys2"]:
        problems.append(f"CDR-L2 이후 첫 [YF]C 위치가 다르다: {first_c} != {anchors['cys2']}")

    return problems


# --------------------------------------------------------------------------
# Case list
# --------------------------------------------------------------------------

def make_cases() -> list[Case]:
    return [
        Case("SYN-VH-001", "H", "정상 VH — 통상 길이 CDR, 주입한 liability 없음", (5, 16, 10)),
        Case("SYN-VH-002", "H", "정상 VH — CDR-H3 가 긴 편(18 잔기)", (5, 16, 18)),
        Case("SYN-VL-001", "L", "정상 VL — kappa 유사 길이 (L1 11, L3 9)", (11, 7, 9)),
        Case("SYN-VL-002", "L", "정상 VL — CDR-L1 이 짧은 편(8 잔기)", (8, 7, 9)),

        Case(
            "SYN-VH-003", "H",
            "CDR-H2 에 NG — deamidation 최고 위험, CDR 내부라 위험도 상향",
            (5, 16, 10),
            cdr_inserts={"H2": [(6, "NG")]},
            notes=["CDR-H2 내 NG"],
        ),
        Case(
            "SYN-VH-004", "H",
            "CDR-H1 에 N-glycosylation 서열(NLS) — CDR 내부 당화 위험",
            (6, 16, 10),
            cdr_inserts={"H1": [(1, "NLS")]},
            notes=["CDR-H1 내 N[^P][ST]"],
        ),
        Case(
            "SYN-VL-003", "L",
            "CDR-L3 에 NG — 경쇄 CDR 내 deamidation",
            (11, 7, 9),
            cdr_inserts={"L3": [(3, "NG")]},
            notes=["CDR-L3 내 NG"],
        ),
        Case(
            "SYN-VH-005", "H",
            "비결합 Cys — CDR-H3 에 Cys 를 하나 더 넣어 도메인 내 Cys 가 3개(홀수)",
            (5, 16, 12),
            cdr_inserts={"H3": [(5, "C")]},
            notes=["도메인 내 Cys 3개 → 홀수"],
        ),
        Case(
            "SYN-VH-006", "H",
            "framework 에만 liability — FR3 에 DG, FR1 에 Met. CDR 내부 건과 위험도 비교용",
            (5, 16, 10),
            fr_inserts={"fr3_rest": [(8, "DG")], "fr1": [(12, "M")]},
            notes=["framework 내 DG", "framework 내 Met"],
        ),
        Case(
            "SYN-VH-007", "H",
            "CDR-H3 에 DG — isomerization 최고 위험",
            (5, 16, 11),
            cdr_inserts={"H3": [(4, "DG")]},
            notes=["CDR-H3 내 DG"],
        ),
        Case(
            "SYN-VH-008", "H",
            "CDR-H3 에 Met 과 Trp — CDR 내 산화 위험",
            (5, 16, 12),
            cdr_inserts={"H3": [(3, "M"), (7, "W")]},
            notes=["CDR-H3 내 Met", "CDR-H3 내 Trp"],
        ),
        Case(
            "SYN-VL-004", "L",
            "CDR-L1 에 DP — fragmentation 위험",
            (11, 7, 9),
            cdr_inserts={"L1": [(4, "DP")]},
            notes=["CDR-L1 내 DP"],
        ),
        Case(
            "SYN-VL-005", "L",
            "CDR 에 Lys/Arg 를 몰아넣은 염기성 서열 — pI/net charge 계산 검증용",
            (11, 7, 9),
            cdr_inserts={"L1": [(2, "KRKRK")], "L3": [(1, "RKRK")]},
            notes=["염기성 CDR → 높은 pI"],
        ),

        Case(
            "SYN-VH-014", "H",
            "CDR-H2 안에 FR3 모티프처럼 보이는 조각(RVTIS)이 있는 VH — 모티프 탐색이 "
            "여기 먼저 걸려 CDR-H2 가 잘린다. 보존 Cys 까지의 거리(29)로 교차검증해 "
            "바로잡아야 한다",
            (5, 16, 10),
            cdr_inserts={"H2": [(10, "RVTIS")]},
            decoy_fr3_motif=True,
            fr3_motif_broken=True,
            notes=["CDR-H2 내 가짜 FR3 모티프"],
        ),
        Case(
            "SYN-VH-013", "H",
            "FR3 모티프가 체세포 과변이로 깨진 VH (EATLT — 첫 자리가 R/K/Q 가 아니다). "
            "모티프 탐색이 실패하므로 보존 Cys 에서 29잔기 역산하는 대체 경로를 타야 한다. "
            "위치는 정확히 복원하되 근거가 약하므로 confident=False 여야 한다",
            (5, 16, 10),
            fr3_motif="EATLT",
            fr3_motif_broken=True,
            notes=["FR3 모티프 깨짐 — Cys 역산 경로"],
        ),
        Case(
            "SYN-VL-007", "L",
            "FR3 안에 여분 Cys 가 있는 VL — 실제 서열에 있다(…LTIS<C>LQSEDFATYY<C>…). "
            "'CDR-L2 이후 첫 Cys' 를 앵커로 쓰면 CDR-L3 가 10잔기 넘게 앞에서 시작한다",
            (11, 7, 9),
            fr_inserts={"fr3_rest": [(6, "C")]},
            notes=["FR3 내 여분 Cys", "도메인 내 Cys 3개 → 홀수"],
        ),
        Case(
            "SYN-VH-011", "H",
            "FR3 모티프가 QVTIS 로 시작하는 VH — IGHV2 계열에서 흔하다. "
            "첫 자리가 R/K 가 아니라 Q 라 좁은 정규식으로는 못 찾는다",
            (5, 16, 10),
            fr3_motif="QVTIS",
            notes=["FR3 모티프 첫 자리 Q"],
        ),
        Case(
            "SYN-VH-012", "H",
            "FR3 모티프가 RITIN 으로 끝나는 VH — 끝 자리가 N 이라 좁은 정규식으로는 못 찾는다",
            (5, 16, 10),
            fr3_motif="RITIN",
            notes=["FR3 모티프 끝 자리 N"],
        ),
        Case(
            "SYN-VH-010", "H",
            "CDR-H2 에 NT — 약한 deamidation 등급. CDR 상향으로 '높음'이 되지만 "
            "최고 위험 모티프가 아니라 선별 규칙에는 걸리지 않아야 한다",
            (5, 16, 10),
            cdr_inserts={"H2": [(5, "NT")]},
            notes=["CDR-H2 내 NT (약한 등급)"],
        ),
        Case(
            "SYN-VH-009", "H",
            "CDR-H1 안에 Trp 이 있는 VH — 실제 항체에서 흔하다(예: SYWMQ). "
            "'CDR1 이후 첫 Trp'를 앵커로 쓰면 CDR-H1 이 잘려 나간다",
            (5, 16, 10),
            cdr_inserts={"H1": [(2, "W")]},
            notes=["CDR-H1 내 Trp"],
        ),
        Case(
            "SYN-VL-006", "L",
            "CDR-L1 안에 Trp 이 있는 VL — 중쇄와 같은 함정",
            (11, 7, 9),
            cdr_inserts={"L1": [(3, "W")]},
            notes=["CDR-L1 내 Trp"],
        ),

        # --- Cases whose anchors are broken so annotation must fail ---
        Case(
            "SYN-VH-901", "H",
            "말단 WGxG 가 깨진 VH — 체인 판별부터 실패해야 한다 (confident=False)",
            (5, 16, 10),
            breakages=[("terminal_motif", "S")],
            anchors_broken=True,
            notes=["말단 WGxG 의 Trp → Ser"],
        ),
        Case(
            "SYN-VH-902", "H",
            "보존 Cys22 가 Ser 로 치환된 VH — CDR-H1 앵커를 못 찾아야 한다 (confident=False)",
            (5, 16, 10),
            breakages=[("cys1", "S")],
            anchors_broken=True,
            notes=["Cys22 → Ser"],
        ),

        # --- Cases with intact anchors but lengths outside the typical range ---
        Case(
            "SYN-VH-903", "H",
            "CDR-H1 이 9 잔기 — 앵커는 찾히지만 통상 범위(5-7) 밖이라 confident=False",
            (9, 16, 10),
            length_outlier=True,
            notes=["CDR-H1 길이 9 (통상 5-7)"],
        ),
    ]


def _apply_breakages(seq: str, case: Case, anchors: dict) -> str:
    buf = list(seq)
    for anchor_name, new_res in case.breakages:
        pos = anchors[anchor_name]
        buf[pos - 1] = new_res
    return "".join(buf)


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

def build_record(case: Case, heavy: HeavyScaffold, light: LightScaffold) -> dict:
    """Turn one case into a record. Redraw CDR contents until the self-check passes."""
    last_problems: list[str] = []
    for attempt in range(20000):
        rng = random.Random(f"{case.seq_id}/{attempt}")
        if case.chain == "H":
            seq, cdrs, anchors = _assemble_heavy(case, heavy, rng)
            problems = verify_heavy(
                seq, cdrs, anchors,
                allow_motif_mismatch=case.decoy_fr3_motif,
            )
        else:
            seq, cdrs, anchors = _assemble_light(case, light, rng)
            problems = verify_light(seq, cdrs, anchors)
        # Redraw if unintended CDR liabilities slipped in.
        # This way each case's CDR liabilities are only those injected via case.cdr_inserts.
        for name, (start, end) in cdrs.items():
            short = name.split("-")[1]          # "CDR-H1" → "H1"
            stray = unintended_cdr_liabilities_in_context(
                seq, start, end, case.cdr_inserts.get(short, [])
            )
            if stray:
                problems.append(
                    f"{name} 에 의도하지 않은 liability: "
                    + ", ".join(f"{lab}@{s}" for s, _, lab in stray)
                )

        if problems:
            last_problems = problems
            continue

        intact_seq = seq
        final_seq = _apply_breakages(seq, case, anchors)

        expected_confident = not (
            case.anchors_broken or case.length_outlier or case.fr3_motif_broken
        )
        record = {
            "id": case.seq_id,
            "chain": case.chain,
            "description": case.description,
            "sequence": final_seq,
            "length": len(final_seq),
            "synthetic": True,
            "notes": case.notes,
            "expected": {
                # Cases with broken anchors define no ground-truth CDR positions.
                # The annotator must not guess; it must drop to confident=False.
                "cdrs": None if case.anchors_broken else {
                    name: {"start": s, "end": e, "seq": intact_seq[s - 1:e]}
                    for name, (s, e) in cdrs.items()
                },
                "chain": None if case.anchors_broken else case.chain,
                "confident": expected_confident,
                "anchors_broken": case.anchors_broken,
                "length_outlier": case.length_outlier,
                "fr3_motif_broken": case.fr3_motif_broken,
            },
            "anchors_intended": None if case.anchors_broken else anchors,
        }
        return record

    raise RuntimeError(
        f"{case.seq_id}: 자기 검증을 통과하는 서열을 만들지 못했다. 마지막 문제: {last_problems}"
    )


def build_corpus() -> dict:
    scaffold_rng = random.Random("ab-seq-selfcheck/scaffold/v1")
    heavy = build_heavy_scaffold(scaffold_rng)
    light = build_light_scaffold(scaffold_rng)

    records = [build_record(c, heavy, light) for c in make_cases()]
    return {
        "schema": "ab-seq-selfcheck/synthetic-corpus/v1",
        "generator": "scripts/make_synthetic_sequences.py",
        "warning": (
            "전부 합성 서열이다. 실제 치료용 항체나 germline DB 레코드를 복사하지 않았으며, "
            "주석 로직 회귀 테스트 용도로만 쓴다. 생물학적 해석 대상이 아니다."
        ),
        "coordinates": "CDR 위치는 1-based, 양끝 포함(start/end inclusive)",
        "typical_cdr_lengths": {k: list(v) for k, v in TYPICAL_LENGTHS.items()},
        "count": len(records),
        "sequences": records,
    }


def write_fasta(corpus: dict, path: Path) -> None:
    lines = []
    for rec in corpus["sequences"]:
        lines.append(f">{rec['id']} chain={rec['chain']} synthetic | {rec['description']}")
        seq = rec["sequence"]
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="합성 항체 서열 코퍼스를 생성한다.")
    parser.add_argument(
        "--out-dir", default="data/synthetic", type=Path,
        help="출력 디렉터리 (기본: data/synthetic)",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="파일을 쓰지 않고, 기존 sequences.json 과 내용이 같은지만 확인한다",
    )
    args = parser.parse_args()

    corpus = build_corpus()
    out_json = args.out_dir / "sequences.json"
    payload = json.dumps(corpus, ensure_ascii=False, indent=2) + "\n"

    if args.check:
        if not out_json.exists():
            print(f"없음: {out_json}")
            return 1
        same = out_json.read_text(encoding="utf-8") == payload
        print("일치" if same else "불일치 — 다시 생성이 필요하다")
        return 0 if same else 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_json.write_text(payload, encoding="utf-8")
    write_fasta(corpus, args.out_dir / "sequences.fasta")

    print(f"{out_json} 에 {corpus['count']} 개 서열을 썼다.")
    for rec in corpus["sequences"]:
        exp = rec["expected"]
        flag = "confident" if exp["confident"] else "NOT-confident"
        print(f"  {rec['id']}  {rec['chain']}  len={rec['length']:3d}  {flag:13s}  {rec['description']}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    from abcheck._streams import configure_streams

    configure_streams()      # Korean output must not crash a Windows pipe
    raise SystemExit(main())
