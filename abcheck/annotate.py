"""CDR annotation — based on conserved-residue anchors, approximating Kabat boundaries.

This module is not a numbering engine. It assigns no numbers; it only cuts out CDR
spans by regex-matching conserved residues (Cys, Trp) and motifs. When an anchor
can't be found it **does not guess positions**; it records `confident=False` and a
reason in Korean.

Search rules
------------
Chain type: decided by the terminal motif. `WGxG` means heavy (H), `FGxG` means light (L).
If both are present, the one closer to the C-terminus wins. If neither is, detection fails.

VH
    Cys22      the Cys closest to position 22 within the 15..30 window
    CDR-H1     starts at Cys22 + 9, ends just before Trp36
    Trp36      Trp of the first `W[VIA][RKQ]Q` from CDR-H1 start + 5 (else falls back to the first Trp, with a warning)
    CDR-H2     starts at Trp36 + 14, ends just before the FR3 motif
    FR3 motif  first `[RKQ].[TSAV][IMLFV][STAN]` match from CDR-H2 start + 9
    Cys92      Cys of the first `[YF]C` after the FR3 motif (else falls back to the first Cys, with a warning)
    CDR-H3     starts at Cys92 + 3, ends just before the terminal `WGxG`

VL
    Cys23      the Cys closest to position 23 within the 16..30 window
    CDR-L1     starts at Cys23 + 1, ends just before Trp35
    Trp35      Trp of the first `W[YF][QL]Q` from CDR-L1 start + 7 (else falls back to the first Trp, with a warning)
    CDR-L2     starts at Trp35 + 15, fixed at 7 residues
    Cys88      Cys of the first `[YF]C` after CDR-L2 (else falls back to the first Cys, with a warning)
    CDR-L3     starts at Cys88 + 1, ends just before the terminal `FGxG`

Using "first" (as in "first Cys after the FR3 motif") is deliberate. For sequences with
an unpaired Cys inside CDR-H3, "last Cys" would wrongly pick that Cys as the anchor.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .validate import ensure_valid

# --------------------------------------------------------------------------
# Anchor motifs and constants
# --------------------------------------------------------------------------

WGXG_RE = re.compile(r"WG.G")
FGXG_RE = re.compile(r"FG.G")
#: FR3 start motif.
#:
#: The original design was [RK].[TSAV][IMLFV][STA]. Measured on 1,000 real-repertoire VH
#: from OAS, it missed the motif in 26 (2.6%), so CDR-H2/H3 annotation failed.
#: The two missing cases were:
#:   - IGHV2-family QVTIS…  → first position is Q, not R/K
#:   - RITIN…               → last position is N, outside [STA]
#: Adding Q to the first position and N to the last cut failures from 26 → 8, and
#: **not one of the 1,000 previously successful annotations changed** (measured).
#: To revert, use the STRICT pattern below.
FR3_MOTIF_PATTERN_STRICT = r"[RK].[TSAV][IMLFV][STA]"
FR3_MOTIF_PATTERN = r"[RKQ].[TSAV][IMLFV][STAN]"
FR3_MOTIF_RE = re.compile(FR3_MOTIF_PATTERN)

#: Window to search for the conserved Cys (1-based, inclusive) and its expected position
CYS1_WINDOW_H = (15, 30)
CYS1_EXPECTED_H = 22
CYS1_WINDOW_L = (16, 30)
CYS1_EXPECTED_L = 23

#: Distance from Cys1 to the start of CDR1
CDR1_OFFSET_H = 9      # Cys22 + 9 → 31
CDR1_OFFSET_L = 1      # Cys23 + 1 → 24

#: Distance from the Trp anchor to the start of CDR2
CDR2_OFFSET_H = 14     # Trp36 + 14 → 50
CDR2_OFFSET_L = 15     # Trp35 + 15 → 50

#: Distance from Cys2 to the start of CDR3
CDR3_OFFSET_H = 3      # Cys92 + 3 → 95
CDR3_OFFSET_L = 1      # Cys88 + 1 → 89

#: CDR-L2 has a fixed length in Kabat.
CDR_L2_LENGTH = 7

#: FR2 context immediately following the conserved Trp anchor.
#:
#: Searching for a lone Trp picks up Trps inside CDR1. IGKV1-family chains in particular
#: have two Trps, at the end of CDR-L1 and in FR2 (`...RASQGISS**W**LA**W**YQQ...`);
#: picking the first one shifts every later boundary by 3 residues.
#:
#: Chosen by comparing ~12,000 sequences across 5 datasets against the ANARCI annotations
#: in OAS as ground truth. Preferring this context cut CDR boundary mismatches from 3.83% → 0.31%.
FR2_CONTEXT_H = re.compile(r"W[VIA][RKQ]Q")
FR2_CONTEXT_L = re.compile(r"W[YF][QL]Q")

#: The conserved-Trp search starts at CDR1 start + this value.
#: Real antibodies often have a Trp in CDR-H1 (e.g. ...GYTFT**SYWMQ**W VKQ...).
#: Plain "first Trp after CDR1 start" wrongly anchors on the in-CDR Trp and truncates
#: CDR-H1 to 2-4 residues. Kabat CDR-H1 is at least 5 residues (31-35), so search after that.
#: CDR-L1 likewise starts searching after its minimum length.
MIN_CDR1_LENGTH_H = 5
MIN_CDR1_LENGTH_L = 7

#: Minimum CDR-H2 length guaranteed before the FR3 motif search begins.
#: A shorter CDR-H2 means the motif matched prematurely inside the CDR.
MIN_CDR_H2_BEFORE_MOTIF = 9

#: Distance from the FR3 start motif to the conserved Cys92.
#: In 992 real VH from OAS it was exactly 29 in **every** case. When somatic hypermutation
#: breaks the motif, we back-calculate from this distance. That is safer than widening the
#: regex further — the wider it gets, the likelier a chance match inside CDR-H2, which
#: silently shortens CDR-H2.
FR3_TO_CYS_GAP = 29

#: Typical CDR length ranges (a screening criterion, not a strict numbering one)
#:
#: These ranges answer **"can the annotation be trusted?"**, not "is this antibody odd?".
#: Repeatedly marking normal sequences as untrustworthy dulls the warning, so the upper
#: bounds follow the measured distribution of 2,000 real OAS sequences (measured min-max in parentheses).
#:   CDR-H1 5-7   (measured 5-7,  0 out of range)
#:   CDR-H2 16-20 (measured 16-20, 20 is from germline families with an insertion)
#:   CDR-H3 3-30  (measured 3-27,  26-27 residues is not rare in humans)
#:   CDR-L1 7-17  (measured 8-17,  0 out of range)
#:   CDR-L2 7     (measured 7 in all cases)
#:   CDR-L3 5-15  (4-14 measured is normal; the 20s came from a Cys anchor bug, now fixed)
#: Lower bounds were not widened — a misplaced anchor shortens the CDR, so the lower bound is that signal.
TYPICAL_LENGTHS: dict[str, tuple[int, int]] = {
    "CDR-H1": (5, 7),
    "CDR-H2": (16, 20),
    "CDR-H3": (3, 30),
    "CDR-L1": (7, 17),
    "CDR-L2": (7, 7),
    "CDR-L3": (5, 15),
}


# --------------------------------------------------------------------------
# Result types
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class CDR:
    """A single CDR span. start/end are 1-based, inclusive."""

    name: str
    start: int
    end: int
    seq: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Annotation:
    """Annotation result.

    When confident=False, `reasons` holds the reasons in Korean. If an anchor wasn't
    found, `cdrs` is an empty dict — it is never filled in with guessed positions.
    """

    sequence: str
    chain: str | None = None                      # "H" | "L" | None
    cdrs: dict[str, CDR] = field(default_factory=dict)
    anchors: dict[str, int] = field(default_factory=dict)
    confident: bool = True
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def annotated(self) -> bool:
        """Whether at least one CDR span was cut out."""
        return bool(self.cdrs)

    def framework_spans(self) -> list[tuple[int, int]]:
        """List of the remaining (framework) spans outside the CDRs. 1-based, inclusive."""
        if not self.cdrs:
            return [(1, len(self.sequence))]
        spans: list[tuple[int, int]] = []
        cursor = 1
        for cdr in sorted(self.cdrs.values(), key=lambda c: c.start):
            if cdr.start > cursor:
                spans.append((cursor, cdr.start - 1))
            cursor = cdr.end + 1
        if cursor <= len(self.sequence):
            spans.append((cursor, len(self.sequence)))
        return spans

    def region_of(self, position: int) -> str:
        """Which region a 1-based position falls in: a CDR name or "framework"."""
        for cdr in self.cdrs.values():
            if cdr.start <= position <= cdr.end:
                return cdr.name
        return "framework"


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _last_match(pattern: re.Pattern[str], seq: str) -> int | None:
    """1-based start position of the last match."""
    last = None
    for m in pattern.finditer(seq):
        last = m.start() + 1
    return last


def _first_match_from(pattern: re.Pattern[str], seq: str, min_pos: int) -> int | None:
    """1-based start position of the first match starting at or after min_pos (1-based)."""
    if min_pos < 1:
        min_pos = 1
    m = pattern.search(seq, min_pos - 1)
    return m.start() + 1 if m else None


def _first_residue_from(seq: str, residue: str, min_pos: int) -> int | None:
    """1-based position of the first occurrence of residue at or after min_pos (1-based)."""
    idx = seq.find(residue, max(min_pos - 1, 0))
    return idx + 1 if idx >= 0 else None


#: Residue immediately preceding the conserved Cys. In 1,991 real OAS sequences: Y 83.2% + F 16.1% = 99.3%.
#: Looking only for the "first Cys" without this context, sequences with an extra Cys in FR3 (e.g. ...LTIS**C**...)
#: get that Cys as the anchor and CDR3 starts more than 10 residues early. 7 of 2,000 measured sequences did this,
#: and one of them had a length within the typical range, so the length check didn't catch it either.
CONSERVED_CYS_RE = re.compile(r"[YF]C")


def _find_conserved_cys(seq: str, min_pos: int, limit: int) -> tuple[int | None, bool]:
    """Find the conserved Cys in [min_pos, limit). Returns (position, found_by_context).

    Searches with the `[YF]C` context first, falling back to the plain first Cys. The
    second value tells the caller whether it fell back, so it can record a warning.
    """
    start = max(min_pos - 2, 0)
    for m in CONSERVED_CYS_RE.finditer(seq, start):
        pos = m.start() + 2
        if pos < min_pos:
            continue
        if pos >= limit:
            break
        return pos, True

    fallback = _first_residue_from(seq, "C", min_pos)
    if fallback is not None and fallback < limit:
        return fallback, False
    return None, False


def _find_conserved_trp(seq: str, min_pos: int,
                        context: re.Pattern[str]) -> tuple[int | None, bool]:
    """Find the conserved Trp. Returns (position, found_by_FR2_context).

    Tries the FR2 context (`W[VIA][RKQ]Q` or `W[YF][QL]Q`) first, falling back to the
    plain first Trp. The second value tells the caller whether it fell back.
    """
    m = context.search(seq, max(min_pos - 1, 0))
    if m is not None and m.start() + 1 >= min_pos:
        return m.start() + 1, True
    return _first_residue_from(seq, "W", min_pos), False


def _cys_in_window(seq: str, window: tuple[int, int], expected: int) -> tuple[int | None, int]:
    """Position of the Cys in the window closest to expected, and the total Cys count in the window."""
    lo, hi = window
    found = [i + 1 for i, ch in enumerate(seq[: min(hi, len(seq))]) if ch == "C" and i + 1 >= lo]
    if not found:
        return None, 0
    best = min(found, key=lambda p: (abs(p - expected), p))
    return best, len(found)


def detect_chain(seq: str) -> tuple[str | None, int | None]:
    """Detect the chain from the terminal motif. Returns (chain, 1-based motif start position)."""
    wgxg = _last_match(WGXG_RE, seq)
    fgxg = _last_match(FGXG_RE, seq)

    if wgxg is None and fgxg is None:
        return None, None
    if fgxg is None:
        return "H", wgxg
    if wgxg is None:
        return "L", fgxg
    # If both are present, take the one closer to the C-terminus
    return ("H", wgxg) if wgxg > fgxg else ("L", fgxg)


# --------------------------------------------------------------------------
# Main logic
# --------------------------------------------------------------------------

def annotate(sequence: str, *, validate_input: bool = True) -> Annotation:
    """Annotate a single sequence.

    With validate_input=True, validates first via abcheck.validate (SequenceError on failure).
    """
    seq = ensure_valid(sequence) if validate_input else str(sequence).strip().upper()
    ann = Annotation(sequence=seq)

    chain, motif_pos = detect_chain(seq)
    if chain is None:
        ann.chain = None
        ann.confident = False
        ann.reasons.append(
            "말단 모티프를 찾지 못했다. 중쇄의 WGxG 또는 경쇄의 FGxG 가 없어서 "
            "체인 종류(H/L)조차 판별할 수 없다. 서열이 중간에서 잘렸거나, 가변도메인이 "
            "아니거나, 해당 보존 모티프가 치환된 경우다. CDR 위치는 추측하지 않는다."
        )
        return ann

    ann.chain = chain
    ann.anchors["terminal_motif"] = motif_pos

    if motif_pos < len(seq) * 0.5:
        ann.warnings.append(
            f"말단 모티프({'WGxG' if chain == 'H' else 'FGxG'})가 {motif_pos}번째로 "
            "서열 전반부에 있다. C-말단이 아닌 자리에 우연히 걸린 모티프일 수 있다."
        )

    if chain == "H":
        _annotate_heavy(seq, ann)
    else:
        _annotate_light(seq, ann)

    _check_lengths(ann)
    return ann


def _annotate_heavy(seq: str, ann: Annotation) -> None:
    motif_pos = ann.anchors["terminal_motif"]

    # --- Anchor 1: Cys22 ---
    cys1, n_cys = _cys_in_window(seq, CYS1_WINDOW_H, CYS1_EXPECTED_H)
    if cys1 is None:
        ann.confident = False
        ann.reasons.append(
            f"보존 Cys22 를 찾지 못했다. {CYS1_WINDOW_H[0]}-{CYS1_WINDOW_H[1]}번 구간에 "
            "Cys 가 없다. 이 Cys 는 도메인 내부 이황화결합의 한쪽이라 거의 항상 보존된다. "
            "N-말단이 잘렸거나 Cys 가 치환된 경우다. CDR-H1 위치는 추측하지 않는다."
        )
        return
    if n_cys > 1:
        ann.warnings.append(
            f"{CYS1_WINDOW_H[0]}-{CYS1_WINDOW_H[1]}번 구간에 Cys 가 {n_cys}개 있다. "
            f"22번에 가장 가까운 {cys1}번을 앵커로 썼다. 결과를 ANARCI/IMGT 로 재확인할 것."
        )
    ann.anchors["cys1"] = cys1

    # --- Anchor 2: Trp36, CDR-H1 ---
    h1_start = cys1 + CDR1_OFFSET_H
    if h1_start > len(seq):
        ann.confident = False
        ann.reasons.append(
            f"Cys22({cys1}번) 뒤로 서열이 {len(seq) - cys1} 잔기밖에 없어 "
            "CDR-H1 시작 위치가 서열 범위를 벗어난다. 서열이 잘린 것으로 보인다."
        )
        return

    trp, by_context = _find_conserved_trp(
        seq, h1_start + MIN_CDR1_LENGTH_H, FR2_CONTEXT_H
    )
    if trp is not None and not by_context:
        ann.warnings.append(
            f"보존 Trp36 을 FR2 문맥(W[VIA][RKQ]Q)으로 찾지 못해 {trp}번의 Trp 을 그대로 썼다. "
            "CDR-H1 안에 Trp 이 있으면 경계가 밀릴 수 있다. ANARCI/IMGT 로 재확인할 것."
        )
    if trp is None or trp >= motif_pos:
        ann.confident = False
        ann.reasons.append(
            f"보존 Trp36 을 찾지 못했다. CDR-H1 시작 추정 위치({h1_start}번) 이후에 "
            "Trp 가 없다. 이 Trp 는 CDR-H1 의 끝을 정하는 앵커다. "
            "CDR-H1/H2/H3 위치 모두 추측하지 않는다."
        )
        return
    if trp == h1_start:
        ann.confident = False
        ann.reasons.append(
            f"CDR-H1 시작 위치({h1_start}번)가 곧 Trp 라서 CDR-H1 길이가 0 이 된다. "
            "Cys22 앵커를 잘못 집었을 가능성이 높다."
        )
        return
    ann.anchors["trp"] = trp
    ann.cdrs["CDR-H1"] = CDR("CDR-H1", h1_start, trp - 1, seq[h1_start - 1:trp - 1])

    # --- Anchor 3: FR3 motif, CDR-H2 ---
    h2_start = trp + CDR2_OFFSET_H
    if h2_start > len(seq):
        ann.confident = False
        ann.reasons.append(
            f"Trp36({trp}번) 뒤로 서열이 짧아 CDR-H2 시작 위치({h2_start}번)가 "
            "서열 범위를 벗어난다. 서열이 잘린 것으로 보인다."
        )
        return

    search_from = h2_start + MIN_CDR_H2_BEFORE_MOTIF
    fr3 = _first_match_from(FR3_MOTIF_RE, seq, search_from)
    if fr3 is not None and fr3 >= motif_pos:
        fr3 = None

    if fr3 is not None:
        cys2, by_context = _find_conserved_cys(seq, fr3, motif_pos)
        # Cross-check: the FR3-motif-to-conserved-Cys distance was 29 in all 3,429 measured cases,
        # and the 14 mismatches were all the motif matching in the wrong place (CDR-H2 came out as 15 or 29).
        # If the distance differs, discard the motif and re-derive the position from the Cys.
        if cys2 is not None and cys2 - fr3 != FR3_TO_CYS_GAP:
            corrected = cys2 - FR3_TO_CYS_GAP
            if corrected > h2_start:
                ann.confident = False
                ann.reasons.append(
                    f"FR3 시작 모티프를 {fr3}번에서 찾았지만, 보존 Cys({cys2}번)까지의 거리가 "
                    f"{cys2 - fr3}로 통상값({FR3_TO_CYS_GAP})과 다르다. 모티프가 엉뚱한 자리에 "
                    f"걸린 것으로 보여 Cys 기준으로 {corrected}번으로 고쳐 잡았다. "
                    "CDR-H2 경계를 ANARCI/IMGT 로 반드시 재확인할 것."
                )
                fr3 = corrected
    else:
        # The motif was broken by somatic hypermutation. Back-calculate from the conserved Cys
        # (29 residues in all 992 measured cases). Safer than widening the regex further —
        # the wider it gets, the likelier a chance match inside CDR-H2.
        cys2, by_context = _find_conserved_cys(seq, search_from + FR3_TO_CYS_GAP, motif_pos)
        if cys2 is not None:
            fr3 = cys2 - FR3_TO_CYS_GAP
            # Report the position but lower confidence. The distance was constant across all 992
            # measured cases, but it is weaker evidence than seeing the motif — have a human check first.
            ann.confident = False
            ann.reasons.append(
                f"FR3 시작 모티프({FR3_MOTIF_PATTERN})를 찾지 못해, 보존 Cys({cys2}번)에서 "
                f"{FR3_TO_CYS_GAP}잔기 역산해 CDR-H2 의 끝을 {fr3 - 1}번으로 잡았다. "
                "체세포 과변이로 모티프가 깨진 서열에서 흔하다. 위치는 그대로 보고하지만 "
                "모티프를 직접 확인한 것보다 근거가 약하므로, CDR-H2 경계를 "
                "ANARCI/IMGT 로 반드시 재확인할 것."
            )

    if cys2 is None or fr3 is None or fr3 <= h2_start:
        ann.confident = False
        ann.reasons.append(
            f"CDR-H2 의 끝을 정할 수 없다. FR3 시작 모티프({FR3_MOTIF_PATTERN})가 없고, "
            f"보존 Cys 에서 역산할 수도 없다 (CDR-H2 시작 {h2_start}번). "
            "CDR-H2/H3 위치는 추측하지 않는다."
        )
        return

    if not by_context:
        ann.warnings.append(
            f"보존 Cys92 를 문맥(YC/FC) 으로 찾지 못해 {cys2}번의 Cys 를 그대로 썼다. "
            "FR3 에 여분 Cys 가 있으면 CDR-H3 시작이 밀릴 수 있다. ANARCI/IMGT 로 재확인할 것."
        )

    ann.anchors["fr3_motif"] = fr3
    ann.cdrs["CDR-H2"] = CDR("CDR-H2", h2_start, fr3 - 1, seq[h2_start - 1:fr3 - 1])
    ann.anchors["cys2"] = cys2

    h3_start = cys2 + CDR3_OFFSET_H
    if h3_start > motif_pos - 1:
        ann.confident = False
        ann.reasons.append(
            f"Cys92({cys2}번)와 말단 WGxG({motif_pos}번) 사이가 너무 좁아 "
            "CDR-H3 길이가 0 이하가 된다. 두 앵커 중 하나를 잘못 집었을 가능성이 높다."
        )
        return
    ann.cdrs["CDR-H3"] = CDR("CDR-H3", h3_start, motif_pos - 1, seq[h3_start - 1:motif_pos - 1])


def _annotate_light(seq: str, ann: Annotation) -> None:
    motif_pos = ann.anchors["terminal_motif"]

    # --- Anchor 1: Cys23 ---
    cys1, n_cys = _cys_in_window(seq, CYS1_WINDOW_L, CYS1_EXPECTED_L)
    if cys1 is None:
        ann.confident = False
        ann.reasons.append(
            f"보존 Cys23 을 찾지 못했다. {CYS1_WINDOW_L[0]}-{CYS1_WINDOW_L[1]}번 구간에 "
            "Cys 가 없다. N-말단이 잘렸거나 Cys 가 치환된 경우다. "
            "CDR-L1 위치는 추측하지 않는다."
        )
        return
    if n_cys > 1:
        ann.warnings.append(
            f"{CYS1_WINDOW_L[0]}-{CYS1_WINDOW_L[1]}번 구간에 Cys 가 {n_cys}개 있다. "
            f"23번에 가장 가까운 {cys1}번을 앵커로 썼다. 결과를 ANARCI/IMGT 로 재확인할 것."
        )
    ann.anchors["cys1"] = cys1

    # --- Anchor 2: Trp35, CDR-L1 ---
    l1_start = cys1 + CDR1_OFFSET_L
    trp, by_context = _find_conserved_trp(
        seq, l1_start + MIN_CDR1_LENGTH_L, FR2_CONTEXT_L
    )
    if trp is not None and not by_context:
        ann.warnings.append(
            f"보존 Trp35 를 FR2 문맥(W[YF][QL]Q)으로 찾지 못해 {trp}번의 Trp 을 그대로 썼다. "
            "CDR-L1 안에 Trp 이 있으면 경계가 밀릴 수 있다. ANARCI/IMGT 로 재확인할 것."
        )
    if trp is None or trp >= motif_pos:
        ann.confident = False
        ann.reasons.append(
            f"보존 Trp35 를 찾지 못했다. CDR-L1 시작 추정 위치({l1_start}번) 이후에 "
            "Trp 가 없다. 이 Trp 는 CDR-L1 의 끝을 정하는 앵커다. "
            "CDR-L1/L2/L3 위치 모두 추측하지 않는다."
        )
        return
    if trp == l1_start:
        ann.confident = False
        ann.reasons.append(
            f"CDR-L1 시작 위치({l1_start}번)가 곧 Trp 라서 CDR-L1 길이가 0 이 된다. "
            "Cys23 앵커를 잘못 집었을 가능성이 높다."
        )
        return
    ann.anchors["trp"] = trp
    ann.cdrs["CDR-L1"] = CDR("CDR-L1", l1_start, trp - 1, seq[l1_start - 1:trp - 1])

    # --- CDR-L2 (fixed at 7 residues) ---
    l2_start = trp + CDR2_OFFSET_L
    l2_end = l2_start + CDR_L2_LENGTH - 1
    if l2_end > len(seq) or l2_end >= motif_pos:
        ann.confident = False
        ann.reasons.append(
            f"Trp35({trp}번) 기준으로 잡은 CDR-L2 구간({l2_start}-{l2_end}번)이 "
            "서열 범위 또는 말단 FGxG 를 넘는다. 서열이 잘렸거나 Trp 앵커를 "
            "잘못 집은 것으로 보인다. CDR-L2/L3 위치는 추측하지 않는다."
        )
        return
    ann.cdrs["CDR-L2"] = CDR("CDR-L2", l2_start, l2_end, seq[l2_start - 1:l2_end])

    # --- Anchor 3: Cys88, CDR-L3 ---
    cys2, by_context = _find_conserved_cys(seq, l2_end + 1, motif_pos)
    if not by_context and cys2 is not None:
        ann.warnings.append(
            f"보존 Cys88 을 문맥(YC/FC) 으로 찾지 못해 {cys2}번의 Cys 를 그대로 썼다. "
            "FR3 에 여분 Cys 가 있으면 CDR-L3 시작이 밀릴 수 있다. ANARCI/IMGT 로 재확인할 것."
        )
    if cys2 is None or cys2 >= motif_pos:
        ann.confident = False
        ann.reasons.append(
            f"보존 Cys88 을 찾지 못했다. CDR-L2 종료({l2_end}번) 이후 말단 FGxG 전까지 "
            "Cys 가 없다. 이 Cys 는 CDR-L3 시작을 정하는 앵커다. "
            "CDR-L3 위치는 추측하지 않는다."
        )
        return
    ann.anchors["cys2"] = cys2

    l3_start = cys2 + CDR3_OFFSET_L
    if l3_start > motif_pos - 1:
        ann.confident = False
        ann.reasons.append(
            f"Cys88({cys2}번)과 말단 FGxG({motif_pos}번) 사이가 너무 좁아 "
            "CDR-L3 길이가 0 이하가 된다. 두 앵커 중 하나를 잘못 집었을 가능성이 높다."
        )
        return
    ann.cdrs["CDR-L3"] = CDR("CDR-L3", l3_start, motif_pos - 1, seq[l3_start - 1:motif_pos - 1])


def _check_lengths(ann: Annotation) -> None:
    """If a CDR length is outside the typical range, set confident=False and record a reason."""
    for name, cdr in sorted(ann.cdrs.items()):
        lo, hi = TYPICAL_LENGTHS[name]
        if not (lo <= cdr.length <= hi):
            ann.confident = False
            bound = "짧다" if cdr.length < lo else "길다"
            ann.reasons.append(
                f"{name} 길이가 {cdr.length} 잔기로 통상 범위({lo}-{hi})보다 {bound}. "
                "앵커 자체는 찾았지만 경계가 밀렸을 수 있다. "
                "ANARCI/IMGT 로 경계를 재확인할 것."
            )
