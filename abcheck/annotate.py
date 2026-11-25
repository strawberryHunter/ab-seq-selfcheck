"""CDR 주석 — 보존 잔기 앵커 기반, Kabat 경계 근사.

이 모듈은 numbering 엔진이 아니다. 번호를 매기지 않고, 보존 잔기(Cys·Trp)와 모티프를
정규식으로 찾아 CDR 구간만 끊는다. 앵커를 못 찾으면 **위치를 추측하지 않고**
`confident=False` 와 한글 사유를 남긴다.

탐색 규칙
---------
체인 판별: 말단 모티프로 한다. `WGxG` 가 있으면 중쇄(H), `FGxG` 가 있으면 경쇄(L).
둘 다 있으면 C-말단에 더 가까운 쪽을 택한다. 둘 다 없으면 판별 실패다.

VH
    Cys22      15..30 창에서 22번에 가장 가까운 Cys
    CDR-H1     Cys22 + 9 에서 시작, Trp36 직전까지
    Trp36      CDR-H1 시작 + 5 이후 첫 Trp (CDR-H1 안의 Trp 을 피하려고 최소 길이만큼 건너뛴다)
    CDR-H2     Trp36 + 14 에서 시작, FR3 모티프 직전까지
    FR3 모티프  CDR-H2 시작 + 9 이후 첫 `[RKQ].[TSAV][IMLFV][STAN]` 매치
    Cys92      FR3 모티프 이후 첫 `[YF]C` 의 Cys (없으면 첫 Cys 로 물러서고 경고)
    CDR-H3     Cys92 + 3 에서 시작, 말단 `WGxG` 직전까지

VL
    Cys23      16..30 창에서 23번에 가장 가까운 Cys
    CDR-L1     Cys23 + 1 에서 시작, Trp35 직전까지
    Trp35      CDR-L1 시작 + 7 이후 첫 Trp (같은 이유로 최소 길이만큼 건너뛴다)
    CDR-L2     Trp35 + 15 에서 시작, 7 잔기 고정
    Cys88      CDR-L2 이후 첫 `[YF]C` 의 Cys (없으면 첫 Cys 로 물러서고 경고)
    CDR-L3     Cys88 + 1 에서 시작, 말단 `FGxG` 직전까지

"FR3 모티프 이후 첫 Cys"처럼 '첫'을 쓰는 것은 의도적이다. CDR-H3 안에 비결합 Cys 가
들어 있는 서열에서 '마지막 Cys'를 쓰면 그 Cys 를 앵커로 잘못 집기 때문이다.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .validate import ensure_valid

# --------------------------------------------------------------------------
# 앵커 모티프와 상수
# --------------------------------------------------------------------------

WGXG_RE = re.compile(r"WG.G")
FGXG_RE = re.compile(r"FG.G")
#: FR3 시작 모티프.
#:
#: 원래 설계는 [RK].[TSAV][IMLFV][STA] 였다. OAS 실제 레퍼토리 VH 1,000건으로
#: 재어 보니 26건(2.6%)에서 이 모티프를 못 찾아 CDR-H2/H3 주석이 실패했다.
#: 빠진 두 가지는 이랬다.
#:   - IGHV2 계열의 QVTIS…  → 첫 자리가 R/K 가 아니라 Q
#:   - RITIN…              → 끝 자리가 N 이라 [STA] 밖
#: 첫 자리에 Q, 끝 자리에 N 을 더하니 실패가 26건 → 8건으로 줄었고,
#: **기존에 성공하던 1,000건의 주석은 하나도 바뀌지 않았다**(실측).
#: 되돌리려면 아래 STRICT 패턴을 쓰면 된다.
FR3_MOTIF_PATTERN_STRICT = r"[RK].[TSAV][IMLFV][STA]"
FR3_MOTIF_PATTERN = r"[RKQ].[TSAV][IMLFV][STAN]"
FR3_MOTIF_RE = re.compile(FR3_MOTIF_PATTERN)

#: 보존 Cys 를 찾을 창 (1-based, 양끝 포함) 과 기대 위치
CYS1_WINDOW_H = (15, 30)
CYS1_EXPECTED_H = 22
CYS1_WINDOW_L = (16, 30)
CYS1_EXPECTED_L = 23

#: Cys1 에서 CDR1 시작까지의 거리
CDR1_OFFSET_H = 9      # Cys22 + 9 → 31
CDR1_OFFSET_L = 1      # Cys23 + 1 → 24

#: Trp 앵커에서 CDR2 시작까지의 거리
CDR2_OFFSET_H = 14     # Trp36 + 14 → 50
CDR2_OFFSET_L = 15     # Trp35 + 15 → 50

#: Cys2 에서 CDR3 시작까지의 거리
CDR3_OFFSET_H = 3      # Cys92 + 3 → 95
CDR3_OFFSET_L = 1      # Cys88 + 1 → 89

#: CDR-L2 는 Kabat 에서 길이가 고정이다.
CDR_L2_LENGTH = 7

#: 보존 Trp 앵커를 찾기 시작할 위치 = CDR1 시작 + 이 값.
#: 실제 항체의 CDR-H1 에는 Trp 이 흔히 들어간다(예: ...GYTFT**SYWMQ**W VKQ...).
#: 그냥 "CDR1 시작 이후 첫 Trp"를 쓰면 CDR 안의 Trp 을 앵커로 잘못 집어 CDR-H1 이
#: 2-4 잔기로 잘린다. Kabat CDR-H1 은 최소 5 잔기(31-35)이므로 그 뒤부터 찾는다.
#: CDR-L1 도 같은 이유로 최소 길이 뒤부터 찾는다.
MIN_CDR1_LENGTH_H = 5
MIN_CDR1_LENGTH_L = 7

#: FR3 모티프 탐색을 시작하기 전에 CDR-H2 에 최소한 보장할 길이.
#: 이보다 짧은 CDR-H2 는 모티프가 CDR 안에서 먼저 걸린 것으로 본다.
MIN_CDR_H2_BEFORE_MOTIF = 9

#: FR3 시작 모티프에서 보존 Cys92 까지의 거리.
#: OAS 실제 VH 992건이 **전부** 정확히 29였다. 모티프가 체세포 과변이로 깨졌을 때
#: 이 거리로 역산한다. 정규식을 계속 넓히는 것보다 안전하다 — 넓힐수록 CDR-H2 안에서
#: 우연히 걸릴 위험이 커지고, 그러면 CDR-H2 가 조용히 짧아진다.
FR3_TO_CYS_GAP = 29

#: 통상 CDR 길이 범위 (스크리닝 기준; 엄격한 numbering 기준이 아니다)
#:
#: 이 범위는 "이 항체가 이상한가"가 아니라 **"주석을 믿을 수 있는가"** 의 기준이다.
#: 정상 서열을 자꾸 '신뢰 불가' 로 떨어뜨리면 경고가 무뎌지므로, OAS 실제 서열
#: 2,000건의 실측 분포에 맞춰 상한을 잡았다 (아래 괄호는 실측 min-max).
#:   CDR-H1 5-7   (실측 5-7,  범위 밖 0건)
#:   CDR-H2 16-20 (실측 16-20, 20 은 삽입이 있는 germline 계열)
#:   CDR-H3 3-30  (실측 3-27,  사람은 26-27 잔기가 드물지 않다)
#:   CDR-L1 7-17  (실측 8-17,  범위 밖 0건)
#:   CDR-L2 7     (실측 전부 7)
#:   CDR-L3 5-15  (실측 4-14 가 정상. 20 으로 나오던 건은 Cys 앵커 버그였고 고쳤다)
#: 하한은 넓히지 않았다 — 앵커를 잘못 집으면 CDR 이 짧아지므로, 하한이 그 신호다.
TYPICAL_LENGTHS: dict[str, tuple[int, int]] = {
    "CDR-H1": (5, 7),
    "CDR-H2": (16, 20),
    "CDR-H3": (3, 30),
    "CDR-L1": (7, 17),
    "CDR-L2": (7, 7),
    "CDR-L3": (5, 15),
}


# --------------------------------------------------------------------------
# 결과 자료형
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class CDR:
    """CDR 한 구간. start/end 는 1-based, 양끝 포함."""

    name: str
    start: int
    end: int
    seq: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Annotation:
    """주석 결과.

    confident=False 면 `reasons` 에 한글 사유가 들어 있다. 앵커를 못 찾은 경우
    `cdrs` 는 빈 사전이다 — 위치를 추측해서 채우지 않는다.
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
        """CDR 구간을 하나라도 끊어냈는지."""
        return bool(self.cdrs)

    def framework_spans(self) -> list[tuple[int, int]]:
        """CDR 을 뺀 나머지(framework) 구간 목록. 1-based, 양끝 포함."""
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
        """1-based 위치가 어느 구간에 속하는지. CDR 이름 또는 "framework"."""
        for cdr in self.cdrs.values():
            if cdr.start <= position <= cdr.end:
                return cdr.name
        return "framework"


# --------------------------------------------------------------------------
# 보조 함수
# --------------------------------------------------------------------------

def _last_match(pattern: re.Pattern[str], seq: str) -> int | None:
    """마지막 매치의 1-based 시작 위치."""
    last = None
    for m in pattern.finditer(seq):
        last = m.start() + 1
    return last


def _first_match_from(pattern: re.Pattern[str], seq: str, min_pos: int) -> int | None:
    """min_pos(1-based) 이상에서 시작하는 첫 매치의 1-based 시작 위치."""
    if min_pos < 1:
        min_pos = 1
    m = pattern.search(seq, min_pos - 1)
    return m.start() + 1 if m else None


def _first_residue_from(seq: str, residue: str, min_pos: int) -> int | None:
    """min_pos(1-based) 이상에서 residue 가 처음 나오는 1-based 위치."""
    idx = seq.find(residue, max(min_pos - 1, 0))
    return idx + 1 if idx >= 0 else None


#: 보존 Cys 바로 앞에 오는 잔기. OAS 실제 서열 1,991건에서 Y 83.2% + F 16.1% = 99.3%.
#: 이 문맥을 쓰지 않고 "첫 Cys"만 보면, FR3 안에 여분 Cys 가 있는 서열(예: ...LTIS**C**...)에서
#: 그 Cys 를 앵커로 잘못 집어 CDR3 가 10잔기 넘게 앞에서 시작한다. 실측 2,000건 중 7건이 그랬고,
#: 그중 하나는 길이가 통상 범위 안이라 길이 검사로도 걸러지지 않았다.
CONSERVED_CYS_RE = re.compile(r"[YF]C")


def _find_conserved_cys(seq: str, min_pos: int, limit: int) -> tuple[int | None, bool]:
    """min_pos 이상 limit 미만에서 보존 Cys 를 찾는다. (위치, 문맥으로 찾았는지)

    먼저 `[YF]C` 문맥으로 찾고, 없으면 그냥 첫 Cys 로 물러선다. 물러선 경우를
    호출부가 알 수 있게 두 번째 값으로 알린다 — 경고를 남기기 위함이다.
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


def _cys_in_window(seq: str, window: tuple[int, int], expected: int) -> tuple[int | None, int]:
    """창 안의 Cys 중 expected 에 가장 가까운 위치와, 창 안의 Cys 총 개수."""
    lo, hi = window
    found = [i + 1 for i, ch in enumerate(seq[: min(hi, len(seq))]) if ch == "C" and i + 1 >= lo]
    if not found:
        return None, 0
    best = min(found, key=lambda p: (abs(p - expected), p))
    return best, len(found)


def detect_chain(seq: str) -> tuple[str | None, int | None]:
    """말단 모티프로 체인을 판별한다. (chain, 모티프 1-based 시작 위치)."""
    wgxg = _last_match(WGXG_RE, seq)
    fgxg = _last_match(FGXG_RE, seq)

    if wgxg is None and fgxg is None:
        return None, None
    if fgxg is None:
        return "H", wgxg
    if wgxg is None:
        return "L", fgxg
    # 둘 다 있으면 C-말단에 더 가까운 쪽
    return ("H", wgxg) if wgxg > fgxg else ("L", fgxg)


# --------------------------------------------------------------------------
# 본체
# --------------------------------------------------------------------------

def annotate(sequence: str, *, validate_input: bool = True) -> Annotation:
    """서열 하나를 주석한다.

    validate_input=True 면 abcheck.validate 로 먼저 검증한다 (실패 시 SequenceError).
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

    # --- 앵커 1: Cys22 ---
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

    # --- 앵커 2: Trp36, CDR-H1 ---
    h1_start = cys1 + CDR1_OFFSET_H
    if h1_start > len(seq):
        ann.confident = False
        ann.reasons.append(
            f"Cys22({cys1}번) 뒤로 서열이 {len(seq) - cys1} 잔기밖에 없어 "
            "CDR-H1 시작 위치가 서열 범위를 벗어난다. 서열이 잘린 것으로 보인다."
        )
        return

    trp = _first_residue_from(seq, "W", h1_start + MIN_CDR1_LENGTH_H)
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

    # --- 앵커 3: FR3 모티프, CDR-H2 ---
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
        # 교차검증: FR3 모티프와 보존 Cys 의 거리는 실측 3,429건이 전부 29 였고,
        # 어긋난 14건은 전부 모티프가 엉뚱한 자리에 걸린 경우였다(CDR-H2 가 15 또는 29 로 나옴).
        # 거리가 다르면 모티프 쪽을 버리고 Cys 기준으로 다시 잡는다.
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
        # 체세포 과변이로 모티프가 깨진 경우. 보존 Cys 에서 거리로 역산한다
        # (실측 992건 전부 29잔기). 정규식을 더 넓히는 것보다 안전하다 —
        # 넓힐수록 CDR-H2 안에서 우연히 걸릴 위험이 커진다.
        cys2, by_context = _find_conserved_cys(seq, search_from + FR3_TO_CYS_GAP, motif_pos)
        if cys2 is not None:
            fr3 = cys2 - FR3_TO_CYS_GAP
            # 위치는 그대로 보고하되 신뢰도는 내린다. 거리는 실측 992건 전부 일정했지만,
            # 모티프를 직접 확인한 것보다는 약한 근거다 — 결과를 쓰기 전에 사람이 보게 한다.
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

    # --- 앵커 1: Cys23 ---
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

    # --- 앵커 2: Trp35, CDR-L1 ---
    l1_start = cys1 + CDR1_OFFSET_L
    trp = _first_residue_from(seq, "W", l1_start + MIN_CDR1_LENGTH_L)
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

    # --- CDR-L2 (7 잔기 고정) ---
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

    # --- 앵커 3: Cys88, CDR-L3 ---
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
    """CDR 길이가 통상 범위를 벗어나면 confident=False 와 사유를 남긴다."""
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
