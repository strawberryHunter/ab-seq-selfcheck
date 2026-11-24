#!/usr/bin/env python3
"""합성 항체 서열 생성기 (표준 라이브러리만 사용).

목적
----
abcheck 의 CDR 주석/liability 검출 로직을 검증할 수 있는 테스트 코퍼스를 만든다.
"germline 유사 framework 스캐폴드 + CDR 삽입" 방식으로 서열을 조립하므로,
삽입 시점에 CDR 의 정답 위치(1-based, 양끝 포함)를 정확히 알 수 있고 그것을 함께 저장한다.

중요 — 데이터 출처에 관한 명시
------------------------------
이 파일이 만드는 서열은 전부 합성이다. 실제 치료용/승인 항체 서열이나 특정 germline
데이터베이스 레코드를 복사해 오지 않는다. framework 부분은 V 도메인에서 흔히 보이는
잔기 조성(소수성/소형 잔기 우세)을 흉내낸 가중 난수에 보존 잔기(Cys, Trp, WGxG/FGxG 등)만
정해진 자리에 고정해 만든 "germline 유사" 스캐폴드일 뿐이며, 어떤 실제 서열과도 대응되지 않는다.
따라서 이 데이터는 주석 로직의 회귀 테스트용이고, 생물학적 의미를 해석할 대상이 아니다.

스캐폴드 좌표 설계 (VH, 1-based 선형 위치; Kabat 경계를 모사)
------------------------------------------------------------
    1..21      FR1 filler
    22         Cys  (앵커 1)
    23..30     filler (8)
    31..       CDR-H1  (= Cys22 + 9 에서 시작)
    H1 직후    Trp     (앵커 2, 통상 Trp36)
    Trp+1..+13 filler (13)
    Trp+14..   CDR-H2
    H2 직후    FR3 모티프 5 잔기  [RK].[TSAV][IMLFV][STA]  (앵커 3)
    모티프+5..  filler (21) → 그 다음 Cys (앵커 4, 통상 Cys92)
    Cys+1..+2  filler (2)
    Cys+3..    CDR-H3
    H3 직후    WGxG (앵커 5)
    그 뒤      FR4 filler (7)

VL
--
    1..22      FR1 filler
    23         Cys  (앵커 1)
    24..       CDR-L1 (= Cys23 + 1)
    L1 직후    Trp    (앵커 2, 통상 Trp35)
    Trp+1..+14 filler (14)
    Trp+15..   CDR-L2 (7 잔기 고정)
    L2 직후    filler (31) → 그 다음 Cys (앵커 3, 통상 Cys88)
    Cys+1..    CDR-L3
    L3 직후    FGxG (앵커 4)
    그 뒤      FR4 filler (7)

생성기는 조립 후 "주석기가 쓸 탐색 규칙"이 의도한 앵커를 정확히 집어내는지 스스로 검증한다
(ASSERTIONS 참고). 통과하지 못하면 시드를 바꿔 다시 만든다. 즉 이 코퍼스는 탐색 규칙이
우연한 중복 모티프에 걸리지 않는, 통제된 입력만 담는다. 앵커가 일부러 깨진 케이스는
expected.confident = False 로 표시해 따로 구분한다.
"""

from __future__ import annotations

import argparse
import json
import random
import re
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# 알파벳 / 모티프
# --------------------------------------------------------------------------

#: framework filler 에 쓰는 잔기와 가중치.
#: Cys/Trp/Met/Asn/Asp 는 일부러 제외한다 — 보존 Cys·Trp 앵커를 유일하게 유지하고,
#: framework 에 우연한 liability 모티프(NG, DG, DP, N[^P][ST], 산화 M/W)가 생기지 않게 한다.
#: framework liability 가 필요한 케이스는 아래에서 명시적으로 주입한다.
FRAMEWORK_WEIGHTS: dict[str, int] = {
    "A": 7, "G": 7, "V": 7, "L": 8, "I": 5, "P": 5, "F": 4,
    "S": 9, "T": 8, "Y": 4, "Q": 6, "E": 6, "K": 6, "R": 5, "H": 2,
}

#: CDR filler 알파벳. Cys/Trp 는 제외(의도적으로 주입할 때만 넣는다).
CDR_WEIGHTS: dict[str, int] = {
    "A": 5, "G": 7, "V": 4, "L": 4, "I": 4, "P": 3, "F": 3,
    "S": 7, "T": 6, "Y": 7, "Q": 3, "E": 3, "K": 3, "R": 4,
    "H": 2, "N": 4, "D": 5,
}

#: CDR 청결성 검사에 쓰는 liability 모티프.
#: 생성기가 "의도한 것 말고는 CDR 안에 liability 가 없다"를 보장하기 위한 용도이며,
#: 실제 위험도 판정 로직(abcheck.liabilities)과는 별개로 유지한다.
CDR_LIABILITY_PATTERNS: list[tuple[str, str]] = [
    ("glycosylation", r"N[^P][ST]"),
    ("deamidation", r"N[GSTNH]"),
    ("isomerization", r"D[GSTD]"),
    ("fragmentation", r"DP"),
    ("oxidation", r"[MW]"),
    ("cysteine", r"C"),
]

FR3_MOTIF_RE = re.compile(r"[RK].[TSAV][IMLFV][STA]")

#: 보존 Trp 앵커를 찾기 시작할 오프셋. abcheck.annotate 와 같은 값을 쓴다.
MIN_CDR1_LENGTH_H = 5
MIN_CDR1_LENGTH_L = 7
WGXG_RE = re.compile(r"WG.G")
FGXG_RE = re.compile(r"FG.G")

# 주석기가 "통상 범위"로 볼 CDR 길이 (스크리닝 휴리스틱 기준; 엄격한 numbering 기준이 아니다)
TYPICAL_LENGTHS = {
    "H1": (5, 7),
    "H2": (16, 19),
    "H3": (3, 25),
    "L1": (7, 17),
    "L2": (7, 7),
    "L3": (5, 13),
}


# --------------------------------------------------------------------------
# 스캐폴드
# --------------------------------------------------------------------------

@dataclass
class HeavyScaffold:
    """VH framework 조각들. 모든 VH 케이스가 같은 스캐폴드를 공유한다."""

    fr1: str          # 21
    post_cys1: str    # 8
    fr2: str          # 13
    fr3_motif: str    # 5
    fr3_rest: str     # 21
    post_cys2: str    # 2
    fr4_rest: str     # 7


@dataclass
class LightScaffold:
    """VL framework 조각들."""

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
        fr2=_draw(rng, FRAMEWORK_WEIGHTS, 13),
        fr3_motif=motif,
        fr3_rest=_draw(rng, FRAMEWORK_WEIGHTS, 21),
        post_cys2=_draw(rng, FRAMEWORK_WEIGHTS, 2),
        fr4_rest=_draw(rng, FRAMEWORK_WEIGHTS, 7),
    )


def build_light_scaffold(rng: random.Random) -> LightScaffold:
    return LightScaffold(
        fr1=_draw(rng, FRAMEWORK_WEIGHTS, 22),
        fr2=_draw(rng, FRAMEWORK_WEIGHTS, 14),
        fr3_rest=_draw(rng, FRAMEWORK_WEIGHTS, 31),
        fr4_rest=_draw(rng, FRAMEWORK_WEIGHTS, 7),
    )


# --------------------------------------------------------------------------
# 케이스 정의
# --------------------------------------------------------------------------

@dataclass
class Case:
    """합성 서열 한 건의 설계도."""

    seq_id: str
    chain: str                       # "H" | "L"
    description: str                 # 한글 설명
    cdr_lengths: tuple[int, int, int]
    #: CDR 안에 강제로 박아 넣을 조각. {"H1": [(offset, "NLS")], ...} (offset 은 0-based)
    cdr_inserts: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    #: framework 조각에 강제로 박아 넣을 조각. {"fr3_rest": [(offset, "DG")], ...}
    fr_inserts: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    #: 조립이 끝난 서열에 가하는 치환 (1-based 위치, 새 잔기). 앵커를 깨는 데 쓴다.
    breakages: list[tuple[int, str]] = field(default_factory=list)
    #: 앵커가 깨져 주석이 실패해야 하는 케이스
    anchors_broken: bool = False
    #: 앵커는 멀쩡하지만 CDR 길이가 통상 범위를 벗어나 confident=False 가 되어야 하는 케이스
    length_outlier: bool = False
    notes: list[str] = field(default_factory=list)


def _apply_inserts(seq: str, inserts: list[tuple[int, str]]) -> str:
    """offset 위치에서 길이를 유지한 채로 조각을 덮어쓴다 (길이 보존이 좌표 보존의 핵심)."""
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
    """CDR 안의 liability 모티프 위치를 (0-based 시작, 끝+1, 라벨) 로 돌려준다."""
    hits: list[tuple[int, int, str]] = []
    for label, pattern in CDR_LIABILITY_PATTERNS:
        for m in re.finditer(pattern, cdr):
            hits.append((m.start(), m.end(), label))
    return hits


#: liability 모티프 최대 길이(N[^P][ST] = 3)에서 1을 뺀 값.
#: CDR 경계를 걸치는 모티프를 보려면 양옆으로 이만큼 넓게 봐야 한다.
_MOTIF_PAD = 2


def unintended_cdr_liabilities_in_context(
    seq: str, cdr_start: int, cdr_end: int, inserts: list[tuple[int, str]]
) -> list[tuple[int, int, str]]:
    """CDR 에 **걸치는** liability 중 주입 조각 안에 들어가지 않는 것들.

    CDR 문자열만 보면 경계를 걸친 모티프를 놓친다. 예를 들어 CDR 이 ...SNA 로 끝나고
    바로 뒤 framework 가 S 로 시작하면 NAS 가 당화 수용 서열이 되는데, CDR 안쪽만
    훑으면 이걸 못 본다. 그래서 양옆 2 잔기까지 넓혀 훑고, CDR 과 겹치는 hit 만 본다.

    좌표는 모두 1-based 이고 양끝을 포함한다.
    """
    lo = max(1, cdr_start - _MOTIF_PAD)
    hi = min(len(seq), cdr_end + _MOTIF_PAD)
    window = seq[lo - 1:hi]

    #: 주입 조각의 전역 좌표 구간
    spans = [(cdr_start + off, cdr_start + off + len(frag) - 1) for off, frag in inserts]

    out: list[tuple[int, int, str]] = []
    for w_start, w_end, label in _cdr_liability_hits(window):
        start = lo + w_start           # 1-based 시작
        end = lo + w_end - 1           # 1-based 끝(포함)
        if end < cdr_start or start > cdr_end:
            continue                   # CDR 과 겹치지 않으면 framework 쪽 일이다
        if any(s_lo <= start and end <= s_hi for s_lo, s_hi in spans):
            continue                   # 의도적으로 주입한 것
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

    parts = [fr1, "C", post_cys1, h1, "W", fr2, h2,
             scaf.fr3_motif, fr3_rest, "C", scaf.post_cys2, h3, "WGQG", fr4_rest]
    seq = "".join(parts)

    # 1-based 좌표 계산
    cys1 = len(fr1) + 1                                   # 22
    h1_start = cys1 + 9                                   # 31
    h1_end = h1_start + n1 - 1
    trp = h1_end + 1
    h2_start = trp + 14
    h2_end = h2_start + n2 - 1
    fr3_start = h2_end + 1
    cys2 = fr3_start + len(scaf.fr3_motif) + len(fr3_rest)
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
# 자기 검증 — 주석기가 쓸 탐색 규칙이 의도한 앵커를 집는지 확인
# --------------------------------------------------------------------------

def _first_index(pattern: re.Pattern[str], seq: str, min_pos1: int) -> int | None:
    """min_pos1(1-based) 이상에서 시작하는 첫 매치의 1-based 시작 위치."""
    m = pattern.search(seq, min_pos1 - 1)
    return m.start() + 1 if m else None


def _last_index(pattern: re.Pattern[str], seq: str) -> int | None:
    last = None
    for m in pattern.finditer(seq):
        last = m.start() + 1
    return last


def verify_heavy(seq: str, cdrs: dict, anchors: dict) -> list[str]:
    """의도한 앵커가 탐색 규칙으로 유일하게 복원되는지 검사. 문제를 문자열 목록으로 돌려준다."""
    problems: list[str] = []

    # 말단 WGxG 는 서열의 마지막 WGxG 매치여야 하고, 그 뒤에 FGxG 가 있으면 체인 판별이 흔들린다
    if _last_index(WGXG_RE, seq) != anchors["terminal_motif"]:
        problems.append("말단 WGxG 가 마지막 매치가 아니다")
    last_fgxg = _last_index(FGXG_RE, seq)
    if last_fgxg is not None and last_fgxg > anchors["terminal_motif"]:
        problems.append("말단 WGxG 뒤에 FGxG 가 있어 체인 판별이 모호하다")

    # Cys22: 15..30 창에서 유일해야 한다
    window = [i + 1 for i, r in enumerate(seq[:30]) if r == "C" and i + 1 >= 15]
    if window != [anchors["cys1"]]:
        problems.append(f"15..30 창의 Cys 가 유일하지 않다: {window}")

    # Trp36: CDR-H1 시작 + 최소 길이 이후 첫 W
    # (CDR-H1 안에 Trp 이 있는 실제 항체를 흉내내려면 이 규칙이 주석기와 같아야 한다)
    h1_start = cdrs["CDR-H1"][0]
    first_w = seq.find("W", h1_start - 1 + MIN_CDR1_LENGTH_H) + 1
    if first_w != anchors["trp"]:
        problems.append(f"CDR-H1 최소 길이 이후 첫 Trp 위치가 다르다: {first_w} != {anchors['trp']}")

    # FR3 모티프: CDR-H2 시작 + 9 이후 첫 매치
    h2_start = cdrs["CDR-H2"][0]
    found = _first_index(FR3_MOTIF_RE, seq, h2_start + 9)
    if found != anchors["fr3_motif"]:
        problems.append(f"FR3 모티프 첫 매치가 CDR-H2 안에 걸렸다: {found} != {anchors['fr3_motif']}")

    # Cys92: FR3 모티프 이후 첫 Cys
    first_c = seq.find("C", anchors["fr3_motif"] - 1) + 1
    if first_c != anchors["cys2"]:
        problems.append(f"FR3 모티프 이후 첫 Cys 위치가 다르다: {first_c} != {anchors['cys2']}")

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
    first_w = seq.find("W", l1_start - 1 + MIN_CDR1_LENGTH_L) + 1
    if first_w != anchors["trp"]:
        problems.append(f"CDR-L1 최소 길이 이후 첫 Trp 위치가 다르다: {first_w} != {anchors['trp']}")

    l2_end = cdrs["CDR-L2"][1]
    first_c = seq.find("C", l2_end) + 1
    if first_c != anchors["cys2"]:
        problems.append(f"CDR-L2 이후 첫 Cys 위치가 다르다: {first_c} != {anchors['cys2']}")

    return problems


# --------------------------------------------------------------------------
# 케이스 목록
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

        # --- 앵커가 깨져 주석이 실패해야 하는 케이스 ---
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

        # --- 앵커는 멀쩡하지만 길이가 통상 범위를 벗어나는 케이스 ---
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
# 빌드
# --------------------------------------------------------------------------

def build_record(case: Case, heavy: HeavyScaffold, light: LightScaffold) -> dict:
    """케이스 하나를 레코드로 만든다. 자기 검증을 통과할 때까지 CDR 내용을 재추첨한다."""
    last_problems: list[str] = []
    for attempt in range(20000):
        rng = random.Random(f"{case.seq_id}/{attempt}")
        if case.chain == "H":
            seq, cdrs, anchors = _assemble_heavy(case, heavy, rng)
            problems = verify_heavy(seq, cdrs, anchors)
        else:
            seq, cdrs, anchors = _assemble_light(case, light, rng)
            problems = verify_light(seq, cdrs, anchors)
        # 의도하지 않은 CDR liability 가 섞였으면 재추첨한다.
        # 이렇게 하면 각 케이스의 CDR liability 는 case.cdr_inserts 로 주입한 것만 남는다.
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

        expected_confident = not (case.anchors_broken or case.length_outlier)
        record = {
            "id": case.seq_id,
            "chain": case.chain,
            "description": case.description,
            "sequence": final_seq,
            "length": len(final_seq),
            "synthetic": True,
            "notes": case.notes,
            "expected": {
                # 앵커가 깨진 케이스는 정답 CDR 위치를 정의하지 않는다.
                # 주석기는 추측하지 말고 confident=False 로 떨어져야 한다.
                "cdrs": None if case.anchors_broken else {
                    name: {"start": s, "end": e, "seq": intact_seq[s - 1:e]}
                    for name, (s, e) in cdrs.items()
                },
                "chain": None if case.anchors_broken else case.chain,
                "confident": expected_confident,
                "anchors_broken": case.anchors_broken,
                "length_outlier": case.length_outlier,
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
    raise SystemExit(main())
