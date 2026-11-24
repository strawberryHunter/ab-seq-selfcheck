"""물성 계산 — pI, GRAVY, 특정 pH 에서의 net charge.

전부 **서열 기반 추정치**다. 실측 pI 는 구조·당화·완충액 조성에 따라 달라지므로
이 값은 상대 비교용으로만 쓴다.

도메인 전체와 CDR별로 각각 계산한다. CDR 처럼 잘라낸 조각의 pI 는 그 조각을
독립된 펩타이드로 보고 계산한다 — 즉 조각의 양 말단에 N-말단 아미노기와 C-말단
카복실기가 있다고 가정한다. 실제 CDR 은 도메인 안에 박혀 있어 자유 말단이 없으므로,
CDR pI 는 "이 구간의 전하 조성"을 보는 지표로만 읽어야 한다 (DECISIONS.md 참고).
"""

from __future__ import annotations

from dataclasses import dataclass

# --------------------------------------------------------------------------
# pKa — EMBOSS (iep) 값
# --------------------------------------------------------------------------

#: 측쇄 pKa. EMBOSS 의 Epk.dat 에서 쓰는 값이다.
SIDE_CHAIN_PKA: dict[str, float] = {
    "D": 3.9,
    "E": 4.1,
    "C": 8.5,
    "Y": 10.1,
    "H": 6.5,
    "K": 10.8,
    "R": 12.5,
}

#: 양전하를 띠는(양성자화 시 +1) 측쇄
POSITIVE_RESIDUES = frozenset("KRH")

#: 음전하를 띠는(탈양성자화 시 -1) 측쇄
NEGATIVE_RESIDUES = frozenset("DECY")

PKA_N_TERMINUS = 8.6
PKA_C_TERMINUS = 3.6

#: Kyte-Doolittle hydropathy 지수
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

#: net charge 를 보고할 기본 pH. 생리적 조건에 가깝다.
DEFAULT_PH = 7.4

#: pI 탐색 범위와 수렴 조건
_PI_LOW, _PI_HIGH = 0.0, 14.0
_PI_TOLERANCE = 1e-6
_PI_MAX_ITERATIONS = 200


# --------------------------------------------------------------------------
# 계산
# --------------------------------------------------------------------------

def net_charge(seq: str, ph: float = DEFAULT_PH, *, include_termini: bool = True) -> float:
    """주어진 pH 에서의 net charge.

    Henderson-Hasselbalch 로 각 이온화 그룹의 평균 전하를 더한다.
    양전하 그룹: 1 / (1 + 10^(pH - pKa))
    음전하 그룹: -1 / (1 + 10^(pKa - pH))
    """
    seq = seq.upper()
    charge = 0.0

    if include_termini and seq:
        charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))
        charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))

    for residue in seq:
        pka = SIDE_CHAIN_PKA.get(residue)
        if pka is None:
            continue
        if residue in POSITIVE_RESIDUES:
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        elif residue in NEGATIVE_RESIDUES:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(seq: str, *, include_termini: bool = True) -> float | None:
    """net charge 가 0 이 되는 pH. bisection 으로 찾는다.

    이온화 그룹이 하나도 없으면(예: 'AAAA' 에서 말단까지 제외한 경우) None 을 돌려준다.
    """
    seq = seq.upper()
    if not seq:
        return None

    has_group = include_termini or any(r in SIDE_CHAIN_PKA for r in seq)
    if not has_group:
        return None

    low, high = _PI_LOW, _PI_HIGH
    charge_low = net_charge(seq, low, include_termini=include_termini)
    charge_high = net_charge(seq, high, include_termini=include_termini)

    # 전 구간에서 부호가 같으면 pI 가 탐색 범위 밖이다 — 경계값을 돌려준다.
    if charge_low * charge_high > 0:
        return low if abs(charge_low) < abs(charge_high) else high

    for _ in range(_PI_MAX_ITERATIONS):
        mid = (low + high) / 2.0
        charge_mid = net_charge(seq, mid, include_termini=include_termini)
        if abs(charge_mid) < _PI_TOLERANCE or (high - low) < _PI_TOLERANCE:
            return mid
        # charge 는 pH 에 대해 단조 감소한다
        if charge_mid > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def gravy(seq: str) -> float | None:
    """GRAVY — Kyte-Doolittle hydropathy 의 산술평균.

    양수면 전체적으로 소수성, 음수면 친수성이다. 표에 없는 문자는 건너뛴다.
    """
    seq = seq.upper()
    values = [KYTE_DOOLITTLE[r] for r in seq if r in KYTE_DOOLITTLE]
    if not values:
        return None
    return sum(values) / len(values)


# --------------------------------------------------------------------------
# 묶음 결과
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Properties:
    """한 구간의 물성."""

    label: str
    sequence: str
    length: int
    pi: float | None
    gravy: float | None
    net_charge_at_ph: float
    ph: float

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "sequence": self.sequence,
            "length": self.length,
            "pI": None if self.pi is None else round(self.pi, 2),
            "GRAVY": None if self.gravy is None else round(self.gravy, 3),
            "net_charge": round(self.net_charge_at_ph, 2),
            "pH": self.ph,
        }


def compute(label: str, seq: str, ph: float = DEFAULT_PH) -> Properties:
    """한 구간의 물성을 계산한다."""
    seq = seq.upper()
    return Properties(
        label=label,
        sequence=seq,
        length=len(seq),
        pi=isoelectric_point(seq),
        gravy=gravy(seq),
        net_charge_at_ph=net_charge(seq, ph),
        ph=ph,
    )


def compute_all(sequence: str, cdrs: dict, ph: float = DEFAULT_PH) -> dict[str, Properties]:
    """도메인 전체와 각 CDR 의 물성을 계산한다.

    cdrs 는 abcheck.annotate.Annotation.cdrs 와 같은 모양이면 된다
    ({이름: start/end/seq 를 가진 객체}). CDR 이 없으면 도메인만 돌려준다.
    """
    result: dict[str, Properties] = {"도메인 전체": compute("도메인 전체", sequence, ph)}
    for name, cdr in sorted(cdrs.items(), key=lambda kv: kv[1].start):
        result[name] = compute(name, cdr.seq, ph)
    return result
