"""developability liability 모티프 검출.

핵심은 **CDR 안에 있는 것과 framework 에 있는 것의 위험도를 구분**하는 것이다.
같은 `NG` 라도 CDR 안이면 결합에 직접 영향을 주므로 더 위험하고, framework 면
상대적으로 덜 위험하다. 그래서 모티프별 기본 등급을 두고 구간에 따라 한 단계
올리거나 내린다.

이 모듈은 **서열 모티프만** 본다. 구조적 용매 접근성을 보지 않으므로, 묻혀 있는
잔기도 그대로 플래그된다. 실제 위험도는 다를 수 있다 — 플래그는 근거이지 판정이 아니다.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .annotate import Annotation

# --------------------------------------------------------------------------
# 위험도 등급
# --------------------------------------------------------------------------

#: 낮은 쪽부터 높은 쪽으로. 등급 조정은 이 순서를 따라 한 칸씩 움직인다.
RISK_LEVELS: tuple[str, ...] = ("낮음", "중간", "높음")


def _shift(level: str, step: int) -> str:
    """등급을 step 칸 움직인다. 양 끝에서는 멈춘다."""
    idx = RISK_LEVELS.index(level)
    return RISK_LEVELS[max(0, min(len(RISK_LEVELS) - 1, idx + step))]


# --------------------------------------------------------------------------
# 모티프 정의
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class MotifSpec:
    """liability 모티프 하나의 정의."""

    kind: str                 # 분류 (예: "deamidation")
    pattern: str              # 정규식
    base_risk: str            # 구간 조정 전 기본 등급
    label: str                # 한글 이름
    rationale: str            # 왜 문제인지 (한글 근거)
    top_within_kind: bool = False   # 같은 분류 안에서 최고 위험인지

    @property
    def regex(self) -> re.Pattern[str]:
        return re.compile(self.pattern)


#: 검출 대상 모티프. 같은 분류 안에서 더 위험한 것을 먼저 둔다 (겹치면 앞선 것이 이긴다).
MOTIF_SPECS: tuple[MotifSpec, ...] = (
    MotifSpec(
        kind="glycosylation",
        pattern=r"N[^P][ST]",
        base_risk="높음",
        label="N-결합 당화 서열",
        rationale=(
            "N-X-S/T (X≠Pro) 는 N-결합 당화 수용 서열이다. 당화가 일어나면 "
            "불균일성(heterogeneity)이 생기고 배치 간 재현성이 떨어진다."
        ),
        top_within_kind=True,
    ),
    MotifSpec(
        kind="deamidation",
        pattern=r"NG",
        base_risk="높음",
        label="deamidation (NG)",
        rationale=(
            "Asn-Gly 는 deamidation 속도가 가장 빠른 조합이다. Asn 이 Asp/isoAsp 로 "
            "바뀌면 전하와 구조가 변해 결합력 손실로 이어질 수 있다."
        ),
        top_within_kind=True,
    ),
    MotifSpec(
        kind="deamidation",
        pattern=r"N[STNH]",
        base_risk="중간",
        label="deamidation (NS/NT/NN/NH)",
        rationale=(
            "Asn 뒤 잔기에 따라 deamidation 이 일어날 수 있다. NG 보다는 느리지만 "
            "장기 보관 안정성에서 확인이 필요하다."
        ),
    ),
    MotifSpec(
        kind="isomerization",
        pattern=r"DG",
        base_risk="높음",
        label="isomerization (DG)",
        rationale=(
            "Asp-Gly 는 isomerization 속도가 가장 빠른 조합이다. Asp 가 isoAsp 로 "
            "바뀌면 주쇄가 한 원자 길어져 국소 구조가 틀어진다."
        ),
        top_within_kind=True,
    ),
    MotifSpec(
        kind="isomerization",
        pattern=r"D[STD]",
        base_risk="중간",
        label="isomerization (DS/DT/DD)",
        rationale=(
            "Asp 뒤 잔기에 따라 isomerization 이 일어날 수 있다. DG 보다는 느리지만 "
            "가속 안정성 시험에서 확인이 필요하다."
        ),
    ),
    MotifSpec(
        kind="fragmentation",
        pattern=r"DP",
        base_risk="중간",
        label="fragmentation (DP)",
        rationale=(
            "Asp-Pro 펩타이드 결합은 산성 조건에서 쉽게 끊어진다. 저장 완충액 pH 가 "
            "낮으면 단편화가 나타날 수 있다."
        ),
        top_within_kind=True,
    ),
    MotifSpec(
        kind="oxidation",
        pattern=r"M",
        base_risk="낮음",
        label="산화 (Met)",
        rationale=(
            "Met 은 산화되기 쉽다. 노출된 Met 이 산화되면 결합력이나 안정성이 떨어질 수 있다."
        ),
    ),
    MotifSpec(
        kind="oxidation",
        pattern=r"W",
        base_risk="낮음",
        label="산화 (Trp)",
        rationale=(
            "Trp 은 빛·금속 이온에 의해 산화될 수 있다. 보존 Trp 앵커(Trp36, WGxG 의 Trp)는 "
            "대개 묻혀 있어 위험이 낮지만, 검출 자체는 숨기지 않는다."
        ),
    ),
)


# --------------------------------------------------------------------------
# 결과 자료형
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Hit:
    """liability 모티프 한 건."""

    kind: str
    label: str
    match: str
    start: int                # 1-based, 양끝 포함
    end: int
    region: str               # "CDR-H2" 또는 "framework" 또는 "구간 미확정"
    in_cdr: bool
    risk: str                 # 구간 조정을 반영한 최종 등급
    base_risk: str
    rationale: str
    note: str                 # 등급 조정 근거 (한글)


@dataclass
class LiabilityReport:
    hits: list[Hit] = field(default_factory=list)
    #: 도메인 단위 플래그 (비결합 Cys 처럼 위치 하나로 특정할 수 없는 것)
    domain_flags: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def cdr_hits(self) -> list[Hit]:
        return [h for h in self.hits if h.in_cdr]

    @property
    def framework_hits(self) -> list[Hit]:
        return [h for h in self.hits if not h.in_cdr]

    def by_risk(self, level: str) -> list[Hit]:
        return [h for h in self.hits if h.risk == level]

    def counts_by_risk(self) -> dict[str, int]:
        return {lvl: len(self.by_risk(lvl)) for lvl in reversed(RISK_LEVELS)}


# --------------------------------------------------------------------------
# 검출
# --------------------------------------------------------------------------

def _region_name(ann: Annotation, start: int, end: int) -> tuple[str, bool]:
    """모티프 구간이 어디에 걸치는지. (구간 이름, CDR 여부)

    모티프가 CDR 경계에 걸쳐 있으면 CDR 쪽으로 본다 — 위험도를 과소평가하지 않기 위한
    보수적 선택이다.
    """
    if not ann.cdrs:
        return "구간 미확정", False
    for cdr in ann.cdrs.values():
        if start <= cdr.end and end >= cdr.start:
            return cdr.name, True
    return "framework", False


def _adjust(base: str, in_cdr: bool, region: str) -> tuple[str, str]:
    """구간에 따라 등급을 조정한다. (최종 등급, 한글 근거)"""
    if region == "구간 미확정":
        return base, (
            "CDR 주석이 실패해 구간을 구분할 수 없다. 기본 등급을 그대로 쓴다 — "
            "CDR 안일 수도 있으니 과소평가하지 말 것."
        )
    if in_cdr:
        adjusted = _shift(base, +1)
        if adjusted == base:
            return adjusted, f"{region} 안이고 기본 등급이 이미 최고({base})라 그대로 둔다."
        return adjusted, f"{region} 안이라 기본 {base} 에서 한 단계 올렸다 (결합에 직접 영향)."
    adjusted = _shift(base, -1)
    if adjusted == base:
        return adjusted, f"framework 이고 기본 등급이 이미 최저({base})라 그대로 둔다."
    return adjusted, f"framework 이라 기본 {base} 에서 한 단계 내렸다 (결합 영향이 간접적)."


def find_liabilities(ann: Annotation) -> LiabilityReport:
    """주석 결과를 받아 liability 를 검출한다.

    CDR 주석이 실패한 서열(`ann.cdrs` 가 빈 경우)에서도 동작한다. 그때는 구간을
    "구간 미확정"으로 표시하고 등급 조정을 하지 않으며 경고를 남긴다.
    """
    seq = ann.sequence
    report = LiabilityReport()

    if not ann.cdrs:
        report.warnings.append(
            "CDR 주석이 실패해 CDR/framework 구간을 구분할 수 없다. "
            "모든 모티프를 '구간 미확정'으로 보고하며 구간별 위험도 조정을 적용하지 않았다."
        )

    #: 같은 분류 안에서 더 위험한 모티프가 이미 덮은 자리는 건너뛴다.
    #: (예: "NG" 를 잡았으면 같은 자리를 "N[STNH]" 로 다시 잡지 않는다)
    claimed: dict[str, set[int]] = {}

    for spec in MOTIF_SPECS:
        taken = claimed.setdefault(spec.kind, set())
        for m in spec.regex.finditer(seq):
            start, end = m.start() + 1, m.end()
            positions = set(range(start, end + 1))
            if positions & taken:
                continue
            taken |= positions

            region, in_cdr = _region_name(ann, start, end)
            risk, note = _adjust(spec.base_risk, in_cdr, region)
            report.hits.append(
                Hit(
                    kind=spec.kind,
                    label=spec.label,
                    match=m.group(),
                    start=start,
                    end=end,
                    region=region,
                    in_cdr=in_cdr,
                    risk=risk,
                    base_risk=spec.base_risk,
                    rationale=spec.rationale,
                    note=note,
                )
            )

    report.hits.sort(key=lambda h: (h.start, h.kind))
    _add_cysteine_flag(seq, ann, report)
    return report


def _add_cysteine_flag(seq: str, ann: Annotation, report: LiabilityReport) -> None:
    """도메인 내 Cys 개수가 홀수면 비결합 Cys 플래그를 단다."""
    positions = [i + 1 for i, ch in enumerate(seq) if ch == "C"]
    count = len(positions)
    if count % 2 == 0:
        return

    regions = [ann.region_of(p) if ann.cdrs else "구간 미확정" for p in positions]
    in_cdr = [p for p, r in zip(positions, regions) if r.startswith("CDR")]
    risk = "높음" if in_cdr else "중간"

    detail = ", ".join(f"{p}번({r})" for p, r in zip(positions, regions))
    rationale = (
        f"도메인 내 Cys 가 {count}개로 홀수다. 가변도메인은 보존 Cys 한 쌍이 "
        "도메인 내부 이황화결합을 이루므로 개수는 보통 짝수다. 홀수면 짝을 못 이룬 "
        "Cys 가 있다는 뜻이고, 이는 응집·이분자 가교·불균일한 이황화 짝짓기를 일으킬 수 있다."
    )
    if in_cdr:
        rationale += (
            f" 특히 CDR 안의 Cys({', '.join(str(p) for p in in_cdr)}번)는 "
            "결합 부위에 있어 위험도가 더 높다."
        )

    report.domain_flags.append(
        {
            "kind": "unpaired_cysteine",
            "label": "비결합 Cys 가능성",
            "risk": risk,
            "cysteine_count": count,
            "positions": positions,
            "regions": regions,
            "rationale": rationale,
        }
    )
