"""사용자가 정한 선별 규칙을 적용해 '검토 필요'를 표시한다.

이 모듈은 도구가 항체의 좋고 나쁨을 판정하는 곳이 **아니다**. 연구자가 미리 정해 둔
규칙을 기계적으로 적용해, 수십 개 후보 중 사람이 먼저 들여다볼 것을 골라낼 뿐이다.
그래서 결과도 '합격/불합격'이 아니라 '검토 필요 / 기준 해당 없음'으로 적는다.

왜 이렇게 나누는가
------------------
점수를 매기려면 "CDR 의 NG 하나"와 "framework 의 DG 둘" 중 무엇이 더 나쁜지 정해야
하는데, 여기에는 합의된 답이 없다. 투여 경로·보관 조건·타깃에 따라 달라진다. 도구가
임의의 가중치로 숫자를 만들면 그 숫자는 근거가 아니라 추측이고, 사람은 근거 대신
숫자만 보게 된다. 그래서 가중치를 만들지 않고, 규칙을 눈에 보이게 적어 두고 어떤
규칙이 걸렸는지 그대로 보여준다.

규칙을 바꾸려면 DEFAULT_RULES 를 고치면 된다. 기준이 코드에 명시적으로 남으므로
나중에 "왜 이 후보를 먼저 봤는지" 추적할 수 있다.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from .report import SequenceReport

#: 선별 결과
STATUS_REVIEW = "검토 필요"
STATUS_CLEAR = "기준 해당 없음"
STATUS_UNKNOWN = "판정 불가"

#: 이 모듈이 하는 일과 하지 않는 일. 출력에 그대로 싣는다.
TRIAGE_NOTICE: str = (
    "'검토 필요'는 아래 선별 규칙에 걸렸다는 뜻이지, 이 항체가 나쁘다는 뜻이 아니다. "
    "'기준 해당 없음'도 안전하다는 뜻이 아니라 이 규칙에 걸리지 않았다는 뜻일 뿐이다. "
    "규칙은 연구자가 정한 것이고, 도구는 그대로 적용만 한다."
)


@dataclass(frozen=True)
class Rule:
    """선별 규칙 하나."""

    key: str
    description: str
    #: report 를 받아 걸린 사유 목록을 돌려준다. 비어 있으면 걸리지 않은 것이다.
    check: Callable[["SequenceReport"], list[str]]


@dataclass
class TriageResult:
    status: str
    reasons: list[str] = field(default_factory=list)
    #: 규칙을 적용할 수 있었는지 (CDR 주석이 실패하면 적용 불가)
    evaluable: bool = True
    note: str = ""

    @property
    def needs_review(self) -> bool:
        return self.status == STATUS_REVIEW

    def as_dict(self) -> dict:
        return {
            "status": self.status,
            "reasons": list(self.reasons),
            "evaluable": self.evaluable,
            "note": self.note,
            "notice": TRIAGE_NOTICE,
        }


# --------------------------------------------------------------------------
# 규칙
# --------------------------------------------------------------------------

def _cdr_high_risk(report: "SequenceReport") -> list[str]:
    """CDR 안에 위험도 '높음' liability 가 있으면 걸린다."""
    if report.liabilities is None:
        return []
    reasons = []
    for hit in report.liabilities.cdr_hits:
        if hit.risk == "높음":
            reasons.append(
                f"{hit.region} {hit.start}번에 {hit.label} ({hit.match}) — 위험도 높음"
            )
    return reasons


def _domain_flag_high_risk(report: "SequenceReport") -> list[str]:
    """도메인 단위 플래그(비결합 Cys)가 '높음'이면 걸린다.

    위치 하나로 특정할 수 없어 CDR hit 목록에는 안 들어가지만, 비결합 Cys 가 CDR 안에
    있으면 위험도가 '높음'으로 매겨진다. CDR 의 높음을 보려는 규칙이 이것을 놓치면
    곤란하므로 같이 본다.
    """
    if report.liabilities is None:
        return []
    return [
        f"{flag['label']} — 위험도 {flag['risk']} (Cys {flag['cysteine_count']}개, 홀수)"
        for flag in report.liabilities.domain_flags
        if flag.get("risk") == "높음"
    ]


#: 기본 선별 규칙. 바꾸려면 여기를 고친다.
DEFAULT_RULES: tuple[Rule, ...] = (
    Rule(
        key="cdr-high",
        description="CDR 안에 위험도 '높음' liability 가 1건 이상 있으면 검토 필요",
        check=_cdr_high_risk,
    ),
    Rule(
        key="domain-high",
        description="비결합 Cys 등 도메인 플래그가 위험도 '높음'이면 검토 필요",
        check=_domain_flag_high_risk,
    ),
)


# --------------------------------------------------------------------------
# 적용
# --------------------------------------------------------------------------

def evaluate(report: "SequenceReport", rules: tuple[Rule, ...] = DEFAULT_RULES) -> TriageResult:
    """선별 규칙을 적용한다.

    CDR 주석이 실패한 서열은 CDR 과 framework 를 가를 수 없으므로 규칙을 적용할 수 없다.
    그때는 '기준 해당 없음'으로 넘기지 않고 '판정 불가'로 표시한다 — 못 본 것을 괜찮다고
    적으면 그게 가장 위험하다.
    """
    if not report.ok:
        return TriageResult(
            status=STATUS_UNKNOWN,
            evaluable=False,
            note="입력이 검증을 통과하지 못해 규칙을 적용할 수 없다.",
        )

    if report.annotation is None or not report.annotation.cdrs:
        return TriageResult(
            status=STATUS_UNKNOWN,
            evaluable=False,
            note=(
                "CDR 주석이 실패해 CDR 과 framework 를 가를 수 없다. "
                "'CDR 안의 높음' 규칙을 적용할 수 없으므로 사람이 직접 확인해야 한다."
            ),
        )

    reasons: list[str] = []
    for rule in rules:
        reasons.extend(rule.check(report))

    if reasons:
        return TriageResult(status=STATUS_REVIEW, reasons=reasons)

    note = ""
    if not report.annotation.confident:
        note = (
            "규칙에는 걸리지 않았지만 주석 신뢰도가 낮다(confident=False). "
            "CDR 경계가 밀렸을 수 있으므로 ANARCI/IMGT 로 재확인할 것."
        )
    return TriageResult(status=STATUS_CLEAR, note=note)


def rule_descriptions(rules: tuple[Rule, ...] = DEFAULT_RULES) -> list[str]:
    """적용 중인 규칙 설명 목록. 출력에 실어 기준을 눈에 보이게 한다."""
    return [rule.description for rule in rules]
