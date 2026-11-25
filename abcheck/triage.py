"""Apply user-defined triage rules and mark sequences as "needs review".

This module is **not** where the tool judges whether an antibody is good or bad. It just
applies rules the researcher set in advance, mechanically, to pick out which of dozens of
candidates a person should look at first. So results read "needs review / no rule matched"
(STATUS_REVIEW / STATUS_CLEAR), not "pass/fail".

Why split it this way
---------------------
Scoring would require deciding which is worse, "one NG in a CDR" or "two DGs in the
framework", and there is no agreed answer: it depends on route of administration, storage
conditions and target. If the tool produced a number from arbitrary weights, that number
would be a guess, not evidence, and people would look at the number instead of the evidence.
So no weights are made up; the rules are written out in plain sight and the output shows
exactly which rule fired.

To change the rules, edit DEFAULT_RULES. The criteria stay explicit in code, so you can
later trace "why was this candidate looked at first".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from .report import SequenceReport

#: Triage outcomes
STATUS_REVIEW = "검토 필요"
STATUS_CLEAR = "기준 해당 없음"
STATUS_UNKNOWN = "판정 불가"

#: What this module does and doesn't do. Included verbatim in the output.
TRIAGE_NOTICE: str = (
    "'검토 필요'는 아래 선별 규칙에 걸렸다는 뜻이지, 이 항체가 나쁘다는 뜻이 아니다. "
    "'기준 해당 없음'도 안전하다는 뜻이 아니라 이 규칙에 걸리지 않았다는 뜻일 뿐이다. "
    "규칙은 연구자가 정한 것이고, 도구는 그대로 적용만 한다."
)


@dataclass(frozen=True)
class Rule:
    """A single triage rule."""

    key: str
    description: str
    #: Takes a report and returns the list of reasons it matched. Empty means no match.
    check: Callable[["SequenceReport"], list[str]]


@dataclass
class TriageResult:
    status: str
    reasons: list[str] = field(default_factory=list)
    #: Whether the rules could be applied (not possible if CDR annotation failed)
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
# Rules
# --------------------------------------------------------------------------

def _cdr_high_risk(report: "SequenceReport") -> list[str]:
    """Matches when a CDR contains a "high"-risk motif that is **the top-risk motif of its category**.

    Why restrict to top-risk motifs
    -------------------------------
    Being inside a CDR raises the risk by one step, so even weak motifs with a "medium"
    base (NS/NT/NN/NH, DS/DT/DD) become "high". On real repertoire sequences, 13 of 25
    "high" hits inside CDRs were NS/NT/NN/NH, and 12 of 16 sequences got flagged for review.
    Flagging that many defeats the purpose of "picking what to look at first".

    So triage only considers the fastest-reacting motif in each category — N-linked
    glycosylation, NG, DG, DP. Weak motifs still show as "high" in the table and details,
    so no information is lost.
    """
    if report.liabilities is None:
        return []
    reasons = []
    for hit in report.liabilities.cdr_hits:
        if hit.risk == "높음" and hit.top_within_kind:
            reasons.append(
                f"{hit.region} {hit.start}번에 {hit.label} ({hit.match}) — 위험도 높음"
            )
    return reasons


def _domain_flag_high_risk(report: "SequenceReport") -> list[str]:
    """Matches when a domain-level flag (unpaired Cys) is "high" risk.

    It can't be pinned to one position, so it isn't in the CDR hit list, but an unpaired
    Cys inside a CDR is rated "high". A rule meant to catch high risk in CDRs shouldn't
    miss it, so it is checked here too.
    """
    if report.liabilities is None:
        return []
    return [
        f"{flag['label']} — 위험도 {flag['risk']} (Cys {flag['cysteine_count']}개, 홀수)"
        for flag in report.liabilities.domain_flags
        if flag.get("risk") == "높음"
    ]


#: Default triage rules. Edit here to change them.
DEFAULT_RULES: tuple[Rule, ...] = (
    Rule(
        key="cdr-high",
        description=(
            "CDR 안에 위험도 '높음'이면서 각 분류의 최고 위험 모티프"
            "(N-결합 당화, NG, DG, DP)가 1건 이상 있으면 검토 필요"
        ),
        check=_cdr_high_risk,
    ),
    Rule(
        key="domain-high",
        description="비결합 Cys 등 도메인 플래그가 위험도 '높음'이면 검토 필요",
        check=_domain_flag_high_risk,
    ),
)


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------

def evaluate(report: "SequenceReport", rules: tuple[Rule, ...] = DEFAULT_RULES) -> TriageResult:
    """Apply the triage rules.

    If CDR annotation failed, CDRs can't be separated from the framework, so the rules
    can't be applied. Such sequences are marked "cannot evaluate" (STATUS_UNKNOWN) rather
    than passed as "no rule matched" — reporting something unseen as fine is the riskiest outcome.
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
    """Descriptions of the active rules. Included in the output so the criteria are visible."""
    return [rule.description for rule in rules]
