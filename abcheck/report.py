"""주석 + liability + 물성을 하나의 결과로 합친다.

여기서 만든 dict 가 터미널 표·HTML·JSON 세 출력의 공통 입력이다. 한계 고지와
"종합 판정 없음" 문구를 결과 안에 넣어 두기 때문에, 어느 출력을 쓰든 고지가
빠질 수 없다 — 이게 이 모듈에서 가장 중요한 설계다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import LIMITATION_NOTICE, NO_VERDICT_NOTICE, __version__
from .annotate import Annotation, annotate
from .liabilities import RISK_LEVELS, LiabilityReport, find_liabilities
from .properties import DEFAULT_PH, Properties, compute_all
from .triage import TRIAGE_NOTICE, TriageResult, evaluate, rule_descriptions
from .validate import SequenceError, validate


@dataclass
class SequenceReport:
    """서열 한 건의 점검 결과."""

    seq_id: str
    sequence: str
    ok: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    annotation: Annotation | None = None
    liabilities: LiabilityReport | None = None
    properties: dict[str, Properties] = field(default_factory=dict)
    triage: TriageResult | None = None

    # ---- 편의 접근자 ----

    @property
    def chain(self) -> str | None:
        return self.annotation.chain if self.annotation else None

    @property
    def confident(self) -> bool:
        return bool(self.annotation and self.annotation.confident)

    def as_dict(self) -> dict:
        """JSON 출력용 dict. 한계 고지를 반드시 포함한다."""
        payload: dict = {
            "id": self.seq_id,
            "tool": {"name": "abcheck", "version": __version__},
            "input_ok": self.ok,
            "errors": self.errors,
            "warnings": list(self.warnings),
            "limitations": list(LIMITATION_NOTICE),
            "no_verdict": NO_VERDICT_NOTICE,
        }
        if self.triage is not None:
            payload["triage"] = self.triage.as_dict()
        if not self.ok:
            return payload

        ann = self.annotation
        assert ann is not None
        payload["sequence"] = self.sequence
        payload["length"] = len(self.sequence)
        payload["annotation"] = {
            "chain": ann.chain,
            "confident": ann.confident,
            "reasons": list(ann.reasons),
            "warnings": list(ann.warnings),
            "anchors": dict(ann.anchors),
            "coordinates": "1-based, 양끝 포함",
            "cdrs": {
                name: {
                    "start": cdr.start,
                    "end": cdr.end,
                    "length": cdr.length,
                    "seq": cdr.seq,
                }
                for name, cdr in sorted(ann.cdrs.items(), key=lambda kv: kv[1].start)
            },
        }

        lia = self.liabilities
        assert lia is not None
        payload["liabilities"] = {
            "warnings": list(lia.warnings),
            "counts_by_risk": lia.counts_by_risk(),
            "hits": [
                {
                    "kind": h.kind,
                    "label": h.label,
                    "match": h.match,
                    "start": h.start,
                    "end": h.end,
                    "region": h.region,
                    "in_cdr": h.in_cdr,
                    "risk": h.risk,
                    "base_risk": h.base_risk,
                    "rationale": h.rationale,
                    "risk_note": h.note,
                }
                for h in lia.hits
            ],
            "domain_flags": list(lia.domain_flags),
        }

        payload["properties"] = {
            name: props.as_dict() for name, props in self.properties.items()
        }
        return payload


def build(seq_id: str, sequence: str, ph: float = DEFAULT_PH) -> SequenceReport:
    """서열 하나를 점검한다. 입력이 잘못돼도 예외를 올리지 않는다."""
    validation = validate(sequence)
    if not validation.ok:
        failed = SequenceReport(
            seq_id=seq_id,
            sequence=validation.cleaned,
            ok=False,
            errors=list(validation.errors),
            warnings=list(validation.warnings),
        )
        failed.triage = evaluate(failed)
        return failed

    cleaned = validation.cleaned
    ann = annotate(cleaned, validate_input=False)
    lia = find_liabilities(ann)
    props = compute_all(cleaned, ann.cdrs, ph)

    report = SequenceReport(
        seq_id=seq_id,
        sequence=cleaned,
        ok=True,
        warnings=list(validation.warnings),
        annotation=ann,
        liabilities=lia,
        properties=props,
    )
    report.triage = evaluate(report)
    return report


def build_many(items: list[tuple[str, str]], ph: float = DEFAULT_PH) -> list[SequenceReport]:
    """(id, 서열) 목록을 한꺼번에 점검한다."""
    return [build(seq_id, seq, ph) for seq_id, seq in items]


# --------------------------------------------------------------------------
# 집계
# --------------------------------------------------------------------------

#: 집계 표의 구간 칸. 중쇄(CDR-H1)와 경쇄(CDR-L1)를 같은 칸에 모으기 위해 번호만 쓴다.
SUMMARY_REGIONS: tuple[str, ...] = ("CDR1", "CDR2", "CDR3", "framework")

SUMMARY_NOTE: str = (
    "아래는 건수 집계일 뿐이며 종합 판정이 아니다. 건수가 적다고 안전한 것도, "
    "많다고 탈락시킬 근거도 되지 않는다. 각 서열의 상세 항목과 근거를 함께 볼 것."
)

#: 주석이 실패해 구간을 가를 수 없는 서열의 모티프가 들어가는 칸
UNASSIGNED_REGION = "framework"

#: 집계 표를 체인별로 나눌 때 쓰는 제목. None 은 체인을 판별하지 못했거나 입력이 잘못된 서열.
CHAIN_GROUP_TITLES: dict[str | None, str] = {
    "H": "중쇄 (VH)",
    "L": "경쇄 (VL)",
    None: "체인 판별 실패 · 입력 오류",
}


def cdr_labels(chain: str | None) -> list[str]:
    """체인별 CDR 열 이름. 표를 나누면 CDR-H1 / CDR-L1 로 정확히 쓸 수 있다."""
    if chain in ("H", "L"):
        return [f"CDR-{chain}{i}" for i in (1, 2, 3)]
    return ["CDR1", "CDR2", "CDR3"]


def _slot_of(region: str) -> str:
    """liability hit 의 구간 이름을 집계 칸으로 바꾼다.

    "CDR-H2" → "CDR2", "CDR-L1" → "CDR1", "framework" → "framework".
    "구간 미확정"(주석 실패)은 framework 칸으로 모으고, 행에 별도 표시를 단다.
    """
    if region.startswith("CDR-") and region[-1].isdigit():
        return f"CDR{region[-1]}"
    if region == "구간 미확정":
        return UNASSIGNED_REGION
    return "framework"


def _empty_counts() -> dict[str, int]:
    return {level: 0 for level in reversed(RISK_LEVELS)}


def summarize(reports: list[SequenceReport]) -> dict:
    """서열 × 구간 × 위험도 건수를 교차 집계한다.

    후보 여러 건을 나란히 비교할 때 쓴다. 점수를 매기거나 순위를 내지 않는다 —
    건수만 센다.
    """
    rows: list[dict] = []
    totals: dict[str, dict[str, int]] = {r: _empty_counts() for r in SUMMARY_REGIONS}
    flag_total = 0

    for report in reports:
        regions: dict[str, dict[str, int]] = {r: _empty_counts() for r in SUMMARY_REGIONS}
        unassigned = False

        if report.ok and report.liabilities is not None:
            for hit in report.liabilities.hits:
                if hit.region == "구간 미확정":
                    unassigned = True
                slot = _slot_of(hit.region)
                regions[slot][hit.risk] += 1
                totals[slot][hit.risk] += 1

        flags = [
            {"label": f["label"], "risk": f["risk"]}
            for f in (report.liabilities.domain_flags if report.ok and report.liabilities else [])
        ]
        flag_total += len(flags)

        rows.append({
            "id": report.seq_id,
            "input_ok": report.ok,
            "chain": report.chain,
            "confident": report.confident,
            "annotated": bool(report.annotation and report.annotation.cdrs),
            "regions_unassigned": unassigned,
            "regions": regions,
            "domain_flags": flags,
            "triage": report.triage.as_dict() if report.triage else None,
        })

    return {
        "note": SUMMARY_NOTE,
        "regions": list(SUMMARY_REGIONS),
        "risk_levels": list(reversed(RISK_LEVELS)),
        "rows": rows,
        "totals": totals,
        "domain_flag_total": flag_total,
        "triage_rules": rule_descriptions(),
        "triage_notice": TRIAGE_NOTICE,
        "review_count": sum(
            1 for r in rows
            if r["triage"] and r["triage"]["status"] == "검토 필요"
        ),
        "groups": _group_by_chain(rows),
    }


def _group_by_chain(rows: list[dict]) -> list[dict]:
    """집계 행을 중쇄 / 경쇄 / 판별 실패로 나눈다.

    표를 나누면 열 이름을 CDR-H1, CDR-L1 처럼 정확히 쓸 수 있고, 중쇄에만 있는
    CDR-H3 길이 특성 같은 것을 경쇄와 섞어 보지 않게 된다. 빈 묶음은 돌려주지 않는다.
    """
    groups: list[dict] = []

    for chain in ("H", "L", None):
        members = [r for r in rows if (r["chain"] if r["input_ok"] else None) == chain]
        if not members:
            continue

        totals = {region: _empty_counts() for region in SUMMARY_REGIONS}
        for row in members:
            for region in SUMMARY_REGIONS:
                for level, n in row["regions"][region].items():
                    totals[region][level] += n

        groups.append({
            "chain": chain,
            "title": CHAIN_GROUP_TITLES[chain],
            "cdr_labels": cdr_labels(chain),
            "rows": members,
            "totals": totals,
            "domain_flag_total": sum(len(r["domain_flags"]) for r in members),
        })
    return groups


__all__ = [
    "SequenceReport",
    "build",
    "build_many",
    "summarize",
    "cdr_labels",
    "SUMMARY_REGIONS",
    "SUMMARY_NOTE",
    "CHAIN_GROUP_TITLES",
    "SequenceError",
]
