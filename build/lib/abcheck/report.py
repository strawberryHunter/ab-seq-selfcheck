"""Combine annotation + liabilities + properties into a single result.

The dict built here is the shared input for all three outputs: terminal table, HTML and
JSON. Because the limitation notice and the "no overall verdict" text live inside the
result itself, no output can drop them — that is the most important design choice in this module.
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
    """Check result for a single sequence."""

    seq_id: str
    sequence: str
    ok: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    annotation: Annotation | None = None
    liabilities: LiabilityReport | None = None
    properties: dict[str, Properties] = field(default_factory=dict)
    triage: TriageResult | None = None

    # ---- convenience accessors ----

    @property
    def chain(self) -> str | None:
        return self.annotation.chain if self.annotation else None

    @property
    def confident(self) -> bool:
        return bool(self.annotation and self.annotation.confident)

    def as_dict(self) -> dict:
        """Dict for JSON output. Always includes the limitation notice."""
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
    """Check a single sequence. Does not raise even on invalid input."""
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
    """Check a list of (id, sequence) in one go."""
    return [build(seq_id, seq, ph) for seq_id, seq in items]


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

#: Region columns of the summary table. Numbers only, so heavy (CDR-H1) and light (CDR-L1) share a column.
SUMMARY_REGIONS: tuple[str, ...] = ("CDR1", "CDR2", "CDR3", "framework")

SUMMARY_NOTE: str = (
    "아래는 건수 집계일 뿐이며 종합 판정이 아니다. 건수가 적다고 안전한 것도, "
    "많다고 탈락시킬 근거도 되지 않는다. 각 서열의 상세 항목과 근거를 함께 볼 것."
)

#: Column for motifs from sequences whose regions can't be split because annotation failed
UNASSIGNED_REGION = "framework"

#: Titles used when splitting the summary table by chain. None means chain detection failed or the input was invalid.
CHAIN_GROUP_TITLES: dict[str | None, str] = {
    "H": "중쇄 (VH)",
    "L": "경쇄 (VL)",
    None: "체인 판별 실패 · 입력 오류",
}


def cdr_labels(chain: str | None) -> list[str]:
    """CDR column names per chain. With split tables they can be named exactly, CDR-H1 / CDR-L1."""
    if chain in ("H", "L"):
        return [f"CDR-{chain}{i}" for i in (1, 2, 3)]
    return ["CDR1", "CDR2", "CDR3"]


def _slot_of(region: str) -> str:
    """Map a liability hit's region name to a summary column.

    "CDR-H2" → "CDR2", "CDR-L1" → "CDR1", "framework" → "framework".
    The undetermined-region marker (annotation failed) goes into the framework column and the row is marked separately.
    """
    if region.startswith("CDR-") and region[-1].isdigit():
        return f"CDR{region[-1]}"
    if region == "구간 미확정":
        return UNASSIGNED_REGION
    return "framework"


def _empty_counts() -> dict[str, int]:
    return {level: 0 for level in reversed(RISK_LEVELS)}


def summarize(reports: list[SequenceReport]) -> dict:
    """Cross-tabulate counts by sequence × region × risk level.

    Used to compare several candidates side by side. It does not score or rank —
    it only counts.
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
    """Split summary rows into heavy / light / detection failed.

    Split tables can name columns exactly (CDR-H1, CDR-L1), and heavy-only traits such as
    CDR-H3 length aren't mixed in with light chains. Empty groups are not returned.
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
