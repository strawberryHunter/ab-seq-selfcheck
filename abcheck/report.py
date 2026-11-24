"""주석 + liability + 물성을 하나의 결과로 합친다.

여기서 만든 dict 가 터미널 표·HTML·JSON 세 출력의 공통 입력이다. 한계 고지와
"종합 판정 없음" 문구를 결과 안에 넣어 두기 때문에, 어느 출력을 쓰든 고지가
빠질 수 없다 — 이게 이 모듈에서 가장 중요한 설계다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import LIMITATION_NOTICE, NO_VERDICT_NOTICE, __version__
from .annotate import Annotation, annotate
from .liabilities import LiabilityReport, find_liabilities
from .properties import DEFAULT_PH, Properties, compute_all
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
        return SequenceReport(
            seq_id=seq_id,
            sequence=validation.cleaned,
            ok=False,
            errors=list(validation.errors),
            warnings=list(validation.warnings),
        )

    cleaned = validation.cleaned
    ann = annotate(cleaned, validate_input=False)
    lia = find_liabilities(ann)
    props = compute_all(cleaned, ann.cdrs, ph)

    return SequenceReport(
        seq_id=seq_id,
        sequence=cleaned,
        ok=True,
        warnings=list(validation.warnings),
        annotation=ann,
        liabilities=lia,
        properties=props,
    )


def build_many(items: list[tuple[str, str]], ph: float = DEFAULT_PH) -> list[SequenceReport]:
    """(id, 서열) 목록을 한꺼번에 점검한다."""
    return [build(seq_id, seq, ph) for seq_id, seq in items]


__all__ = ["SequenceReport", "build", "build_many", "SequenceError"]
