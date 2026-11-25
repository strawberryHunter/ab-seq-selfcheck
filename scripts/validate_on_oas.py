#!/usr/bin/env python3
"""Measure the annotator's accuracy on public OAS data.

Why this script exists
----------------------
With synthetic data the generation rules are the ground truth, so it is only right within "assumptions I made".
Real sequences have patterns never put into the synthetic data, and in fact all four bugs in this repository surfaced that way.
So the procedure for measuring on public data is kept as code -- others must be able to re-check the numbers.

What it measures
----------------
OAS paired data has VH/VL paired up, and which one is the heavy chain is known. So
**chain detection becomes validation against ground truth**. CDR positions have no ground truth, so these stand in:
  - annotation success rate -- were anchors found and CDRs cut out?
  - confident rate -- fraction it judged 'trustworthy' itself
  - failure reason distribution -- when it failed, what reason does it give?
  - FR3<->conserved Cys distance -- is a structurally constant value actually constant? (internal consistency)

Usage
-----
    # 1) Download a paired data unit (.csv.gz) from OAS
    #    https://opig.stats.ox.ac.uk/webapps/oas/  ->  Paired Sequences -> Search
    # 2) Convert to FASTA
    python3 scripts/oas_to_fasta.py downloaded.csv.gz candidates.fasta --count 500
    # 3) Measure
    python3 scripts/validate_on_oas.py candidates.fasta [other.fasta ...]

Do not commit downloaded sequence files to the repository. This repository holds synthetic data only.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from abcheck.annotate import FR3_TO_CYS_GAP, annotate  # noqa: E402
from abcheck.validate import validate  # noqa: E402


def read_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    name, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(chunks)))
            name, chunks = line[1:].split()[0], []
        elif line:
            chunks.append(line)
    if name is not None:
        entries.append((name, "".join(chunks)))
    return entries


def truth_chain(seq_id: str) -> str | None:
    """Read the true chain from the name oas_to_fasta.py assigned (...-VH / ...-VL)."""
    if seq_id.endswith("-VH"):
        return "H"
    if seq_id.endswith("-VL"):
        return "L"
    return None


def measure(entries: list[tuple[str, str]]) -> dict:
    stats = collections.Counter()
    reasons: collections.Counter = collections.Counter()
    errors: collections.Counter = collections.Counter()
    gaps: collections.Counter = collections.Counter()
    lengths: dict[str, list[int]] = collections.defaultdict(list)

    for seq_id, raw in entries:
        stats["total"] += 1
        checked = validate(raw)
        if not checked.ok:
            for err in checked.errors:
                errors[err.split("(")[0].strip()[:50]] += 1
            continue
        stats["valid"] += 1

        ann = annotate(checked.cleaned, validate_input=False)
        truth = truth_chain(seq_id)
        if truth is not None:
            if ann.chain == truth:
                stats["chain_right"] += 1
            elif ann.chain is None:
                stats["chain_unknown"] += 1
            else:
                stats["chain_wrong"] += 1

        if ann.cdrs:
            stats["annotated"] += 1
            for name, cdr in ann.cdrs.items():
                lengths[name].append(cdr.length)
        if ann.confident:
            stats["confident"] += 1
        for reason in ann.reasons:
            reasons[reason.split(".")[0][:60]] += 1
        if "fr3_motif" in ann.anchors and "cys2" in ann.anchors:
            gaps[ann.anchors["cys2"] - ann.anchors["fr3_motif"]] += 1

    return {
        "stats": stats,
        "reasons": reasons,
        "errors": errors,
        "gaps": gaps,
        "lengths": lengths,
    }


def _pct(n: int, d: int) -> str:
    return f"{n * 100 / d:.1f}%" if d else "-"


def report(name: str, result: dict) -> None:
    s = result["stats"]
    total, valid = s["total"], s["valid"]
    print(f"\n{'=' * 72}\n{name} — {total} 서열\n{'=' * 72}")
    print(f"  입력 검증 통과   {valid:5d} / {total}  {_pct(valid, total)}")
    print(f"  CDR 주석 성공    {s['annotated']:5d} / {total}  {_pct(s['annotated'], total)}")
    print(f"  신뢰 가능        {s['confident']:5d} / {total}  {_pct(s['confident'], total)}")
    graded = s["chain_right"] + s["chain_wrong"] + s["chain_unknown"]
    if graded:
        print(f"  체인 판별 정확   {s['chain_right']:5d} / {graded}  "
              f"{_pct(s['chain_right'], graded)}"
              f"   (오판 {s['chain_wrong']}, 판별 실패 {s['chain_unknown']})")

    gaps = result["gaps"]
    if gaps:
        tot = sum(gaps.values())
        same = gaps.get(FR3_TO_CYS_GAP, 0)
        print(f"  FR3↔보존Cys {FR3_TO_CYS_GAP}잔기  {same}/{tot}  {_pct(same, tot)}"
              + (f"   그 외: {dict(sorted((k, v) for k, v in gaps.items() if k != FR3_TO_CYS_GAP))}"
                 if same != tot else ""))

    if result["errors"]:
        print("\n  [입력 거부 사유]")
        for text, n in result["errors"].most_common(5):
            print(f"    {n:4d}  {text}")

    if result["reasons"]:
        print("\n  [신뢰 불가 사유]")
        for text, n in result["reasons"].most_common(8):
            print(f"    {n:4d}  {text}")

    if result["lengths"]:
        print("\n  [CDR 길이]  최소-최대 (중앙)")
        for cdr in sorted(result["lengths"]):
            vals = sorted(result["lengths"][cdr])
            print(f"    {cdr}  {vals[0]:2d}-{vals[-1]:2d}  "
                  f"(중앙 {vals[len(vals) // 2]})  n={len(vals)}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="OAS FASTA 로 주석기 정확도를 측정한다.",
        epilog="먼저 scripts/oas_to_fasta.py 로 OAS 데이터 단위를 FASTA 로 바꿀 것.",
    )
    parser.add_argument("fasta", nargs="+", type=Path, help="측정할 FASTA 파일들")
    parser.add_argument("--no-detail", action="store_true", help="파일별 상세를 생략한다")
    args = parser.parse_args(argv)

    merged: list[tuple[str, str]] = []
    for path in args.fasta:
        if not path.exists():
            print(f"파일을 찾을 수 없다: {path}", file=sys.stderr)
            return 2
        entries = read_fasta(path)
        merged += entries
        if not args.no_detail:
            report(path.stem, measure(entries))

    if len(args.fasta) > 1:
        report(f"전체 합계 ({len(args.fasta)}개 파일)", measure(merged))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
