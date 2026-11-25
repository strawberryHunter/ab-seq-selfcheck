#!/usr/bin/env python3
"""OAS 공개 데이터로 주석기의 정확도를 측정한다.

왜 이 스크립트가 있나
--------------------
합성 데이터는 생성 규칙이 곧 정답이라 "내가 만든 가정" 안에서만 맞는다. 실제 서열에는
합성 데이터에 넣지 않은 패턴이 있고, 실제로 이 저장소의 버그 네 개가 전부 그렇게 드러났다.
그래서 공개 데이터로 재는 절차를 코드로 남긴다 — 숫자를 남이 다시 확인할 수 있어야 한다.

무엇을 재나
----------
OAS paired 데이터는 VH/VL 이 짝지어져 있고 어느 쪽이 중쇄인지 알려져 있다. 그래서
**체인 판별은 정답이 있는 검증**이 된다. CDR 위치는 정답이 없으므로 아래로 대신한다.
  - 주석 성공률 — 앵커를 찾아 CDR 을 끊어냈는가
  - 신뢰 가능률 — 스스로 '믿을 만하다'고 한 비율
  - 실패 사유 분포 — 못 한 경우 왜 못 했다고 말하는가
  - FR3↔보존Cys 거리 — 구조적으로 일정해야 하는 값이 실제로 일정한가 (내부 일관성)

쓰는 법
------
    # 1) OAS 에서 paired 데이터 단위(.csv.gz)를 받는다
    #    https://opig.stats.ox.ac.uk/webapps/oas/  →  Paired Sequences → Search
    # 2) FASTA 로 바꾼다
    python3 scripts/oas_to_fasta.py 받은파일.csv.gz 후보.fasta --count 500
    # 3) 측정한다
    python3 scripts/validate_on_oas.py 후보.fasta [다른.fasta ...]

받은 서열 파일은 저장소에 커밋하지 말 것. 이 저장소는 합성 데이터만 담는다.
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
    """oas_to_fasta.py 가 붙인 이름에서 정답 체인을 읽는다 (…-VH / …-VL)."""
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
