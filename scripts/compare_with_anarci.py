#!/usr/bin/env python3
"""Compare this tool's CDR boundaries with the ANARCI annotations in OAS data.

Why 'offsets' rather than an 'agreement rate'
---------------------------------------------
OAS's cdr*_aa_* columns were numbered by ANARCI in **IMGT**, while this tool approximates
**Kabat** boundaries. The two schemes define boundaries differently, so the strings cannot be
identical, nor should they be. Forcing an 'agreement rate' gives a low number, but that does not mean it is wrong.

So we look at it this way:
  1. Find where ANARCI's CDR string sits in the sequence to get the reference coordinates.
  2. Compute the offset (my start - reference start, my end - reference end).
  3. **If the offset is constant across sequences**, it is just a definitional difference and my
     boundaries are systematically right. **The inconsistent ones** are the real errors.

In other words, "the fraction deviating from the modal offset" is this tool's CDR boundary error rate.

Usage
-----
    python3 scripts/compare_with_anarci.py downloaded.csv.gz [--count 500]

Do not commit downloaded sequence files to the repository.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from abcheck.annotate import annotate  # noqa: E402
from abcheck.validate import validate  # noqa: E402

#: This tool's CDR names <-> OAS column suffixes
PAIRS_HEAVY = [("CDR-H1", "cdr1_aa_heavy"), ("CDR-H2", "cdr2_aa_heavy"),
               ("CDR-H3", "cdr3_aa_heavy")]
PAIRS_LIGHT = [("CDR-L1", "cdr1_aa_light"), ("CDR-L2", "cdr2_aa_light"),
               ("CDR-L3", "cdr3_aa_light")]


def _family(v_call: str | None) -> str:
    """IGHV4-39*01 -> IGHV4. The Kabat<->IMGT relationship differs per germline family."""
    if not v_call:
        return "?"
    return v_call.split(",")[0].split("-")[0].split("*")[0]


def _clean(seq: str) -> str:
    return seq.replace(".", "").replace("-", "").strip().upper()


def read_rows(path: Path, count: int):
    """Yield (sequence, {CDR name: ANARCI string}) pairs for each of the heavy and light chains."""
    out = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()                       # Metadata line
        for row in csv.DictReader(fh):
            if row.get("productive_heavy") != "T" or row.get("productive_light") != "T":
                continue
            vh = _clean(row.get("sequence_alignment_aa_heavy") or "")
            vl = _clean(row.get("sequence_alignment_aa_light") or "")
            if not vh or not vl or "*" in vh or "*" in vl:
                continue
            truth_h = {n: _clean(row.get(c) or "") for n, c in PAIRS_HEAVY}
            truth_l = {n: _clean(row.get(c) or "") for n, c in PAIRS_LIGHT}
            fam_h = _family(row.get("v_call_heavy"))
            fam_l = _family(row.get("v_call_light"))
            out.append((vh, truth_h, fam_h, vl, truth_l, fam_l))
            if len(out) >= count:
                break
    return out


def locate(seq: str, fragment: str) -> tuple[int, int] | None:
    """1-based coordinates of the ANARCI string. Used only when it occurs exactly once in the sequence."""
    if not fragment:
        return None
    first = seq.find(fragment)
    if first < 0 or seq.find(fragment, first + 1) >= 0:
        return None                          # If absent or occurring multiple times, the coordinates can't be trusted
    return first + 1, first + len(fragment)


def compare(rows) -> dict:
    offsets: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    by_len: dict[str, dict[int, collections.Counter]] = collections.defaultdict(
        lambda: collections.defaultdict(collections.Counter)
    )
    by_group: dict[str, dict[tuple, collections.Counter]] = collections.defaultdict(
        lambda: collections.defaultdict(collections.Counter)
    )
    compared: collections.Counter = collections.Counter()
    skipped: collections.Counter = collections.Counter()
    outliers: dict[str, list] = collections.defaultdict(list)

    for vh, truth_h, fam_h, vl, truth_l, fam_l in rows:
        for seq, truth, family in ((vh, truth_h, fam_h), (vl, truth_l, fam_l)):
            checked = validate(seq)
            if not checked.ok:
                skipped["입력 검증 실패"] += 1
                continue
            ann = annotate(checked.cleaned, validate_input=False)
            if not ann.cdrs:
                skipped["주석 실패"] += 1
                continue
            for name, fragment in truth.items():
                if name not in ann.cdrs:
                    continue
                spot = locate(checked.cleaned, fragment)
                if spot is None:
                    skipped[f"{name}: 정답 위치 특정 불가"] += 1
                    continue
                mine = ann.cdrs[name]
                pair = (mine.start - spot[0], mine.end - spot[1])
                offsets[name][pair] += 1
                # IMGT CDR lengths vary, but Kabat CDR-H1/L2 have fixed lengths.
                # So it is normal for the offset to vary **with the reference length**.
                # Only what still varies with length held fixed is a real boundary error.
                by_len[name][len(fragment)][pair] += 1
                by_group[name][(family, len(fragment))][pair] += 1
                compared[name] += 1
                outliers[name].append((*pair, mine.seq, fragment))
    return {"offsets": offsets, "compared": compared, "by_len": by_len,
            "by_group": by_group, "skipped": skipped, "outliers": outliers}


def report(result: dict, show_examples: int = 3) -> None:
    print(f"\n{'CDR':8s} {'비교':>6s} {'최빈 오프셋':>12s} {'일치':>8s} {'벗어남':>8s}")
    print("-" * 52)
    total_c = total_off = 0
    for name in ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]:
        counts = result["offsets"].get(name)
        if not counts:
            continue
        n = sum(counts.values())
        (ds, de), top = counts.most_common(1)[0]
        off = n - top
        total_c += n
        total_off += off
        print(f"{name:8s} {n:6d} {f'({ds:+d}, {de:+d})':>12s} "
              f"{f'{top / n * 100:.1f}%':>8s} {off:8d}")
    if total_c:
        print("-" * 52)
        print(f"{'합계':8s} {total_c:6d} {'':>12s} "
              f"{f'{(total_c - total_off) / total_c * 100:.2f}%':>8s} {total_off:8d}")
        print(f"\n→ 최빈 오프셋에서 벗어난 비율 = CDR 경계 오류율: "
              f"{total_off * 100 / total_c:.2f}%")

    print("\n[오프셋 분포]  (내 시작 - 정답 시작, 내 끝 - 정답 끝)")
    for name in ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]:
        counts = result["offsets"].get(name)
        if not counts:
            continue
        top = counts.most_common(4)
        shown = "  ".join(f"({a:+d},{b:+d})×{v}" for (a, b), v in top)
        more = f"  …외 {len(counts) - len(top)}종" if len(counts) > len(top) else ""
        print(f"  {name}  {shown}{more}")

    # Recount with the reference length held fixed -- this is the real boundary error rate
    print("\n[정답 CDR 길이를 고정했을 때]  IMGT 길이는 가변, Kabat 은 고정이라 "
          "길이가 다르면 오프셋도 달라지는 것이 정상이다")
    g_tot = g_bad = 0
    for name in ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]:
        groups = result["by_len"].get(name)
        if not groups:
            continue
        n = bad = 0
        for length, counts in groups.items():
            sub = sum(counts.values())
            n += sub
            bad += sub - counts.most_common(1)[0][1]
        g_tot += n
        g_bad += bad
        mark = "" if bad == 0 else f"   ← {bad}건"
        print(f"  {name}  길이 {len(groups)}종  일치 {(n - bad) * 100 / n:.1f}%{mark}")
    if g_tot:
        print(f"\n→ 길이를 고정했을 때 경계 오류율: {g_bad * 100 / g_tot:.2f}% "
              f"({g_bad}/{g_tot})")

    # Also hold the germline family fixed -- insertion positions differ per family, so the Kabat<->IMGT relationship differs
    print("\n[germline 계열 + 정답 길이를 모두 고정했을 때]")
    t_tot = t_bad = 0
    for name in ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]:
        groups = result["by_group"].get(name)
        if not groups:
            continue
        n = bad = 0
        for counts in groups.values():
            sub = sum(counts.values())
            n += sub
            bad += sub - counts.most_common(1)[0][1]
        t_tot += n
        t_bad += bad
        mark = "" if bad == 0 else f"   ← {bad}건"
        print(f"  {name}  묶음 {len(groups)}개  일치 {(n - bad) * 100 / n:.1f}%{mark}")
    if t_tot:
        print(f"\n→ 계열·길이를 고정했을 때 경계 불일치: {t_bad * 100 / t_tot:.2f}% "
              f"({t_bad}/{t_tot})")
        print("  이 값이 0 에 가까우면, 앞의 오프셋 차이는 Kabat↔IMGT 체계 차이일 뿐이고")
        print("  경계 자체는 germline 계열 안에서 일정하다는 뜻이다.")

    if result["skipped"]:
        print("\n[비교 제외]")
        for reason, n in result["skipped"].most_common(8):
            print(f"  {n:5d}  {reason}")

    if show_examples:
        print("\n[최빈에서 벗어난 예시]")
        for name in ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]:
            counts = result["offsets"].get(name)
            if not counts:
                continue
            common = counts.most_common(1)[0][0]
            odd = [r for r in result["outliers"][name] if (r[0], r[1]) != common]
            for ds, de, mine, truth in odd[:show_examples]:
                print(f"  {name}  오프셋({ds:+d},{de:+d})  내것={mine!r}  ANARCI={truth!r}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="OAS 안의 ANARCI 주석과 이 도구의 CDR 경계를 비교한다.",
    )
    parser.add_argument("source", type=Path, help="OAS paired 데이터 단위 (.csv.gz)")
    parser.add_argument("--count", type=int, default=500, help="비교할 VH/VL 쌍 수")
    parser.add_argument("--examples", type=int, default=3, help="보여줄 이상치 예시 수")
    args = parser.parse_args(argv)

    if not args.source.exists():
        print(f"파일을 찾을 수 없다: {args.source}", file=sys.stderr)
        return 2

    rows = read_rows(args.source, args.count)
    if not rows:
        print("비교할 쌍이 없다. paired 데이터 단위가 맞는지 확인하라.", file=sys.stderr)
        return 1
    print(f"{args.source.name}: {len(rows)} 쌍 ({len(rows) * 2} 서열) 비교")
    report(compare(rows), args.examples)
    return 0


if __name__ == "__main__":
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    from abcheck._streams import configure_streams

    configure_streams()      # Korean output must not crash a Windows pipe
    raise SystemExit(main())
