#!/usr/bin/env python3
"""Convert an OAS (Observed Antibody Space) paired data unit into FASTA that abcheck reads.

OAS is a public antibody repertoire database (CC-BY 4.0).
    https://opig.stats.ox.ac.uk/webapps/oas/

How to download
---------------
1. On the site above, pick criteria under Paired Sequences -> Search and run the search.
2. The script obtained from 'A shell-script with the commands to download...' on the results page
   contains a list of data unit URLs. Pick one and download it.
3. Feed the downloaded .csv.gz to this script.

    python3 scripts/oas_to_fasta.py downloaded.csv.gz candidates.fasta --count 20
    python3 -m abcheck.cli candidates.fasta --html

Notes
-----
- OAS's `*_alignment_aa` columns are IMGT alignments and contain '.' gaps. Passed as-is they
  would be rejected as non-amino-acid characters, so the gaps are removed.
- Do not commit downloaded sequence files to the repository. This repository holds synthetic data only.
- The first line of the file is metadata JSON; the CSV header starts on the second line.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

#: OAS's IMGT alignment gap character
GAP_CHARS = ".-"


def read_metadata(path: Path) -> dict:
    """Read the metadata JSON on the first line of a data unit."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        line = fh.readline().strip()
    if line.startswith('"') and line.endswith('"'):
        line = line[1:-1].replace('""', '"')
    try:
        return json.loads(line)
    except json.JSONDecodeError:
        return {}


def _clean(seq: str) -> str:
    for ch in GAP_CHARS:
        seq = seq.replace(ch, "")
    return seq.strip().upper()


def extract(path: Path, count: int, productive_only: bool = True) -> list[dict]:
    """Extract VH/VL pairs."""
    pairs: list[dict] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()                       # Skip the metadata line
        for row in csv.DictReader(fh):
            if productive_only and (
                row.get("productive_heavy") != "T" or row.get("productive_light") != "T"
            ):
                continue

            vh = _clean(row.get("sequence_alignment_aa_heavy") or "")
            vl = _clean(row.get("sequence_alignment_aa_light") or "")
            if not vh or not vl:
                continue
            if "*" in vh or "*" in vl:      # Discard sequences containing stop codons
                continue

            pairs.append({
                "vh": vh,
                "vl": vl,
                "v_heavy": (row.get("v_call_heavy") or "").split(",")[0],
                "v_light": (row.get("v_call_light") or "").split(",")[0],
                "locus_light": row.get("locus_light") or "",
            })
            if len(pairs) >= count:
                break
    return pairs


def write_fasta(pairs: list[dict], out: Path, prefix: str) -> None:
    lines: list[str] = []
    for i, p in enumerate(pairs, start=1):
        lines.append(f">{prefix}-{i:02d}-VH len={len(p['vh'])} {p['v_heavy']}")
        lines.append(p["vh"])
        lines.append(f">{prefix}-{i:02d}-VL len={len(p['vl'])} {p['locus_light']} {p['v_light']}")
        lines.append(p["vl"])
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="OAS paired 데이터 단위(.csv.gz)를 abcheck 용 FASTA 로 바꾼다.",
    )
    parser.add_argument("source", type=Path, help="OAS 에서 받은 .csv.gz 파일")
    parser.add_argument("output", type=Path, help="만들 FASTA 파일 경로")
    parser.add_argument("--count", type=int, default=10, help="뽑을 VH/VL 쌍 수 (기본 10)")
    parser.add_argument("--prefix", default="OAS", help="서열 이름 접두사 (기본 OAS)")
    parser.add_argument(
        "--all", action="store_true",
        help="productive 가 아닌 서열도 포함한다 (기본은 productive 만)",
    )
    args = parser.parse_args(argv)

    if not args.source.exists():
        print(f"파일을 찾을 수 없다: {args.source}", file=sys.stderr)
        return 2

    meta = read_metadata(args.source)
    if meta:
        fields = ["Species", "BSource", "BType", "Disease", "Vaccine", "Age"]
        summary = ", ".join(f"{k}={meta[k]}" for k in fields if meta.get(k))
        print(f"데이터 단위: {summary}")
        if meta.get("Link"):
            print(f"출처: {meta['Link']}")

    pairs = extract(args.source, args.count, productive_only=not args.all)
    if not pairs:
        print("뽑아낼 쌍이 없다. paired 데이터 단위가 맞는지 확인하라.", file=sys.stderr)
        return 1

    write_fasta(pairs, args.output, args.prefix)
    print(f"{len(pairs)} 쌍 ({len(pairs) * 2} 서열) → {args.output}")

    vh_len = [len(p["vh"]) for p in pairs]
    vl_len = [len(p["vl"]) for p in pairs]
    print(f"  VH 길이 {min(vh_len)}-{max(vh_len)} / VL 길이 {min(vl_len)}-{max(vl_len)}")
    print(f"\n다음 단계:  python3 -m abcheck.cli {args.output} --html")
    return 0


if __name__ == "__main__":
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    from abcheck._streams import configure_streams

    configure_streams()      # Korean output must not crash a Windows pipe
    raise SystemExit(main())
