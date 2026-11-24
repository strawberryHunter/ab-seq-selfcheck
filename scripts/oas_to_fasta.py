#!/usr/bin/env python3
"""OAS(Observed Antibody Space) paired 데이터 단위를 abcheck 가 읽는 FASTA 로 바꾼다.

OAS 는 공개 항체 레퍼토리 데이터베이스다 (CC-BY 4.0).
    https://opig.stats.ox.ac.uk/webapps/oas/

받는 방법
--------
1. 위 사이트에서 Paired Sequences → Search 로 조건을 고르고 검색한다.
2. 결과 페이지의 'A shell-script with the commands to download...' 로 받은 스크립트 안에
   데이터 단위 URL 목록이 들어 있다. 그중 하나를 골라 받는다.
3. 받은 .csv.gz 를 이 스크립트에 넣는다.

    python3 scripts/oas_to_fasta.py 받은파일.csv.gz 후보.fasta --count 20
    python3 -m abcheck.cli 후보.fasta --html

주의
----
- OAS 의 `*_alignment_aa` 컬럼은 IMGT 정렬이라 '.' 갭이 들어 있다. 그대로 넣으면
  아미노산이 아닌 문자로 거부당하므로 갭을 제거한다.
- 받은 서열 파일은 저장소에 커밋하지 말 것. 이 저장소는 합성 데이터만 담는다.
- 파일 첫 줄은 메타데이터 JSON 이고 둘째 줄부터 CSV 헤더다.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

#: OAS 의 IMGT 정렬 갭 문자
GAP_CHARS = ".-"


def read_metadata(path: Path) -> dict:
    """데이터 단위 첫 줄의 메타데이터 JSON 을 읽는다."""
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
    """VH/VL 쌍을 뽑아낸다."""
    pairs: list[dict] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()                       # 메타데이터 줄 건너뛰기
        for row in csv.DictReader(fh):
            if productive_only and (
                row.get("productive_heavy") != "T" or row.get("productive_light") != "T"
            ):
                continue

            vh = _clean(row.get("sequence_alignment_aa_heavy") or "")
            vl = _clean(row.get("sequence_alignment_aa_light") or "")
            if not vh or not vl:
                continue
            if "*" in vh or "*" in vl:      # 종결코돈이 섞인 것은 버린다
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
    raise SystemExit(main())
