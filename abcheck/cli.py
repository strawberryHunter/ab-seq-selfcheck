"""명령줄 인터페이스.

    python -m abcheck.cli <입력> [--html] [--json] [--out-dir DIR]

입력은 세 가지를 자동 판별한다.
  - `.json`  — 합성 코퍼스 형식 ({"sequences": [{"id", "sequence"}, ...]})
  - `.fasta` / `.fa` / `.faa` — FASTA 파일
  - 그 외    — 아미노산 서열 문자열 그대로 (붙여넣기용)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .properties import DEFAULT_PH
from .render import as_html, as_json, as_text
from .report import SequenceReport, build_many

JSON_SUFFIXES = {".json"}
FASTA_SUFFIXES = {".fasta", ".fa", ".faa", ".fas"}


# --------------------------------------------------------------------------
# 입력 읽기
# --------------------------------------------------------------------------

def read_fasta(text: str) -> list[tuple[str, str]]:
    """FASTA 텍스트를 (id, 서열) 목록으로 바꾼다."""
    entries: list[tuple[str, str]] = []
    seq_id: str | None = None
    chunks: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id is not None:
                entries.append((seq_id, "".join(chunks)))
            seq_id = line[1:].split()[0] if len(line) > 1 else f"seq{len(entries) + 1}"
            chunks = []
        else:
            chunks.append(line)
    if seq_id is not None:
        entries.append((seq_id, "".join(chunks)))
    return entries


def read_corpus_json(text: str) -> list[tuple[str, str]]:
    """코퍼스 JSON 을 (id, 서열) 목록으로 바꾼다.

    {"sequences": [...]} 또는 최상위 리스트를 받는다. 각 항목은 "sequence" 키가
    있어야 하고, "id" 가 없으면 순번으로 이름을 붙인다.
    """
    data = json.loads(text)
    items = data.get("sequences", data) if isinstance(data, dict) else data
    if not isinstance(items, list):
        raise ValueError(
            "JSON 형식을 알 수 없다. {\"sequences\": [{\"id\": ..., \"sequence\": ...}]} "
            "또는 같은 모양의 리스트여야 한다."
        )

    entries: list[tuple[str, str]] = []
    for i, item in enumerate(items, start=1):
        if isinstance(item, str):
            entries.append((f"seq{i}", item))
            continue
        if not isinstance(item, dict) or "sequence" not in item:
            raise ValueError(f"{i}번째 항목에 'sequence' 키가 없다.")
        entries.append((str(item.get("id", f"seq{i}")), str(item["sequence"])))
    return entries


def load_input(raw: str) -> tuple[list[tuple[str, str]], str]:
    """입력을 (id, 서열) 목록과 입력 종류 설명으로 바꾼다."""
    path = Path(raw)
    looks_like_path = len(raw) < 4096 and ("/" in raw or "." in raw or path.exists())

    if looks_like_path and path.exists():
        text = path.read_text(encoding="utf-8")
        suffix = path.suffix.lower()
        if suffix in JSON_SUFFIXES:
            return read_corpus_json(text), f"코퍼스 JSON ({path})"
        if suffix in FASTA_SUFFIXES or text.lstrip().startswith(">"):
            return read_fasta(text), f"FASTA ({path})"
        # 확장자를 모르면 서열 파일로 본다
        return [(path.stem, text)], f"서열 파일 ({path})"

    if looks_like_path and not path.exists() and path.suffix:
        raise FileNotFoundError(f"파일을 찾을 수 없다: {raw}")

    if raw.lstrip().startswith(">"):
        return read_fasta(raw), "붙여넣은 FASTA"
    return [("입력 서열", raw)], "붙여넣은 서열"


# --------------------------------------------------------------------------
# 실행
# --------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m abcheck.cli",
        description=(
            "항체 가변도메인(VH/VL) 서열의 CDR 위치와 developability 플래그를 점검한다. "
            "스크리닝 휴리스틱이며 numbering 엔진이 아니다."
        ),
        epilog=(
            "예시:\n"
            "  python -m abcheck.cli QVQL...SS\n"
            "  python -m abcheck.cli data/synthetic/sequences.json --html\n"
            "  python -m abcheck.cli my_seqs.fasta --json > result.json\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="아미노산 서열, FASTA 파일, 또는 코퍼스 JSON 파일",
    )
    parser.add_argument("--html", action="store_true", help="HTML 리포트를 파일로 만든다")
    parser.add_argument("--json", action="store_true", help="JSON 을 표준출력으로 낸다")
    parser.add_argument(
        "--out-dir", default="sample_output", type=Path,
        help="HTML 리포트를 쓸 디렉터리 (기본: sample_output)",
    )
    parser.add_argument(
        "--out-name", default="report.html",
        help="HTML 파일 이름 (기본: report.html)",
    )
    parser.add_argument(
        "--ph", default=DEFAULT_PH, type=float,
        help=f"net charge 를 계산할 pH (기본: {DEFAULT_PH})",
    )
    parser.add_argument("--version", action="version", version=f"abcheck {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    try:
        entries, source = load_input(args.input)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
        print(f"입력을 읽을 수 없다: {exc}", file=sys.stderr)
        return 2

    if not entries:
        print("처리할 서열이 없다.", file=sys.stderr)
        return 2

    reports: list[SequenceReport] = build_many(entries, args.ph)

    if args.json:
        sys.stdout.write(as_json(reports))
    elif not args.html:
        print(as_text(reports))

    if args.html:
        out_dir: Path = args.out_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / args.out_name
        out_path.write_text(as_html(reports), encoding="utf-8")
        print(f"HTML 리포트를 썼다: {out_path}", file=sys.stderr)

    rejected = [r for r in reports if not r.ok]
    if rejected:
        print("", file=sys.stderr)
        print(f"입력 검증에서 걸러진 서열 {len(rejected)}건:", file=sys.stderr)
        for report in rejected:
            print(f"  [{report.seq_id}]", file=sys.stderr)
            for err in report.errors:
                print(f"    - {err}", file=sys.stderr)
        return 1

    not_confident = [r for r in reports if not r.confident]
    if not_confident:
        print("", file=sys.stderr)
        print(
            f"주석을 신뢰할 수 없는 서열 {len(not_confident)}건 "
            f"(ANARCI/IMGT 재확인 필요): "
            + ", ".join(r.seq_id for r in not_confident),
            file=sys.stderr,
        )

    print(f"처리 완료: {source} — 서열 {len(reports)}건", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
