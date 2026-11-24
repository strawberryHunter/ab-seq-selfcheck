"""출력 — 터미널 표 / 단일 HTML 리포트 / JSON.

세 출력 모두 `LIMITATION_NOTICE` 와 `NO_VERDICT_NOTICE` 를 싣는다. 이 도구가
numbering 엔진이 아니라는 사실과, 종합 판정을 내리지 않는다는 사실이 결과만
보는 사람에게도 전달돼야 하기 때문이다.

HTML 은 의도적으로 수수하다 — Arial, 검정 헤더, 회색 지브라, 컬러 강조 없음.
CDR 구간은 굵게만 표시한다. 위험도는 색이 아니라 '높음/중간/낮음' 글자로 읽는다.
색으로 위험도를 칠하면 "빨간 게 몇 개인지"로 종합 판정을 하게 되는데, 그건
이 도구가 하지 않기로 한 일이다.
"""

from __future__ import annotations

import html
import json
import unicodedata
from datetime import date

from . import LIMITATION_NOTICE, NO_VERDICT_NOTICE, __version__
from .report import SUMMARY_REGIONS, SequenceReport, summarize

#: 서열을 보여줄 때 한 줄에 넣을 잔기 수
SEQUENCE_LINE_WIDTH = 60

#: 집계 표에서 건수가 0 인 칸에 쓸 표시
EMPTY_CELL = "-"

#: 집계 표 열 제목
SUMMARY_HEADERS: dict[str, str] = {
    "CDR1": "CDR1",
    "CDR2": "CDR2",
    "CDR3": "CDR3",
    "framework": "framework",
}


def _counts_cell(counts: dict[str, int]) -> str:
    """{'높음': 1, '중간': 0, '낮음': 2} → '높음 1, 낮음 2'."""
    parts = [f"{level} {n}" for level, n in counts.items() if n]
    return ", ".join(parts) if parts else EMPTY_CELL


def _summary_row_cells(row: dict) -> tuple[str, str, str, list[str], str]:
    """집계 표 한 행에 들어갈 칸들을 만든다."""
    if not row["input_ok"]:
        return row["id"], "-", "입력 오류", [EMPTY_CELL] * len(SUMMARY_REGIONS), EMPTY_CELL

    chain = {"H": "VH", "L": "VL", None: "판별 실패"}[row["chain"]]
    confidence = "신뢰 가능" if row["confident"] else "신뢰 불가"

    cells: list[str] = []
    for region in SUMMARY_REGIONS:
        if region.startswith("CDR") and not row["annotated"]:
            cells.append("주석 실패")
        else:
            cells.append(_counts_cell(row["regions"][region]))

    flags = row["domain_flags"]
    flag_cell = (
        ", ".join(f"{f['label']} ({f['risk']})" for f in flags) if flags else EMPTY_CELL
    )
    return row["id"], chain, confidence, cells, flag_cell


# --------------------------------------------------------------------------
# 터미널 폭 계산 — 한글은 두 칸을 차지한다
# --------------------------------------------------------------------------

def _display_width(text: str) -> int:
    return sum(2 if unicodedata.east_asian_width(ch) in "WF" else 1 for ch in text)


def _pad(text: str, width: int, align: str = "left") -> str:
    gap = max(0, width - _display_width(text))
    if align == "right":
        return " " * gap + text
    return text + " " * gap


def _table(headers: list[str], rows: list[list[str]], aligns: list[str] | None = None) -> list[str]:
    """고정폭 표를 줄 목록으로 만든다. 컬러는 쓰지 않는다."""
    if aligns is None:
        aligns = ["left"] * len(headers)
    widths = [
        max([_display_width(h)] + [_display_width(r[i]) for r in rows])
        for i, h in enumerate(headers)
    ]
    sep = "-+-".join("-" * w for w in widths)
    out = [
        " | ".join(_pad(h, widths[i], aligns[i]) for i, h in enumerate(headers)),
        sep,
    ]
    for row in rows:
        out.append(" | ".join(_pad(c, widths[i], aligns[i]) for i, c in enumerate(row)))
    return out


def _wrap(text: str, width: int, first: str = "", rest: str | None = None) -> list[str]:
    """한글 폭을 고려해 줄바꿈한다.

    first 는 첫 줄 앞에 붙는 문자열(글머리표 등), rest 는 이어지는 줄의 들여쓰기다.
    rest 를 생략하면 first 와 같은 폭의 공백을 쓴다 — 글머리표가 매 줄 반복되지 않는다.
    """
    if rest is None:
        rest = " " * _display_width(first)

    lines: list[str] = []
    current = ""
    prefix = first
    for word in text.split():
        candidate = f"{current} {word}".strip()
        if current and _display_width(candidate) + _display_width(prefix) > width:
            lines.append(prefix + current)
            prefix = rest
            current = word
        else:
            current = candidate
    if current:
        lines.append(prefix + current)
    return lines


# --------------------------------------------------------------------------
# 서열 표시
# --------------------------------------------------------------------------

def _cdr_position_map(report: SequenceReport) -> dict[int, str]:
    """1-based 위치 → CDR 이름."""
    mapping: dict[int, str] = {}
    if report.annotation:
        for name, cdr in report.annotation.cdrs.items():
            for pos in range(cdr.start, cdr.end + 1):
                mapping[pos] = name
    return mapping


def _sequence_block_text(report: SequenceReport) -> list[str]:
    """터미널용 서열 블록. CDR 구간을 아래 줄에 표시한다 (컬러를 쓰지 않으므로)."""
    seq = report.sequence
    cdr_at = _cdr_position_map(report)
    lines: list[str] = []

    for offset in range(0, len(seq), SEQUENCE_LINE_WIDTH):
        chunk = seq[offset:offset + SEQUENCE_LINE_WIDTH]
        start = offset + 1
        lines.append(f"  {start:>4}  {chunk}")
        marker = "".join("^" if (offset + i + 1) in cdr_at else " " for i in range(len(chunk)))
        if marker.strip():
            names: list[str] = []
            for i in range(len(chunk)):
                name = cdr_at.get(offset + i + 1)
                if name and name not in names:
                    names.append(name)
            lines.append(f"        {marker}  {', '.join(names)}")
    return lines


def _sequence_block_html(report: SequenceReport) -> str:
    """HTML 용 서열 블록. CDR 구간만 <b> 로 굵게 한다."""
    seq = report.sequence
    cdr_at = _cdr_position_map(report)
    rows: list[str] = []

    for offset in range(0, len(seq), SEQUENCE_LINE_WIDTH):
        chunk = seq[offset:offset + SEQUENCE_LINE_WIDTH]
        pieces: list[str] = []
        buffer = ""
        buffer_in_cdr = False

        def flush() -> None:
            if not buffer:
                return
            escaped = html.escape(buffer)
            pieces.append(f"<b>{escaped}</b>" if buffer_in_cdr else escaped)

        for i, residue in enumerate(chunk):
            in_cdr = (offset + i + 1) in cdr_at
            if in_cdr != buffer_in_cdr:
                flush()
                buffer = ""
                buffer_in_cdr = in_cdr
            buffer += residue
        flush()

        labels = []
        for i in range(len(chunk)):
            name = cdr_at.get(offset + i + 1)
            if name and name not in labels:
                labels.append(name)
        rows.append(
            "<tr>"
            f'<td class="num">{offset + 1}</td>'
            f'<td class="seq">{"".join(pieces)}</td>'
            f'<td class="note">{html.escape(", ".join(labels))}</td>'
            "</tr>"
        )
    return "<table class=\"seqblock\">" + "".join(rows) + "</table>"


# --------------------------------------------------------------------------
# 터미널 출력
# --------------------------------------------------------------------------

def as_text(reports: list[SequenceReport], *, width: int = 100) -> str:
    """터미널 표 출력."""
    out: list[str] = []
    out.append("=" * width)
    out.append(f"abcheck {__version__} — 항체 가변도메인 서열 셀프 점검")
    out.append("=" * width)
    out.append("")
    out.append("[이 도구의 한계]")
    for note in LIMITATION_NOTICE:
        out.extend(_wrap(note, width - 4, first="  - "))
    out.append("")
    out.extend(_wrap(NO_VERDICT_NOTICE, width - 4, first="  "))
    out.append("")

    out.extend(_summary_text(reports, width))
    out.append("")

    for report in reports:
        out.extend(_render_one_text(report, width))
        out.append("")

    out.append("=" * width)
    out.append(f"서열 {len(reports)}건 처리. 생성 {date.today().isoformat()}")
    return "\n".join(out)


def _summary_text(reports: list[SequenceReport], width: int) -> list[str]:
    """CDR별 liability 집계 표 (터미널용)."""
    summary = summarize(reports)
    out: list[str] = ["-" * width, "[CDR별 liability 집계]", "-" * width]
    out.extend(_wrap(summary["note"], width - 4, first="  "))
    out.append("")

    headers = ["서열", "체인", "주석"] + [SUMMARY_HEADERS[r] for r in SUMMARY_REGIONS] \
        + ["도메인 플래그"]
    rows: list[list[str]] = []

    for row in summary["rows"]:
        seq_id, chain, confidence, cells, flag_cell = _summary_row_cells(row)
        rows.append([seq_id, chain, confidence] + cells + [flag_cell])

    totals = summary["totals"]
    rows.append(
        ["합계", "", ""]
        + [_counts_cell(totals[r]) for r in SUMMARY_REGIONS]
        + [f"{summary['domain_flag_total']}건" if summary["domain_flag_total"] else EMPTY_CELL]
    )

    out.extend("  " + line for line in _table(headers, rows))

    if any(r["regions_unassigned"] for r in summary["rows"]):
        out.append("")
        out.extend(_wrap(
            "'주석 실패' 로 표시된 서열은 CDR 구간을 가를 수 없어 모티프를 framework 칸에 "
            "모았다. 실제로는 CDR 안일 수 있으니 과소평가하지 말 것.",
            width - 4, first="  ",
        ))
    return out


def _render_one_text(report: SequenceReport, width: int) -> list[str]:
    out: list[str] = ["-" * width, f"[{report.seq_id}]", "-" * width]

    if not report.ok:
        out.append("  입력을 쓸 수 없다:")
        for err in report.errors:
            out.extend(_wrap(err, width - 6, first="    * "))
        return out

    ann = report.annotation
    assert ann is not None

    chain_label = {"H": "중쇄 (VH)", "L": "경쇄 (VL)", None: "판별 실패"}[ann.chain]
    verdict = "신뢰 가능" if ann.confident else "신뢰 불가 (confident=False)"
    out.append(f"  체인: {chain_label}    길이: {len(report.sequence)} 잔기    주석 신뢰도: {verdict}")
    out.append("")

    if ann.reasons:
        out.append("  [주석을 신뢰할 수 없는 이유]")
        for reason in ann.reasons:
            out.extend(_wrap(reason, width - 6, first="    * "))
        out.append("")

    for note in ann.warnings + report.warnings + (report.liabilities.warnings if report.liabilities else []):
        out.extend(_wrap(f"참고: {note}", width - 4, first="  ", rest="        "))
    if ann.warnings or report.warnings or (report.liabilities and report.liabilities.warnings):
        out.append("")

    # --- CDR ---
    out.append("  [CDR 주석]  좌표는 1-based, 양끝 포함")
    if ann.cdrs:
        rows = [
            [cdr.name, f"{cdr.start}-{cdr.end}", str(cdr.length), cdr.seq]
            for cdr in sorted(ann.cdrs.values(), key=lambda c: c.start)
        ]
        out.extend(
            "    " + line
            for line in _table(["CDR", "위치", "길이", "서열"], rows,
                               ["left", "right", "right", "left"])
        )
    else:
        out.append("    앵커를 찾지 못해 CDR 을 끊지 못했다. 위치를 추측하지 않는다.")
    out.append("")

    if ann.anchors:
        anchor_labels = {
            "cys1": "보존 Cys(1)", "trp": "보존 Trp", "fr3_motif": "FR3 모티프",
            "cys2": "보존 Cys(2)", "terminal_motif": "말단 모티프",
        }
        found = ", ".join(
            f"{anchor_labels.get(k, k)} {v}번" for k, v in ann.anchors.items()
        )
        out.extend(_wrap(f"앵커: {found}", width - 4, first="  ", rest="        "))
        out.append("")

    # --- 서열 ---
    out.append("  [서열]  ^ 표시가 CDR 구간")
    out.extend(_sequence_block_text(report))
    out.append("")

    # --- liability ---
    lia = report.liabilities
    assert lia is not None
    out.append("  [liability 플래그]  CDR 안은 위험도를 한 단계 올리고, framework 는 한 단계 내린다")
    if lia.hits:
        counts = ", ".join(f"{lvl} {n}건" for lvl, n in lia.counts_by_risk().items())
        out.append(f"    등급 분포: {counts}")
        rows = [
            [h.risk, h.label, h.match, f"{h.start}-{h.end}" if h.end > h.start else str(h.start),
             h.region]
            for h in lia.hits
        ]
        out.extend(
            "    " + line
            for line in _table(["위험도", "항목", "매치", "위치", "구간"], rows,
                               ["left", "left", "left", "right", "left"])
        )
    else:
        out.append("    검출된 모티프 없음.")

    for flag in lia.domain_flags:
        out.append("")
        out.append(f"    [도메인 플래그] {flag['label']} — 위험도 {flag['risk']}")
        out.extend(_wrap(flag["rationale"], width - 8, first="      "))
    out.append("")

    if lia.hits:
        out.append("  [플래그 근거]")
        seen: set[str] = set()
        for hit in lia.hits:
            if hit.label in seen:
                continue
            seen.add(hit.label)
            out.append(f"    - {hit.label}")
            out.extend(_wrap(hit.rationale, width - 10, first="        "))
        out.append("")

    # --- 물성 ---
    out.append(f"  [물성]  pI 는 EMBOSS pKa 기준, GRAVY 는 Kyte-Doolittle, net charge 는 pH 7.4")
    rows = []
    for name, props in report.properties.items():
        d = props.as_dict()
        rows.append([
            name,
            str(d["length"]),
            "-" if d["pI"] is None else f"{d['pI']:.2f}",
            "-" if d["GRAVY"] is None else f"{d['GRAVY']:+.3f}",
            f"{d['net_charge']:+.2f}",
        ])
    out.extend(
        "    " + line
        for line in _table(["구간", "길이", "pI", "GRAVY", "net charge"], rows,
                           ["left", "right", "right", "right", "right"])
    )
    return out


# --------------------------------------------------------------------------
# JSON 출력
# --------------------------------------------------------------------------

def as_json(reports: list[SequenceReport], *, indent: int = 2) -> str:
    payload = {
        "tool": {"name": "abcheck", "version": __version__},
        "generated": date.today().isoformat(),
        "coordinates": "CDR 위치는 1-based, 양끝 포함",
        "limitations": list(LIMITATION_NOTICE),
        "no_verdict": NO_VERDICT_NOTICE,
        "count": len(reports),
        "summary": summarize(reports),
        "sequences": [r.as_dict() for r in reports],
    }
    return json.dumps(payload, ensure_ascii=False, indent=indent) + "\n"


# --------------------------------------------------------------------------
# HTML 출력
# --------------------------------------------------------------------------

_CSS = """
:root { color-scheme: light; }
body {
  font-family: Arial, Helvetica, sans-serif;
  font-size: 13px;
  line-height: 1.55;
  color: #000000;
  background: #ffffff;
  margin: 0;
  padding: 24px;
}
.wrap { max-width: 960px; margin: 0 auto; }
h1 { font-size: 20px; margin: 0 0 4px; }
h2 { font-size: 16px; margin: 28px 0 8px; padding-bottom: 4px; border-bottom: 1px solid #000000; }
h3 { font-size: 13px; margin: 18px 0 6px; }
p { margin: 6px 0; }
.meta { font-size: 12px; color: #555555; margin-bottom: 18px; }
.notice {
  border: 1px solid #000000;
  padding: 12px 14px;
  margin: 16px 0 24px;
  background: #f5f5f5;
}
.notice ul { margin: 6px 0 0; padding-left: 20px; }
.notice li { margin: 3px 0; }
.notice strong { display: block; margin-bottom: 4px; }
table { border-collapse: collapse; width: 100%; margin: 8px 0 4px; }
th {
  background: #000000;
  color: #ffffff;
  text-align: left;
  font-weight: bold;
  padding: 6px 8px;
  border: 1px solid #000000;
}
td { padding: 5px 8px; border: 1px solid #cccccc; vertical-align: top; }
tbody tr:nth-child(odd) { background: #ffffff; }
tbody tr:nth-child(even) { background: #f0f0f0; }
td.num, td.right { text-align: right; white-space: nowrap; }
.seqblock { font-family: "Courier New", Courier, monospace; font-size: 13px; }
.seqblock td { border: none; padding: 1px 8px 1px 0; }
.seqblock tr { background: #ffffff !important; }
.seqblock td.seq { letter-spacing: 1px; white-space: pre; }
.seqblock td.num { color: #555555; }
.seqblock td.note { font-family: Arial, Helvetica, sans-serif; color: #555555; }
.seq b { font-weight: bold; }
.domain-flag { border: 1px solid #000000; padding: 10px 12px; margin: 10px 0; background: #f5f5f5; }
.failed { border: 1px solid #000000; padding: 10px 12px; background: #f5f5f5; }
.rationale { font-size: 12px; color: #333333; margin: 2px 0 10px; }
.seq-section { margin-bottom: 36px; padding-bottom: 8px; border-bottom: 2px solid #000000; }
.summary { margin-bottom: 40px; padding-bottom: 12px; border-bottom: 2px solid #000000; }
.summary td:first-child { white-space: nowrap; }
.summary a { color: #000000; }
.summary tbody tr:last-child { border-top: 2px solid #000000; }
.kv { font-size: 13px; margin: 2px 0; }
.kv span { display: inline-block; min-width: 120px; color: #555555; }
ul.reasons { margin: 6px 0; padding-left: 20px; }
ul.reasons li { margin: 4px 0; }
footer { margin-top: 32px; font-size: 12px; color: #555555; }
"""


def _html_table(headers: list[str], rows: list[list[str]],
                right_cols: set[int] | None = None) -> str:
    right_cols = right_cols or set()
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = []
    for row in rows:
        cells = []
        for i, cell in enumerate(row):
            cls = ' class="right"' if i in right_cols else ""
            cells.append(f"<td{cls}>{cell}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def as_html(reports: list[SequenceReport], *, title: str = "항체 서열 셀프 점검 리포트") -> str:
    """단일 HTML 리포트. 컬러 강조 없음, CDR 은 굵게만."""
    parts: list[str] = [
        "<!doctype html>",
        '<html lang="ko"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{html.escape(title)}</title>",
        f"<style>{_CSS}</style>",
        "</head><body><div class=\"wrap\">",
        f"<h1>{html.escape(title)}</h1>",
        f'<p class="meta">abcheck {html.escape(__version__)} · '
        f"생성 {date.today().isoformat()} · 서열 {len(reports)}건 · "
        "CDR 좌표는 1-based, 양끝 포함</p>",
        '<div class="notice"><strong>이 도구의 한계</strong><ul>'
        + "".join(f"<li>{html.escape(n)}</li>" for n in LIMITATION_NOTICE)
        + f"</ul><p>{html.escape(NO_VERDICT_NOTICE)}</p></div>",
        _summary_html(reports),
    ]

    for report in reports:
        parts.append(_render_one_html(report))

    parts.append(
        "<footer>합성 데이터로 검증한 스크리닝 휴리스틱이다. "
        "결과를 근거로 후보를 탈락시키기 전에 ANARCI/IMGT 로 재확인할 것.</footer>"
    )
    parts.append("</div></body></html>")
    return "\n".join(parts)


def _summary_html(reports: list[SequenceReport]) -> str:
    """CDR별 liability 집계 표 (HTML용). 리포트 맨 앞에 놓는다."""
    summary = summarize(reports)

    headers = ["서열", "체인", "주석"] + [SUMMARY_HEADERS[r] for r in SUMMARY_REGIONS] \
        + ["도메인 플래그"]
    rows: list[list[str]] = []

    for row in summary["rows"]:
        seq_id, chain, confidence, cells, flag_cell = _summary_row_cells(row)
        link = f'<a href="#{html.escape(_anchor(seq_id))}">{html.escape(seq_id)}</a>'
        rows.append(
            [link, html.escape(chain), html.escape(confidence)]
            + [html.escape(c) for c in cells]
            + [html.escape(flag_cell)]
        )

    totals = summary["totals"]
    rows.append(
        ["<b>합계</b>", "", ""]
        + [f"<b>{html.escape(_counts_cell(totals[r]))}</b>" for r in SUMMARY_REGIONS]
        + [f"<b>{summary['domain_flag_total']}건</b>" if summary["domain_flag_total"]
           else EMPTY_CELL]
    )

    out = [
        '<div class="summary">',
        "<h2>CDR별 liability 집계</h2>",
        f'<p class="rationale">{html.escape(summary["note"])}</p>',
        _html_table(headers, rows),
    ]
    if any(r["regions_unassigned"] for r in summary["rows"]):
        out.append(
            '<p class="rationale">&lsquo;주석 실패&rsquo; 로 표시된 서열은 CDR 구간을 '
            "가를 수 없어 모티프를 framework 칸에 모았다. 실제로는 CDR 안일 수 있으니 "
            "과소평가하지 말 것.</p>"
        )
    out.append("</div>")
    return "\n".join(out)


def _anchor(seq_id: str) -> str:
    """서열 ID 를 HTML 앵커로 쓸 수 있게 다듬는다."""
    return "seq-" + "".join(ch if ch.isalnum() else "-" for ch in seq_id)


def _render_one_html(report: SequenceReport) -> str:
    out: list[str] = [
        f'<div class="seq-section" id="{html.escape(_anchor(report.seq_id))}">',
        f"<h2>{html.escape(report.seq_id)}</h2>",
    ]

    if not report.ok:
        out.append('<div class="failed"><strong>입력을 쓸 수 없다</strong><ul class="reasons">')
        out.extend(f"<li>{html.escape(e)}</li>" for e in report.errors)
        out.append("</ul></div></div>")
        return "\n".join(out)

    ann = report.annotation
    assert ann is not None
    lia = report.liabilities
    assert lia is not None

    chain_label = {"H": "중쇄 (VH)", "L": "경쇄 (VL)", None: "판별 실패"}[ann.chain]
    verdict = "신뢰 가능" if ann.confident else "신뢰 불가 (confident=False)"
    out.append(
        f'<p class="kv"><span>체인</span>{html.escape(chain_label)}</p>'
        f'<p class="kv"><span>길이</span>{len(report.sequence)} 잔기</p>'
        f'<p class="kv"><span>주석 신뢰도</span>{html.escape(verdict)}</p>'
    )

    if ann.reasons:
        out.append("<h3>주석을 신뢰할 수 없는 이유</h3><ul class=\"reasons\">")
        out.extend(f"<li>{html.escape(r)}</li>" for r in ann.reasons)
        out.append("</ul>")

    notes = list(ann.warnings) + list(report.warnings) + list(lia.warnings)
    if notes:
        out.append("<h3>참고</h3><ul class=\"reasons\">")
        out.extend(f"<li>{html.escape(n)}</li>" for n in notes)
        out.append("</ul>")

    # --- CDR ---
    out.append("<h3>CDR 주석</h3>")
    if ann.cdrs:
        rows = [
            [html.escape(c.name), f"{c.start}-{c.end}", str(c.length),
             f"<b>{html.escape(c.seq)}</b>"]
            for c in sorted(ann.cdrs.values(), key=lambda c: c.start)
        ]
        out.append(_html_table(["CDR", "위치", "길이", "서열"], rows, right_cols={1, 2}))
    else:
        out.append(
            '<p class="failed">앵커를 찾지 못해 CDR 을 끊지 못했다. 위치를 추측하지 않는다.</p>'
        )

    # --- 서열 ---
    out.append("<h3>서열 (CDR 구간은 굵게)</h3>")
    out.append(_sequence_block_html(report))

    # --- liability ---
    out.append("<h3>liability 플래그</h3>")
    out.append(
        '<p class="rationale">CDR 안은 기본 위험도에서 한 단계 올리고, '
        "framework 는 한 단계 내린다. 위험도는 색이 아니라 글자로 표시한다.</p>"
    )
    if lia.hits:
        rows = [
            [html.escape(h.risk), html.escape(h.label), f"<code>{html.escape(h.match)}</code>",
             f"{h.start}-{h.end}" if h.end > h.start else str(h.start),
             html.escape(h.region), html.escape(h.note)]
            for h in lia.hits
        ]
        out.append(_html_table(
            ["위험도", "항목", "매치", "위치", "구간", "등급 조정 근거"], rows, right_cols={3}
        ))
    else:
        out.append("<p>검출된 모티프 없음.</p>")

    for flag in lia.domain_flags:
        out.append(
            f'<div class="domain-flag"><strong>{html.escape(flag["label"])} '
            f"— 위험도 {html.escape(flag['risk'])}</strong>"
            f'<p class="rationale">{html.escape(flag["rationale"])}</p></div>'
        )

    if lia.hits:
        out.append("<h3>플래그 근거</h3>")
        seen: set[str] = set()
        for hit in lia.hits:
            if hit.label in seen:
                continue
            seen.add(hit.label)
            out.append(
                f"<p class=\"kv\"><b>{html.escape(hit.label)}</b></p>"
                f'<p class="rationale">{html.escape(hit.rationale)}</p>'
            )

    # --- 물성 ---
    out.append("<h3>물성</h3>")
    out.append(
        '<p class="rationale">pI 는 EMBOSS pKa 기준 bisection, GRAVY 는 Kyte-Doolittle 평균, '
        "net charge 는 pH 7.4. 모두 서열 기반 추정치다. CDR 조각의 pI 는 그 조각을 "
        "독립 펩타이드로 보고 계산한 값이므로 구간 간 상대 비교에만 쓸 것.</p>"
    )
    rows = []
    for name, props in report.properties.items():
        d = props.as_dict()
        rows.append([
            html.escape(name),
            str(d["length"]),
            "-" if d["pI"] is None else f"{d['pI']:.2f}",
            "-" if d["GRAVY"] is None else f"{d['GRAVY']:+.3f}",
            f"{d['net_charge']:+.2f}",
        ])
    out.append(_html_table(
        ["구간", "길이", "pI", "GRAVY", "net charge (pH 7.4)"], rows,
        right_cols={1, 2, 3, 4},
    ))

    out.append("</div>")
    return "\n".join(out)
