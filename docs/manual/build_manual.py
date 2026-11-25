#!/usr/bin/env python3
"""사용자 매뉴얼(HTML)을 현재 코드의 실제 출력으로 생성한다.

매뉴얼에 붙이는 터미널 출력·오류 메시지·JSON 은 손으로 옮겨 적지 않고 전부 여기서
abcheck 를 직접 호출해 만든다. 그래야 코드를 바꿨을 때 매뉴얼이 조용히 낡지 않는다
(tests/test_manual.py 가 커밋된 manual.html 이 지금 코드의 출력과 같은지 확인한다).

PDF 로 만드는 전체 절차는 docs/manual/build.sh 참고.

    python3 docs/manual/build_manual.py           # manual.html 생성
    python3 docs/manual/build_manual.py --stamp   # PDF 를 만든 뒤 원본 해시 기록
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from abcheck import LIMITATION_NOTICE, NO_VERDICT_NOTICE, __version__  # noqa: E402
from abcheck.annotate import TYPICAL_LENGTHS  # noqa: E402
from abcheck.liabilities import MOTIF_SPECS  # noqa: E402
from abcheck.render import as_text  # noqa: E402
from abcheck.report import build, build_many  # noqa: E402
from abcheck.triage import rule_descriptions  # noqa: E402
from abcheck.validate import validate  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT_HTML = HERE / "manual.html"
OUT_PDF = HERE / "abcheck_manual.pdf"
STAMP = HERE / "abcheck_manual.pdf.source-sha256"

#: 터미널 출력 폭. PDF 한 줄에 들어가도록 CLI 기본값(100)보다 조금 좁힌다.
TEXT_WIDTH = 92

CORPUS = json.loads((ROOT / "data/synthetic/sequences.json").read_text(encoding="utf-8"))
BY_ID = {r["id"]: r for r in CORPUS["sequences"]}


# --------------------------------------------------------------------------
# 실제 출력 만들기
# --------------------------------------------------------------------------

def seq_of(seq_id: str) -> str:
    return BY_ID[seq_id]["sequence"]


def section_of(text: str, seq_label: str) -> str:
    """as_text 전체 출력에서 한 서열의 섹션만 잘라낸다."""
    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.strip() == f"[{seq_label}]")
    end = len(lines)
    for i in range(start + 2, len(lines)):
        if lines[i].startswith("-" * 20) and i + 1 < len(lines) and lines[i + 1].startswith("["):
            end = i
            break
        if lines[i].startswith("=" * 20):
            end = i
            break
    return "\n".join(lines[start - 1:end]).rstrip()


def terminal_for(seq_id: str, label: str = "입력 서열") -> str:
    report = build(label, seq_of(seq_id))
    return section_of(as_text([report], width=TEXT_WIDTH), label)


def terminal_header() -> str:
    """터미널 출력 맨 앞(한계 고지)."""
    text = as_text([build("입력 서열", seq_of("SYN-VH-001"))], width=TEXT_WIDTH)
    lines = text.splitlines()
    end = next(i for i, ln in enumerate(lines) if ln.startswith("[CDR별 liability 집계]"))
    return "\n".join(lines[:end - 1]).rstrip()


def terminal_summary() -> str:
    """코퍼스 전체를 돌렸을 때의 집계 표 부분."""
    reports = build_many([(r["id"], r["sequence"]) for r in CORPUS["sequences"]])
    lines = as_text(reports, width=TEXT_WIDTH).splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("[CDR별 liability 집계]"))
    end = next(i for i in range(start, len(lines)) if lines[i].startswith("[SYN-"))
    return "\n".join(lines[start - 1:end - 1]).rstrip()


def json_excerpt(seq_id: str) -> str:
    """JSON 출력에서 핵심 필드만 발췌한다."""
    full = build(seq_id, seq_of(seq_id)).as_dict()
    hits = full["liabilities"]["hits"]
    excerpt = {
        "id": full["id"],
        "annotation": {
            "chain": full["annotation"]["chain"],
            "confident": full["annotation"]["confident"],
            "cdrs": full["annotation"]["cdrs"],
        },
        "liabilities": {
            "counts_by_risk": full["liabilities"]["counts_by_risk"],
            "hits": [
                {k: h[k] for k in ("label", "match", "start", "end", "region", "risk")}
                for h in hits if h["in_cdr"]
            ][:2],
        },
        "triage": {k: full["triage"][k] for k in ("status", "reasons")},
        "properties": {"도메인 전체": {
            k: v for k, v in full["properties"]["도메인 전체"].items() if k != "sequence"
        }},
    }
    return json.dumps(excerpt, ensure_ascii=False, indent=2)


def error_examples() -> list[tuple[str, str, str]]:
    """(상황, 넣은 것, 실제 오류 메시지) 목록."""
    vh = seq_of("SYN-VH-001")
    vl = seq_of("SYN-VL-001")
    cases = [
        ("서열이 너무 짧다", "ACDEFGHIKL", "ACDEFGHIKL"),
        ("FASTA 헤더가 섞였다", ">후보1\n" + vh, ">후보1 (줄바꿈) QVQL…"),
        ("아미노산이 아닌 문자", "S" * 50 + "123" + "S" * 50, "…SSS123SSS…"),
        ("모호한 아미노산 코드", vh[:40] + "X" + vh[41:], "…X…"),
        ("DNA 서열을 넣었다", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 3, "ATGGCCATTG…"),
        ("scFv 를 통째로 넣었다", vh + "GGGGSGGGGSGGGGS" + vl,
         f"VH + 링커 + VL ({len(vh) + 15 + len(vl)} 잔기)"),
    ]
    out = []
    for title, raw, shown in cases:
        result = validate(raw)
        assert not result.ok, title
        out.append((title, shown, result.errors[0]))
    return out


def failure_examples() -> list[tuple[str, str, str]]:
    """(서열 ID, 상황, 실제 신뢰 불가 사유)."""
    picks = [
        ("SYN-VH-901", "말단 모티프(WGxG)가 깨진 서열"),
        ("SYN-VH-902", "보존 Cys22 가 치환된 서열 (N-말단이 잘린 경우와 같은 증상)"),
        ("SYN-VH-903", "CDR-H1 이 비정상적으로 긴 서열"),
        ("SYN-VH-013", "FR3 모티프가 체세포 과변이로 깨진 서열"),
    ]
    out = []
    for seq_id, title in picks:
        report = build(seq_id, seq_of(seq_id))
        out.append((seq_id, title, report.annotation.reasons[0]))
    return out


# --------------------------------------------------------------------------
# HTML 조립
# --------------------------------------------------------------------------

def esc(text: str) -> str:
    return html.escape(text, quote=False)


def pre(text: str, cls: str = "term") -> str:
    return f'<pre class="{cls}">{esc(text)}</pre>'


def table(headers: list[str], rows: list[list[str]], raw: bool = False) -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{c if raw else esc(c)}</td>" for c in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def note(text: str, kind: str = "note") -> str:
    return f'<div class="{kind}">{text}</div>'


CSS = """
@page { size: A4; margin: 20mm 17mm 20mm 17mm; }
:root { color-scheme: light; }
* { box-sizing: border-box; }
body {
  font-family: "WenQuanYi Zen Hei", "Noto Sans KR", "Apple SD Gothic Neo",
               "Malgun Gothic", Arial, sans-serif;
  font-size: 10pt; line-height: 1.62; color: #111111; background: #ffffff; margin: 0;
  font-synthesis: weight style;
}
b, strong, th, h1, h2, h3 { font-weight: 700; }
h1 {
  font-size: 18pt; margin: 0 0 4mm; padding-bottom: 2mm;
  border-bottom: 2px solid #111111; page-break-before: always; break-before: page;
}
h1.first { page-break-before: avoid; break-before: avoid; }
h2 { font-size: 12.5pt; margin: 7mm 0 2mm; }
h3 { font-size: 10.5pt; margin: 5mm 0 1.5mm; }
p { margin: 1.5mm 0; }
ul, ol { margin: 1.5mm 0; padding-left: 6mm; }
li { margin: 0.8mm 0; }
code {
  font-family: "WenQuanYi Zen Hei Mono", "DejaVu Sans Mono", monospace;
  font-size: 8.8pt; background: #f0f0f0; padding: 0 1mm; border-radius: 1mm;
}
pre {
  font-family: "WenQuanYi Zen Hei Mono", "DejaVu Sans Mono", monospace;
  white-space: pre; overflow: hidden; margin: 2mm 0 3mm;
  page-break-inside: avoid; break-inside: avoid;
}
pre.allowbreak { page-break-inside: auto; break-inside: auto; }
pre.wrap { white-space: pre-wrap; word-break: break-all; font-size: 8pt; line-height: 1.5; }
pre.term.small { font-size: 6.2pt; line-height: 1.35; }
pre.term {
  font-size: 7.4pt; line-height: 1.38; background: #f6f6f6;
  border: 1px solid #cfcfcf; border-left: 3px solid #111111; padding: 2.5mm 3mm;
}
pre.cmd {
  font-size: 8.6pt; background: #111111; color: #ffffff; padding: 2mm 3mm; border-radius: 1mm;
}
pre.cmd::before { content: "$ "; color: #9a9a9a; }
pre.json { font-size: 7.4pt; line-height: 1.38; background: #f6f6f6; border: 1px solid #cfcfcf; padding: 2.5mm 3mm; }
table { border-collapse: collapse; width: 100%; margin: 2mm 0 3mm; font-size: 9pt;
        page-break-inside: avoid; break-inside: avoid; }
th { background: #111111; color: #ffffff; text-align: left; padding: 1.4mm 2mm; font-weight: bold; }
td { border: 1px solid #cfcfcf; padding: 1.3mm 2mm; vertical-align: top; }
tbody tr:nth-child(even) td { background: #f3f3f3; }
.note, .warn {
  border: 1px solid #111111; padding: 2.5mm 3.5mm; margin: 3mm 0; background: #f6f6f6;
  page-break-inside: avoid; break-inside: avoid;
}
.warn { border-width: 2px; }
.note b:first-child, .warn b:first-child {
  display: block; margin-bottom: 1.2mm; font-size: 10.5pt;
  padding-bottom: 1mm; border-bottom: 1px solid #bdbdbd;
}
figure { margin: 3mm 0; page-break-inside: avoid; break-inside: avoid; }
figure img {
  max-width: 100%; max-height: 200mm; width: auto; height: auto;
  border: 1px solid #cfcfcf; display: block; margin: 0 auto;
}
figcaption { font-size: 8.5pt; color: #555555; margin-top: 1.2mm; }
.cover { height: 247mm; display: flex; flex-direction: column; justify-content: center; }
.cover .title { font-size: 30pt; font-weight: bold; letter-spacing: -0.5pt; margin: 0; }
.cover .sub { font-size: 14pt; margin: 3mm 0 10mm; color: #333333; }
.cover .meta { font-size: 9.5pt; color: #555555; line-height: 1.8; }
.cover .box { border-top: 2px solid #111111; border-bottom: 2px solid #111111;
              padding: 5mm 0; margin: 8mm 0; font-size: 10pt; }
.toc ol { padding-left: 6mm; }
.toc li { margin: 1.5mm 0; font-size: 10.5pt; }
.step { display: flex; gap: 3mm; margin: 2.5mm 0; page-break-inside: avoid; break-inside: avoid; }
.step .n { flex: 0 0 7mm; height: 7mm; border-radius: 50%; background: #111111; color: #ffffff;
           text-align: center; line-height: 7mm; font-weight: bold; font-size: 9pt; }
.step .b { flex: 1; }
.keep { page-break-inside: avoid; break-inside: avoid; }
.kv { display: grid; grid-template-columns: 34mm 1fr; gap: 1mm 4mm; margin: 2mm 0; }
.kv div:nth-child(odd) { font-weight: bold; }
"""


def build_html() -> str:
    vh001 = seq_of("SYN-VH-001")
    vh004 = seq_of("SYN-VH-004")

    parts: list[str] = []
    add = parts.append

    # ---------------- 표지 ----------------
    add(f"""
<section class="cover">
  <p class="title">abcheck 사용 설명서</p>
  <p class="sub">항체 가변도메인(VH/VL) 서열 셀프 점검 도구</p>
  <div class="box">
    서열을 붙여넣으면 <b>CDR 위치</b>, <b>developability liability</b>, <b>물성</b>을 보여준다.<br>
    후보 선별 단계에서 "이 서열 일단 괜찮은지" 30초 만에 보는 용도.
  </div>
  <div class="meta">
    버전 {esc(__version__)} · 파이썬 3.10 이상 · 외부 라이브러리 없음 · 네트워크를 쓰지 않음<br>
    저장소 https://github.com/strawberryHunter/ab-seq-selfcheck<br>
    이 문서의 모든 출력 예시는 해당 버전의 도구를 실제로 실행해 얻은 것이다.<br>
    예시 서열은 전부 합성 데이터다 (실제 항체 서열이 아니다).
  </div>
</section>
""")

    # ---------------- 목차 ----------------
    add("""
<h1 class="first" style="page-break-before: always; break-before: page;">목차</h1>
<div class="toc"><ol>
  <li>이 도구는 무엇이고, 무엇을 하지 않는가</li>
  <li>준비 — 설치와 확인</li>
  <li>입력 — 무엇을 어떤 형식으로 넣는가</li>
  <li>실행 — 명령 네 가지</li>
  <li>결과 ① 터미널 출력 읽는 법</li>
  <li>결과 ② HTML 리포트</li>
  <li>결과 ③ JSON (엑셀·스크립트로 넘길 때)</li>
  <li>판정 기준 — 위험도 등급과 선별 규칙</li>
  <li>'신뢰 불가'가 나왔을 때</li>
  <li>오류 메시지와 대처</li>
  <li>실전 — 후보 여러 개를 선별하는 흐름</li>
  <li>검증 결과와 한계</li>
  <li>부록 — 탐색 규칙과 기준표</li>
</ol></div>
""")

    # ---------------- 1. 소개 ----------------
    add("""
<h1>1. 이 도구는 무엇이고, 무엇을 하지 않는가</h1>
<p>항체의 중쇄 가변도메인(VH) 또는 경쇄 가변도메인(VL) 아미노산 서열 하나를 받아 다음 세 가지를 보여준다.</p>
<div class="kv">
  <div>CDR 위치</div><div>항원에 붙는 고리 세 개(CDR1/2/3)가 서열 어디부터 어디까지인지. Kabat 경계를 근사한다.</div>
  <div>liability</div><div>개발 과정에서 문제를 일으킬 수 있는 서열 모티프(당화, deamidation, isomerization, fragmentation, 산화, 비결합 Cys). CDR 안이면 위험도를 한 단계 올린다.</div>
  <div>물성</div><div>pI, GRAVY(소수성), pH 7.4 에서의 net charge. 도메인 전체와 CDR별로.</div>
</div>
<p>여러 서열을 한 번에 넣으면 리포트 맨 앞에 <b>후보 비교표</b>가 나오고, 연구자가 정한 규칙에 따라 <b>'검토 필요'</b>를 표시한다.</p>
<h2>하지 않는 것</h2>
""")
    add("<ul>" + "".join(f"<li>{esc(n)}</li>" for n in LIMITATION_NOTICE) + "</ul>")
    add(note(f"<b>종합 판정을 내리지 않는다</b>{esc(NO_VERDICT_NOTICE)}"))
    add("""
<p>구조 예측, 도킹, 결합력 예측은 범위 밖이다. 외부 데이터베이스나 API 를 호출하지 않으므로
<b>사내 후보 서열을 넣어도 밖으로 나가지 않는다.</b></p>
""")

    # ---------------- 2. 준비 ----------------
    add("""
<h1>2. 준비 — 설치와 확인</h1>
<p>터미널을 연다. VS Code 를 쓴다면 <code>Ctrl</code> + <code>`</code> (맥은 <code>Cmd</code> + <code>`</code>).</p>
<div class="step"><div class="n">1</div><div class="b"><b>파이썬 확인</b> — 3.10 이상이면 된다. 없으면 python.org 에서 설치한다.
<pre class="cmd">python3 --version</pre></div></div>
<div class="step"><div class="n">2</div><div class="b"><b>코드 받기</b> — 설치할 라이브러리는 없다.
<pre class="cmd">git clone https://github.com/strawberryHunter/ab-seq-selfcheck.git</pre>
<pre class="cmd">cd ab-seq-selfcheck</pre>
git 이 없으면 GitHub 화면의 초록색 <b>Code → Download ZIP</b> 으로 받아 압축을 풀어도 된다.</div></div>
<div class="step"><div class="n">3</div><div class="b"><b>동작 확인</b> — 도움말이 나오면 준비 끝.
<pre class="cmd">python3 -m abcheck.cli --help</pre></div></div>
<div class="note"><b>모든 명령은 저장소 폴더 안에서 실행한다.</b>
<code>No module named abcheck</code> 가 나오면 <code>cd ab-seq-selfcheck</code> 를 먼저 했는지 확인한다.</div>
""")

    # ---------------- 3. 입력 ----------------
    add("""
<h1>3. 입력 — 무엇을 어떤 형식으로 넣는가</h1>
<h2>넣는 것은 하나: 가변도메인 하나의 아미노산 서열</h2>
<p><b>VH 하나 또는 VL 하나</b>를 넣는다. 중쇄인지 경쇄인지는 도구가 말단 모티프(<code>WGxG</code> / <code>FGxG</code>)로 스스로 판별하므로 알려주지 않아도 된다.</p>
""")
    add(table(
        ["형식", "언제 쓰나", "명령"],
        [
            ["서열 문자열", "서열 하나를 바로 볼 때", "python3 -m abcheck.cli \"서열\""],
            ["FASTA 파일 (.fasta .fa .faa)", "후보 여러 개", "python3 -m abcheck.cli 후보.fasta"],
            ["JSON 파일 (.json)", "다른 스크립트가 만든 목록", "python3 -m abcheck.cli 목록.json"],
        ],
    ))
    add("<h3>FASTA 파일 만드는 법</h3><p>메모장이나 VS Code 로 아래처럼 쓰고 <code>후보.fasta</code> 로 저장한다. <code>&gt;</code> 로 시작하는 줄이 이름, 그 아래가 서열이다. 서열은 여러 줄로 나눠 써도 된다.</p>")
    add(pre(f">후보-01\n{vh001[:60]}\n{vh001[60:]}\n>후보-02\n{vh004}", "term"))
    add("<h3>JSON 파일 형식</h3>")
    add(pre('{"sequences": [\n  {"id": "후보-01", "sequence": "QVQL..."},\n  {"id": "후보-02", "sequence": "EVQL..."}\n]}', "json"))
    add("<h2>입력 규칙</h2>")
    add(table(
        ["항목", "규칙"],
        [
            ["허용 문자", "표준 아미노산 20종 (ACDEFGHIKLMNPQRSTVWY). 대소문자 무관"],
            ["자동 정리", "공백·줄바꿈은 알아서 지운다. 소문자는 대문자로 바꾼다"],
            ["길이", "70 ~ 200 잔기 (VH 는 보통 115~130, VL 은 105~115)"],
            ["거부", "숫자·구두점, 모호 코드(X B Z J U O), FASTA 헤더가 섞인 문자열, DNA 서열"],
        ],
    ))
    add("<h2>긴 서열에서 가변도메인만 잘라내기</h2><p>전장 중쇄, Fab, scFv 는 그대로 넣으면 '너무 길다'로 거부된다. 가변도메인만 잘라낸다.</p>")
    add(table(
        ["", "시작 (대개)", "끝 — 여기까지 자른다", "그 뒤는"],
        [
            ["VH", "EVQL… / QVQL…", "WGxG 뒤 …VTVSS", "ASTKGPSV… (CH1)"],
            ["VL kappa", "DIQM… / EIVL…", "FGxG 뒤 …KVEIK", "RTVAAPSV… (CL)"],
            ["VL lambda", "QSVL… / SYEL…", "FGxG 뒤 …TVL", "GQPKAAPS… (CL)"],
        ],
    ))
    add(note("<b>신호서열(signal peptide)을 꼭 떼어낸다</b>앞에 <code>MGWSWIFLFLLSG…</code> 같은 19잔기 안팎의 신호서열이 붙어 있으면 보존 Cys 가 탐색 범위를 벗어나 주석이 실패한다. 발현 태그(<code>MA…</code>) 정도의 짧은 꼬리는 괜찮다.", "warn"))

    # ---------------- 4. 실행 ----------------
    add("""
<h1>4. 실행 — 명령 네 가지</h1>
<h2>① 서열 하나를 터미널에서 보기</h2>
<p>따옴표 안에 서열을 붙여넣는다.</p>
""")
    add(pre(f'python3 -m abcheck.cli "{vh004}"', "cmd"))
    add("""
<h2>② 후보 여러 개를 HTML 리포트로</h2>
<p><code>sample_output/report.html</code> 이 생긴다. 더블클릭하면 브라우저에서 열린다.</p>
<pre class="cmd">python3 -m abcheck.cli 후보.fasta --html</pre>
<p>다른 폴더·이름으로 저장하려면:</p>
<pre class="cmd">python3 -m abcheck.cli 후보.fasta --html --out-dir 결과 --out-name 1차선별.html</pre>
<h2>③ JSON 으로 받기 (엑셀·스크립트용)</h2>
<pre class="cmd">python3 -m abcheck.cli 후보.fasta --json &gt; 결과.json</pre>
<h2>④ net charge 를 다른 pH 에서</h2>
<p>제형 pH 를 검토할 때 쓴다. 기본은 7.4.</p>
<pre class="cmd">python3 -m abcheck.cli 후보.fasta --ph 6.0</pre>
""")
    add(table(
        ["옵션", "뜻", "기본값"],
        [
            ["--html", "HTML 리포트 파일을 만든다", "꺼짐 (터미널 표)"],
            ["--json", "JSON 을 화면으로 출력한다 (> 로 파일에 저장)", "꺼짐"],
            ["--out-dir", "HTML 을 저장할 폴더", "sample_output"],
            ["--out-name", "HTML 파일 이름", "report.html"],
            ["--ph", "net charge 를 계산할 pH", "7.4"],
        ],
    ))
    add(table(
        ["종료 코드", "뜻"],
        [
            ["0", "정상 처리 (주석을 신뢰할 수 없는 서열이 있어도 0 — 그건 실패가 아니라 결과다)"],
            ["1", "입력 검증에서 걸러진 서열이 있다 (나머지는 처리된다)"],
            ["2", "입력 파일을 읽을 수 없다 (파일이 없거나 형식이 틀림)"],
        ],
    ))

    # ---------------- 5. 터미널 출력 ----------------
    add("""
<h1>5. 결과 ① 터미널 출력 읽는 법</h1>
<p>4장의 ① 명령을 실행했을 때 실제로 나오는 출력이다. 덩어리별로 읽는 법을 적었다.</p>
<h2>맨 앞 — 한계 고지 (매번 나온다)</h2>
<p>결과만 떼어 남에게 보여줄 때도 이 고지가 따라가도록 일부러 매번 붙인다.</p>
""")
    add(pre(terminal_header()))
    add("<h2>서열별 결과 — 실제 출력 전체</h2><p>예시는 CDR-H1 안에 N-결합 당화 서열이 있는 합성 서열이다.</p>")
    add(pre(terminal_for("SYN-VH-004"), "term allowbreak"))
    add("<h2>덩어리별로 읽기</h2>")
    add(table(
        ["덩어리", "읽는 법"],
        [
            ["체인 · 주석 신뢰도",
             "<b>신뢰 가능</b>이면 CDR 위치를 그대로 써도 된다. <b>신뢰 불가</b>면 바로 아래 '이유'를 읽는다 (9장)."],
            ["선별 분류",
             "8장의 규칙에 걸렸는지. <b>검토 필요</b>면 → 다음 줄에 어디에 무엇이 걸렸는지가 나온다."],
            ["[CDR 주석]",
             "좌표는 <b>1-based, 양끝 포함</b>. <code>31-36</code> 은 31번 잔기부터 36번 잔기까지."],
            ["앵커",
             "CDR 경계를 정하는 데 쓴 보존 잔기의 위치. 이상한 결과가 나왔을 때 원인을 찾는 단서."],
            ["[서열]", "<code>^</code> 표시가 CDR 구간이다. 오른쪽에 어느 CDR 인지 적힌다."],
            ["[liability 플래그]",
             "위험도는 <b>높음/중간/낮음</b> 글자로만 표시한다. 같은 모티프라도 CDR 안이면 한 단계 높다."],
            ["[플래그 근거]", "각 항목이 왜 문제인지 한 문장씩. 회의 자료에 그대로 옮겨도 된다."],
            ["[물성]",
             "pI 는 등전점, GRAVY 는 양수면 소수성, net charge 는 pH 7.4 의 알짜 전하. CDR 별 pI 는 구간 간 상대 비교에만 쓴다."],
        ],
        raw=True,
    ))
    add(note("<b>보존 Trp 의 '산화 (Trp)' 는 대개 무시해도 된다</b>Trp36 과 <code>WGxG</code> 의 Trp 은 거의 모든 항체에 있는 보존 잔기라 매번 framework 산화(낮음)로 잡힌다. 단백질 안쪽에 묻혀 있어 실제 위험은 낮지만, 검출 자체를 숨기지는 않는다."))
    add("<h2>도메인 플래그 — 비결합 Cys</h2><p>Cys 개수가 홀수면 짝을 못 이룬 Cys 가 있다는 뜻이라 서열 전체에 대한 플래그로 따로 나온다.</p>")
    unpaired = terminal_for("SYN-VH-005")
    block = [ln for ln in unpaired.splitlines()]
    start = next(i for i, ln in enumerate(block) if "[도메인 플래그]" in ln)
    add(pre("\n".join(block[start:start + 6]).rstrip()))

    # ---------------- 6. HTML 리포트 ----------------
    add("""
<h1>6. 결과 ② HTML 리포트</h1>
<p><code>--html</code> 로 만든 리포트다. 아래 화면은 저장소의 합성 예시 24개를 넣어 만든 <code>sample_output/report.html</code> 이다.</p>
<div class="keep"><h2>맨 앞 — 후보 비교표</h2>
<p>중쇄·경쇄·체인 판별 실패가 따로 나온다. <b>분류</b> 열로 먼저 볼 후보를 고르고, 서열 이름을 누르면 상세로 이동한다.</p>
<figure><img src="img/report_summary.png" alt="HTML 리포트 비교표">
<figcaption>그림 1. 리포트 맨 앞의 비교표. 위험도는 색이 아니라 글자로 표시한다.</figcaption></figure></div>
<div class="keep"><h2>서열별 상세</h2>
<figure><img src="img/report_detail.png" alt="HTML 리포트 서열 상세">
<figcaption>그림 2. 서열 상세. CDR 구간은 굵게만 표시한다.</figcaption></figure></div>
<div class="keep"><h2>주석이 실패한 서열</h2>
<p>위치를 지어내지 않고, 왜 못 했는지를 적는다.</p>
<figure><img src="img/report_failed.png" alt="HTML 리포트 실패 사례">
<figcaption>그림 3. 말단 모티프가 깨진 서열. CDR 을 끊지 않고 이유를 적는다.</figcaption></figure></div>
<div class="note"><b>왜 색을 쓰지 않나</b>위험도를 빨강·노랑으로 칠하면 "빨간 게 몇 개인지"로 후보를 고르게 되는데, 그건 이 도구가 하지 않기로 한 종합 판정이다. 대신 무엇이 왜 걸렸는지를 읽게 한다.</div>
""")

    # ---------------- 7. JSON ----------------
    add("""
<h1>7. 결과 ③ JSON (엑셀·스크립트로 넘길 때)</h1>
<p>터미널 표에 안 나오는 것까지 전부 들어 있다. 아래는 4장 ① 서열의 실제 JSON 출력에서 핵심 필드만 발췌한 것이다.</p>
""")
    add(pre(json_excerpt("SYN-VH-004"), "json"))
    add(table(
        ["키", "내용"],
        [
            ["limitations / no_verdict", "한계 고지. 파일 하나만 떼어 넘겨도 따라간다"],
            ["summary", "후보 비교표와 같은 집계 (중쇄/경쇄 묶음, 합계)"],
            ["sequences[].annotation", "체인, 신뢰도, 사유, 앵커 위치, CDR 좌표"],
            ["sequences[].liabilities", "모티프별 위치·구간·위험도·근거, 도메인 플래그"],
            ["sequences[].properties", "도메인 전체와 CDR 별 pI · GRAVY · net charge"],
            ["sequences[].triage", "선별 분류와 걸린 사유"],
        ],
    ))

    # ---------------- 8. 기준 ----------------
    add('<h1>8. 판정 기준 — 위험도 등급과 선별 규칙</h1><h2>liability 위험도</h2><p>모티프마다 기본 등급이 있고, <b>CDR 안이면 한 단계 올리고 framework 면 한 단계 내린다</b> (높음·낮음에서는 멈춘다).</p>')
    rows = []
    for spec in MOTIF_SPECS:
        cdr = {"낮음": "중간", "중간": "높음", "높음": "높음"}[spec.base_risk]
        fw = {"낮음": "낮음", "중간": "낮음", "높음": "중간"}[spec.base_risk]
        rows.append([spec.label, spec.pattern, spec.base_risk, cdr, fw])
    rows.append(["비결합 Cys", "도메인 내 Cys 홀수", "-", "높음 (Cys 가 CDR 안)", "중간"])
    add(table(["항목", "모티프", "기본", "CDR 안", "framework"], rows))
    add("<h2>선별 규칙 ('검토 필요')</h2><p>기본 규칙은 아래 두 개다. 하나라도 걸리면 <b>검토 필요</b>.</p>")
    add("<ul>" + "".join(f"<li>{esc(d)}</li>" for d in rule_descriptions()) + "</ul>")
    add(table(
        ["분류", "뜻"],
        [
            ["검토 필요", "규칙에 걸렸다. 항체가 나쁘다는 뜻이 아니다 — 먼저 들여다볼 후보라는 뜻"],
            ["기준 해당 없음", "규칙에 걸리지 않았다. 안전하다는 뜻이 아니다"],
            ["판정 불가", "CDR 주석이 실패해 규칙을 적용할 수 없다. 사람이 직접 봐야 한다"],
        ],
    ))
    add(note("<b>규칙은 바꿀 수 있다</b><code>abcheck/triage.py</code> 의 <code>DEFAULT_RULES</code> 에 글로 적혀 있다. 기준이 코드에 남으므로 나중에 \"왜 이 후보를 먼저 봤지\"를 추적할 수 있다. 약한 등급(NS/NT/DS 등)이 CDR 상향으로 '높음'이 되어도 선별에는 걸리지 않게 한 것은, 실제 레퍼토리에서 그런 모티프가 너무 흔해 모든 후보가 걸리기 때문이다."))

    # ---------------- 9. 신뢰 불가 ----------------
    add("<h1>9. '신뢰 불가'가 나왔을 때</h1><p>도구는 앵커를 못 찾으면 <b>위치를 지어내지 않고</b> 이유를 적는다. 아래는 실제 출력되는 사유들이다.</p>")
    for seq_id, title, reason in failure_examples():
        add(f"<h3>{esc(title)}</h3>")
        add(pre(reason, "term wrap"))
    add(table(
        ["사유에 이런 말이 있으면", "의심할 것", "할 일"],
        [
            ["보존 Cys22 / Cys23 를 찾지 못했다", "신호서열이 붙어 있다 / N-말단이 잘렸다", "서열 앞부분 확인, 신호서열 제거"],
            ["말단 모티프를 찾지 못했다", "C-말단이 잘렸다 / 가변도메인이 아니다", "서열 끝이 …WGxG… 또는 …FGxG… 인지 확인"],
            ["CDR-? 길이가 통상 범위보다", "앵커는 찾았지만 경계가 밀렸을 수 있다", "ANARCI/IMGT 로 경계 재확인"],
            ["… 역산해 …", "체세포 과변이로 FR3 모티프가 깨졌다", "위치는 대개 맞다. ANARCI 로 확인 후 사용"],
            ["… 통상값과 다르다 … 고쳐 잡았다", "모티프가 엉뚱한 자리에 걸렸다", "고친 위치를 ANARCI 로 확인"],
        ],
    ))
    add(note("<b>결과를 근거로 후보를 탈락시키기 전에는 반드시 ANARCI/IMGT 로 재확인한다</b>이 도구는 빠른 1차 스크리닝용이다. 신뢰 가능으로 나왔더라도, 의사결정에 쓰기 전에는 정식 numbering 도구로 한 번 더 본다.", "warn"))

    # ---------------- 10. 오류 ----------------
    add("<h1>10. 오류 메시지와 대처</h1><p>입력이 잘못되면 처리하지 않고 한글로 이유를 알려준다. 아래는 실제 메시지다.</p>")
    for title, shown, message in error_examples():
        add(f"<h3>{esc(title)}</h3><p>넣은 것: <code>{esc(shown)}</code></p>")
        add(pre(message, "term wrap"))
    add(table(
        ["메시지", "원인", "고치는 법"],
        [
            ["파일을 찾을 수 없다", "경로가 틀렸다", "파일 이름·확장자 확인, 저장소 폴더 안에서 실행"],
            ["'sequence' 키가 없다", "JSON 형식이 다르다", "3장의 JSON 형식으로 맞춘다"],
            ["No module named abcheck", "저장소 폴더 밖에서 실행했다", "cd ab-seq-selfcheck"],
        ],
    ))

    # ---------------- 11. 실전 ----------------
    add("<h1>11. 실전 — 후보 여러 개를 선별하는 흐름</h1><p>파지 디스플레이나 B 세포 클로닝에서 후보가 수십 개 나왔을 때의 순서다.</p>")
    add("""
<div class="step"><div class="n">1</div><div class="b"><b>가변도메인만 FASTA 로 모은다</b> — 3장의 자르는 기준대로. 신호서열을 뗐는지 확인.</div></div>
<div class="step"><div class="n">2</div><div class="b"><b>한 번에 돌려 리포트를 만든다</b>
<pre class="cmd">python3 -m abcheck.cli 후보.fasta --html --out-name 1차선별.html</pre></div></div>
<div class="step"><div class="n">3</div><div class="b"><b>터미널 마지막 줄을 본다</b> — 입력에서 걸러진 서열과 '신뢰 불가' 서열 목록이 나온다. 이것부터 처리한다 (9·10장).</div></div>
<div class="step"><div class="n">4</div><div class="b"><b>리포트 맨 앞 비교표에서 '검토 필요'를 본다</b> — 어느 CDR 에 무엇이 걸렸는지 한 줄로 보인다.</div></div>
<div class="step"><div class="n">5</div><div class="b"><b>서열 이름을 눌러 상세와 근거를 읽는다</b> — 그 모티프가 결합에 필수적인 자리인지는 도구가 모른다. 그 판단은 사람이 한다.</div></div>
<div class="step"><div class="n">6</div><div class="b"><b>진행할 후보를 ANARCI/IMGT 로 재확인한다</b> — 탈락·진행 결정 전에.</div></div>
<h2>터미널에서 본 비교표 (합성 예시 24개)</h2>
""")
    add(pre(terminal_summary(), "term allowbreak small"))

    # ---------------- 12. 검증 ----------------
    add("""
<h1>12. 검증 결과와 한계</h1>
<p>공개 항체 레퍼토리 데이터베이스(OAS)의 실제 서열 <b>7,000개</b>로 측정했다. 측정 절차는 저장소의
<code>scripts/validate_on_oas.py</code>, <code>scripts/compare_with_anarci.py</code> 로 누구나 다시 돌려볼 수 있다.</p>
""")
    add(table(
        ["항목", "결과"],
        [
            ["체인(H/L) 판별", "오판 0건 / 6,984"],
            ["CDR 주석 성공", "98.3%"],
            ["CDR 경계 (ANARCI 정답 대비, 16,427건)", "불일치 0.35%"],
            ["사람 naive B 세포 데이터", "신뢰 가능 99.8%"],
            ["과변이가 심한 memory/plasma B 세포", "신뢰 가능 93~95%"],
        ],
    ))
    add("""
<h2>이런 서열에서는 약하다</h2>
<ul>
  <li><b>체세포 과변이가 심한 서열</b> — 보존 모티프가 깨져 '신뢰 불가'가 늘어난다. 위치는 대개 복원하지만 보수적으로 표시한다.</li>
  <li><b>잘린 리드</b> — 앵커가 없으면 주석하지 않는다.</li>
  <li><b>보존 Trp 이 치환된 서열</b> — 드물지만(약 0.03%) 앵커를 못 찾는다.</li>
  <li><b>나노바디(VHH)</b> — 주석은 나오지만 정답과 비교해 검증하지 못했다. FR2 문맥 경고가 항상 뜬다 (VHH 특유의 V37F 치환 때문).</li>
  <li><b>구조</b> — liability 는 서열 모티프만 본다. 묻혀 있는 자리도 플래그된다.</li>
</ul>
<p>자세한 측정 방법, 실제 데이터로 찾아 고친 버그, 고치지 않기로 한 것은 저장소의 <code>VALIDATION.md</code> 에 있다.</p>
""")

    # ---------------- 13. 부록 ----------------
    add("<h1>13. 부록 — 탐색 규칙과 기준표</h1><h2>CDR 을 찾는 순서</h2>")
    add(table(
        ["체인", "순서"],
        [
            ["VH", "Cys22 → CDR-H1 (Cys22+9 부터) → Trp36 (뒤에 [VIA][RKQ]Q) → CDR-H2 (Trp36+14 부터) → "
                   "FR3 모티프 [RKQ].[TSAV][IMLFV][STAN] → Cys92 (앞에 Y/F) → CDR-H3 (Cys92+3 부터) → WGxG"],
            ["VL", "Cys23 → CDR-L1 (Cys23+1 부터) → Trp35 (뒤에 [YF][QL]Q) → CDR-L2 (Trp35+15 부터, 7잔기) → "
                   "Cys88 (앞에 Y/F) → CDR-L3 (Cys88+1 부터) → FGxG"],
        ],
    ))
    add("<p>FR3 모티프를 못 찾으면 보존 Cys 에서 29잔기 역산한다 (실제 VH 3,443건이 전부 29 였다). 이 경우 '신뢰 불가'로 표시한다.</p>")
    add("<h2>CDR 통상 길이</h2><p>이 범위를 벗어나면 '신뢰 불가'가 된다. 이 범위는 \"항체가 이상한가\"가 아니라 \"주석을 믿을 수 있는가\"의 기준이다.</p>")
    add(table(
        ["CDR", "통상 길이 (잔기)"],
        [[name, str(lo) if lo == hi else f"{lo}-{hi}"] for name, (lo, hi) in TYPICAL_LENGTHS.items()],
    ))
    add("<h2>물성 계산</h2>")
    add(table(
        ["값", "방법"],
        [
            ["pI", "EMBOSS pKa 값 (N말단 8.6, C말단 3.6, D 3.9, E 4.1, C 8.5, Y 10.1, H 6.5, K 10.8, R 12.5), bisection"],
            ["GRAVY", "Kyte-Doolittle hydropathy 평균. 양수 = 소수성"],
            ["net charge", "Henderson-Hasselbalch. 기본 pH 7.4"],
        ],
    ))
    add('<div class="keep"><h2>더 읽을 것 (저장소)</h2>')
    add(table(
        ["파일", "내용"],
        [
            ["README.md / README.en.md", "개요와 검증 요약 (한글 / 영문)"],
            ["TUTORIAL.md", "이 설명서의 짧은 판"],
            ["VALIDATION.md", "측정 방법, 결과, 알려진 한계"],
            ["DECISIONS.md", "설계 판단과 그 이유"],
        ],
    ))
    add("</div>")

    body = "\n".join(parts)
    return f"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8">
<title>abcheck 사용 설명서</title>
<style>{CSS}</style>
</head><body>
{body}
</body></html>
"""


def write_stamp() -> None:
    """PDF 를 만든 원본 HTML 의 해시를 남긴다 — PDF 가 낡았는지 테스트가 확인한다."""
    digest = hashlib.sha256(OUT_HTML.read_bytes()).hexdigest()
    STAMP.write_text(digest + "\n", encoding="utf-8")
    print(f"원본 해시 기록: {STAMP.relative_to(ROOT)}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="사용 설명서 HTML 을 생성한다.")
    parser.add_argument("--stamp", action="store_true",
                        help="PDF 를 만든 뒤 원본 HTML 해시를 기록한다")
    args = parser.parse_args(argv)
    if args.stamp:
        write_stamp()
        return 0
    OUT_HTML.write_text(build_html(), encoding="utf-8")
    print(f"생성: {OUT_HTML.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
