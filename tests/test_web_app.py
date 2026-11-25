"""The browser app (web/abcheck.html) must be built from the current template and engine.

The engine itself is checked against Python in tests/test_js_parity.py.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "web" / "abcheck.html"


def test_committed_app_is_current():
    result = subprocess.run(
        [sys.executable, "web/build.py", "--check"],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr + " — run: python3 web/build.py"


def test_app_loads_nothing_but_fonts_from_outside():
    """Sequences must never leave the browser: no fetch/XHR, no external scripts."""
    html = APP.read_text(encoding="utf-8")
    hosts = set(re.findall(r'(?:src|href)="(https?://[^/"]+)', html))
    assert hosts <= {"https://fonts.googleapis.com", "https://fonts.gstatic.com", "https://github.com"}, hosts
    for api in ["fetch(", "XMLHttpRequest", "WebSocket", "sendBeacon"]:
        assert api not in html, api


def test_app_carries_the_limitations_and_no_verdict_notice():
    html = APP.read_text(encoding="utf-8")
    assert "<title>AbCheck 항체 서열 점검</title>" in html
    assert "numbering 엔진이 아니다" in html
    assert "종합 판정을 내리지 않는다" in html


def test_example_uses_only_synthetic_sequences():
    html = APP.read_text(encoding="utf-8")
    example = re.search(r"const EXAMPLE = (\".*?\");", html).group(1)
    assert ">SYN-" in example and "OAS" not in example


# --------------------------------------------------------------------------
# web/app-lib.js -- the data helpers behind the app screens, run under Node
# --------------------------------------------------------------------------

import json  # noqa: E402
import shutil  # noqa: E402

import pytest  # noqa: E402

from abcheck.liabilities import MOTIF_SPECS  # noqa: E402

NODE = shutil.which("node")
needs_node = pytest.mark.skipif(NODE is None, reason="node is not installed; app-lib.js runs in the browser")


def _node(script: str):
    """Run script with ABCHECK (engine) and ABUI (app-lib) loaded; it prints one JSON value."""
    prelude = (
        "const ABCHECK = require(process.argv[1]);"
        "const ABUI = require(process.argv[2]);"
        "const fs = require('fs');"
        "const corpus = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));"
        "const entries = corpus.sequences.map(r => [r.id, r.sequence]);"
    )
    result = subprocess.run(
        [NODE, "-e", prelude + script, str(ROOT / "web" / "abcheck-engine.js"),
         str(ROOT / "web" / "app-lib.js"), str(ROOT / "data" / "synthetic" / "sequences.json")],
        capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


@needs_node
def test_build_rows_joins_reports_and_summary_rows():
    out = _node("""
      const p = ABCHECK.analyze(entries.concat([["bad", "ACDEFGHIK"]]), 7.4);
      const rows = ABUI.buildRows(p);
      console.log(JSON.stringify({
        n: rows.length, ids: rows.map(r => r.id), review_count: p.summary.review_count,
        tally: ABUI.tally(rows),
        vh003: rows.find(r => r.id === "SYN-VH-003"),
        bad: rows.find(r => r.id === "bad"),
      }));""")
    assert out["n"] == 25
    assert out["ids"][-1] == "bad"
    assert out["tally"]["review"] == out["review_count"]
    assert out["tally"]["errors"] == 1
    assert out["tally"]["all"] == 25
    vh = out["vh003"]
    assert vh["chain"] == "H" and vh["status"] == "검토 필요" and vh["cats"] == ["review"]
    assert [c["name"] for c in vh["cdrs"]] == ["CDR-H1", "CDR-H2", "CDR-H3"]
    assert vh["cdr3"]["seq"] == "YNKQQQFYLK"
    assert vh["cdrCounts"] == {"높음": 1, "중간": 0, "낮음": 0}
    assert vh["pI"] == 9.65 and vh["gravy"] == -0.258
    bad = out["bad"]
    assert bad["inputOk"] is False and bad["cats"] == ["errors"] and bad["pI"] is None
    assert "너무 짧다" in bad["errors"][0]


@needs_node
def test_unconfident_sequences_land_in_recheck():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries, 7.4));
      console.log(JSON.stringify(rows.filter(r => r.cats.includes("recheck")).map(r => r.id).sort()));""")
    assert out == ["SYN-VH-013", "SYN-VH-014", "SYN-VH-901", "SYN-VH-902", "SYN-VH-903"]


@needs_node
def test_filter_and_sort():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries, 7.4));
      const ids = (rs) => rs.map(r => r.id);
      const byPi = ABUI.sortRows(rows, "pi-desc");
      console.log(JSON.stringify({
        vl: ids(ABUI.filterRows(rows, {view: "all", chain: "L"})),
        name: ids(ABUI.filterRows(rows, {query: "vh-00"})).length,
        residues: ids(ABUI.filterRows(rows, {query: "ngg"})),
        twoLetters: ids(ABUI.filterRows(rows, {query: "NG"})).length,
        input: ids(ABUI.sortRows(rows, "input")).join() === ids(rows).join(),
        byName: ids(ABUI.sortRows(rows, "name")).slice(0, 3),
        piDesc: byPi.filter(r => r.pI !== null).every((r, i, a) => i === 0 || a[i - 1].pI >= r.pI),
        piLast: byPi[byPi.length - 1].pI,
      }));""")
    assert out["vl"] == ["SYN-VL-001", "SYN-VL-002", "SYN-VL-003", "SYN-VL-004", "SYN-VL-005", "SYN-VL-007", "SYN-VL-006"]
    assert out["name"] == 9                     # SYN-VH-001 … SYN-VH-009
    assert "SYN-VH-003" in out["residues"]     # LSTAARNGGQQQ…
    assert out["twoLetters"] == 0              # 2 letters match names only, never residues
    assert out["input"] is True
    assert out["byName"] == ["SYN-VH-001", "SYN-VH-002", "SYN-VH-003"]
    assert out["piDesc"] is True
    assert out["piLast"] is not None            # every corpus sequence passes validation


@needs_node
def test_shared_cdrs_finds_identical_cdrs():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze([["a", entries[0][1]], ["b", entries[0][1]], ["c", entries[2][1]]], 7.4));
      console.log(JSON.stringify(ABUI.sharedCdrs(rows)));""")
    assert [g["name"] for g in out] == ["CDR-H1", "CDR-H2", "CDR-H3"]
    assert all(g["ids"] == ["a", "b"] for g in out)


@needs_node
def test_csv_and_tsv_exports():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze([["a, \\"x\\"", entries[0][1]], ["bad", "ACDEF"]], 7.4));
      console.log(JSON.stringify({csv: ABUI.toCsv(rows), tsv: ABUI.toTsv(rows), fasta: ABUI.toFasta(rows), head: ABUI.TABLE_HEAD}));""")
    csv = out["csv"]
    assert csv.startswith("﻿순번,서열 이름,")           # BOM so Excel reads Korean
    lines = csv.rstrip("\r\n").split("\r\n")
    assert len(lines) == 3
    assert lines[1].startswith('1,"a, ""x""",VH,기준 해당 없음,신뢰 가능,RAHSI,')
    assert "서열이 너무 짧다" in lines[2]
    tsv = out["tsv"].split("\n")
    assert all(len(line.split("\t")) == len(out["head"]) for line in tsv)
    fasta = out["fasta"].splitlines()
    assert fasta[0] == '>a, "x"' and max(len(x) for x in fasta[1:]) == 60 and len(fasta) == 3


@needs_node
def test_nice_ticks_cover_the_range():
    out = _node("""
      console.log(JSON.stringify([[6.2, 11.2, 7], [-0.6, 0.15, 6], [9.65, 9.65, 5]].map(a => ABUI.niceTicks(...a))));""")
    for t, (lo, hi) in zip(out, [(6.2, 11.2), (-0.6, 0.15), (9.65, 9.65)]):
        assert t["ticks"][0] <= lo and t["ticks"][-1] >= hi
        assert t["ticks"] == sorted(t["ticks"]) and len(t["ticks"]) >= 2
    assert out[0]["ticks"] == [6, 7, 8, 9, 10, 11, 12]


@needs_node
def test_motif_table_matches_python():
    out = _node("console.log(JSON.stringify(ABUI.MOTIFS));")
    assert [(m["pattern"], m["base"], m["top"]) for m in out] == [
        (s.pattern, s.base_risk, s.top_within_kind) for s in MOTIF_SPECS
    ]


@needs_node
def test_validation_figures_match_the_documents():
    """The 검증 기록 screen repeats README/VALIDATION numbers; they must not drift."""
    out = _node("console.log(JSON.stringify(ABUI.VALIDATION));")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    validation = (ROOT / "VALIDATION.md").read_text(encoding="utf-8")
    for key in ["sequences", "passedInput", "annotation", "boundary", "fr3Gap"]:
        assert out[key] in readme, (key, out[key])
    assert out["boundaryOf"] in validation and out["confident"] in validation
    assert out["chainCorrect"].replace(",", "") in validation
    for name, _cells, annot, conf, chain, mismatch in out["rows"]:
        lines = [line for line in validation.splitlines() if line.startswith("| " + name)]
        assert any(all(v in line for v in (annot, conf, chain)) for line in lines), (name, annot, conf, chain)
        assert any(mismatch in line for line in lines), (name, mismatch)


def test_app_screens_are_all_present():
    html = APP.read_text(encoding="utf-8")
    for view in ["new", "all", "review", "recheck", "errors", "compare", "properties", "rules", "validation", "guide"]:
        assert f'"{view}"' in html, view
    for label in ["+ 새 분석", "후보 전체", "CDR 비교", "물성 분포", "선별 규칙", "검증 기록", "사용법 · 한계"]:
        assert label in html, label
