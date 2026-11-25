"""The browser app (web/abcheck.html) and the helpers behind it (web/app-lib.js).

The engine itself is checked against Python in tests/test_js_parity.py. Here:
- the committed app is built from the current sources and makes no outside requests;
- app-lib.js, run under Node against real engine output, does what the screens rely on --
  in particular, the default triage rule set reproduces the engine's triage exactly.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from abcheck.liabilities import MOTIF_SPECS

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "web" / "abcheck.html"


def test_committed_app_is_current():
    result = subprocess.run(
        [sys.executable, "web/build.py", "--check"],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr + " — run: python3 web/build.py"


def test_app_makes_no_outside_requests():
    """Sequences must never leave the browser: no fetch/XHR, no external scripts, styles or fonts."""
    html = APP.read_text(encoding="utf-8")
    loaded = re.findall(r'<(?:script|link|img|iframe)[^>]+(?:src|href)="(https?://[^"]+)"', html)
    assert loaded == [], loaded
    assert "@import" not in html and "url(http" not in html
    for api in ["fetch(", "XMLHttpRequest", "WebSocket", "sendBeacon"]:
        assert api not in html, api


def test_fonts_are_bundled():
    html = APP.read_text(encoding="utf-8")
    assert html.count("font/woff2;base64,") == 2          # Pretendard 400 and 700
    assert (ROOT / "web" / "fonts" / "Pretendard-LICENSE.txt").exists()


def test_app_carries_the_limitation_marks():
    html = APP.read_text(encoding="utf-8")
    assert "<title>AbCheck 항체 서열 점검</title>" in html
    assert "numbering 엔진이 아니다" in html               # engine notice, shown on the guide page
    assert "종합 판정을 내리지 않는다" in html
    for mark in ["로컬 처리", "Kabat 근사"]:
        assert mark in html, mark


def test_example_uses_only_synthetic_sequences():
    html = APP.read_text(encoding="utf-8")
    example = re.search(r"const EXAMPLE = (\".*?\");", html).group(1)
    assert ">SYN-" in example and "OAS" not in example


def test_app_screens_are_all_present():
    html = APP.read_text(encoding="utf-8")
    for view in ["new", "all", "review", "unknown", "errors", "picked", "compare", "properties",
                 "onepage", "rules", "validation", "guide"]:
        assert f'"{view}"' in html, view
    for label in ["+ 새 분석", "검토 필요", "판정 불가", "CDR 비교", "물성 분포", "선별 규칙", "검증 기록", "사용법",
                  "개인 프로젝트 · 셀트리온과 무관"]:
        assert label in html, label


def test_design_tokens_follow_the_series():
    """Shared design system: green brand, sizes 13/15/17/19/24/32, weights 400/700 only."""
    css = re.search(r"<style>(.*?)</style>", (ROOT / "web" / "app.html").read_text(encoding="utf-8"), re.S).group(1)
    assert "--brand: #00733C" in css and "--brand: #5FE65A" in css
    sizes = set(re.findall(r"font-size:\s*([0-9.]+)px", css))
    assert sizes <= {"13", "15", "17", "19", "24", "32"}, sizes
    weights = set(re.findall(r"font-weight:\s*([0-9]+)", css))
    assert weights <= {"400", "700"}, weights
    assert "box-shadow" not in css
    for risky in ["#D92D20", "#C2410C"]:                      # status colours stay unused: risk is text only
        assert risky.lower() not in css.lower(), risky


# --------------------------------------------------------------------------
# web/app-lib.js under Node
# --------------------------------------------------------------------------

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


EDGE = """
  const s0 = entries[0][1], l0 = entries[2][1];
  const edge = [
    ["bad-short", "ACDEFGHIK"], ["bad-dna", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"],
    ["gapped", s0.slice(0, 30) + ".." + s0.slice(30)],
    ["vh-j-tail", s0.replace(/WG.G[A-Z]*$/, "WAQGTLVTVSS")],
    ["vl-j-tail", l0.replace(/FG.G[A-Z]*$/, "LGGGTKLTVL")],
  ];
  const all = entries.concat(edge);
"""


@needs_node
def test_default_rule_set_reproduces_the_engine_triage():
    """applyRules(DEFAULT_RULESET) must equal abcheck/triage.py for every sequence."""
    out = _node(EDGE + """
      const diffs = [];
      for (const ph of [7.4, 5.0]) {
        const p = ABCHECK.analyze(all, ph);
        p.sequences.forEach((s) => {
          const t = ABUI.applyRules(s, ABUI.DEFAULT_RULESET);
          const e = s.triage;
          if (t.status !== e.status || JSON.stringify(t.reasons) !== JSON.stringify(e.reasons) ||
              t.evaluable !== e.evaluable || t.note !== e.note) diffs.push([s.id, t, e]);
        });
      }
      console.log(JSON.stringify({n: all.length, diffs}));""")
    assert out["n"] == 29
    assert out["diffs"] == []


@needs_node
def test_custom_rule_sets_change_only_what_they_say():
    out = _node("""
      const p = ABCHECK.analyze(entries, 7.4);
      const by = (rs) => Object.fromEntries(ABUI.buildRows(p, rs).map(r => [r.id, [r.status, r.reasons]]));
      const base = ABUI.DEFAULT_RULESET;
      console.log(JSON.stringify({
        base: by(base),
        weak: by(Object.assign({}, base, {motifs: base.motifs.concat(["N[STNH]"])})),
        cdr3: by(Object.assign({}, base, {where: "cdr3"})),
        noCys: by(Object.assign({}, base, {cys: "off"})),
        anyCys: by(Object.assign({}, base, {cys: "any", motifs: []})),
      }));""")
    base, weak, cdr3 = out["base"], out["weak"], out["cdr3"]
    assert base["SYN-VH-010"][0] == "기준 해당 없음"           # NT in CDR-H2: high, but not top of kind
    assert weak["SYN-VH-010"][0] == "검토 필요"
    assert "deamidation (NS/NT/NN/NH)" in weak["SYN-VH-010"][1][0] or "(NT)" in weak["SYN-VH-010"][1][0]
    assert cdr3["SYN-VH-003"][0] == "기준 해당 없음"           # NG sits in CDR-H2, not CDR-H3
    assert cdr3["SYN-VH-007"][0] == "검토 필요"                # DG in CDR-H3
    assert out["noCys"]["SYN-VH-005"][0] == "기준 해당 없음"
    assert out["anyCys"]["SYN-VH-902"][0] == "판정 불가"        # annotation failed: never evaluated
    assert out["anyCys"]["SYN-VH-005"][1] == ["비결합 Cys 가능성 — 위험도 높음 (Cys 3개, 홀수)"]
    # an unconfident sequence never becomes "기준 해당 없음" under any rule set
    for rs in out.values():
        assert rs["SYN-VH-903"][0] in ("판정 불가", "검토 필요")


@needs_node
def test_json_export_carries_the_active_rule_set():
    out = _node("""
      const payload = ABCHECK.analyze(entries, 7.4);
      const before = ABCHECK.toJSON(payload);
      const team = {id: "team", name: "팀 규칙", motifs: ["NG", "D[STD]", "W"], where: "cdr", minRisk: "중간", cys: "off"};
      const same = ABUI.withRuleResults(payload, ABUI.DEFAULT_RULESET);
      const added = ABUI.withRuleResults(payload, team);
      const text = ABCHECK.toJSON(added);
      const rows = ABUI.buildRows(payload, team);
      console.log(JSON.stringify({
        defaultIsEngineOutput: same === payload && ABCHECK.toJSON(same) === before,
        payloadUntouched: ABCHECK.toJSON(payload) === before && !("rule_set" in payload) && !("rule_result" in payload.sequences[0]),
        topKeys: Object.keys(added), seqKeys: Object.keys(added.sequences[0]),
        ruleSet: added.rule_set,
        engineTriageKept: added.sequences.every((s, i) => s.triage === payload.sequences[i].triage),
        agreesWithScreen: added.sequences.every((s, i) => s.rule_result.status === rows[i].status
          && JSON.stringify(s.rule_result.reasons) === JSON.stringify(rows[i].reasons)),
        differsFromEngine: added.sequences.filter(s => s.rule_result.status !== s.triage.status).length,
        rereadIds: ABCHECK.readInput(text, "result.json")[0].length === entries.length,
        validJson: JSON.parse(text).rule_set.name}));""")
    assert out["defaultIsEngineOutput"] and out["payloadUntouched"]
    assert out["topKeys"] == ["tool", "generated", "coordinates", "limitations", "no_verdict", "count", "summary",
                              "rule_set", "sequences"]
    assert out["seqKeys"][out["seqKeys"].index("triage") + 1] == "rule_result"
    assert out["ruleSet"] == {"name": "팀 규칙", "motifs": ["NG", "D[STD]", "W"], "where": "cdr", "min_risk": "중간",
                              "cys": "off", "description": out["ruleSet"]["description"]}
    assert out["ruleSet"]["description"]
    assert out["engineTriageKept"] and out["agreesWithScreen"]
    assert out["differsFromEngine"] > 0             # the case that used to be exported wrongly
    assert out["rereadIds"] and out["validJson"] == "팀 규칙"


@needs_node
def test_rule_sets_are_validated_and_described():
    out = _node("""
      const errors = [{}, {name: "x", motifs: ["XX"], where: "cdr", minRisk: "높음", cys: "cdr"},
        {name: "x", motifs: [], where: "fw", minRisk: "높음", cys: "cdr"},
        {name: "", motifs: [], where: "cdr", minRisk: "높음", cys: "cdr"}].map((r) => {
          try { ABUI.normalizeRuleSet(r); return null; } catch (e) { return e.message; } });
      const ok = ABUI.normalizeRuleSet({name: " 팀 ", motifs: ["DP", "NG"], where: "cdr3", minRisk: "중간", cys: "off"});
      console.log(JSON.stringify({errors, ok, d0: ABUI.describeRuleSet(ABUI.DEFAULT_RULESET), d1: ABUI.describeRuleSet(ok)}));""")
    assert all(out["errors"]), out["errors"]
    assert out["ok"]["name"] == "팀" and out["ok"]["motifs"] == ["NG", "DP"]
    assert out["d0"] == ["CDR 안 · N-결합 당화 · NG · DG · DP · 위험도 높음", "비결합 Cys · CDR 안 Cys"]
    assert out["d1"] == ["CDR3 안 · NG · DP · 위험도 중간 이상", "비결합 Cys · 끔"]


@needs_node
def test_build_rows_joins_reports_and_summary_rows():
    out = _node("""
      const p = ABCHECK.analyze(entries.concat([["bad", "ACDEFGHIK"]]), 7.4);
      const rows = ABUI.buildRows(p);
      console.log(JSON.stringify({
        n: rows.length, ids: rows.map(r => r.id), review_count: p.summary.review_count,
        tally: ABUI.tally(rows, new Set([0, 2])),
        vh003: rows.find(r => r.id === "SYN-VH-003"),
        bad: rows.find(r => r.id === "bad"),
      }));""")
    t = out["tally"]
    assert out["n"] == 25 and out["ids"][-1] == "bad"
    assert t["review"] == out["review_count"]
    assert t["review"] + t["unknown"] + t["clear"] + t["errors"] == t["all"] == 25
    assert t["errors"] == 1 and t["picked"] == 2
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
def test_unknown_list_holds_the_unconfident_and_failed():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries, 7.4));
      console.log(JSON.stringify(rows.filter(r => r.cats[0] === "unknown").map(r => r.id).sort()));""")
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
        cdrResidues: ids(ABUI.filterRows(rows, {query: "ng"})),
        seqResidues: ids(ABUI.filterRows(rows, {query: "YITTIE"})).length,
        motifName: ids(ABUI.filterRows(rows, {query: "fragmentation"})),
        lowConf: ids(ABUI.filterRows(rows, {conf: "no"})).length,
        picked: ids(ABUI.filterRows(rows, {view: "picked", picked: new Set([1, 3])})),
        input: ids(ABUI.sortRows(rows, "input")).join() === ids(rows).join(),
        byName: ids(ABUI.sortRows(rows, "name")).slice(0, 3),
        piDesc: byPi.filter(r => r.pI !== null).every((r, i, a) => i === 0 || a[i - 1].pI >= r.pI),
      }));""")
    assert out["vl"] == ["SYN-VL-001", "SYN-VL-002", "SYN-VL-003", "SYN-VL-004", "SYN-VL-005", "SYN-VL-007", "SYN-VL-006"]
    assert out["name"] == 9                             # SYN-VH-001 … SYN-VH-009
    assert "SYN-VH-003" in out["cdrResidues"]          # NG inside CDR-H2
    assert out["seqResidues"] == 17                     # every VH-named record shares the scaffold start
    assert out["motifName"] == ["SYN-VL-004"]          # the only DP
    assert out["lowConf"] == 5                          # 013, 014, 901, 902, 903
    assert out["picked"] == ["SYN-VH-002", "SYN-VL-002"]
    assert out["input"] is True
    assert out["byName"] == ["SYN-VH-001", "SYN-VH-002", "SYN-VH-003"]
    assert out["piDesc"] is True


@needs_node
def test_combine_entries_from_several_files():
    out = _node("""
      const a = ABCHECK.readInput(">a\\n" + entries[0][1] + "\\n", "x.fasta");
      const b = ABCHECK.readInput("id,sequence\\nb," + entries[2][1] + "\\n", "y.csv");
      const one = ABUI.combineEntries([{name: "x.fasta", entries: a[0], source: a[1]}]);
      const two = ABUI.combineEntries([{name: "x.fasta", entries: a[0], source: a[1]}, {name: "y.csv", entries: b[0], source: b[1]}]);
      let err = null; try { ABUI.combineEntries([]); } catch (e) { err = e.message; }
      console.log(JSON.stringify({one, two: {ids: two.entries.map(e => e[0]), source: two.source}, err}));""")
    assert out["one"]["source"] == "FASTA (x.fasta)"
    assert out["two"] == {"ids": ["a", "b"], "source": "x.fasta 외 1개"}
    assert out["err"]


@needs_node
def test_pasted_tables_are_recognised():
    """Rows copied from Excel used to be read as one sequence and fail with a misleading error."""
    out = _node("""
      const s = entries.map(e => e[1]);
      const kinds = {
        fasta: ABUI.pastedKind(">a\\n" + s[0]),
        single: ABUI.pastedKind(s[0]),
        wrapped: ABUI.pastedKind(s[0].slice(0, 60) + "\\n" + s[0].slice(60)),
        rows: ABUI.pastedKind("mAb1_VH\\t" + s[0] + "\\nmAb1_VL\\t" + s[2]),
        header: ABUI.pastedKind("이름\\t서열\\r\\nmAb1_VH\\t" + s[0] + "\\r\\n"),
        column: ABUI.pastedKind(s[0] + "\\n" + s[2] + "\\n"),
        csv: ABUI.pastedKind("id,sequence\\nx," + s[0]),
        json: ABUI.pastedKind('  [{"id": "x"}]'),
        empty: ABUI.pastedKind("  "),
      };
      const read = (t) => ABCHECK.readInput(t, "pasted" + ABUI.pastedKind(t))[0].map(e => e[0]);
      console.log(JSON.stringify({kinds,
        rows: read("mAb1_VH\\t" + s[0] + "\\nmAb1_VL\\t" + s[2]),
        header: read("이름\\t서열\\nmAb1_VH\\t" + s[0]),
        column: read(s[0] + "\\n" + s[2])}));""")
    assert out["kinds"] == {"fasta": "", "single": "", "wrapped": "", "rows": ".tsv", "header": ".tsv",
                            "column": ".tsv", "csv": ".csv", "json": ".json", "empty": ""}
    assert out["rows"] == ["mAb1_VH", "mAb1_VL"]
    assert out["header"] == ["mAb1_VH"]
    assert out["column"] == ["seq1", "seq2"]


@needs_node
def test_position_numbers_and_closing_stop_are_removed_on_request():
    out = _node("""
      const seq = entries[0][1];
      const numbered = ">mAb12 lot 2024\\n" + seq.match(/.{1,10}/g).map((c, i) => String(i * 10 + 1).padStart(6) + " " + c.toLowerCase()).join("\\n") + "*\\n";
      const a = ABUI.stripNumbering(numbered);
      const again = ABCHECK.analyze(ABCHECK.readInput(a.text, "pasted.fasta")[0], 7.4).sequences[0];
      const clean = ABCHECK.analyze([["x", seq]], 7.4).sequences[0];
      const mid = ABUI.stripNumbering(">a\\nACDEF*GHIK\\n>b\\nLMNPQ*  \\n");
      const table = "mAb01\\t" + seq + "\\nmAb02\\t" + seq;
      console.log(JSON.stringify({
        header: a.text.split("\\n")[0], digits: a.digits, stops: a.stops,
        rejectedBefore: ABCHECK.analyze(ABCHECK.readInput(numbered, "pasted.fasta")[0], 7.4).sequences[0].input_ok,
        okAfter: again.input_ok, sameSequence: again.sequence === clean.sequence,
        sameCdrs: JSON.stringify(again.annotation.cdrs) === JSON.stringify(clean.annotation.cdrs),
        mid: mid, table: ABUI.stripNumbering(table), tableSame: ABUI.stripNumbering(table).text === table,
        nothing: ABUI.stripNumbering(">a\\n" + seq)}));""")
    assert out["header"] == ">mAb12 lot 2024"                    # digits in a name are not touched
    assert out["digits"] > 0 and out["stops"] == 1
    assert out["rejectedBefore"] is False and out["okAfter"] is True
    assert out["sameSequence"] and out["sameCdrs"]
    assert out["mid"] == {"text": ">a\nACDEF*GHIK\n>b\nLMNPQ  \n", "digits": 0, "stops": 1}   # a stop in the middle stays
    assert out["tableSame"] and out["table"]["digits"] == 0      # names in a pasted table keep their digits
    assert out["nothing"]["digits"] == 0 and out["nothing"]["stops"] == 0


@needs_node
def test_partner_chains_are_found_by_name():
    out = _node("""
      const s = Object.fromEntries(entries);
      const named = [["mAb01_VH", s["SYN-VH-001"]], ["mAb02_VH", s["SYN-VH-003"]], ["mAb01_VL", s["SYN-VL-001"]],
                     ["mAb02_VL", s["SYN-VL-004"]], ["lonely_VH", s["SYN-VH-007"]], ["dup_VH", s["SYN-VH-001"]],
                     ["dup_VH", s["SYN-VH-003"]], ["dup_VL", s["SYN-VL-001"]], ["plain", s["SYN-VH-002"]]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const pairs = ABUI.pairMap(rows);
      const picked = rows.filter(r => ["mAb02_VL", "mAb01_VL", "mAb01_VH", "plain"].includes(r.id));
      console.log(JSON.stringify({
        keys: ["mAb01_VH", "SYN-VL-003", "Ab7 heavy", "mAb01VH", "plain", "VH"].map(ABUI.pairKey),
        pairs, order: ABUI.withPartnersAdjacent(picked, pairs).map(r => r.id)}));""")
    assert out["keys"] == ["mab01", "syn-003", "ab7", "mab01", None, None]
    assert out["pairs"] == {"0": 2, "2": 0, "1": 3, "3": 1}      # lonely, ambiguous and plain names stay unpaired
    assert out["order"] == ["mAb01_VH", "mAb01_VL", "mAb02_VL", "plain"]


@needs_node
def test_antibody_units_keep_both_chains_and_their_own_results():
    out = _node("""
      const s = Object.fromEntries(entries);
      const named = [["mAb01_VL", s["SYN-VL-001"]], ["mAb02_VH", s["SYN-VH-003"]], ["mAb01_VH", s["SYN-VH-001"]],
                     ["mAb02_VL", s["SYN-VL-005"]], ["solo-VH", s["SYN-VH-007"]]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const pairs = ABUI.pairMap(rows);
      const show = (units) => units.map(u => u.h ? [u.name, u.h.id, u.h.status, u.l.id, u.l.status] : [u.row.id]);
      console.log(JSON.stringify({
        names: ["mAb02_VH", "SYN-VL-003", "Ab7 heavy", "mAb01VH", "plain"].map(ABUI.pairName),
        all: show(ABUI.antibodyUnits(rows, rows, pairs)),
        review: show(ABUI.antibodyUnits(rows.filter(r => r.cats[0] === "review"), rows, pairs))}));""")
    assert out["names"] == ["mAb02", "SYN-003", "Ab7", "mAb01", None]
    assert out["all"] == [["mAb01", "mAb01_VH", "기준 해당 없음", "mAb01_VL", "기준 해당 없음"],
                          ["mAb02", "mAb02_VH", "검토 필요", "mAb02_VL", "기준 해당 없음"], ["solo-VH"]]
    # filtered to 검토 필요: the flagged chain brings its partner, and each keeps its own status
    assert out["review"] == [["mAb02", "mAb02_VH", "검토 필요", "mAb02_VL", "기준 해당 없음"], ["solo-VH"]]


@needs_node
def test_shared_cdrs_finds_identical_cdrs():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze([["a", entries[0][1]], ["b", entries[0][1]], ["c", entries[2][1]]], 7.4));
      console.log(JSON.stringify(ABUI.sharedCdrs(rows)));""")
    assert [g["name"] for g in out] == ["CDR-H1", "CDR-H2", "CDR-H3"]
    assert all(g["ids"] == ["a", "b"] for g in out)


@needs_node
def test_csv_tsv_fasta_exports():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze([["a, \\"x\\"", entries[0][1]], ["bad", "ACDEF"]], 7.4));
      console.log(JSON.stringify({csv: ABUI.toCsv(rows, "기본"), tsv: ABUI.toTsv(rows, "기본"), fasta: ABUI.toFasta(rows), head: ABUI.TABLE_HEAD}));""")
    csv = out["csv"]
    assert csv.startswith("﻿순번,서열 이름,")           # BOM so Excel reads Korean
    lines = csv.rstrip("\r\n").split("\r\n")
    assert len(lines) == 3
    assert lines[1].startswith('1,"a, ""x""",VH,기준 해당 없음,신뢰 가능,RAHSI,')
    assert ",기본," in lines[1] and "서열이 너무 짧다" in lines[2]
    tsv = out["tsv"].split("\n")
    assert all(len(line.split("\t")) == len(out["head"]) for line in tsv)
    fasta = out["fasta"].splitlines()
    assert fasta[0] == '>a, "x"' and max(len(x) for x in fasta[1:]) == 60 and len(fasta) == 3


@needs_node
def test_mail_text_lists_rule_hits_without_verdicts():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries.concat([["bad", "ACDEF"]]), 7.4));
      console.log(JSON.stringify(ABUI.toMailText(rows, {ruleName: "기본", when: "2026-11-25 18:40"})));""")
    lines = out.split("\n")
    assert lines[0] == "AbCheck · 후보 25건 · 규칙 기본 · 2026-11-25 18:40"
    assert lines[1] == "검토 필요 6 · 판정 불가 5 · 기준 해당 없음 13 · 입력 오류 1"
    assert "- SYN-VH-003 (VH): CDR-H2 56번에 deamidation (NG) — 위험도 높음" in lines
    assert "- SYN-VH-902 (VH): CDR 주석 실패" in lines and "- SYN-VH-013 (VH): 주석 신뢰 불가" in lines
    assert any(line.startswith("- bad: 서열이 너무 짧다 (5 잔기).") for line in lines)
    for word in ["통과", "불합격", "안전", "pass", "fail", "점수"]:
        assert word not in out, word


@needs_node
def test_mail_text_keeps_an_antibodys_chains_together():
    out = _node("""
      const s = Object.fromEntries(entries);
      // as a sheet gives them: heavy chains first, light chains after
      const named = [["mAb01_VH", s["SYN-VH-003"]], ["mAb02_VH", s["SYN-VH-007"]], ["mAb01_VL", s["SYN-VL-004"]], ["mAb02_VL", s["SYN-VL-001"]]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const ids = (text) => text.split("\\n").filter(l => l.startsWith("- ")).map(l => l.slice(2).split(" ")[0]);
      console.log(JSON.stringify({
        plain: ids(ABUI.toMailText(rows, {ruleName: "기본", when: "t"})),
        paired: ids(ABUI.toMailText(rows, {ruleName: "기본", when: "t", pairs: ABUI.pairMap(rows)}))}));""")
    assert out["plain"] == ["mAb01_VH", "mAb02_VH", "mAb01_VL"]        # input order
    assert out["paired"] == ["mAb01_VH", "mAb01_VL", "mAb02_VH"]       # partner moved next to its heavy chain


@needs_node
def test_one_page_sheet():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries, 7.4)).filter(r => ["SYN-VH-003", "SYN-VH-902"].includes(r.id));
      console.log(JSON.stringify({frag: ABUI.onePageFragment(rows, {ruleName: "기본", when: "11/25"}),
        doc: ABUI.onePageDoc(rows, {ruleName: "<b>", when: "11/25"})}));""")
    frag, doc = out["frag"], out["doc"]
    assert "후보 2건" in frag and "규칙 기본" in frag
    assert "LSTAAR<u " in frag                          # NG underlined inside CDR-H2
    assert "CDR 주석 실패" in frag and "검토 필요" not in frag
    assert doc.startswith("<!doctype html>") and "@page{size:A4 landscape" in doc
    assert "&lt;b&gt;" in doc and "<b></span>" not in doc  # rule names are escaped


@needs_node
def test_scatter_answers_for_the_nearest_point():
    out = _node("""
      const pts = [{x: 100, y: 100}, {x: 108, y: 104}, {x: 300, y: 50}];
      console.log(JSON.stringify([
        ABUI.nearestPoint(pts, 100, 100, 14),      // on the first mark, although the second is within reach
        ABUI.nearestPoint(pts, 107, 104, 14),
        ABUI.nearestPoint(pts, 104, 102, 14),      // exactly between: the earlier point keeps it
        ABUI.nearestPoint(pts, 200, 200, 14),      // nothing near
        ABUI.nearestPoint(pts, 300, 64, 14),       // at the edge of the reach
        ABUI.nearestPoint(pts, 300, 65, 14),
        ABUI.nearestPoint([], 1, 1, 14)]));""")
    assert out == [0, 1, 0, -1, 2, -1, -1]


@needs_node
def test_scatter_names_never_collide():
    out = _node("""
      const opts = {charW: 8, h: 16, gap: 10, r: 8, left: 50, right: 700, top: 10, bottom: 400};
      const hits = (a, b) => a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;
      const check = (pts) => {
        const boxes = ABUI.placeLabels(pts, opts);
        const marks = pts.map(p => ({x0: p.x - 8, x1: p.x + 8, y0: p.y - 8, y1: p.y + 8}));
        return {
          n: boxes.length, sides: boxes.map(b => b.i + b.side[0]).join(" "),
          labelOnLabel: boxes.some((a, i) => boxes.some((b, j) => i < j && hits(a, b))),
          labelOnMark: boxes.some(b => marks.some((m, j) => j !== b.i && hits(b, m))),
          outside: boxes.some(b => b.x0 < opts.left || b.x1 > opts.right || b.y0 < opts.top || b.y1 > opts.bottom)};
      };
      // the real case: every synthetic sequence at its pI / GRAVY position
      const rows = ABUI.buildRows(ABCHECK.analyze(entries, 7.4)).filter(r => r.pI !== null);
      const real = rows.map(r => ({x: 50 + (r.pI - 4) * 80, y: 400 - (r.gravy + 1.2) * 250, text: r.id}));
      let seed = 7; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
      const dense = Array.from({length: 300}, (_, i) => ({x: 300 + rnd() * 120, y: 150 + rnd() * 80, text: "cand-" + i}));
      console.log(JSON.stringify({
        apart: check([{x: 100, y: 100, text: "aaaa"}, {x: 400, y: 300, text: "bbbb"}]),
        edge: check([{x: 690, y: 100, text: "near-right-edge"}]),
        stacked: check([{x: 300, y: 200, text: "one"}, {x: 303, y: 201, text: "two"}, {x: 301, y: 199, text: "three"}]),
        // the same chain under three names sits on one spot: one name, the first one's
        same: check([{x: 300, y: 200, text: "first"}, {x: 300, y: 200, text: "second"}, {x: 300.4, y: 199.8, text: "third"}]),
        // a point with a neighbour close on each side: its name goes above, the next one's below
        fenced: check([{x: 300, y: 200, text: "middle-one"}, {x: 280, y: 200, text: ""}, {x: 320, y: 200, text: ""},
                       {x: 300, y: 172, text: ""}]),
        real: check(real), realPoints: real.length,
        dense: check(dense), empty: ABUI.placeLabels([], opts)}));""")
    assert out["apart"]["sides"] == "0r 1r"
    assert out["edge"]["sides"] == "0l"                      # no room on the right: the name goes left
    stacked = out["stacked"]["sides"].split()                # three points on one spot: each name takes another side
    assert len(stacked) == 3 and len({side[-1] for side in stacked}) == 3
    assert out["same"]["sides"] == "0r"
    assert "0b" in out["fenced"]["sides"].split()            # sides and the top are taken: the name goes below
    for case in ["apart", "edge", "stacked", "fenced", "real", "dense"]:
        assert not out[case]["labelOnLabel"] and not out[case]["labelOnMark"] and not out[case]["outside"], case
    assert 0 < out["real"]["n"] <= out["realPoints"]
    assert out["dense"]["n"] < 300 and out["empty"] == []


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
    assert [(m["pattern"], m["label"], m["base"], m["top"]) for m in out] == [
        (s.pattern, s.label, s.base_risk, s.top_within_kind) for s in MOTIF_SPECS
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
