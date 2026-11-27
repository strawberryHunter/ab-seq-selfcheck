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
    for view in ["new", "all", "review", "unknown", "errors", "picked", "compare", "properties", "motifs",
                 "onepage", "rules", "validation", "guide"]:
        assert f'"{view}"' in html, view
    for label in ['${icon("plus", 18)}새 분석', "검토 필요", "판정 불가", "CDR 비교", "물성 분포", "모티프 분포", "선별 규칙", "검증 기록", "사용법",
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


def test_shell_follows_the_series():
    """Values shared by the three tools (portfolio-prompts/4-시리즈-맞추기.md); the measuring is done by series-check.mjs."""
    src = (ROOT / "web" / "app.html").read_text(encoding="utf-8")
    css = re.search(r"<style>(.*?)</style>", src, re.S).group(1)
    assert "--edge: #86908C" in css                           # outer border of inputs, dropdowns, segments (3:1 on white)
    assert "background: var(--bright)" not in css             # no green dot on the status chip: green is brand / action only
    assert 'stroke-width="2"' in src and 'stroke-width="1.8"' not in src
    # the series line: the two other tools are links the viewer may follow, this tool is not a link
    links = re.findall(r'\["([^"]+)", "(https://claude\.ai/artifact/[A-Za-z0-9]+)"\]', src)
    assert [name for name, _ in links] == ["RegRadar", "Deviation 데스크"]
    assert '["AbCheck", ""]' in src and 'target="_blank" rel="noopener noreferrer"' in src
    # every address in the built page: these links and the validation record, nothing else
    addresses = set(re.findall(r"https?://[^\s\"'<>)]+", APP.read_text(encoding="utf-8")))
    assert addresses == {url for _, url in links} | {"https://github.com/strawberryHunter/ab-seq-selfcheck/blob/main/VALIDATION.md"}


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
def test_another_ph_changes_net_charge_and_nothing_else():
    """The pH field on 물성 분포 re-runs the list and keeps picks and notes: nothing they rest on may move."""
    out = _node(EDGE + r"""
      const strip = (v) => JSON.stringify(v, (k, x) => (k === "net_charge" || k === "pH" || k === "charge" || k === "generated") ? undefined : x);
      const a = ABCHECK.analyze(all, 7.4), b = ABCHECK.analyze(all, 5.5);
      const ra = ABUI.buildRows(a), rb = ABUI.buildRows(b);
      console.log(JSON.stringify({
        count: all.length,
        payloadSame: strip(a) === strip(b),
        rowsSame: strip(ra) === strip(rb),
        phs: [...new Set(rb.filter(r => r.pH !== null).map(r => r.pH))],
        chargeMoved: ra.filter((r, i) => r.charge !== null && r.charge !== rb[i].charge).length,
        withCharge: ra.filter(r => r.charge !== null).length,
        lowerPhMorePositive: ra.every((r, i) => r.charge === null || rb[i].charge >= r.charge)}));""")
    assert out["payloadSame"] and out["rowsSame"]              # classification, CDRs, motifs, pI, GRAVY: unchanged
    assert out["phs"] == [5.5]
    assert out["chargeMoved"] == out["withCharge"] > 0 and out["lowerPhMorePositive"]


@needs_node
def test_json_result_carries_the_marks_and_gives_them_back():
    out = _node(r"""
      const payload = ABCHECK.analyze(entries, 6.0);
      const before = ABCHECK.toJSON(payload);
      const team = {id: "team", name: "팀 규칙", motifs: ["NG", "D[STD]", "W"], where: "cdr", minRisk: "중간", cys: "off"};
      const ruled = ABUI.withRuleResults(payload, team);
      const marked = ABUI.withReview(ruled, new Set([0, 2]), {2: "발현량 높음", 4: "재확인"});
      const text = ABCHECK.toJSON(marked);
      const reread = ABCHECK.readInput(text, "result.json")[0];
      const back = ABUI.reviewFromJson(text, reread.length);
      // a file edited by hand: marks of another shape are left out, a note becomes one line of 200
      const edited = JSON.parse(text);
      edited.sequences[0].review = {picked: "yes", note: 5};
      edited.sequences[1].review = {picked: true, note: " 두 줄\n메모 " + "가".repeat(300)};
      edited.sequences[2].review = "picked";
      edited.rule_set = {name: "깨진 규칙", motifs: ["XX"], where: "cdr", min_risk: "중간", cys: "off"};
      const hand = ABUI.reviewFromJson(JSON.stringify(edited), reread.length);
      let mismatch = null; try { ABUI.reviewFromJson(text, reread.length + 1); } catch (e) { mismatch = e.message; }
      console.log(JSON.stringify({
        noMarksSameObject: ABUI.withReview(payload, new Set(), {}) === payload && ABUI.withReview(payload, [], undefined) === payload,
        payloadUntouched: ABCHECK.toJSON(payload) === before && payload.sequences.every(s => !("review" in s)),
        withReview: marked.sequences.map((s, i) => "review" in s ? i : -1).filter(i => i >= 0),
        lastKey: Object.keys(marked.sequences[2]).pop(),
        reviews: [0, 2, 4].map(i => marked.sequences[i].review),
        unmarkedIsTheSameObject: marked.sequences[1] === ruled.sequences[1],
        engineKeysKept: Object.keys(payload.sequences[2]).every(k => marked.sequences[2][k] === ruled.sequences[2][k]),
        listAsPicked: JSON.stringify(ABUI.withReview(payload, [0, 2], {})) === JSON.stringify(ABUI.withReview(payload, new Set([0, 2]), {})),
        rereadSame: JSON.stringify(reread) === JSON.stringify(entries),
        back, sameRules: ABUI.sameRules(back.ruleSet, team), defaultRules: ABUI.sameRules(back.ruleSet, ABUI.DEFAULT_RULESET),
        renamedIsSame: ABUI.sameRules(Object.assign({}, team, {name: "다른 이름", motifs: ["W", "NG", "D[STD]"]}), team),
        hand: {picked: hand.picked, notes: hand.notes, ruleSet: hand.ruleSet, problem: hand.ruleProblem},
        cli: ABUI.reviewFromJson(ABCHECK.toJSON(ABCHECK.analyze(entries, 7.4)), entries.length),
        corpus: ABUI.reviewFromJson(JSON.stringify({sequences: entries.map(e => ({id: e[0], sequence: e[1]}))}), entries.length),
        list: ABUI.reviewFromJson(JSON.stringify(entries.map(e => e[1])), entries.length),
        nan: ABUI.reviewFromJson('{"sequences": [], "x": NaN}', 0), mismatch,
        bom: ABUI.reviewFromJson("﻿" + text, reread.length).picked}));""")
    assert out["noMarksSameObject"] and out["payloadUntouched"]          # without marks: the file as it was before
    assert out["withReview"] == [0, 2, 4] and out["lastKey"] == "review"
    assert out["reviews"] == [{"picked": True, "note": ""}, {"picked": True, "note": "발현량 높음"},
                              {"picked": False, "note": "재확인"}]
    assert out["unmarkedIsTheSameObject"] and out["engineKeysKept"] and out["listAsPicked"]
    assert out["rereadSame"]                                             # the engine's reader takes the file as before
    back = out["back"]
    assert back["picked"] == [0, 2] and back["notes"] == {"2": "발현량 높음", "4": "재확인"}
    assert back["ph"] == 6.0 and back["ruleProblem"] == "" and back["unread"] == ""
    assert back["ruleSet"] == {"name": "팀 규칙", "motifs": ["NG", "D[STD]", "W"], "where": "cdr", "minRisk": "중간", "cys": "off"}
    assert out["sameRules"] and not out["defaultRules"] and out["renamedIsSame"]
    hand = out["hand"]
    assert hand["picked"] == [1] and list(hand["notes"]) == ["1", "4"]
    assert hand["notes"]["1"].startswith("두 줄 메모 가") and len(hand["notes"]["1"]) == 200 and "\n" not in hand["notes"]["1"]
    assert hand["ruleSet"] is None and "motifs must be a list" in hand["problem"]
    assert out["cli"] == {"picked": [], "notes": {}, "refs": [], "ph": 7.4, "ruleSet": None, "ruleProblem": "", "unread": ""}
    assert out["corpus"]["ph"] is None and out["corpus"]["picked"] == []
    assert out["list"] is None                                           # a bare list of sequences has nothing to give back
    assert out["nan"]["unread"] and out["nan"]["picked"] == []
    assert "entries read" in out["mismatch"]
    assert out["bom"] == [0, 2]


@needs_node
def test_the_comparison_reference_is_kept_with_the_analysis():
    out = _node(r"""
      const s = Object.fromEntries(entries);
      const named = [["h1", s["SYN-VH-003"]], ["l1", s["SYN-VL-003"]], ["h2", s["SYN-VH-001"]],
                     ["noCdr", ABUI.applySubstitutions(s["SYN-VH-003"], [{pos: 36, from: "W", to: "A"}])], ["bad", "ACDEFGHIK"], ["l2", s["SYN-VL-001"]]];
      const payload = ABCHECK.analyze(named, 7.4), rows = ABUI.buildRows(payload);
      const R = (v) => ABUI.refOfRows(rows, v);
      const saved = ABUI.withReview(payload, [], {}, {H: 2, L: 1});
      const text = ABCHECK.toJSON(saved);
      const back = ABUI.reviewFromJson(text, named.length);
      const stored = (ref) => ABUI.cleanRecent([{id: "r", entries: named, ph: 7.4, at: 1764000000000, ref}]).list[0].ref;
      console.log(JSON.stringify({
        chains: rows.map(r => r.id + ":" + r.chain + ":" + r.cdrs.length),
        ok: R({H: 2, L: 1}), none: [R(null), R(undefined), R({}), R("x"), R({H: null, L: null})],
        // a place of the other chain, a row without CDRs, an input error, a place outside the list: not a reference
        wrong: [R({H: 1, L: 0}), R({H: 3, L: 4}), R({H: 99, L: -1}), R({H: "2", L: 1.5})],
        fromPlaces: [R([2, 1]), R([5, 1, 0, 2]), R([3, 4, 99, "x"]), R([])],
        noRefSameObject: ABUI.withReview(payload, [], {}, {H: null, L: null}) === payload && ABUI.withReview(payload, [], {}) === payload,
        reviews: saved.sequences.map(q => q.review || null),
        back: {picked: back.picked, refs: back.refs}, restored: R(back.refs),
        withPick: ABUI.withReview(payload, [2], {}, {H: 2, L: null}).sequences[2].review,
        stored: [stored({H: 2, L: 1}), stored(undefined), stored({H: 99, L: "1"}), stored([2, 1]), stored("x")]}));""")
    assert out["chains"] == ["h1:H:3", "l1:L:3", "h2:H:3", "noCdr:H:0", "bad:null:0", "l2:L:3"]
    none = {"H": None, "L": None}
    assert out["ok"] == {"H": 2, "L": 1} and out["none"] == [none] * 5
    assert out["wrong"] == [none] * 4
    assert out["fromPlaces"] == [{"H": 2, "L": 1}, {"H": 0, "L": 5}, none, none]      # the first of a chain wins
    assert out["noRefSameObject"]                                    # no reference: the file as it was
    assert out["reviews"] == [None, {"picked": False, "note": "", "reference": True}, {"picked": False, "note": "", "reference": True},
                              None, None, None]
    assert out["back"] == {"picked": [], "refs": [1, 2]} and out["restored"] == {"H": 2, "L": 1}
    assert out["withPick"] == {"picked": True, "note": "", "reference": True}
    assert out["stored"] == [{"H": 2, "L": 1}, none, none, none, none]  # shape only; the rows are checked after the analysis


@needs_node
def test_a_result_with_rejected_inputs_can_be_opened_again():
    """The engine's JSON reader wants a sequence in every item; a result has none for a rejected input."""
    out = _node(r"""
      const named = entries.slice(0, 3).concat([["bad", "ACDEFGHIK"], ["empty", ""], ["dna", "atggcc attgta\nATGGGC"], ["dup", entries[0][1]], [entries[0][0], entries[1][1]]]);
      const payload = ABCHECK.analyze(named, 6.0);
      const before = ABCHECK.toJSON(payload);
      let engineSays = null; try { ABCHECK.readInput(before, "result.json"); } catch (e) { engineSays = e.message; }
      const saved = ABUI.withReview(payload, [3, 6], {4: "서열 다시 받기"}, {H: 0, L: null}, named);
      const text = ABCHECK.toJSON(saved);
      let stillRefused = null; try { ABCHECK.readInput(text, "result.json"); } catch (e) { stillRefused = e.message; }
      const own = ABUI.entriesFromResult(text), old = ABUI.entriesFromResult(before);
      const again = ABCHECK.analyze(own.entries, 6.0);
      const rv = ABUI.reviewFromJson(text, own.entries.length);
      const clean = entries.slice(0, 3), cleanPayload = ABCHECK.analyze(clean, 7.4);
      let wrongCount = null; try { ABUI.withReview(payload, [], {}, null, named.slice(1)); } catch (e) { wrongCount = e.message; }
      console.log(JSON.stringify({
        engineSays, stillRefused,
        payloadUntouched: ABCHECK.toJSON(payload) === before,
        reviews: saved.sequences.map(q => q.review || null),
        entriesBack: JSON.stringify(own.entries) === JSON.stringify(named.map(([id, x], i) => [id, payload.sequences[i].input_ok ? payload.sequences[i].sequence : x])),
        lost: own.lost, sameAnalysis: ABCHECK.toJSON(again) === before,
        marks: {picked: rv.picked, notes: rv.notes, refs: rv.refs, ph: rv.ph},
        oldFile: {lost: old.lost, entries: old.entries.map(e => e[0] + ":" + e[1].length)},
        noRejectedSameObject: ABUI.withReview(cleanPayload, [], {}, null, clean) === cleanPayload,
        notAResult: [ABUI.entriesFromResult("{"), ABUI.entriesFromResult("[1, 2]"), ABUI.entriesFromResult('{"sequences": []}'),
                     ABUI.entriesFromResult('{"sequences": [{"id": "a", "sequence": "ACD"}]}'), ABUI.entriesFromResult('{"sequences": ["ACD"]}')],
        wrongCount}));""")
    assert "'sequence' 키가 없다" in out["engineSays"] and out["stillRefused"] == out["engineSays"]      # the engine's reader is as it was
    assert out["payloadUntouched"]
    assert out["reviews"] == [
        {"picked": False, "note": "", "reference": True}, None, None,
        {"picked": True, "note": "", "input": "ACDEFGHIK"}, {"picked": False, "note": "서열 다시 받기", "input": ""},
        {"picked": False, "note": "", "input": "atggcc attgta\nATGGGC"}, {"picked": True, "note": ""}, None]
    assert out["entriesBack"] and out["lost"] == 0 and out["sameAnalysis"]       # same list, same result, same errors
    assert out["marks"] == {"picked": [3, 6], "notes": {"4": "서열 다시 받기"}, "refs": [0], "ph": 6.0}
    # a file saved before the text was kept: the entry stays in its place with an empty sequence, and is counted
    assert out["oldFile"]["lost"] == 3 and out["oldFile"]["entries"][3:6] == ["bad:0", "empty:0", "dna:0"]
    assert out["noRejectedSameObject"]                                          # no rejected input, no mark: the file as it was
    assert out["notAResult"] == [None] * 5                                      # a corpus is the engine reader's business
    assert "8 sequences" in out["wrongCount"]


@needs_node
def test_marks_of_files_opened_together_are_renumbered():
    out = _node(r"""
      const a = {picked: [0, 2], notes: {1: "가"}, refs: [3], ph: 6.0, ruleSet: {name: "팀"}, ruleProblem: "", unread: ""};
      const b = {picked: [1], notes: {0: "나", 3: "다"}, ph: 7.4, ruleSet: null, ruleProblem: "", unread: "Unexpected token"};
      console.log(JSON.stringify({
        one: ABUI.joinReviews([{count: 5, review: a}]),
        three: ABUI.joinReviews([{count: 5, review: a}, {count: 2, review: null}, {count: 4, review: b}]),
        fastaOnly: ABUI.joinReviews([{count: 3, review: null}, {count: 2, review: null}]),
        untouched: JSON.stringify(a.picked) === "[0,2]"}));""")
    assert out["one"] == {"picked": [0, 2], "notes": {"1": "가"}, "refs": [3], "ph": 6.0, "ruleSet": {"name": "팀"}, "ruleProblem": "", "unread": ""}
    # the second file starts at place 5, the third at place 7; pH and rules of several files are not taken over
    assert out["three"] == {"picked": [0, 2, 8], "notes": {"1": "가", "7": "나", "10": "다"}, "refs": [3], "ph": None, "ruleSet": None,
                            "ruleProblem": "", "unread": "Unexpected token"}
    assert out["fastaOnly"] is None and out["untouched"]


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
def test_a_rule_set_travels_as_text():
    out = _node("""
      const team = {id: "s123", name: "DS 포함", motifs: ["NG", "D[STD]", "W"], where: "cdr3", minRisk: "중간", cys: "off"};
      const line = ABUI.ruleSetToText(team);
      const back = ABUI.ruleSetFromText(line);
      const payload = ABCHECK.analyze(entries, 7.4);
      const savedJson = ABCHECK.toJSON(ABUI.withRuleResults(payload, team));
      const fails = (text) => { try { ABUI.ruleSetFromText(text); return null; } catch (e) { return e.message; } };
      const same = (a, b) => JSON.stringify(ABUI.buildRows(payload, Object.assign({id: "x"}, a)).map(r => [r.status, r.reasons]))
                          === JSON.stringify(ABUI.buildRows(payload, Object.assign({id: "y"}, b)).map(r => [r.status, r.reasons]));
      console.log(JSON.stringify({
        line, oneLine: !line.includes("\\n"), back,
        sameResults: same(team, back),
        inMail: ABUI.ruleSetFromText("민준님, 규칙 보냅니다 > " + line + " 확인 부탁"),
        fromSavedJson: ABUI.ruleSetFromText(savedJson),
        defaultSet: ABUI.ruleSetFromText(ABUI.ruleSetToText(ABUI.DEFAULT_RULESET)),
        refused: [fails(""), fails("그냥 글"), fails('{"name":"x"'), fails('{"name":"x","motifs":["ZZZ"],"where":"cdr","minRisk":"높음","cys":"cdr"}'),
                  fails('{"name":"","motifs":["NG"],"where":"cdr","minRisk":"높음","cys":"cdr"}'), fails("[1,2]"), fails('{"a":1}')]}));""")
    assert out["oneLine"] and '"abcheck":"rule-set"' in out["line"] and "s123" not in out["line"]     # the id stays behind
    expected = {"name": "DS 포함", "motifs": ["NG", "D[STD]", "W"], "where": "cdr3", "minRisk": "중간", "cys": "off"}
    assert out["back"] == expected and out["inMail"] == expected and out["fromSavedJson"] == expected
    assert out["sameResults"]                                  # the colleague gets the same classification
    assert out["defaultSet"]["name"] == "기본" and out["defaultSet"]["motifs"] == ["N[^P][ST]", "NG", "DG", "DP"]
    assert all(out["refused"]), out["refused"]                 # every bad text is refused with a reason


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
def test_partner_chains_under_the_naming_schemes_people_use():
    out = _node("""
      const s = Object.fromEntries(entries);
      const VH = s["SYN-VH-001"], VL = s["SYN-VL-001"], VH2 = s["SYN-VH-003"], VL2 = s["SYN-VL-004"];
      // every scheme is tried with two antibodies at once, so pairing across antibodies would show
      const second = (name) => name.replace("mAb01", "mAb02").replace("1A2", "3B4").replace("clone 7", "clone 8");
      const tryScheme = ([h, l]) => {
        const rows = ABUI.buildRows(ABCHECK.analyze([[h, VH], [second(h), VH2], [l, VL], [second(l), VL2]], 7.4));
        const p = ABUI.pairMap(rows);
        if (p[0] === 2 && p[2] === 0 && p[1] === 3 && p[3] === 1) return ABUI.pairName(h) + " + " + ABUI.pairName(second(l));
        return Object.keys(p).length ? "WRONG " + JSON.stringify(p) : null;
      };
      const paired = [["mAb01_VH", "mAb01_VL"], ["mAb01-HC", "mAb01-LC"], ["mAb01 heavy", "mAb01 light"], ["mAb01.H", "mAb01.L"],
        ["mAb01_VH", "mAb01_VK"], ["mAb01_VH", "mAb01_Vκ"], ["mAb01_VH", "mAb01_Vλ"], ["mAb01_VH", "mAb01_VLambda"],
        ["mAb01H", "mAb01L"], ["1A2H", "1A2L"], ["1A2_h", "1A2_l"], ["mAb01VH", "mAb01VL"], ["mAb01HC", "mAb01LC"],
        ["mAb01_H_chain", "mAb01_L_chain"], ["mAb01_heavy_chain", "mAb01_light_chain"], ["mAb01|VH", "mAb01|VL"],
        ["mAb01/VH", "mAb01/VL"], ["mAb01:VH", "mAb01:VL"], ["VH_mAb01", "VL_mAb01"], ["HC mAb01", "LC mAb01"],
        ["mAb01(VH)", "mAb01(VL)"], ["mAb01 (VH)", "mAb01 (VL)"], ["mAb01 중쇄", "mAb01 경쇄"], ["mAb01_IGH", "mAb01_IGK"],
        ["mAb01_IGH", "mAb01_IGL"], ["clone 7 VH", "clone 7 VL"], ["IGH_mAb01_VH", "IGK_mAb01_VL"]];
      const unpaired = [["mAb01_VH1", "mAb01_VL1"], ["mAb01_HV", "mAb01_LV"], ["mAb01_IGHV3", "mAb01_IGKV1"],
        ["mAb01-CTRH", "mAb01-CTRL"], ["mAb01aH", "mAb01aL"], ["mAb01", "mAb01"]];
      // three chains under one name are ambiguous and stay apart; a mark alone is not a name
      const three = ABUI.pairMap(ABUI.buildRows(ABCHECK.analyze([["x_VH", VH], ["x_VL", VL], ["x_VK", VL2]], 7.4)));
      console.log(JSON.stringify({
        paired: paired.map(tryScheme), unpaired: unpaired.map(tryScheme), three,
        names: ["mAb01 (clone 7) VH", "SYN-VL-003", "mAb01H", "CTRL", "VH", "H_chain", "mAb01_heavy_chain"].map(ABUI.pairName),
        keys: ["mAb01|VH", "1A2H", "mAb01 중쇄", "CTRL", "mAb01_VH1"].map(ABUI.pairKey)}));""")
    assert all(name is not None and not name.startswith("WRONG") for name in out["paired"]), out["paired"]
    assert out["paired"][0] == "mAb01 + mAb02" and out["paired"][9] == "1A2 + 3B4" and out["paired"][25] == "clone 7 + clone 8"
    assert out["unpaired"] == [None] * 6
    assert out["three"] == {}
    assert out["names"] == ["mAb01 clone 7", "SYN-003", "mAb01", None, None, None, "mAb01"]
    assert out["keys"] == ["mab01", "1a2", "mab01", None, None]


@needs_node
def test_default_rule_set_reproduces_the_engine_triage_on_every_parity_case(tmp_path):
    """The same comparison as above over everything scripts/js_parity.py analyses: corpus and all edge cases, three pH values."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("js_parity", ROOT / "scripts" / "js_parity.py")
    parity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parity)
    cases = [{"label": label, "entries": [list(e) for e in entries], "ph": ph} for label, entries, ph in parity.standard_cases()]
    path = tmp_path / "cases.json"
    path.write_text(json.dumps(cases, ensure_ascii=False), encoding="utf-8")
    out = _node("const cases = JSON.parse(fs.readFileSync(" + json.dumps(str(path)) + """, "utf8"));
      let n = 0; const diffs = [], seen = {};
      for (const c of cases) {
        const p = ABCHECK.analyze(c.entries, c.ph);
        const rows = ABUI.buildRows(p, ABUI.DEFAULT_RULESET);
        p.sequences.forEach((s, i) => {
          n++;
          const t = ABUI.applyRules(s, ABUI.DEFAULT_RULESET), e = s.triage;
          seen[e.status] = (seen[e.status] || 0) + 1;
          if (t.status !== e.status || JSON.stringify(t.reasons) !== JSON.stringify(e.reasons) ||
              t.evaluable !== e.evaluable || t.note !== e.note) diffs.push([c.label, s.id, t, e]);
          if (s.input_ok && rows[i].status !== e.status) diffs.push([c.label, s.id, "row", rows[i].status, e.status]);
        });
      }
      console.log(JSON.stringify({n, diffs: diffs.slice(0, 5), statuses: Object.keys(seen).length}));""")
    assert out["n"] == sum(len(c["entries"]) for c in cases) and out["n"] > 500
    assert out["statuses"] == 3                      # the cases reach all three outcomes
    assert out["diffs"] == []


@needs_node
def test_stored_analyses_are_checked_when_read_back():
    out = _node("""
      const seq = entries[0][1];
      const good = {id: "r1", entries: [["a", seq], ["b", seq]], source: "붙여넣은 FASTA", ph: 5.5, text: ">a", at: 1764000000000, picked: [1], count: 2};
      const older = {id: "r0", entries: [["old", seq]], source: "예전 분석", ph: 7.4, at: 1763000000000};     // no picked, no text
      const c = ABUI.cleanRecent;
      console.log(JSON.stringify({
        nothing: [c(null), c(undefined)],
        notAList: [c({}), c("x"), c(42)],
        good: c([good]), older: c([older]),
        junk: c([null, 42, {entries: "oops"}, {entries: []}, good]),
        rows: c([Object.assign({}, good, {entries: [["a"], "b", [1, 2]]})]),
        noPh: c([Object.assign({}, good, {ph: "x"}), Object.assign({}, good, {ph: 99})]),
        noTime: c([Object.assign({}, good, {at: "yesterday"})]),
        soft: c([Object.assign({}, good, {id: 7, source: 123, text: 5, picked: [0, 7, -1, "x", 1.5], count: "many"})])}));""")
    empty = {"list": [], "dropped": 0}
    assert out["nothing"] == [empty, empty]
    assert out["notAList"] == [{"list": [], "dropped": 1}] * 3
    assert out["good"]["dropped"] == 0 and out["good"]["list"][0]["picked"] == [1] and out["good"]["list"][0]["ph"] == 5.5
    assert out["older"]["list"][0] == {"id": "r0", "entries": out["older"]["list"][0]["entries"], "ph": 7.4,
                                       "at": 1763000000000, "count": 1, "source": "예전 분석", "text": "", "picked": [], "notes": {},
                                       "ref": {"H": None, "L": None}}
    assert out["junk"]["dropped"] == 4 and [r["id"] for r in out["junk"]["list"]] == ["r1"]
    # an analysis whose rows, pH or time cannot be trusted is left out, not restored with made-up values
    assert out["rows"] == {"list": [], "dropped": 1}
    assert out["noPh"] == {"list": [], "dropped": 2} and out["noTime"] == {"list": [], "dropped": 1}
    soft = out["soft"]["list"][0]
    assert out["soft"]["dropped"] == 0
    assert (soft["id"], soft["source"], soft["text"], soft["picked"], soft["count"]) == ("r1764000000000", "저장된 분석", "", [0], 2)


@needs_node
def test_stored_settings_and_rule_sets_are_checked_when_read_back():
    out = _node("""
      const views = ["new", "all", "review", "rules"];
      const set = (id, extra) => Object.assign({id, name: id, motifs: ["NG"], where: "all", minRisk: "낮음", cys: "off"}, extra || {});
      console.log(JSON.stringify({
        uiNothing: [ABUI.cleanUi(null, views), ABUI.cleanUi([1, 2], views), ABUI.cleanUi("x", views)],
        uiBad: ABUI.cleanUi({chain: "X", conf: "maybe", sort: "liab", view: "zzz", unitChoice: 5}, views),
        uiGood: ABUI.cleanUi({chain: "H", conf: "no", sort: "pi-desc", view: "rules", unitChoice: "seq"}, views),
        uiLegacy: [ABUI.cleanUi({unit: "ab"}, views).unitChoice, ABUI.cleanUi({unit: "seq"}, views).unitChoice],
        rulesNothing: [ABUI.cleanRules(null), ABUI.cleanRules([]), ABUI.cleanRules("x"), ABUI.cleanRules({})],
        rulesNotAList: ABUI.cleanRules({active: "x", sets: "nope"}),
        rulesGood: ABUI.cleanRules({active: "team", sets: [set("team"), set("other")]}),
        rulesMixed: ABUI.cleanRules({active: "gone", sets: [set("team"), set("old", {motifs: ["ZZZ"]}), set("team"), set("default")]})}));""")
    default_ui = {"chain": "all", "conf": "all", "sort": "input", "view": "all", "unitChoice": None}
    assert out["uiNothing"] == [default_ui] * 3 and out["uiBad"] == default_ui
    assert out["uiGood"] == {"chain": "H", "conf": "no", "sort": "pi-desc", "view": "rules", "unitChoice": "seq"}
    assert out["uiLegacy"] == ["ab", None]
    nothing = {"active": "default", "sets": [], "dropped": 0, "errors": []}
    assert out["rulesNothing"] == [nothing] * 4
    assert out["rulesNotAList"]["sets"] == [] and out["rulesNotAList"]["dropped"] == 1
    assert out["rulesGood"]["active"] == "team" and [s["id"] for s in out["rulesGood"]["sets"]] == ["team", "other"]
    mixed = out["rulesMixed"]
    # unknown motif, a repeated id and the built-in id are left out; the active set falls back to the built-in one
    assert [s["id"] for s in mixed["sets"]] == ["team"] and mixed["dropped"] == 3 and mixed["active"] == "default"
    assert len(mixed["errors"]) == 3


@needs_node
def test_list_starts_in_the_antibody_view_when_most_chains_are_paired():
    out = _node("""
      const s = Object.fromEntries(entries);
      const unit = (named) => { const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4)); return ABUI.preferredUnit(rows, ABUI.pairMap(rows)); };
      const vh = s["SYN-VH-001"], vl = s["SYN-VL-001"];
      console.log(JSON.stringify({
        allPaired: unit([["a_VH", vh], ["a_VL", vl]]),
        mostPaired: unit([["a_VH", vh], ["a_VL", vl], ["b_VH", vh], ["b_VL", vl], ["solo", vh]]),
        exactlyHalf: unit([["a_VH", vh], ["a_VL", vl], ["x", vh], ["y", vl]]),
        fewPaired: unit([["a_VH", vh], ["a_VL", vl], ["x", vh], ["y", vh], ["z", vl]]),
        nonePaired: unit([["x", vh], ["y", vl]]),
        // a rejected entry has no chain and does not count against the pairs
        withRejected: unit([["a_VH", vh], ["a_VL", vl], ["bad", "ACDEF"], ["bad2", "ACDEF"], ["bad3", "ACDEF"]]),
        empty: ABUI.preferredUnit([], {})}));""")
    assert out == {"allPaired": "ab", "mostPaired": "ab", "exactlyHalf": "seq", "fewPaired": "seq", "nonePaired": "seq",
                   "withRejected": "ab", "empty": "seq"}


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
def test_cdr_written_against_a_reference():
    out = _node("""
      const ref = {name: "CDR-H2", start: 50, end: 65, length: 16, seq: "LSTAARNGGQQQGFDL"};
      const fixed = Object.assign({}, ref, {seq: "LSTAARQGGQQQGFDL"});
      const longer = {name: "CDR-H2", start: 50, end: 66, length: 17, seq: "LSTAARNGGQQQGFDLA"};
      const hit = {start: 56, end: 57, label: "deamidation (NG)", risk: "높음"};
      const bad = (() => { try { ABUI.markedCdr(ref, [], "short"); return null; } catch (e) { return e.message; } })();
      console.log(JSON.stringify({
        same: ABUI.diffCdr(ref, ref), fixed: ABUI.diffCdr(ref, fixed), longer: ABUI.diffCdr(ref, longer),
        noRef: ABUI.diffCdr(null, fixed), noOther: ABUI.diffCdr(ref, undefined),
        // the variant that still carries NG keeps its underline on the dots; the fixed one has none
        dotsWithMotif: ABUI.markedCdr(ref, [hit], ABUI.diffCdr(ref, ref).text),
        fixedShown: ABUI.markedCdr(fixed, [], ABUI.diffCdr(ref, fixed).text),
        plain: ABUI.markedCdr(ref, [hit]), bad}));""")
    assert out["same"] == {"comparable": True, "text": "·" * 16, "count": 0}
    assert out["fixed"] == {"comparable": True, "text": "······Q·········", "count": 1}
    assert out["longer"] == {"comparable": False, "text": "LSTAARNGGQQQGFDLA", "count": None}     # no alignment is made
    assert out["noRef"] == {"comparable": False, "text": "LSTAARQGGQQQGFDL", "count": None}
    assert out["noOther"] == {"comparable": False, "text": "", "count": None}
    assert out["dotsWithMotif"].startswith("······<u ") and out["dotsWithMotif"].endswith(">··</u>········")
    assert out["fixedShown"] == "······Q·········"
    assert "LSTAAR<u " in out["plain"] and ">NG</u>" in out["plain"]            # unchanged without a display text
    assert out["bad"] and "as long as the CDR" in out["bad"]


@needs_node
def test_motif_table_counts_candidates_like_the_python_engine():
    """Every cell of 모티프 분포 is recounted here with the Python package, and each cell's list filter returns that many."""
    from collections import Counter

    from abcheck.annotate import annotate
    from abcheck.liabilities import find_liabilities

    corpus = json.loads((ROOT / "data" / "synthetic" / "sequences.json").read_text(encoding="utf-8"))["sequences"]
    cells, totals, flags = Counter(), Counter(), Counter()
    for rec in corpus:
        rep = find_liabilities(annotate(rec["sequence"]))
        for key in {(h.label, h.region) for h in rep.hits}:
            cells[key] += 1
        for label in {h.label for h in rep.hits}:
            totals[label] += 1
        for label in {f["label"] if isinstance(f, dict) else f.label for f in rep.domain_flags}:
            flags[label] += 1

    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries.concat([["bad", "ACDEFGHIK"]]), 7.4));
      const t = ABUI.motifTable(rows);
      const n = (motif) => ABUI.filterRows(rows, {view: "all", motif}).length;
      // two NG in one candidate (CDR-H2 and CDR-H3): one candidate in each cell, one in the total
      const s = Object.fromEntries(entries)["SYN-VH-003"];
      const twice = ABUI.motifTable(ABUI.buildRows(ABCHECK.analyze([["v", s.slice(0, 99) + "G" + s.slice(100)]], 7.4))).rows.find(r => r.label === "deamidation (NG)");
      console.log(JSON.stringify({
        evaluated: t.evaluated, regions: t.regions, labels: t.rows.map(r => r.label), table: t.rows, flags: t.flags,
        listed: t.rows.flatMap(r => Object.keys(r.cells).map(reg => [r.label, reg, n({label: r.label, region: reg})])),
        listedTotals: t.rows.map(r => [r.label, n({label: r.label, region: null})]),
        listedFlags: t.flags.map(f => [f.label, n({label: f.label, region: null, flag: true})]),
        twice, motifOrder: ABUI.MOTIFS.map(m => m.label), empty: ABUI.motifTable([])}));""")
    assert out["evaluated"] == len(corpus)                          # the input error is not counted
    assert out["labels"] == out["motifOrder"]                       # every motif of the table, found or not
    order = ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3", "framework", "구간 미확정"]
    assert out["regions"] == [r for r in order if any(k[1] == r for k in cells)]
    got = {(row["label"], reg): n for row in out["table"] for reg, n in row["cells"].items()}
    assert got == dict(cells)
    assert {row["label"]: row["total"] for row in out["table"] if row["total"]} == dict(totals)
    assert {f["label"]: f["total"] for f in out["flags"]} == dict(flags) and flags      # the corpus has an unpaired Cys case
    # the list filter of a cell holds exactly as many candidates as the cell says
    assert {(label, reg): n for label, reg, n in out["listed"]} == dict(cells)
    assert {label: n for label, n in out["listedTotals"] if n} == dict(totals)
    assert {label: n for label, n in out["listedFlags"]} == dict(flags)
    assert out["twice"] == {"label": "deamidation (NG)", "cells": {"CDR-H2": 1, "CDR-H3": 1}, "total": 1}
    assert out["empty"]["evaluated"] == 0 and out["empty"]["regions"] == [] and out["empty"]["flags"] == []


@needs_node
def test_taking_an_entry_out_keeps_picks_and_references_on_the_same_candidates():
    out = _node("""
      const list = [["a", "AAA"], ["b", "CCC"], ["c", "DDD"], ["d", "EEE"], ["e", "FFF"]];
      const names = (x) => ({entries: x.entries.map(e => e[0]), picked: x.picked.map(i => x.entries[i][0]), ref: Object.fromEntries(Object.entries(x.ref).map(([k, v]) => [k, v === null ? null : x.entries[v][0]]))});
      const refused = [-1, 5, 1.5, "1"].map(i => { try { ABUI.withoutEntry(list, i, [], {}); return null; } catch (e) { return e.message; } });
      const untouched = JSON.stringify(list);
      console.log(JSON.stringify({
        middle: names(ABUI.withoutEntry(list, 1, [0, 1, 3, 4], {H: 3, L: 0})),
        ofPickedRef: names(ABUI.withoutEntry(list, 3, [3], {H: 3, L: null})),
        first: names(ABUI.withoutEntry(list, 0, [4], {H: 1, L: 4})),
        last: names(ABUI.withoutEntry(list, 4, [0, 4], {H: 4, L: 2})),
        bare: ABUI.withoutEntry(list, 2).picked, inputUnchanged: JSON.stringify(list) === untouched, refused}));""")
    assert out["middle"] == {"entries": ["a", "c", "d", "e"], "picked": ["a", "d", "e"], "ref": {"H": "d", "L": "a"}}
    assert out["ofPickedRef"] == {"entries": ["a", "b", "c", "e"], "picked": [], "ref": {"H": None, "L": None}}
    assert out["first"] == {"entries": ["b", "c", "d", "e"], "picked": ["e"], "ref": {"H": "b", "L": "e"}}
    assert out["last"] == {"entries": ["a", "b", "c", "d"], "picked": ["a"], "ref": {"H": None, "L": "c"}}
    assert out["bare"] == [] and out["inputUnchanged"] is True
    assert all(m and "withoutEntry: no entry" in m for m in out["refused"])


@needs_node
def test_notes_are_one_line_and_checked_when_read_back():
    out = _node("""
      const rec = (notes) => ({id: "r1", entries: [["a", "AAA"], ["b", "CCC"], ["c", "DDD"]], ph: 7.4, at: 1, notes});
      const back = (notes) => ABUI.cleanRecent([rec(notes)]).list[0].notes;
      console.log(JSON.stringify({
        clean: ["  발현량  높음\\n다음 라운드 ", "", null, undefined, 7, "가".repeat(250)].map(ABUI.cleanNote).map(x => x.length > 40 ? x.length : x),
        kept: back({0: "대표 후보", 2: " H3 NG \\t 제거 필요 "}),
        bad: back({"-1": "x", 3: "밖", 1.5: "x", a: "x", 1: 5, 2: "   ", 0: "남음"}),
        notObject: [back(["x"]), back("x"), back(null), back(undefined)],
        moved: ABUI.withoutEntry([["a", "A"], ["b", "C"], ["c", "D"], ["d", "E"]], 1, [], {}, {0: "첫째", 1: "빠짐", 3: "넷째"}).notes,
        noNotes: ABUI.withoutEntry([["a", "A"], ["b", "C"]], 0, [], {}).notes}));""")
    assert out["clean"] == ["발현량 높음 다음 라운드", "", "", "", "7", 200]
    assert out["kept"] == {"0": "대표 후보", "2": "H3 NG 제거 필요"}
    assert out["bad"] == {"0": "남음"}                         # outside the list, not a place, not text, blank: left out
    assert out["notObject"] == [{}, {}, {}, {}]
    assert out["moved"] == {"0": "첫째", "2": "넷째"}             # the note of the removed entry goes, the later one moves up
    assert out["noNotes"] == {}


@needs_node
def test_picks_and_notes_travel_with_the_exports():
    out = _node("""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries.concat([["bad", "ACDEFGHIK"]]), 7.4));
      const by = Object.fromEntries(rows.map(r => [r.id, r.index]));
      const meta = {ruleName: "기본", when: "2026-11-26 10:00"};
      const notes = {[by["SYN-VH-003"]]: "N56Q 변이체 제작", [by["SYN-VH-001"]]: '대표 <b>후보</b>, "1순위"', [by["bad"]]: "서열 다시 받기"};
      const picked = new Set([by["SYN-VH-001"], by["SYN-VL-001"]]);
      const plain = ABUI.toMailText(rows, meta), marked = ABUI.toMailText(rows, Object.assign({picked, notes}, meta));
      const body = plain.split("\\n\\nKabat")[0];
      const some = rows.slice(0, 5);
      const csv = ABUI.toCsv(some, "기본", notes).split("\\r\\n"), tsv = ABUI.toTsv(some, "기본", notes).split("\\n");
      const page = ABUI.onePageFragment(some, Object.assign({notes}, meta));
      console.log(JSON.stringify({
        ids: some.map(r => r.id),
        plainHasMarks: /고른 후보|메모/.test(plain), startsSame: marked.startsWith(body), tail: marked.slice(body.length),
        emptyMarks: ABUI.toMailText(rows, Object.assign({picked: new Set(), notes: {}}, meta)) === plain,
        csvSame: ABUI.toCsv(some, "기본", {}) === ABUI.toCsv(some, "기본") && ABUI.toCsv(some, "기본", {[by["bad"]]: "목록에 없음"}) === ABUI.toCsv(some, "기본"),
        csvHead: csv[0].split(",").slice(-2), csvRow1: csv[1].slice(-28), csvRow2: csv[2].split(",").slice(-1),
        tsvHead: tsv[0].split("\\t").slice(-1), tsvRow1: tsv[1].split("\\t").slice(-1), tsvRow5: tsv[5].split("\\t").slice(-1), tsvWidths: [...new Set(tsv.map(l => l.split("\\t").length))],
        pageNote: (page.match(/<tr class="op-note">.*?<\\/tr>/g) || []), pagePlain: ABUI.onePageFragment(some, meta).includes("op-note")}));""")
    assert out["ids"][0] == "SYN-VH-001" and out["ids"][4] == "SYN-VH-003"
    assert out["plainHasMarks"] is False and out["emptyMarks"] is True       # nothing marked: the text is what it was
    assert out["startsSame"] is True                                         # the other sections are untouched
    assert out["tail"] == (
        '\n\n고른 후보 2\n- SYN-VH-001 (VH): 대표 <b>후보</b>, "1순위"\n- SYN-VL-001 (VL)'
        "\n\n메모 2\n- SYN-VH-003 (VH): N56Q 변이체 제작\n- bad: 서열 다시 받기"
        "\n\nKabat 근사 · ANARCI/IMGT 아님 · 종합 판정 없음")
    assert out["csvSame"] is True                                            # no note on a listed row: the file is what it was
    assert out["csvHead"] == ["입력 오류", "메모"]
    assert out["csvRow1"].endswith(',"대표 <b>후보</b>, ""1순위"""')          # quoted for Excel
    assert out["csvRow2"] == [""]
    assert out["tsvHead"] == ["메모"] and out["tsvRow1"] == ['대표 <b>후보</b>, "1순위"'] and out["tsvRow5"] == ["N56Q 변이체 제작"]
    assert len(out["tsvWidths"]) == 1                                        # every line has the same number of columns
    assert out["pageNote"] == ['<tr class="op-note"><td colspan="7"><small>메모</small> 대표 &lt;b&gt;후보&lt;/b&gt;, &quot;1순위&quot;</td></tr>',
                               '<tr class="op-note"><td colspan="7"><small>메모</small> N56Q 변이체 제작</td></tr>']
    assert out["pagePlain"] is False


@needs_node
def test_substitutions_are_checked_against_the_sequence():
    out = _node("""
      const seq = Object.fromEntries(entries)["SYN-VH-003"];      // N56 G57 in CDR-H2, W36 and C22 are anchors
      const p = (t) => ABUI.parseSubstitutions(t, seq);
      const two = p("g57a, n56q");
      const refused = (() => { try { ABUI.applySubstitutions(seq, [{pos: 56, from: "G", to: "A"}]); return null; } catch (e) { return e.message; } })();
      const made = ABUI.applySubstitutions(seq, two.subs);
      console.log(JSON.stringify({
        at: seq.slice(54, 58), len: seq.length, one: p("N56Q"), two, separators: p("N56Q/G57A + W36F").subs.length,
        made: made.slice(54, 58), sameElsewhere: made.slice(0, 55) + made.slice(57) === seq.slice(0, 55) + seq.slice(57),
        name: ABUI.variantName("mAb03_VH", two.subs), refused,
        errors: ["N56", "56Q", "N57Q", "N560Q", "N56N", "N56X", "N56Q N56S", "N56Q G57", ""].map(t => p(t).errors),
        partlyGood: p("N56Q A58G").subs}));""")
    assert out["at"] == "RNGG"
    assert out["one"] == {"subs": [{"pos": 56, "from": "N", "to": "Q"}], "errors": []}
    # any case, any order: sorted by place
    assert out["two"] == {"subs": [{"pos": 56, "from": "N", "to": "Q"}, {"pos": 57, "from": "G", "to": "A"}], "errors": []}
    assert out["separators"] == 3
    assert out["made"] == "RQAG" and out["sameElsewhere"] is True
    assert out["name"] == "mAb03_VH_N56Q_G57A"
    assert "place 56 holds N, not G" in out["refused"]
    e = out["errors"]
    assert e[0] == ["N56 뒤에 바꿀 잔기 한 글자 필요"]
    assert e[1] == ["56Q: 치환 표기 아님 (예 N56Q)"]
    assert e[2] == ["N57Q: 57번은 G"]                        # an off-by-one is refused, not corrected
    assert e[3] == [f"N560Q: 서열은 1–{out['len']}번"]
    assert e[4] == ["N56N: 같은 잔기"]
    assert e[5] == ["N56X: X 는 표준 아미노산 아님"]
    assert e[6] == ["N56S: 56번을 두 번 적음"]
    assert e[7] == ["G57 뒤에 바꿀 잔기 한 글자 필요"]
    assert e[8] == []
    assert out["partlyGood"] == [{"pos": 56, "from": "N", "to": "Q"}]


@needs_node
def test_variant_hits_against_the_parent():
    out = _node("""
      const seq = Object.fromEntries(entries)["SYN-VH-003"];
      const row = (name, s) => ABUI.buildRows(ABCHECK.analyze([[name, s]], 7.4))[0];
      const parent = row("p", seq);
      const variant = (text) => { const v = row("v", ABUI.applySubstitutions(seq, ABUI.parseSubstitutions(text, seq).subs)); return [v, ABUI.diffHits(parent.hits, v.hits)]; };
      const show = (h) => h.label + " " + h.start + "-" + h.end + " " + h.region + " " + h.risk;
      const brief = ([v, d]) => ({status: v.status, confident: v.confident, removed: d.removed.map(show), added: d.added.map(show),
        moved: d.moved.map(m => show(m.before) + " -> " + show(m.after)), kept: d.kept, total: v.hits.length});
      console.log(JSON.stringify({
        parent: {status: parent.status, hits: parent.hits.length, k100: seq[99], w36: seq[35]},
        fixed: brief(variant("N56Q")), made: brief(variant("K100G")), anchor: brief(variant("W36A")),
        regions: [31, 36, 56, 100, 120].map(p => ABUI.regionAt(parent, p)), noCdrs: ABUI.regionAt(variant("W36A")[0], 56),
        same: ABUI.diffHits(parent.hits, parent.hits)}));""")
    assert out["parent"] == {"status": "검토 필요", "hits": 3, "k100": "K", "w36": "W"}
    # the substitution that removes the motif: one hit gone, nothing new, the classification follows
    assert out["fixed"] == {"status": "기준 해당 없음", "confident": True, "removed": ["deamidation (NG) 56-57 CDR-H2 높음"],
                            "added": [], "moved": [], "kept": 2, "total": 2}
    # a substitution that makes a new NG next to N99: reported as added
    assert out["made"]["added"] == ["deamidation (NG) 99-100 CDR-H3 높음"] and out["made"]["removed"] == []
    assert out["made"]["kept"] == 3
    # an anchor residue replaced: the Trp hit is gone and whatever the annotation now says is reported, not hidden
    anchor = out["anchor"]
    assert "산화 (Trp) 36-36 framework 낮음" in anchor["removed"]
    assert anchor["status"] == "판정 불가" and anchor["confident"] is False and out["noCdrs"] == "구간 미확정"
    assert anchor["moved"] == ["deamidation (NG) 56-57 CDR-H2 높음 -> deamidation (NG) 56-57 구간 미확정 높음",
                               "산화 (Trp) 108-108 framework 낮음 -> 산화 (Trp) 108-108 구간 미확정 낮음"]
    assert anchor["kept"] + len(anchor["moved"]) + len(anchor["added"]) == anchor["total"]
    assert out["regions"] == ["CDR-H1", "framework", "CDR-H2", "CDR-H3", "framework"]
    assert out["same"] == {"removed": [], "added": [], "moved": [], "kept": 3}


@needs_node
def test_a_sequence_is_cut_at_its_cdrs():
    out = _node(EDGE + r"""
      const rows = ABUI.buildRows(ABCHECK.analyze(all, 7.4));
      const cut = rows.map(r => ABUI.regionsOf(r));
      const ok = rows.filter((r, i) => cut[i]);
      console.log(JSON.stringify({
        cutRows: ok.length, confidentWithThree: rows.filter(r => r.inputOk && r.confident && r.cdrs.length === 3).length,
        notCut: rows.filter((r, i) => !cut[i]).map(r => [r.inputOk, r.confident, r.cdrs.length]),
        names: [...new Set(ok.map(r => ABUI.regionsOf(r).map(p => p.name).join(" ")))],
        joined: ok.every(r => ABUI.regionsOf(r).map(p => p.seq).join("") === r.sequence),
        contiguous: ok.every(r => { const ps = ABUI.regionsOf(r); return ps[0].start === 1 && ps[6].end === r.sequence.length && ps.every((p, i) => i === 0 || p.start === ps[i - 1].end + 1); }),
        cdrPiecesAreTheEngines: ok.every(r => r.cdrs.every(c => ABUI.regionsOf(r).some(p => p.name === c.name && p.start === c.start && p.end === c.end && p.seq === c.seq))),
        first: ABUI.regionsOf(rows.find(r => r.id === "SYN-VH-003")).map(p => p.name + ":" + p.start + "-" + p.end)}));""")
    assert out["cutRows"] == out["confidentWithThree"] > 10
    assert all(not (ok and conf and n == 3) for ok, conf, n in out["notCut"])      # the rest: input error, not confident, or no CDRs
    assert sorted(out["names"]) == ["FR1 CDR-H1 FR2 CDR-H2 FR3 CDR-H3 FR4", "FR1 CDR-L1 FR2 CDR-L2 FR3 CDR-L3 FR4"]
    assert out["joined"] and out["contiguous"] and out["cdrPiecesAreTheEngines"]
    assert out["first"] == ["FR1:1-30", "CDR-H1:31-35", "FR2:36-49", "CDR-H2:50-65", "FR3:66-97", "CDR-H3:98-107", "FR4:108-118"]


@needs_node
def test_a_candidate_against_the_reference_piece_by_piece():
    out = _node(r"""
      const s = Object.fromEntries(entries);
      const vh = s["SYN-VH-003"];                                  // FR1 1-30, CDR-H1 31-35, CDR-H2 50-65, CDR-H3 98-107, 118 aa
      const sub = (seq, text) => ABUI.applySubstitutions(seq, ABUI.parseSubstitutions(text, seq).subs);
      const typed = "F10Y N56Q H80R K100G S115T";
      const named = [
        ["ref", vh], ["same", vh], ["point", sub(vh, typed)],
        ["longer3", sub(vh.slice(0, 100) + "A" + vh.slice(100), "F10Y N56Q")],      // CDR-H3 one residue longer
        ["shortFr1", sub(vh.slice(1), "N55Q")],                      // FR1 one residue shorter: every later place moves by one
        ["noAnchor", sub(vh, "W36A N56Q")],                          // annotation fails, same length
        ["noAnchorLonger", sub(vh, "W36A") + "A"], ["light", s["SYN-VL-003"]], ["bad", "ACDEFGHIK"]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const r = Object.fromEntries(rows.map(x => [x.id, x]));
      const show = (d) => ({mode: d.mode, count: d.count, samePlaces: d.samePlaces,
        subs: d.subs.map(x => x.from + x.pos + x.to + "@" + x.region).join(" "),
        regions: d.regions.map(g => g.name + ":" + (g.comparable ? g.count : g.refLength + "!=" + g.length)).join(" ")});
      console.log(JSON.stringify({
        conf: rows.map(x => x.id + ":" + x.confident + ":" + x.cdrs.length),
        same: show(ABUI.diffRegions(r.ref, r.same)), point: show(ABUI.diffRegions(r.ref, r.point)),
        back: show(ABUI.diffRegions(r.point, r.ref)),
        longer3: show(ABUI.diffRegions(r.ref, r.longer3)), shortFr1: show(ABUI.diffRegions(r.ref, r.shortFr1)),
        noAnchor: show(ABUI.diffRegions(r.ref, r.noAnchor)), noAnchorLonger: show(ABUI.diffRegions(r.ref, r.noAnchorLonger)),
        bad: show(ABUI.diffRegions(r.ref, r.bad)), nothing: show(ABUI.diffRegions(null, r.ref)),
        // what was typed is what is found
        typedBack: ABUI.diffRegions(r.ref, r.point).subs.map(x => x.from + x.pos + x.to).join(" ") === typed,
        // the row is the reference with these substitutions, wherever places correspond
        rebuilds: ["point", "noAnchor"].every(id => { const d = ABUI.diffRegions(r.ref, r[id]); return d.samePlaces && ABUI.applySubstitutions(r.ref.sequence, d.subs) === r[id].sequence; })}));""")
    assert out["conf"][:5] == ["ref:true:3", "same:true:3", "point:true:3", "longer3:true:3", "shortFr1:true:3"]
    assert out["conf"][5] == "noAnchor:false:0"
    assert out["same"] == {"mode": "regions", "count": 0, "samePlaces": True, "subs": "",
                           "regions": "FR1:0 CDR-H1:0 FR2:0 CDR-H2:0 FR3:0 CDR-H3:0 FR4:0"}
    assert out["point"] == {"mode": "regions", "count": 5, "samePlaces": True,
                            "subs": "F10Y@FR1 N56Q@CDR-H2 H80R@FR3 K100G@CDR-H3 S115T@FR4",
                            "regions": "FR1:1 CDR-H1:0 FR2:0 CDR-H2:1 FR3:1 CDR-H3:1 FR4:1"}
    assert out["back"]["subs"] == "Y10F@FR1 Q56N@CDR-H2 R80H@FR3 G100K@CDR-H3 T115S@FR4"
    # a longer CDR-H3 is left out, the pieces around it are still compared; places after it no longer correspond
    assert out["longer3"] == {"mode": "regions", "count": 2, "samePlaces": False, "subs": "F10Y@FR1 N56Q@CDR-H2",
                              "regions": "FR1:1 CDR-H1:0 FR2:0 CDR-H2:1 FR3:0 CDR-H3:10!=11 FR4:0"}
    # a shorter FR1 is left out; N56Q of the reference is at place 55 of this candidate
    assert out["shortFr1"] == {"mode": "regions", "count": 1, "samePlaces": False, "subs": "N55Q@CDR-H2",
                               "regions": "FR1:30!=29 CDR-H1:0 FR2:0 CDR-H2:1 FR3:0 CDR-H3:0 FR4:0"}
    # no CDR borders on one side: same length, so place by place, regions in the engine's words for that row
    assert out["noAnchor"] == {"mode": "whole", "count": 2, "samePlaces": True, "subs": "W36A@구간 미확정 N56Q@구간 미확정",
                               "regions": "전체 서열:2"}
    none = {"mode": "none", "count": None, "samePlaces": False, "subs": "", "regions": ""}
    assert out["noAnchorLonger"] == none and out["bad"] == none and out["nothing"] == none
    assert out["typedBack"] and out["rebuilds"]


@needs_node
def test_differences_from_the_reference_travel_with_the_exports():
    out = _node(r"""
      const s = Object.fromEntries(entries);
      const vh = s["SYN-VH-003"];
      const sub = (seq, text) => ABUI.applySubstitutions(seq, ABUI.parseSubstitutions(text, seq).subs);
      const named = [
        ["parent", vh], ["point", sub(vh, "F10Y N56Q H80R")], ["twin", vh],
        ["longer3", sub(vh.slice(0, 100) + "A" + vh.slice(100), "F10Y")], ["noAnchorLonger", sub(vh, "W36A") + "A"],
        ["light", s["SYN-VL-003"]], ["bad", "ACDEFGHIK"], ["<i>", sub(vh, "K100G")]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const by = Object.fromEntries(rows.map(r => [r.id, r]));
      const refs = ABUI.refNotes(rows, {H: by.parent.index, L: null});
      const meta = {ruleName: "기본", when: "2026-11-26 10:00"};
      const notes = {[by.point.index]: "발현량 높음"};
      const picked = new Set([by.parent.index, by.point.index, by.longer3.index]);
      const plainTsv = ABUI.toTsv(rows, "기본"), tsv = ABUI.toTsv(rows, "기본", notes, refs).split("\n").map(l => l.split("\t"));
      const csv = ABUI.toCsv(rows, "기본", {}, refs).split("\r\n");
      const mail = ABUI.toMailText(rows, Object.assign({picked, notes, refs}, meta));
      const page = ABUI.onePageFragment(rows, Object.assign({notes, refs}, meta));
      let refused = null; try { ABUI.refLine(by.parent, by.point, 0); } catch (e) { refused = e.message; }
      console.log(JSON.stringify({
        lines: Object.fromEntries(Object.keys(refs).map(k => [rows[k].id, refs[k].id + "|" + refs[k].text])),
        capped: ABUI.refLine(by.parent, by.point, 2), refused,
        noRef: [ABUI.refNotes(rows, {H: null, L: null}), ABUI.refNotes(rows, null), ABUI.refNotes(rows, {H: 99, L: null}), ABUI.refNotes(rows, {H: by.bad.index, L: null})],
        lightRef: Object.keys(ABUI.refNotes(rows, {H: null, L: by.light.index})),
        unchanged: ABUI.toTsv(rows, "기본", {}, {}) === plainTsv && ABUI.toTsv(rows, "기본", undefined, ABUI.refNotes(rows, {H: null, L: null})) === plainTsv,
        tsvHead: tsv[0].slice(-2), tsvCells: tsv.slice(1).map(l => l.slice(-2).join(" / ")), widths: [...new Set(tsv.map(l => l.length))],
        csvHead: csv[0].split(",").pop(), csvPoint: csv[2].split(",").pop(),
        mailTail: mail.slice(mail.indexOf("고른 후보")),
        pageRows: (page.match(/<tr class="op-note(?: op-noted)?">.*?<\/tr>/g) || []),
        pagePlain: ABUI.onePageFragment(rows, meta).includes("op-note")}));""")
    assert out["lines"] == {
        "point": "parent|F10Y N56Q H80R", "twin": "parent|차이 없음", "longer3": "parent|F10Y · CDR-H3 길이 다름",
        "noAnchorLonger": "parent|길이 다름 · 견주지 않음", "<i>": "parent|K100G"}       # no line: the reference, the light chain, the input error
    assert out["capped"] == "F10Y N56Q 외 1" and "limit must be a whole number" in out["refused"]
    assert out["noRef"] == [{}, {}, {}, {}] and out["lightRef"] == []
    assert out["unchanged"]                                                   # no reference: the table is what it was
    assert out["tsvHead"] == ["기준 대비", "메모"] and len(out["widths"]) == 1
    assert out["tsvCells"] == [" / ", "parent: F10Y N56Q H80R / 발현량 높음", "parent: 차이 없음 / ", "parent: F10Y · CDR-H3 길이 다름 / ",
                               "parent: 길이 다름 · 견주지 않음 / ", " / ", " / ", "parent: K100G / "]
    assert out["csvHead"] == "기준 대비" and out["csvPoint"] == "parent: F10Y N56Q H80R"
    assert out["mailTail"] == (
        "고른 후보 3\n- parent (VH)\n- point (VH): 발현량 높음\n  기준 parent 대비: F10Y N56Q H80R"
        "\n- longer3 (VH)\n  기준 parent 대비: F10Y · CDR-H3 길이 다름"
        "\n\nKabat 근사 · ANARCI/IMGT 아님 · 종합 판정 없음")
    assert out["pageRows"][0] == '<tr class="op-note op-noted"><td colspan="7"><small>기준 parent 대비</small> <span class="op-subs">F10Y N56Q H80R</span></td></tr>'
    assert out["pageRows"][1] == '<tr class="op-note"><td colspan="7"><small>메모</small> 발현량 높음</td></tr>'
    assert out["pageRows"][2] == '<tr class="op-note"><td colspan="7"><small>기준 parent 대비</small> 차이 없음</td></tr>'
    assert out["pageRows"][3] == '<tr class="op-note"><td colspan="7"><small>기준 parent 대비</small> <span class="op-subs">F10Y</span> · CDR-H3 길이 다름</td></tr>'
    assert len(out["pageRows"]) == 6 and out["pagePlain"] is False


@needs_node
def test_one_record_per_motif_as_the_python_engine_finds_them():
    """모티프 분포 › Excel용 복사: every record is recounted here with the Python package, in the same order."""
    from abcheck.annotate import annotate
    from abcheck.liabilities import find_liabilities

    corpus = json.loads((ROOT / "data" / "synthetic" / "sequences.json").read_text(encoding="utf-8"))["sequences"]
    expected = []
    for rec in corpus:
        rep = find_liabilities(annotate(rec["sequence"]))
        expected += [[rec["id"], h.label, h.match, h.start, h.end, h.region, h.risk, h.base_risk] for h in rep.hits]
        expected += [[rec["id"], f["label"], "", "", "", "도메인 전체", f["risk"], ""] for f in rep.domain_flags]

    out = _node(r"""
      const rows = ABUI.buildRows(ABCHECK.analyze(entries.concat([["bad", "ACDEFGHIK"], ["tab\tname", entries[0][1]]]), 7.4));
      const t = ABUI.hitsTable(rows), tsv = ABUI.toHitsTsv(rows).split("\n").map(l => l.split("\t"));
      console.log(JSON.stringify({
        head: t[0], records: t.slice(1), statuses: Object.fromEntries(rows.map(r => [r.id, r.status])),
        chains: Object.fromEntries(rows.map(r => [r.id, r.chain === "H" ? "VH" : r.chain === "L" ? "VL" : ""])),
        tsvLines: tsv.length, widths: [...new Set(tsv.map(l => l.length))], tsvHead: tsv[0],
        tsvSame: tsv.slice(1).every((l, i) => l.join("|") === t[i + 1].map(v => String(v).replace(/\t/g, " ")).join("|")),
        tabName: tsv.filter(l => l[0] === "tab name").length, empty: ABUI.hitsTable([]), onlyBad: ABUI.hitsTable(rows.filter(r => !r.inputOk))}));""")
    head = ["서열 이름", "체인", "분류", "항목", "모티프", "시작", "끝", "구간", "위험도", "기본 위험도"]
    assert out["head"] == head and out["tsvHead"] == head
    own = [r for r in out["records"] if r[0] != "tab\tname"]          # the extra entry repeats the first sequence under another name
    assert [[r[0]] + r[3:] for r in own] == expected and len(expected) > 40
    assert all(r[1] == out["chains"][r[0]] and r[2] == out["statuses"][r[0]] for r in out["records"])
    assert not any(r[0] == "bad" for r in out["records"])             # an input error has no record
    assert out["tsvLines"] == len(out["records"]) + 1 and out["widths"] == [10] and out["tsvSame"]
    assert out["tabName"] == sum(1 for r in expected if r[0] == corpus[0]["id"]) > 0      # a tab in a name does not shift the columns
    assert out["empty"] == [head] and out["onlyBad"] == [head]


@needs_node
def test_similar_cdr3_groups_by_length_and_a_few_residues():
    out = _node("""
      const s = Object.fromEntries(entries);
      const vh = s["SYN-VH-003"], vl = s["SYN-VL-003"];           // CDR-H3 YNKQQQFYLK (98-107), CDR-L3 LRANGEGAD (89-97)
      const sub = (seq, text) => ABUI.applySubstitutions(seq, ABUI.parseSubstitutions(text, seq).subs);
      const named = [
        ["A", vh], ["other", s["SYN-VH-001"]], ["B", sub(vh, "K100G")], ["C", sub(vh, "K100G Q101E")], ["D", sub(vh, "K100G Q101E Q102E Q103E")],
        ["same", sub(vh, "N56Q")],                                 // CDR-H2 changed, CDR-H3 identical to A
        ["longer", vh.slice(0, 107) + "A" + vh.slice(107)],
        ["La", vl], ["Lb", sub(vl, "A91S")], ["bad", "ACDEFGHIK"]];
      const rows = ABUI.buildRows(ABCHECK.analyze(named, 7.4));
      const show = (k) => ABUI.similarCdr3(rows, k).map(g => [g.name, g.length, g.members.map(m => m.row.id + ":" + m.diff).join(" ")]);
      const refused = [-1, 1.5, "2", undefined].map(k => { try { ABUI.similarCdr3(rows, k); return null; } catch (e) { return e.message; } });
      console.log(JSON.stringify({cdr3: rows.filter(r => r.cdr3).map(r => r.id + "=" + r.cdr3.seq), k0: show(0), k1: show(1), k2: show(2), k3: show(3), none: ABUI.similarCdr3([], 2), refused}));""")
    assert "A=YNKQQQFYLK" in out["cdr3"] and "D=YNGEEEFYLK" in out["cdr3"] and "Lb=LRSNGEGAD" in out["cdr3"]
    assert out["k0"] == [["CDR-H3", 10, "A:0 same:0"]]                                  # identical only
    # one residue: A-B, B-C link into one group although A and C differ in two (each member shows its own count)
    assert out["k1"] == [["CDR-H3", 10, "A:0 B:1 C:2 same:0"], ["CDR-L3", 9, "La:0 Lb:1"]]
    assert out["k2"] == [["CDR-H3", 10, "A:0 B:1 C:2 D:4 same:0"], ["CDR-L3", 9, "La:0 Lb:1"]]    # D joins through C (2 apart), 4 from A
    assert out["k3"] == out["k2"]
    assert all("longer" not in g[2] and "other" not in g[2] for g in out["k3"])        # another length is never grouped; no alignment
    assert out["none"] == []
    assert all(m and "whole number of residues" in m for m in out["refused"])


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
