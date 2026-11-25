/*
 * app-lib.js -- pure helpers behind the AbCheck browser app (no DOM, no dependencies).
 *
 * The app (web/app.html) renders what these functions return; tests/test_web_app.py runs
 * them under Node against the real engine output. Defines the global ABUI, and
 * module.exports under Node.
 */
(function (root, factory) {
  "use strict";
  var api = factory();
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ABUI = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  var RISKS = ["높음", "중간", "낮음"];
  var CDR_REGIONS = ["CDR1", "CDR2", "CDR3"];

  /* The liability motifs as abcheck/liabilities.py MOTIF_SPECS defines them, in the same
     order. Shown on the rules screen; a test keeps this list equal to the Python one. */
  var MOTIFS = [
    { pattern: "N[^P][ST]", name: "N-결합 당화", base: "높음", top: true },
    { pattern: "NG", name: "deamidation", base: "높음", top: true },
    { pattern: "N[STNH]", name: "deamidation", base: "중간", top: false },
    { pattern: "DG", name: "isomerization", base: "높음", top: true },
    { pattern: "D[STD]", name: "isomerization", base: "중간", top: false },
    { pattern: "DP", name: "fragmentation", base: "중간", top: true },
    { pattern: "M", name: "산화 (Met)", base: "낮음", top: false },
    { pattern: "W", name: "산화 (Trp)", base: "낮음", top: false }
  ];

  /* Headline validation figures (README.md / VALIDATION.md). A test checks every one of
     these strings still appears in those documents. */
  var VALIDATION = {
    sequences: "7,000",
    datasets: "7",
    passedInput: "6,984",
    chainCorrect: "6,980",
    chainWrong: "0",
    annotation: "98.7%",
    confident: "96.0%",
    boundary: "0.35%",
    boundaryOf: "58/16,498",
    fr3Gap: "3,429 / 3,443",
    rows: [
      ["Jaffe 2022", "사람 · naive B", "100%", "99.9%", "1000/1000", "0.00%"],
      ["Alsoiussi 2020", "마우스 · plasmablast", "100%", "99.2%", "1000/1000", "0.00%"],
      ["King 2020", "사람 · memory B", "99.7%", "95.3%", "1000/1000", "0.46%"],
      ["Ota 2024", "사람 · memory B", "99.9%", "95.8%", "1000/1000", "0.63%"],
      ["Phad 2022", "사람 · plasma B", "99.9%", "93.9%", "999/1000", "0.84%"],
      ["Setliff 2019", "사람 · unsorted B (HIV)", "100%", "97.8%", "1000/1000", "0.25%"],
      ["Sokal 2021", "사람 · memory B (잘린 read 포함)", "91.6%", "89.8%", "981/984", "0.28%"]
    ]
  };

  function zero() { return { "높음": 0, "중간": 0, "낮음": 0 }; }

  function addCounts(into, c) {
    RISKS.forEach(function (k) { into[k] += c[k] || 0; });
    return into;
  }

  /** "높음 1 · 낮음 2", or "-" when nothing was found. */
  function countText(c) {
    var parts = RISKS.filter(function (k) { return c[k]; }).map(function (k) { return k + " " + c[k]; });
    return parts.length ? parts.join(" · ") : "-";
  }

  /** Which sidebar lists a sequence belongs to. A sequence can be in two (needs review AND
      not confident). Keys double as the app's view names. */
  function categories(row) {
    var c = [];
    if (!row.inputOk) c.push("errors");
    if (row.status === "검토 필요") c.push("review");
    if (row.inputOk && !row.confident) c.push("recheck");
    return c;
  }

  /** One flat record per input sequence, joining the engine's per-sequence report with its
      summary row (same index: the engine keeps input order in both). */
  function buildRows(payload) {
    var srows = payload.summary.rows;
    if (srows.length !== payload.sequences.length) {
      throw new Error("summary rows (" + srows.length + ") and sequences (" +
        payload.sequences.length + ") differ in length");
    }
    return payload.sequences.map(function (s, i) {
      var r = srows[i];
      if (r.id !== s.id) throw new Error("summary row " + i + " is '" + r.id + "', sequence is '" + s.id + "'");
      var a = s.annotation || null;
      var cdrs = a ? Object.keys(a.cdrs).map(function (k) {
        var c = a.cdrs[k];
        return { name: k, start: c.start, end: c.end, length: c.length, seq: c.seq };
      }) : [];
      var whole = s.properties ? s.properties["도메인 전체"] : null;
      var cdrCounts = zero();
      if (r.input_ok && !r.regions_unassigned) {
        CDR_REGIONS.forEach(function (reg) { addCounts(cdrCounts, r.regions[reg]); });
      }
      var row = {
        index: i,
        id: s.id,
        inputOk: s.input_ok,
        chain: a ? a.chain : null,
        confident: !!(a && a.confident),
        annotated: !!r.annotated,
        regionsUnassigned: !!r.regions_unassigned,
        status: s.triage.status,
        reasons: s.triage.reasons,
        note: s.triage.note,
        cdrs: cdrs,
        cdr3: cdrs.filter(function (c) { return /3$/.test(c.name); })[0] || null,
        cdrCounts: cdrCounts,
        fwCounts: r.input_ok ? addCounts(zero(), r.regions.framework) : zero(),
        flags: r.domain_flags || [],
        length: s.input_ok ? s.length : null,
        pI: whole ? whole.pI : null,
        gravy: whole ? whole.GRAVY : null,
        charge: whole ? whole.net_charge : null,
        pH: whole ? whole.pH : null,
        errors: s.errors || [],
        sequence: s.sequence || ""
      };
      row.cats = categories(row);
      return row;
    });
  }

  /** Counts for the sidebar and the stat line. */
  function tally(rows) {
    var t = { all: rows.length, review: 0, recheck: 0, errors: 0, clear: 0, H: 0, L: 0 };
    rows.forEach(function (r) {
      r.cats.forEach(function (c) { t[c] += 1; });
      if (r.status === "기준 해당 없음") t.clear += 1;
      if (r.chain === "H") t.H += 1;
      if (r.chain === "L") t.L += 1;
    });
    return t;
  }

  /** view: all|review|recheck|errors · chain: all|H|L · query: name, or residues in a CDR
      or in the whole sequence (3+ letters). */
  function filterRows(rows, opts) {
    var view = opts.view || "all", chain = opts.chain || "all";
    var q = String(opts.query || "").trim();
    var qName = q.toLowerCase();
    var qSeq = /^[A-Za-z]{3,}$/.test(q) ? q.toUpperCase() : null;
    return rows.filter(function (r) {
      if (view !== "all" && r.cats.indexOf(view) < 0) return false;
      if (chain !== "all" && r.chain !== chain) return false;
      if (!q) return true;
      if (r.id.toLowerCase().indexOf(qName) >= 0) return true;
      return !!(qSeq && r.sequence.indexOf(qSeq) >= 0);
    });
  }

  /** input (default) · name · pi-desc · pi-asc. Sorting never ranks by liabilities: the
      tool gives no overall verdict. Rows without a pI go last. */
  function sortRows(rows, key) {
    var out = rows.slice();
    var byIndex = function (a, b) { return a.index - b.index; };
    if (key === "name") {
      out.sort(function (a, b) {
        return a.id.localeCompare(b.id, "ko", { numeric: true }) || byIndex(a, b);
      });
    } else if (key === "pi-desc" || key === "pi-asc") {
      var dir = key === "pi-desc" ? -1 : 1;
      out.sort(function (a, b) {
        if (a.pI === null && b.pI === null) return byIndex(a, b);
        if (a.pI === null) return 1;
        if (b.pI === null) return -1;
        return (a.pI - b.pI) * dir || byIndex(a, b);
      });
    } else {
      out.sort(byIndex);
    }
    return out;
  }

  /** Identical CDR sequences shared by two or more candidates: [{name, seq, ids}]. */
  function sharedCdrs(rows) {
    var groups = {}, order = [];
    rows.forEach(function (r) {
      r.cdrs.forEach(function (c) {
        var key = c.name + "|" + c.seq;
        if (!groups[key]) { groups[key] = { name: c.name, seq: c.seq, ids: [] }; order.push(key); }
        if (groups[key].ids.indexOf(r.id) < 0) groups[key].ids.push(r.id);
      });
    });
    return order.map(function (k) { return groups[k]; }).filter(function (g) { return g.ids.length > 1; });
  }

  /** Round axis ticks covering [min, max]: {lo, hi, step, ticks}. */
  function niceTicks(min, max, count) {
    count = count || 5;
    if (!(isFinite(min) && isFinite(max))) throw new Error("niceTicks needs finite bounds, got " + min + ", " + max);
    if (min === max) { min -= 0.5; max += 0.5; }
    if (min > max) { var t = min; min = max; max = t; }
    var raw = (max - min) / Math.max(1, count - 1);
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    var f = raw / mag;
    var step = (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * mag;
    var lo = Math.floor(min / step) * step, hi = Math.ceil(max / step) * step;
    var ticks = [], decimals = Math.max(0, -Math.floor(Math.log10(step)) + (step / mag === 2.5 ? 1 : 0));
    for (var v = lo; v <= hi + step / 2; v += step) ticks.push(Number(v.toFixed(decimals)));
    return { lo: lo, hi: hi, step: step, ticks: ticks };
  }

  var TABLE_HEAD = ["순번", "서열 이름", "체인", "분류", "주석", "CDR1", "CDR2", "CDR3",
    "CDR 플래그", "framework 플래그", "도메인 플래그", "길이", "pI", "GRAVY", "net charge", "pH",
    "선별 근거", "입력 오류"];

  function tableRecords(rows) {
    return rows.map(function (r) {
      var cdr = function (n) { return r.cdrs[n] ? r.cdrs[n].seq : ""; };
      return [
        r.index + 1, r.id, r.chain === "H" ? "VH" : r.chain === "L" ? "VL" : "",
        r.status, r.inputOk ? (r.confident ? "신뢰 가능" : "신뢰 불가") : "",
        cdr(0), cdr(1), cdr(2),
        r.inputOk && !r.regionsUnassigned ? countText(r.cdrCounts) : (r.inputOk ? "주석 실패" : ""),
        r.inputOk ? countText(r.fwCounts) : "",
        r.flags.map(function (f) { return f.label + " (" + f.risk + ")"; }).join(", "),
        r.length === null ? "" : r.length,
        r.pI === null ? "" : r.pI, r.gravy === null ? "" : r.gravy,
        r.charge === null ? "" : r.charge, r.pH === null ? "" : r.pH,
        r.reasons.join(" / "), r.errors.join(" / ")
      ];
    });
  }

  /** CSV for Excel: UTF-8 with BOM (so Korean opens correctly), CRLF, RFC 4180 quoting. */
  function toCsv(rows) {
    var q = function (v) {
      var s = String(v);
      return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    var lines = [TABLE_HEAD].concat(tableRecords(rows)).map(function (rec) { return rec.map(q).join(","); });
    return "﻿" + lines.join("\r\n") + "\r\n";
  }

  /** Tab-separated text for pasting into Excel or Sheets. */
  function toTsv(rows) {
    var clean = function (v) { return String(v).replace(/[\t\r\n]+/g, " "); };
    return [TABLE_HEAD].concat(tableRecords(rows)).map(function (rec) { return rec.map(clean).join("\t"); }).join("\n");
  }

  /** Sequences as FASTA text (valid inputs only). */
  function toFasta(rows) {
    return rows.filter(function (r) { return r.inputOk; }).map(function (r) {
      return ">" + r.id + "\n" + r.sequence.replace(/(.{60})/g, "$1\n").replace(/\n$/, "") + "\n";
    }).join("");
  }

  return {
    RISKS: RISKS, MOTIFS: MOTIFS, VALIDATION: VALIDATION, TABLE_HEAD: TABLE_HEAD,
    countText: countText, categories: categories, buildRows: buildRows, tally: tally,
    filterRows: filterRows, sortRows: sortRows, sharedCdrs: sharedCdrs, niceTicks: niceTicks,
    toCsv: toCsv, toTsv: toTsv, toFasta: toFasta
  };
});
