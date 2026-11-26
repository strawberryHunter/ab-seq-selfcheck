/*
 * app-lib.js -- pure helpers behind the AbCheck browser app (no DOM, no dependencies).
 *
 * The app (web/app.html) renders what these functions return; tests/test_web_app.py runs
 * them under Node against the real engine output. Defines the global ABUI, and
 * module.exports under Node.
 *
 * Triage rule sets live here, not in the engine: the engine always applies its default
 * rules (abcheck/triage.py DEFAULT_RULES), and applyRules() re-applies a researcher's rule
 * set to the engine's liability hits. With DEFAULT_RULESET it reproduces the engine's
 * triage exactly (a test checks this on the whole synthetic corpus).
 */
(function (root, factory) {
  "use strict";
  var api = factory();
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ABUI = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  var RISKS = ["높음", "중간", "낮음"];
  var RISK_RANK = { "낮음": 0, "중간": 1, "높음": 2 };
  var CDR_REGIONS = ["CDR1", "CDR2", "CDR3"];
  var STATUS = { REVIEW: "검토 필요", CLEAR: "기준 해당 없음", UNKNOWN: "판정 불가" };

  /* The liability motifs exactly as abcheck/liabilities.py MOTIF_SPECS defines them, in the
     same order (pattern, label, base risk, top-of-kind). A test keeps this equal to Python. */
  var MOTIFS = [
    { pattern: "N[^P][ST]", label: "N-결합 당화 서열", short: "N-결합 당화", base: "높음", top: true },
    { pattern: "NG", label: "deamidation (NG)", short: "NG", base: "높음", top: true },
    { pattern: "N[STNH]", label: "deamidation (NS/NT/NN/NH)", short: "NS/NT/NN/NH", base: "중간", top: false },
    { pattern: "DG", label: "isomerization (DG)", short: "DG", base: "높음", top: true },
    { pattern: "D[STD]", label: "isomerization (DS/DT/DD)", short: "DS/DT/DD", base: "중간", top: false },
    { pattern: "DP", label: "fragmentation (DP)", short: "DP", base: "중간", top: true },
    { pattern: "M", label: "산화 (Met)", short: "Met 산화", base: "낮음", top: false },
    { pattern: "W", label: "산화 (Trp)", short: "Trp 산화", base: "낮음", top: false }
  ];
  var MOTIF_BY_LABEL = {};
  MOTIFS.forEach(function (m) { MOTIF_BY_LABEL[m.label] = m; });

  /* The engine's default rules (abcheck/triage.py DEFAULT_RULES) as a rule set. */
  var DEFAULT_RULESET = {
    id: "default", name: "기본",
    motifs: ["N[^P][ST]", "NG", "DG", "DP"], where: "cdr", minRisk: "높음", cys: "cdr"
  };
  var WHERE = { cdr: "CDR 안", cdr3: "CDR3 안", all: "전체 구간" };
  var CYS = { cdr: "CDR 안 Cys", any: "항상", off: "끔" };

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

  function esc(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

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

  /* ------------------------------------------------------------------ rule sets */

  /** Validate a stored rule set; throws Error naming the bad field (never guesses). */
  function normalizeRuleSet(raw) {
    if (!raw || typeof raw !== "object") throw new Error("rule set is not an object");
    var patterns = MOTIFS.map(function (m) { return m.pattern; });
    if (!Array.isArray(raw.motifs) || raw.motifs.some(function (p) { return patterns.indexOf(p) < 0; })) {
      throw new Error("rule set '" + raw.name + "': motifs must be a list of " + patterns.join(", "));
    }
    if (!WHERE[raw.where]) throw new Error("rule set '" + raw.name + "': where must be cdr, cdr3 or all");
    if (!(raw.minRisk in RISK_RANK)) throw new Error("rule set '" + raw.name + "': minRisk must be 높음, 중간 or 낮음");
    if (!CYS[raw.cys]) throw new Error("rule set '" + raw.name + "': cys must be cdr, any or off");
    var name = String(raw.name || "").trim();
    if (!name) throw new Error("rule set without a name");
    return {
      id: String(raw.id || name), name: name,
      // keep MOTIFS order so two equal sets compare equal
      motifs: patterns.filter(function (p) { return raw.motifs.indexOf(p) >= 0; }),
      where: raw.where, minRisk: raw.minRisk, cys: raw.cys
    };
  }

  /** Short labels for a rule set: [motif line, Cys line]. */
  function describeRuleSet(rs) {
    var motifs = MOTIFS.filter(function (m) { return rs.motifs.indexOf(m.pattern) >= 0; })
      .map(function (m) { return m.short; });
    var risk = rs.minRisk === "높음" ? "위험도 높음" : "위험도 " + rs.minRisk + " 이상";
    return [
      motifs.length ? WHERE[rs.where] + " · " + motifs.join(" · ") + " · " + risk : "모티프 규칙 없음",
      "비결합 Cys · " + CYS[rs.cys]
    ];
  }

  function hitLabel(h) {
    return h.label.indexOf("(" + h.match + ")") >= 0 ? h.label : h.label + " (" + h.match + ")";
  }

  /**
   * Re-apply a rule set to one engine report (an element of payload.sequences).
   * Mirrors abcheck/triage.py evaluate(): input error or no CDRs -> 판정 불가 (not evaluable);
   * any rule hit -> 검토 필요; otherwise an unconfident annotation -> 판정 불가; else
   * 기준 해당 없음. Reason strings use the engine's wording.
   */
  function applyRules(s, rs) {
    if (!s.input_ok) return { status: STATUS.UNKNOWN, reasons: [], evaluable: false, note: s.triage.note };
    var a = s.annotation;
    if (!a || !Object.keys(a.cdrs).length) {
      return { status: STATUS.UNKNOWN, reasons: [], evaluable: false, note: s.triage.note };
    }
    var reasons = [];
    s.liabilities.hits.forEach(function (h) {
      var m = MOTIF_BY_LABEL[h.label];
      if (!m) throw new Error("unknown liability label from the engine: " + h.label);
      if (rs.motifs.indexOf(m.pattern) < 0) return;
      if (rs.where === "cdr" && !h.in_cdr) return;
      if (rs.where === "cdr3" && !(h.in_cdr && /3$/.test(h.region))) return;
      if (RISK_RANK[h.risk] < RISK_RANK[rs.minRisk]) return;
      reasons.push(h.region + " " + h.start + "번에 " + hitLabel(h) + " — 위험도 " + h.risk);
    });
    s.liabilities.domain_flags.forEach(function (f) {
      if (rs.cys === "off" || (rs.cys === "cdr" && f.risk !== "높음")) return;
      reasons.push(f.label + " — 위험도 " + f.risk + " (Cys " + f.cysteine_count + "개, 홀수)");
    });
    if (reasons.length) return { status: STATUS.REVIEW, reasons: reasons, evaluable: true, note: "" };
    if (!a.confident) {
      return {
        status: STATUS.UNKNOWN, reasons: [], evaluable: true,
        note: s.triage.status === STATUS.UNKNOWN ? s.triage.note : "주석 신뢰 불가 · CDR 경계를 ANARCI/IMGT 로 확인할 것."
      };
    }
    return { status: STATUS.CLEAR, reasons: [], evaluable: true, note: "" };
  }

  /* ------------------------------------------------------------------ rows */

  /** Which lists a row belongs to; keys double as the app's view names. */
  function categories(row) {
    if (!row.inputOk) return ["errors"];
    if (row.status === STATUS.REVIEW) return ["review"];
    if (row.status === STATUS.UNKNOWN) return ["unknown"];
    return ["clear"];
  }

  /** One flat record per input sequence: the engine's report and summary row (same index),
      with triage re-applied under the given rule set. */
  function buildRows(payload, ruleSet) {
    var rs = ruleSet || DEFAULT_RULESET;
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
      var t = applyRules(s, rs);
      var row = {
        index: i,
        id: s.id,
        inputOk: s.input_ok,
        chain: a ? a.chain : null,
        confident: !!(a && a.confident),
        annotated: !!r.annotated,
        regionsUnassigned: !!r.regions_unassigned,
        status: t.status,
        reasons: t.reasons,
        note: t.note,
        cdrs: cdrs,
        cdr3: cdrs.filter(function (c) { return /3$/.test(c.name); })[0] || null,
        cdrCounts: cdrCounts,
        fwCounts: r.input_ok ? addCounts(zero(), r.regions.framework) : zero(),
        flags: r.domain_flags || [],
        hits: s.liabilities ? s.liabilities.hits : [],
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

  /** Counts for the sidebar, the stat line and the hero. */
  function tally(rows, picked) {
    var t = { all: rows.length, review: 0, unknown: 0, clear: 0, errors: 0, picked: 0, H: 0, L: 0, low: 0 };
    rows.forEach(function (r) {
      t[r.cats[0]] += 1;
      if (r.chain === "H") t.H += 1;
      if (r.chain === "L") t.L += 1;
      if (r.inputOk && !r.confident) t.low += 1;
      if (picked && picked.has(r.index)) t.picked += 1;
    });
    return t;
  }

  /** view: all|review|unknown|clear|errors|picked · chain: all|H|L · conf: all|yes|no ·
      query: a name, a motif name, residues in a CDR (2+ letters) or in the sequence (3+). */
  function filterRows(rows, opts) {
    var view = opts.view || "all", chain = opts.chain || "all", conf = opts.conf || "all";
    var picked = opts.picked, motif = opts.motif || null;
    var q = String(opts.query || "").trim();
    var qName = q.toLowerCase();
    var letters = /^[A-Za-z]+$/.test(q) ? q.toUpperCase() : null;
    return rows.filter(function (r) {
      if (view === "picked") { if (!picked || !picked.has(r.index)) return false; }
      else if (view !== "all" && r.cats.indexOf(view) < 0) return false;
      if (chain !== "all" && r.chain !== chain) return false;
      if (conf === "yes" && !(r.inputOk && r.confident)) return false;
      if (conf === "no" && !(r.inputOk && !r.confident)) return false;
      if (motif) {
        if (motif.flag) { if (!r.flags.some(function (f) { return f.label === motif.label; })) return false; }
        else if (!r.hits.some(function (h) { return h.label === motif.label && (!motif.region || h.region === motif.region); })) return false;
      }
      if (!q) return true;
      if (r.id.toLowerCase().indexOf(qName) >= 0) return true;
      if (r.hits.some(function (h) { return h.label.toLowerCase().indexOf(qName) >= 0; })) return true;
      if (letters && letters.length >= 2 && r.cdrs.some(function (c) { return c.seq.indexOf(letters) >= 0; })) return true;
      return !!(letters && letters.length >= 3 && r.sequence.indexOf(letters) >= 0);
    });
  }

  /** input (default) · name · pi-desc · pi-asc. Never by liabilities: no verdicts.
      Rows without a pI go last. */
  function sortRows(rows, key) {
    var out = rows.slice();
    var byIndex = function (a, b) { return a.index - b.index; };
    if (key === "name") {
      out.sort(function (a, b) { return a.id.localeCompare(b.id, "ko", { numeric: true }) || byIndex(a, b); });
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

  /**
   * One CDR written against a reference CDR: "·" where the residue is the same, the residue
   * itself where it differs ("··Q·············"). Two CDRs of different length cannot be laid
   * over each other without an alignment, and this tool makes none: then comparable is false
   * and text is the sequence as it is. Returns {comparable, text, count}.
   */
  function diffCdr(ref, other) {
    if (!other) return { comparable: false, text: "", count: null };
    if (!ref || ref.seq.length !== other.seq.length) return { comparable: false, text: other.seq, count: null };
    var text = "", count = 0;
    for (var i = 0; i < other.seq.length; i++) {
      if (other.seq[i] === ref.seq[i]) text += "·"; else { text += other.seq[i]; count += 1; }
    }
    return { comparable: true, text: text, count: count };
  }

  /**
   * The list without one entry. Picks and comparison references are kept by place in the list, so
   * they are renumbered to stay on the same candidates; those of the removed entry are dropped.
   * Returns {entries, picked, ref}.
   */
  function withoutEntry(entries, index, picked, ref) {
    if (!(Number.isInteger(index) && index >= 0 && index < entries.length)) {
      throw new Error("withoutEntry: no entry " + index + " in a list of " + entries.length);
    }
    var shift = function (i) { return i > index ? i - 1 : i; };
    var newRef = {};
    Object.keys(ref || {}).forEach(function (k) {
      var v = ref[k];
      newRef[k] = v === null || v === undefined || v === index ? null : shift(v);
    });
    return {
      entries: entries.slice(0, index).concat(entries.slice(index + 1)),
      picked: (picked || []).filter(function (i) { return i !== index; }).map(shift),
      ref: newRef
    };
  }

  /* ---------- motifs over the panel ---------- */
  var REGION_ORDER = ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3", "framework", "구간 미확정"];

  /**
   * How the detected motifs spread over the panel. For each motif (the engine's table order, every
   * motif listed even when nothing was found) the number of candidates with at least one such hit
   * in each region, and in any region. It counts candidates, not hits: two NG in one CDR-H2 count
   * once. Nobody is ranked by it. Regions are the engine's own names; only those that occur are
   * columns. Returns {evaluated, regions, rows: [{label, cells: {region: n}, total}],
   * flags: [{label, total}]} -- flags are the whole-domain ones (unpaired Cys).
   */
  function motifTable(rows) {
    var ok = rows.filter(function (r) { return r.inputOk; });
    var order = MOTIFS.map(function (m) { return m.label; });
    var by = {}, regionSeen = {}, flagCount = {}, flagOrder = [];
    order.forEach(function (label) { by[label] = { label: label, cells: {}, total: 0 }; });
    ok.forEach(function (r) {
      var mine = {};
      r.hits.forEach(function (h) { (mine[h.label] = mine[h.label] || {})[h.region] = true; });
      Object.keys(mine).forEach(function (label) {
        if (!by[label]) { by[label] = { label: label, cells: {}, total: 0 }; order.push(label); }
        by[label].total += 1;
        Object.keys(mine[label]).forEach(function (region) {
          by[label].cells[region] = (by[label].cells[region] || 0) + 1;
          regionSeen[region] = true;
        });
      });
      var flagged = {};
      r.flags.forEach(function (f) { flagged[f.label] = true; });
      Object.keys(flagged).forEach(function (label) {
        if (!(label in flagCount)) { flagCount[label] = 0; flagOrder.push(label); }
        flagCount[label] += 1;
      });
    });
    return {
      evaluated: ok.length,
      regions: REGION_ORDER.filter(function (x) { return regionSeen[x]; })
        .concat(Object.keys(regionSeen).filter(function (x) { return REGION_ORDER.indexOf(x) < 0; })),
      rows: order.map(function (label) { return by[label]; }),
      flags: flagOrder.map(function (label) { return { label: label, total: flagCount[label] }; })
    };
  }

  /* ---------- substitution trial ---------- */
  var AA20 = "ACDEFGHIKLMNPQRSTVWY";

  /**
   * Substitutions written the way people write them ("N56Q"; several separated by spaces, commas,
   * "/" or "+"), checked against the sequence. The number is the 1-based place in the sequence as
   * entered, the same numbers as the CDR table (not a Kabat number). The original residue is
   * required and must be the one found there: a typo or an off-by-one is refused, never corrected.
   * Returns {subs: [{pos, from, to}] sorted by place, errors: [one line each]}.
   */
  function parseSubstitutions(text, sequence) {
    var subs = [], errors = [], seen = {};
    String(text || "").split(/[\s,;\/+]+/).filter(Boolean).forEach(function (raw) {
      var tok = raw.toUpperCase();
      var m = /^([A-Z])([0-9]+)([A-Z])$/.exec(tok);
      if (!m) {
        errors.push(/^[A-Z][0-9]+$/.test(tok) ? tok + " 뒤에 바꿀 잔기 한 글자 필요" : tok + ": 치환 표기 아님 (예 N56Q)");
        return;
      }
      var from = m[1], pos = Number(m[2]), to = m[3];
      if (pos < 1 || pos > sequence.length) { errors.push(tok + ": 서열은 1–" + sequence.length + "번"); return; }
      if (sequence[pos - 1] !== from) { errors.push(tok + ": " + pos + "번은 " + sequence[pos - 1]); return; }
      if (AA20.indexOf(to) < 0) { errors.push(tok + ": " + to + " 는 표준 아미노산 아님"); return; }
      if (to === from) { errors.push(tok + ": 같은 잔기"); return; }
      if (seen[pos]) { errors.push(tok + ": " + pos + "번을 두 번 적음"); return; }
      seen[pos] = true;
      subs.push({ pos: pos, from: from, to: to });
    });
    subs.sort(function (a, b) { return a.pos - b.pos; });
    return { subs: subs, errors: errors };
  }

  /** The sequence with the substitutions made. Refuses one whose original residue is not there. */
  function applySubstitutions(sequence, subs) {
    var chars = String(sequence).split("");
    subs.forEach(function (s) {
      if (chars[s.pos - 1] !== s.from) {
        throw new Error("applySubstitutions: place " + s.pos + " holds " + chars[s.pos - 1] + ", not " + s.from);
      }
      chars[s.pos - 1] = s.to;
    });
    return chars.join("");
  }

  var subText = function (s) { return s.from + s.pos + s.to; };
  /** Name of a variant: the parent's name and its substitutions ("mAb03_VH_N56Q_K100G"). */
  function variantName(id, subs) { return id + "_" + subs.map(subText).join("_"); }

  /** Where a place lies in a row's own annotation, in the engine's words. */
  function regionAt(row, pos) {
    if (!row.cdrs.length) return "구간 미확정";
    for (var i = 0; i < row.cdrs.length; i++) {
      if (row.cdrs[i].start <= pos && pos <= row.cdrs[i].end) return row.cdrs[i].name;
    }
    return "framework";
  }

  /**
   * The hits of a variant against the hits of its parent. A hit is the same hit when it has the
   * same label at the same place (substitutions keep every place where it was).
   * - removed: in the parent only; added: in the variant only;
   * - moved: in both, but its region or risk differs (the substitution changed the annotation);
   * - kept: how many are in both, unchanged.
   */
  function diffHits(before, after) {
    var key = function (h) { return h.label + "|" + h.start + "|" + h.end; };
    var b = {}, a = {};
    before.forEach(function (h) { b[key(h)] = h; });
    after.forEach(function (h) { a[key(h)] = h; });
    var moved = [];
    after.forEach(function (h) {
      var was = b[key(h)];
      if (was && (was.region !== h.region || was.risk !== h.risk)) moved.push({ before: was, after: h });
    });
    return {
      removed: before.filter(function (h) { return !a[key(h)]; }),
      added: after.filter(function (h) { return !b[key(h)]; }),
      moved: moved,
      kept: after.filter(function (h) { return b[key(h)]; }).length - moved.length
    };
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

  /** Merge several parsed inputs (one per file: {name, entries, source}) into one analysis
      input. The source names the first file so a past analysis can be recognised later. */
  function combineEntries(parts) {
    if (!parts.length) throw new Error("combineEntries needs at least one input");
    var entries = [];
    parts.forEach(function (p) { entries = entries.concat(p.entries); });
    var source = parts.length === 1 ? parts[0].source
      : (parts[0].name || "파일") + " 외 " + (parts.length - 1) + "개";
    return { entries: entries, source: source };
  }

  // Marks for the heavy or the light chain in a name. A mark is a whole token, so the "L" in
  // "CTRL" is not one. Numbered marks ("VH1") and gene names ("IGHV3") are left alone: they say
  // more than which chain it is, and two names differing in them need not be one antibody.
  var CHAIN_TOKEN = /^(vh|vl|vk|vκ|vλ|vkappa|vlambda|hc|lc|h|l|κ|λ|heavy|light|kappa|lambda|igh|igk|igl|중쇄|경쇄)$/i;
  var NAME_SEP = /([-_. |\/:;,()\[\]]+)/;      // separators are kept by split(), at the odd positions

  /**
   * A name cut into tokens and separators, with the positions of its chain marks:
   * {parts, drop: [even indexes]} or null when there is no mark. The word "chain" (or "사슬")
   * right after a mark belongs to it ("H_chain", "heavy chain"). Names written without a
   * separator are read too: upper-case VH / VL / VK / HC / LC at the end ("mAb01VH"), or a single
   * H / L right after a digit ("mAb01H", but not "CTRL").
   */
  function chainMarks(id) {
    id = String(id);
    var parts = id.split(NAME_SEP), drop = [];
    for (var i = 0; i < parts.length; i += 2) {
      if (CHAIN_TOKEN.test(parts[i])) {
        drop.push(i);
        if (i + 2 < parts.length && /^(chain|사슬)$/i.test(parts[i + 2])) { drop.push(i + 2); i += 2; }
      }
    }
    if (drop.length) return { parts: parts, drop: drop };
    var m = /^(.*[0-9a-z])(VH|VL|VK|HC|LC)$/.exec(id) || /^(.*[0-9])(H|L)$/.exec(id);
    return m ? { parts: [m[1], "", m[2]], drop: [2] } : null;
  }

  /** A candidate's name without its chain marks, lower-cased: "mAb01_VH" -> "mab01",
      "SYN-VL-003" -> "syn-003", "Ab7 heavy" -> "ab7", "mAb01|VH" -> "mab01", "1A2H" -> "1a2".
      null when the name has no chain mark, or nothing besides it. */
  function pairKey(id) {
    var f = chainMarks(id);
    if (!f) return null;
    var rest = f.parts.filter(function (t, i) { return i % 2 === 0 && t && f.drop.indexOf(i) < 0; });
    return rest.length ? rest.join("-").toLowerCase() : null;
  }

  /**
   * The point nearest to (x, y), as an index into points, or -1 when none lies within maxDist.
   * Used for hover and click on the scatter: with one hit area per point, the area drawn last
   * covered its neighbours and a click on one mark opened another candidate.
   */
  function nearestPoint(points, x, y, maxDist) {
    var best = -1, bestD = maxDist * maxDist;
    points.forEach(function (p, i) {
      var d = (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
      if (d <= bestD && (best < 0 || d < bestD)) { best = i; bestD = d; }
    });
    return best;
  }

  /**
   * Which points of a scatter plot get their name written next to them. A name goes to the
   * right of its point, else to the left, else above, else below, and only when its box stays
   * inside the plot and overlaps no other point and no name placed before it. Points far
   * from the middle are tried first, so outliers keep their names; the others are named on
   * hover and in the table. Points on the same spot (the same chain under several names) share
   * one name, the first one's: four names around a single mark would read as four marks.
   *   points: [{x, y, text}] in pixels
   *   opts: {charW, h, gap, r, left, right, top, bottom} -- text metrics, mark radius, plot box
   * Returns [{i, side: "right" | "left" | "above" | "below", x0, x1, y0, y1}] in point order.
   */
  function placeLabels(points, opts) {
    var n = points.length;
    if (!n) return [];
    var cx = 0, cy = 0;
    points.forEach(function (p) { cx += p.x / n; cy += p.y / n; });
    var order = points.map(function (_, i) { return i; }).sort(function (a, b) {
      var da = Math.pow(points[a].x - cx, 2) + Math.pow(points[a].y - cy, 2);
      var db = Math.pow(points[b].x - cx, 2) + Math.pow(points[b].y - cy, 2);
      return db - da || a - b;
    });
    var hits = function (a, b) { return a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1; };
    var marks = points.map(function (p) { return { x0: p.x - opts.r, x1: p.x + opts.r, y0: p.y - opts.r, y1: p.y + opts.r }; });
    var placed = [];
    order.forEach(function (i) {
      var p = points[i], w = String(p.text).length * opts.charW;
      if (points.some(function (q, j) { return j < i && Math.abs(q.x - p.x) < 1 && Math.abs(q.y - p.y) < 1; })) return;
      ["right", "left", "above", "below"].some(function (side) {
        var beside = side === "right" || side === "left";
        var x0 = side === "right" ? p.x + opts.gap : side === "left" ? p.x - opts.gap - w : p.x - w / 2;
        var yc = beside ? p.y : side === "above" ? p.y - opts.gap - opts.h / 2 : p.y + opts.gap + opts.h / 2;
        var box = { i: i, side: side, x0: x0, x1: x0 + w, y0: yc - opts.h / 2, y1: yc + opts.h / 2 };
        if (box.x0 < opts.left || box.x1 > opts.right || box.y0 < opts.top || box.y1 > opts.bottom) return false;
        if (marks.some(function (m, j) { return j !== i && hits(box, m); })) return false;
        if (placed.some(function (q) { return hits(box, q); })) return false;
        placed.push(box);
        return true;
      });
    });
    return placed.sort(function (a, b) { return a.i - b.i; });
  }

  /**
   * The object to save as JSON. With the default rule set it is the engine's payload itself,
   * so the file equals the Python CLI's --json output. With another rule set the engine keys
   * stay as they are and two are added: "rule_set" (which rules) at the top and "rule_result"
   * (that set's status and reasons) next to each sequence's own "triage". The payload is not
   * modified.
   */
  function withRuleResults(payload, ruleSet) {
    var rs = normalizeRuleSet(ruleSet);
    if (rs.id === DEFAULT_RULESET.id) return payload;
    var out = {};
    Object.keys(payload).forEach(function (k) {
      if (k === "sequences") {
        out.rule_set = { name: rs.name, motifs: rs.motifs.slice(), where: rs.where, min_risk: rs.minRisk, cys: rs.cys,
                         description: describeRuleSet(rs) };
        out.sequences = payload.sequences.map(function (s) {
          var copy = {};
          Object.keys(s).forEach(function (key) {
            copy[key] = s[key];
            if (key === "triage") copy.rule_result = applyRules(s, rs);
          });
          return copy;
        });
      } else {
        out[k] = payload[k];
      }
    });
    return out;
  }

  /**
   * Position numbers and a closing stop mark that come along when a sequence is copied from a
   * database page or a PDF ("   61 gkglewvsai ..." or a final "*"). Returns {text, digits, stops}:
   * the text without digits on sequence lines and without a "*" that ends a record, and how
   * many of each were removed. Header lines (">name") stay as they are. A "*" inside a sequence
   * stays too: that is a stop in the middle, and the engine should say so. Pasted tables and
   * JSON are returned unchanged -- their names may hold digits.
   */
  function stripNumbering(text) {
    text = String(text);
    if (pastedKind(text)) return { text: text, digits: 0, stops: 0 };
    var lines = text.split("\n"), digits = 0, stops = 0, last = -1;
    var closeRecord = function () {
      if (last < 0) return;
      var m = /^(.*)\*(\s*)$/.exec(lines[last]);
      if (m) { lines[last] = m[1] + m[2]; stops++; }
      last = -1;
    };
    for (var i = 0; i < lines.length; i++) {
      if (/^\s*>/.test(lines[i])) { closeRecord(); continue; }
      var cleaned = lines[i].replace(/[0-9]/g, "");
      digits += lines[i].length - cleaned.length;
      lines[i] = cleaned;
      if (cleaned.trim()) last = i;
    }
    closeRecord();
    return { text: lines.join("\n"), digits: digits, stops: stops };
  }

  /** An antibody's display name from one chain's name: the name without its chain marks
      ("mAb02_VH" -> "mAb02", "SYN-VL-003" -> "SYN-003", "mAb01 (VH)" -> "mAb01"). null when
      there is no chain mark. */
  function pairName(id) {
    var f = chainMarks(id);
    if (!f) return null;
    var keep = f.parts.slice();
    f.drop.forEach(function (i) {            // each mark goes with one separator next to it
      keep[i] = null;
      if (i > 0 && keep[i - 1] !== null) keep[i - 1] = null;
      else if (i + 1 < keep.length) keep[i + 1] = null;
    });
    var name = keep.filter(function (x) { return x !== null; }).join("")
      .replace(/^[-_. |\/:;,)\]]+|[-_. |\/:;,(\[]+$/g, "");
    var opened = (name.match(/[(\[]/g) || []).length, closed = (name.match(/[)\]]/g) || []).length;
    if (opened !== closed) name = name.replace(/[()\[\]]/g, "").trim();     // a bracket left without its partner
    return name || null;
  }

  // ---- a rule set as text, to hand to a colleague -------------------------------------

  /** One line of JSON that ruleSetFromText() reads back. The id stays behind: it belongs to this browser. */
  function ruleSetToText(ruleSet) {
    var rs = normalizeRuleSet(ruleSet);
    return JSON.stringify({ abcheck: "rule-set", name: rs.name, motifs: rs.motifs, where: rs.where, minRisk: rs.minRisk, cys: rs.cys });
  }

  /**
   * The rule set in a pasted text: the line ruleSetToText() wrote, or a JSON result saved from
   * this page (its "rule_set" block is used). Returns {name, motifs, where, minRisk, cys} checked
   * by normalizeRuleSet; throws an Error that says what is wrong. Text around the JSON (a mail
   * quote, a chat prefix) is ignored.
   */
  function ruleSetFromText(text) {
    text = String(text);
    var from = text.indexOf("{"), to = text.lastIndexOf("}");
    if (from < 0 || to <= from) throw new Error("no rule set in the text");
    var value;
    try { value = JSON.parse(text.slice(from, to + 1)); }
    catch (e) { throw new Error("the text is not complete JSON: " + e.message); }
    if (value && typeof value === "object" && value.rule_set && typeof value.rule_set === "object") value = value.rule_set;
    if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("no rule set in the text");
    var rs = normalizeRuleSet({
      id: "imported", name: value.name, motifs: value.motifs, where: value.where,
      minRisk: value.minRisk !== undefined ? value.minRisk : value.min_risk, cys: value.cys
    });
    return { name: rs.name, motifs: rs.motifs, where: rs.where, minRisk: rs.minRisk, cys: rs.cys };
  }

  // ---- values read back from browser storage ------------------------------------------
  // They were written by this page, possibly by an older version of it, possibly cut short.
  // Each reader keeps what has the expected shape and counts what it had to leave out, so
  // the caller can tell the user. Nothing is repaired by guessing: an analysis without a
  // usable pH or time is left out rather than restored with a made-up one.

  /** Stored analyses -> {list, dropped}. */
  function cleanRecent(value) {
    if (value === null || value === undefined) return { list: [], dropped: 0 };
    if (!Array.isArray(value)) return { list: [], dropped: 1 };
    var list = [], dropped = 0;
    value.forEach(function (r) {
      var ok = r && typeof r === "object" && Array.isArray(r.entries) && r.entries.length > 0 &&
        r.entries.every(function (e) { return Array.isArray(e) && e.length === 2 && typeof e[0] === "string" && typeof e[1] === "string"; }) &&
        typeof r.ph === "number" && r.ph >= 0 && r.ph <= 14 && typeof r.at === "number" && isFinite(r.at);
      if (!ok) { dropped += 1; return; }
      var n = r.entries.length;
      list.push({
        id: typeof r.id === "string" && r.id ? r.id : "r" + r.at,
        entries: r.entries, ph: r.ph, at: r.at, count: n,
        source: typeof r.source === "string" && r.source ? r.source : "저장된 분석",
        text: typeof r.text === "string" ? r.text : "",
        picked: Array.isArray(r.picked) ? r.picked.filter(function (i) { return Number.isInteger(i) && i >= 0 && i < n; }) : []
      });
    });
    return { list: list, dropped: dropped };
  }

  /** Stored list settings -> the same fields with allowed values only (views: the screen names). */
  function cleanUi(value, views) {
    var v = value && typeof value === "object" && !Array.isArray(value) ? value : {};
    var pick = function (x, allowed, fallback) { return allowed.indexOf(x) >= 0 ? x : fallback; };
    return {
      chain: pick(v.chain, ["all", "H", "L"], "all"),
      conf: pick(v.conf, ["all", "yes", "no"], "all"),
      sort: pick(v.sort, ["input", "name", "pi-desc", "pi-asc"], "input"),
      view: pick(v.view, views, "all"),
      // the unit the user picked, if any. Up to v15 "seq" was stored without a click, "ab" only after one.
      unitChoice: pick(v.unitChoice, ["seq", "ab"], v.unit === "ab" ? "ab" : null)
    };
  }

  /**
   * Stored rule sets -> {active, sets, dropped, errors}. A set that fails normalizeRuleSet, or
   * that carries the built-in set's id or an id already seen, is left out: under the id
   * "default" it would replace the built-in rules while still being shown as 기본.
   */
  function cleanRules(value) {
    var out = { active: DEFAULT_RULESET.id, sets: [], dropped: 0, errors: [] };
    if (!value || typeof value !== "object" || Array.isArray(value) || value.sets === undefined) return out;
    if (!Array.isArray(value.sets)) { out.dropped = 1; out.errors.push("sets is not a list"); return out; }
    var seen = {};
    seen[DEFAULT_RULESET.id] = true;
    value.sets.forEach(function (s) {
      var rs;
      try { rs = normalizeRuleSet(s); }
      catch (e) { out.dropped += 1; out.errors.push(e.message); return; }
      if (seen[rs.id]) { out.dropped += 1; out.errors.push("rule set id used twice: " + rs.id); return; }
      seen[rs.id] = true;
      out.sets.push(rs);
    });
    if (typeof value.active === "string" && seen[value.active]) out.active = value.active;
    return out;
  }

  /**
   * The list unit to start with while the user has not picked one: "ab" (one card per antibody)
   * when more than half of the chains the engine identified have a partner, else "seq".
   */
  function preferredUnit(rows, pairs) {
    var chains = rows.filter(function (r) { return r.chain === "H" || r.chain === "L"; }).length;
    return chains > 0 && Object.keys(pairs).length * 2 > chains ? "ab" : "seq";
  }

  /** Partner chains by name: {index: partnerIndex} for names that differ only in the chain
      token and whose chains the engine called H and L. Ambiguous groups are left unpaired. */
  function pairMap(rows) {
    var groups = {};
    rows.forEach(function (r) {
      if (!r.inputOk || !r.chain) return;
      var k = pairKey(r.id);
      if (k === null) return;
      (groups[k] = groups[k] || []).push(r);
    });
    var out = {};
    Object.keys(groups).forEach(function (k) {
      var g = groups[k];
      if (g.length === 2 && g[0].chain !== g[1].chain) { out[g[0].index] = g[1].index; out[g[1].index] = g[0].index; }
    });
    return out;
  }

  /** Rows in input order, but with each row's partner (when it is in the list too) right
      after it -- so a sheet reads antibody by antibody. */
  function withPartnersAdjacent(rows, pairs) {
    var byIndex = {}, done = {}, out = [];
    rows.forEach(function (r) { byIndex[r.index] = r; });
    rows.forEach(function (r) {
      if (done[r.index]) return;
      var partner = byIndex[pairs[r.index]];
      var first = partner && !done[partner.index] && partner.chain === "H" ? partner : r;
      [first, first === r ? partner : r].forEach(function (x) { if (x && !done[x.index]) { done[x.index] = true; out.push(x); } });
    });
    return out;
  }

  /**
   * Cards for the antibody view. A listed row that has a partner becomes one unit
   * {name, h, l} -- the partner is attached even when the current filter left it out, so an
   * antibody always shows both chains. Other rows stay single: {row}. Order of first
   * appearance. The two chains keep their own triage results; nothing is merged.
   */
  function antibodyUnits(shown, all, pairs) {
    var done = {}, out = [];
    shown.forEach(function (r) {
      if (done[r.index]) return;
      var p = pairs[r.index];
      done[r.index] = true;
      if (p === undefined) { out.push({ row: r }); return; }
      done[p] = true;
      var partner = all[p];
      var h = r.chain === "H" ? r : partner, l = r.chain === "H" ? partner : r;
      out.push({ name: pairName(h.id) || h.id, h: h, l: l });
    });
    return out;
  }

  /**
   * What pasted text is, when it is not FASTA or a single sequence: ".tsv" (rows copied from
   * Excel: name<TAB>sequence, or one sequence per line), ".csv", ".json", or "" (let the
   * engine treat it as FASTA / one raw sequence). The app then asks the engine to read it
   * as a file of that kind, so the engine's own table and JSON readers do the work.
   */
  function pastedKind(text) {
    var t = String(text).trim();
    if (!t || t.charAt(0) === ">") return "";
    if (t.charAt(0) === "{" || t.charAt(0) === "[") return ".json";
    var lines = t.split(/\r?\n/).filter(function (l) { return l.trim(); });
    if (lines.every(function (l) { return l.indexOf("\t") >= 0; })) return ".tsv";
    if (lines.every(function (l) { return l.indexOf(",") >= 0; })) return ".csv";
    // One domain per line (a single column copied from Excel). A wrapped single sequence
    // has short lines (60-80 residues with a shorter last line), so require every line to
    // be a whole domain by itself.
    if (lines.length >= 2 && lines.every(function (l) { var n = l.replace(/\s+/g, "").length; return n >= 70 && n <= 200; })) return ".tsv";
    return "";
  }

  /* ------------------------------------------------------------------ exports */

  var TABLE_HEAD = ["순번", "서열 이름", "체인", "분류", "주석", "CDR1", "CDR2", "CDR3",
    "CDR 플래그", "framework 플래그", "도메인 플래그", "길이", "pI", "GRAVY", "net charge", "pH",
    "규칙 해당", "규칙 세트", "입력 오류"];

  function chainName(r) { return r.chain === "H" ? "VH" : r.chain === "L" ? "VL" : ""; }

  function tableRecords(rows, ruleName) {
    return rows.map(function (r) {
      var cdr = function (n) { return r.cdrs[n] ? r.cdrs[n].seq : ""; };
      return [
        r.index + 1, r.id, chainName(r),
        r.status, r.inputOk ? (r.confident ? "신뢰 가능" : "신뢰 불가") : "",
        cdr(0), cdr(1), cdr(2),
        r.inputOk && !r.regionsUnassigned ? countText(r.cdrCounts) : (r.inputOk ? "주석 실패" : ""),
        r.inputOk ? countText(r.fwCounts) : "",
        r.flags.map(function (f) { return f.label + " (" + f.risk + ")"; }).join(", "),
        r.length === null ? "" : r.length,
        r.pI === null ? "" : r.pI, r.gravy === null ? "" : r.gravy,
        r.charge === null ? "" : r.charge, r.pH === null ? "" : r.pH,
        r.reasons.join(" / "), ruleName || "", r.errors.join(" / ")
      ];
    });
  }

  /** CSV for Excel: UTF-8 with BOM (so Korean opens correctly), CRLF, RFC 4180 quoting. */
  function toCsv(rows, ruleName) {
    var q = function (v) {
      var s = String(v);
      return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    var lines = [TABLE_HEAD].concat(tableRecords(rows, ruleName)).map(function (rec) { return rec.map(q).join(","); });
    return "﻿" + lines.join("\r\n") + "\r\n";
  }

  /** Tab-separated text for pasting into Excel or Sheets. */
  function toTsv(rows, ruleName) {
    var clean = function (v) { return String(v).replace(/[\t\r\n]+/g, " "); };
    return [TABLE_HEAD].concat(tableRecords(rows, ruleName)).map(function (rec) { return rec.map(clean).join("\t"); }).join("\n");
  }

  /** Sequences as FASTA text (valid inputs only), 60 residues per line. */
  function toFasta(rows) {
    return rows.filter(function (r) { return r.inputOk; }).map(function (r) {
      return ">" + r.id + "\n" + r.sequence.replace(/(.{60})/g, "$1\n").replace(/\n$/, "") + "\n";
    }).join("");
  }

  function firstSentence(s) {
    var i = s.indexOf(". ");
    return i > 0 ? s.slice(0, i + 1) : s;
  }

  /** Plain text for an e-mail: counts, then each listed candidate with its rule hits.
      Facts only -- rule hits, not judgements. */
  function toMailText(rows, meta) {
    var t = tally(rows);
    var out = [
      "AbCheck · 후보 " + t.all + "건 · 규칙 " + meta.ruleName + " · " + meta.when,
      "검토 필요 " + t.review + " · 판정 불가 " + t.unknown + " · 기준 해당 없음 " + t.clear + " · 입력 오류 " + t.errors
    ];
    var section = function (title, list, line) {
      if (!list.length) return;
      if (meta.pairs) list = withPartnersAdjacent(list, meta.pairs);   // an antibody's two chains on adjacent lines
      out.push("", title);
      list.forEach(function (r) { out.push("- " + r.id + (r.chain ? " (" + chainName(r) + ")" : "") + ": " + line(r)); });
    };
    section("검토 필요", rows.filter(function (r) { return r.cats[0] === "review"; }), function (r) { return r.reasons.join("; "); });
    section("판정 불가", rows.filter(function (r) { return r.cats[0] === "unknown"; }), function (r) {
      return r.cdrs.length ? "주석 신뢰 불가" : "CDR 주석 실패";
    });
    section("입력 오류", rows.filter(function (r) { return r.cats[0] === "errors"; }), function (r) { return firstSentence(r.errors[0] || ""); });
    out.push("", "Kabat 근사 · ANARCI/IMGT 아님 · 종합 판정 없음");
    return out.join("\n");
  }

  /** CDR residues with liability positions underlined (positions are 1-based in the domain). */
  function markedCdr(cdr, hits, display) {
    var shown = display === undefined ? cdr.seq : String(display);
    if (shown.length !== cdr.seq.length) throw new Error("markedCdr: display text must be as long as the CDR (" + cdr.name + ")");
    var marks = new Array(cdr.seq.length).fill(null);
    (hits || []).forEach(function (x) {
      for (var p = Math.max(x.start, cdr.start); p <= Math.min(x.end, cdr.end); p++) {
        var k = p - cdr.start;
        marks[k] = (marks[k] ? marks[k] + " / " : "") + x.label + " · " + x.risk;
      }
    });
    var out = "", buf = "", cur = null;
    var flush = function () {
      if (!buf) return;
      out += cur ? '<u title="' + esc(cur) + '">' + esc(buf) + "</u>" : esc(buf);
      buf = "";
    };
    for (var k = 0; k < cdr.seq.length; k++) {
      if (marks[k] !== cur) { flush(); cur = marks[k]; }
      buf += shown[k];
    }
    flush();
    return out;
  }

  /** The one-page sheet for picked candidates: an HTML fragment (no verdict wording). */
  function onePageFragment(rows, meta) {
    var head = '<header class="op-head"><b>AbCheck · 후보 ' + rows.length + "건</b><span>규칙 " + esc(meta.ruleName) +
      " · " + esc(meta.when) + "</span></header>";
    var body = rows.map(function (r) {
      var cdrCells = [0, 1, 2].map(function (n) {
        var c = r.cdrs[n];
        return "<td class=\"op-seq\">" + (c ? "<small>" + esc(c.name) + "</small>" + markedCdr(c, r.hits) : "-") + "</td>";
      }).join("");
      var hits = r.reasons.length ? r.reasons.map(esc).join("<br>") : "-";
      var anno = !r.inputOk ? "입력 오류" : !r.cdrs.length ? "CDR 주석 실패" : r.confident ? "신뢰 가능" : "신뢰 불가";
      return "<tr><td class=\"op-id\"><b>" + esc(r.id) + "</b><small>" + esc(chainName(r) || "-") + " · " + (r.length || "-") +
        " aa</small></td>" + cdrCells + "<td>" + hits + "</td><td>" + anno + "</td><td class=\"op-num\">" +
        (r.pI === null ? "-" : r.pI.toFixed(2)) + "<br>" + (r.gravy === null ? "-" : r.gravy.toFixed(3)) + "</td></tr>";
    }).join("");
    return head + "<table class=\"op-table\"><thead><tr><th>서열</th><th>CDR1</th><th>CDR2</th><th>CDR3</th><th>규칙 해당</th>" +
      "<th>주석</th><th class=\"op-num\">pI<br>GRAVY</th></tr></thead><tbody>" + body + "</tbody></table>" +
      "<footer class=\"op-foot\">밑줄 = liability 모티프 · Kabat 근사 · ANARCI/IMGT 아님 · 종합 판정 없음 · 로컬 처리</footer>";
  }

  var ONE_PAGE_CSS =
    "@page{size:A4 landscape;margin:12mm}" +
    ".op{font-family:Pretendard,'Apple SD Gothic Neo','Malgun Gothic',system-ui,sans-serif;color:#1A1F1C;font-size:13px}" +
    ".op-head{display:flex;justify-content:space-between;align-items:baseline;border-bottom:2px solid #1A1F1C;padding-bottom:6px;margin-bottom:8px}" +
    ".op-head b{font-size:17px}.op-head span{color:#5A6469}" +
    ".op-table{width:100%;border-collapse:collapse}" +
    ".op-table th{text-align:left;font-size:13px;color:#5A6469;border-bottom:1px solid #1A1F1C;padding:6px}" +
    ".op-table td{border-bottom:1px solid #E3E8E6;padding:6px;vertical-align:top}" +
    ".op-id small,.op-seq small{display:block;color:#5A6469;font-size:13px;font-family:inherit}" +
    ".op-seq{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;word-break:break-all}" +
    ".op-seq u{text-decoration:underline 2px;text-underline-offset:3px}" +
    ".op-num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}" +
    ".op-foot{margin-top:10px;color:#5A6469}" +
    "tr{break-inside:avoid}";

  /** The same sheet as a standalone HTML document to open and print. */
  function onePageDoc(rows, meta) {
    return "<!doctype html><html lang=\"ko\"><head><meta charset=\"utf-8\"><title>AbCheck 후보 " + rows.length +
      "건</title><style>body{margin:0;padding:12mm;background:#fff}" + ONE_PAGE_CSS + "</style></head><body><div class=\"op\">" +
      onePageFragment(rows, meta) + "</div></body></html>";
  }

  return {
    RISKS: RISKS, STATUS: STATUS, MOTIFS: MOTIFS, DEFAULT_RULESET: DEFAULT_RULESET, WHERE: WHERE, CYS: CYS,
    VALIDATION: VALIDATION, TABLE_HEAD: TABLE_HEAD, ONE_PAGE_CSS: ONE_PAGE_CSS,
    esc: esc, countText: countText, normalizeRuleSet: normalizeRuleSet, describeRuleSet: describeRuleSet,
    applyRules: applyRules, withRuleResults: withRuleResults, categories: categories, buildRows: buildRows, tally: tally,
    filterRows: filterRows, sortRows: sortRows, sharedCdrs: sharedCdrs, diffCdr: diffCdr, parseSubstitutions: parseSubstitutions, applySubstitutions: applySubstitutions,
    variantName: variantName, regionAt: regionAt, diffHits: diffHits, motifTable: motifTable, withoutEntry: withoutEntry, niceTicks: niceTicks, placeLabels: placeLabels, nearestPoint: nearestPoint,
    combineEntries: combineEntries, pastedKind: pastedKind, stripNumbering: stripNumbering, pairKey: pairKey, pairName: pairName, pairMap: pairMap, antibodyUnits: antibodyUnits, preferredUnit: preferredUnit,
    cleanRecent: cleanRecent, cleanUi: cleanUi, cleanRules: cleanRules, ruleSetToText: ruleSetToText, ruleSetFromText: ruleSetFromText,
    withPartnersAdjacent: withPartnersAdjacent, toCsv: toCsv, toTsv: toTsv, toFasta: toFasta, toMailText: toMailText,
    markedCdr: markedCdr, onePageFragment: onePageFragment, onePageDoc: onePageDoc
  };
});
