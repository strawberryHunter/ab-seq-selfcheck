"""Check that the browser engine (web/abcheck-engine.js) gives exactly Python's output.

    python3 scripts/js_parity.py                  # synthetic corpus + edge cases + helpers
    python3 scripts/js_parity.py extra.fasta ...  # ... plus every record of these FASTA files

Runs the Python package and the JavaScript engine (through `node`) on the same inputs and
compares the JSON text that `abcheck.render.as_json` writes, byte for byte, with only the
`generated` date ignored. When the text differs, both sides are parsed and the differing
values are printed with their JSON path; int and float count as different types (7 is not
7.0). Also compares `ABCHECK.readInput` with `abcheck.cli.load_input` and the Python
number helpers (round, float repr) on many values.

Exit status: 0 when there are no differences, 1 otherwise, 2 when node is missing.
Needs only the standard library and a `node` executable (any version from 12 on).
"""

from __future__ import annotations

import csv
import json
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from abcheck import cli  # noqa: E402
from abcheck.annotate import annotate  # noqa: E402
from abcheck.render import as_json  # noqa: E402
from abcheck.report import build_many  # noqa: E402

ENGINE = ROOT / "web" / "abcheck-engine.js"
CORPUS = ROOT / "data" / "synthetic" / "sequences.json"

#: FASTA files are compared in chunks so a 1,000-record file does not become one 100 MB JSON.
CHUNK = 250
#: How many differences to print before stopping.
MAX_SHOWN = 20

# Driver run by node: reads requests from stdin, answers with one JSON document.
NODE_DRIVER = r"""
const A = require(process.argv[1]);
let input = "";
process.stdin.setEncoding("utf8");
process.stdin.on("data", (d) => { input += d; });
process.stdin.on("end", () => {
  const out = JSON.parse(input).map((req) => {
    try {
      if (req.op === "analyze") return { ok: true, text: A.analyzeJSON(req.entries, req.ph) };
      if (req.op === "readInput") {
        const r = A.readInput(req.text, req.filename);
        return { ok: true, entries: r[0], source: r[1] };
      }
      if (req.op === "round") return { ok: true, values: req.values.map((v) => A._py.floatRepr(A._py.round(v[0], v[1]))) };
      if (req.op === "repr") return { ok: true, values: req.values.map((v) => A._py.floatRepr(v)) };
      if (req.op === "reprStr") return { ok: true, values: req.values.map((v) => A._py.reprStr(v)) };
      if (req.op === "sum") return { ok: true, values: req.values.map((v) => A._py.floatRepr(A._py.sumFloats(v))) };
      return { ok: false, error: "unknown op " + req.op };
    } catch (e) {
      return { ok: false, error: e.message };
    }
  });
  process.stdout.write(JSON.stringify(out));
});
"""


def node_executable() -> str | None:
    return shutil.which("node")


def run_node(requests: list[dict]) -> list[dict]:
    """Send requests to the JS engine and return its answers (one per request)."""
    node = node_executable()
    if node is None:
        raise RuntimeError("node 가 PATH 에 없다")
    proc = subprocess.run(
        [node, "-e", NODE_DRIVER, str(ENGINE)],
        input=json.dumps(requests, ensure_ascii=True).encode("ascii"),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("node failed:\n" + proc.stderr.decode("utf-8", "replace"))
    return json.loads(proc.stdout.decode("utf-8"))


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

def corpus_entries() -> list[tuple[str, str]]:
    data = json.loads(CORPUS.read_text(encoding="utf-8"))
    return [(r["id"], r["sequence"]) for r in data["sequences"]]


def _rescue(seq: str, motif: str, replacement: str, tail: str) -> str:
    """Replace the last `motif` and everything after it with `replacement` + `tail`."""
    idx = seq.rfind(motif)
    assert idx > 0, motif
    return seq[:idx] + replacement + tail


def _anchor_cases(vh: str, vl: str) -> list[tuple[str, str]]:
    """Sequences cut or edited around the anchors Python finds, one per annotation branch."""
    h = annotate(vh).anchors
    l = annotate(vl).anchors
    hc2, hm = h["cys2"], h["terminal_motif"]
    lt, lc2, lm = l["trp"], l["cys2"], l["terminal_motif"]
    after_fr2 = vh[h["trp"] + 3:]
    return [
        # heavy
        ("H-trp-without-fr2", vh[:h["trp"]] + "E" + vh[h["trp"] + 1:]),
        ("H-no-trp", vh[:hm - 1].replace("W", "F") + vh[hm - 1:]),
        ("H-h2-off-end", "A" * 21 + "C" + "A" * 60 + "W" + "AAAA" + "WGQG" + "AA"),
        ("H-no-cys92", vh[:30] + vh[30:].replace("C", "S")),
        ("H-cys92-without-context", vh[:hc2 - 2] + "A" + vh[hc2 - 1:]),
        ("H-h3-empty", vh[:hc2] + vh[hm - 1:]),
        ("H-h3-short", vh[:hc2 + 4] + vh[hm - 1:]),
        ("H-h3-long", vh[:hc2 + 3] + "GSYAG" * 7 + vh[hm - 1:]),
        ("H-fr2-variant", vh[:h["trp"] + 3] + "K" + after_fr2[1:]),
        # light
        ("L-trp-without-fr2", vl[:lt + 1] + "E" + vl[lt + 2:]),
        ("L-l2-past-motif", vl[:lt + 5] + vl[lm - 1:] + "A" * 40),
        ("L-cys88-without-context", vl[:lc2 - 2] + "A" + vl[lc2 - 1:]),
        ("L-l3-empty", vl[:lc2] + vl[lm - 1:]),
        ("L-l3-short", vl[:lc2 + 2] + vl[lm - 1:]),
        ("L-l1-long", vl[:l["cys1"] + 3] + "SGSGS" * 3 + vl[l["cys1"] + 3:]),
    ]


def edge_cases() -> list[tuple[str, str]]:
    """Inputs that exercise the validation messages and the annotation fallbacks."""
    seqs = dict(corpus_entries())
    vh, vl = seqs["SYN-VH-001"], seqs["SYN-VL-001"]
    vh2 = seqs["SYN-VH-002"]
    cases = [
        ("empty", ""),
        ("blank", "  \n\t "),
        ("unicode-blank", "　\x1c\x85"),
        ("too-short", vh[:40]),
        ("too-long", vh + vh),
        ("dna", "ATGGCCGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGG"),
        ("rna-short", "ACGUACGUACGU"),
        ("poly-a", "A" * 80),
        ("ambiguous", vh[:50] + "XBZ" + vh[53:]),
        ("ambiguous-and-bad", vh[:50] + "J1#" + vh[53:]),
        ("fasta-header", ">SYN-VH-001 some description\n" + vh),
        ("fasta-header-long", "  >" + "H" * 60 + "\n" + vh),
        ("stop", vh[:60] + "*" + vh[60:]),
        ("stop-end", vh + "*"),
        ("digits", "1 " + vh + " 99"),
        ("bom", "﻿" + vh),
        ("nbsp", vh[:30] + "\xa0" + vh[30:]),
        ("quote", vh[:30] + "'" + vh[30:]),
        ("dquote-quote", vh[:30] + "'\"\\" + vh[30:]),
        ("emoji", vh[:30] + "\U0001F600" + vh[30:]),
        ("tag-char", vh[:30] + "\U000E0001" + vh[30:]),
        ("many-bad", "".join(c + "%" for c in vh[:12]) + vh),
        ("gaps-dots", vh[:26] + "...." + vh[26:60] + ".." + vh[60:]),
        ("gaps-dash", vh[:26] + "--" + vh[26:]),
        ("mostly-gaps", "." * 150 + vh[:20]),
        ("all-gaps", "-" * 30),
        ("lower-ws", "  " + "\n".join(vl[i:i + 10].lower() for i in range(0, len(vl), 10)) + "\r\n"),
        ("crlf-tabs", vh[:50] + "\r\n\t" + vh[50:] + "\x0b\x0c"),
        ("sharp-s", vh[:40] + "ß" + vh[40:]),
        ("vh-j-tail", _rescue(vh, "WGQG", "WAQG", "TLVTVSS")),
        ("vh-j-tail-wsqg", _rescue(vh2, "WGQG", "WSQG", "TLVTVSS")),
        ("vl-j-tail-lambda", _rescue(vl, "FGQG", "LGGG", "TKLTVL")),
        ("vl-j-tail-kappa", _rescue(vl, "FGQG", "FAQG", "TRLEIK")),
        ("no-motif", vh.replace("WGQG", "WAQA")),
        ("vh-no-cys22", vh[:21] + "S" + vh[22:]),
        ("vh-truncated-n", vh[30:]),
        ("vh-truncated-c", vh[:100]),
        ("vh-fr3-broken", vh.replace("RV", "GG", 1)),
        ("vh-extra-cys", vh[:60] + "C" + vh[61:]),
        ("vh-cdr-cys", vh[:100] + "C" + vh[101:]),
        ("vh-two-cys-window", vh[:17] + "C" + vh[18:]),
        ("vl-no-trp", vl.replace("W", "F")),
        ("vl-no-cys88", vl[:85] + vl[85:].replace("C", "S", 1)),
        ("vh-motif-early", vh[:20] + "WGQG" + vh[24:].replace("WGQG", "WAQA")),
        ("vh-and-vl", vh + vl[:40]),
        ("glyco", vh[:52] + "NAS" + vh[55:]),
        ("deamid", vh[:52] + "NGNS" + vh[56:]),
        ("isomer", vh[:52] + "DGDSDP" + vh[58:]),
    ]
    cases += _anchor_cases(vh, vl)
    cases += _anchor_cases(vh2, seqs["SYN-VL-002"])
    # Every corpus record with a single residue substitution near the middle and the ends,
    # so the anchor searches see shifted contexts.
    rng = random.Random(20240601)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    for seq_id, seq in corpus_entries():
        for _ in range(3):
            pos = rng.randrange(len(seq))
            cases.append((f"{seq_id}-mut{pos}", seq[:pos] + rng.choice(aa) + seq[pos + 1:]))
    return cases


def fasta_entries(path: Path) -> list[tuple[str, str]]:
    return cli.read_fasta(cli._read_text(path))


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------

_GENERATED = re.compile(r'^  "generated": "[^"]*",$', re.MULTILINE)


def _normalise(text: str) -> str:
    return _GENERATED.sub('  "generated": "",', text, count=1)


def _type_name(v) -> str:
    return type(v).__name__


def diff_values(a, b, path: str = "$") -> list[str]:
    """Differences between two parsed JSON values, int and float being distinct types."""
    if type(a) is not type(b):
        return [f"{path}: python={a!r} ({_type_name(a)}) js={b!r} ({_type_name(b)})"]
    if isinstance(a, dict):
        out = []
        if list(a) != list(b):
            out.append(f"{path}: key order/set python={list(a)} js={list(b)}")
        for k in a:
            if k in b:
                out.extend(diff_values(a[k], b[k], f"{path}[{json.dumps(k, ensure_ascii=False)}]"))
        return out
    if isinstance(a, list):
        out = []
        if len(a) != len(b):
            out.append(f"{path}: length python={len(a)} js={len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            out.extend(diff_values(x, y, f"{path}[{i}]"))
        return out
    if a != b and not (a != a and b != b):  # NaN == NaN here
        return [f"{path}: python={a!r} js={b!r}"]
    return []


def compare_text(py_text: str, js_text: str, label: str) -> list[str]:
    py_n, js_n = _normalise(py_text), _normalise(js_text)
    if py_n == js_n:
        return []
    diffs = [f"[{label}] {d}" for d in diff_values(json.loads(py_n), json.loads(js_n))]
    if not diffs:  # same values, different text: show the first differing line
        for i, (x, y) in enumerate(zip(py_n.splitlines(), js_n.splitlines()), start=1):
            if x != y:
                diffs.append(f"[{label}] line {i}: python={x!r} js={y!r}")
                break
        else:
            diffs.append(f"[{label}] text length python={len(py_n)} js={len(js_n)}")
    return diffs


def analyze_cases(cases: list[tuple[str, list[tuple[str, str]], float]]) -> tuple[int, list[str]]:
    """Run (label, entries, ph) cases on both sides. Returns (#sequences, differences)."""
    answers = run_node([{"op": "analyze", "entries": e, "ph": ph} for _, e, ph in cases])
    diffs: list[str] = []
    n = 0
    for (label, entries, ph), ans in zip(cases, answers):
        n += len(entries)
        py_text = as_json(build_many(entries, ph))
        if not ans["ok"]:
            diffs.append(f"[{label}] JS raised: {ans['error']}")
            continue
        diffs.extend(compare_text(py_text, ans["text"], label))
    return n, diffs


# --------------------------------------------------------------------------
# readInput
# --------------------------------------------------------------------------

def py_read_input(text: str, filename: str | None = None) -> dict:
    """What abcheck.cli.load_input returns for this text (saved as `filename` if given)."""
    try:
        if filename:
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / filename
                path.write_bytes(text.encode("utf-8"))
                entries, source = cli.load_input(str(path))
                source = source.replace(str(path), filename)
        else:
            # The pasted-text branch of load_input (its "is this a missing file?" guess
            # needs a file system and does not apply in a browser).
            if text.lstrip().startswith(">"):
                entries, source = cli.read_fasta(text), "붙여넣은 FASTA"
            else:
                entries, source = [("입력 서열", text)], "붙여넣은 서열"
        return {"ok": True, "entries": [list(e) for e in entries], "source": source}
    except (ValueError, csv.Error) as exc:
        return {"ok": False, "error": str(exc)}


def read_input_cases() -> list[tuple[str, str | None]]:
    vh, vl = dict(corpus_entries())["SYN-VH-001"], dict(corpus_entries())["SYN-VL-001"]
    fasta = f">SYN-VH-001 heavy\n{vh[:60]}\n{vh[60:]}\n\n>SYN-VL-001\r\n{vl}\r\n>\n{vl}\n"
    cases: list[tuple[str, str | None]] = [
        (fasta, "후보.fasta"),
        (fasta, "x.FA"),
        (fasta, None),
        ("  " + fasta, None),
        ("﻿" + fasta, "bom.fasta"),
        (fasta, "noext"),
        (" >a b " + vh + "\x85>b\x1c" + vl, "u.fas"),
        (vh, None),
        (vh, "seq.txt"),
        (vh, "plain"),
        (f"이름,서열\nA1,{vh}\nA2,{vl}\n", "표.csv"),
        (f"﻿후보,아미노산 서열,비고\r\n  B1 ,{vh},x\r\n,{vl},\r\n\r\n", "k.CSV"),
        (f"id\tsequence\nT1\t{vh}\n", "t.tsv"),
        (f"Sequence;ID\n{vh};X\n", "semi.csv"),
        (f"a1,{vh}\na2,{vl}\n", "two.csv"),
        (f"{vh}\n{vl}\n", "one.csv"),
        (f"a,b,c\n1,{vh},3\n", "three.csv"),
        (f'id,sequence\n"q,1","{vh}"\n"say ""hi""",{vl}\n"multi\nline",{vh}\n', "quoted.csv"),
        (f'id,sequence\nx,"{vh}"tail\n', "quote-tail.csv"),
        (f'id,sequence\nx,"{vh}', "unterminated.csv"),
        ("id,sequence\rx,AAAA\r", "cr.csv"),
        ("", "empty.csv"),
        ("\n \n", "blank.csv"),
        (json.dumps({"sequences": [{"id": "J1", "sequence": vh}, {"sequence": vl}]}), "c.json"),
        (json.dumps([vh, {"id": 7, "sequence": vl}, {"id": None, "sequence": 12}]), "list.json"),
        ('[{"id": 1.0, "sequence": "A"}, {"id": 1e-7, "sequence": true},'
         ' {"id": [1, "x", {"k": null}], "sequence": -0}]', "types.json"),
        ('{"sequences": [{"id": "a", "id": "b", "sequence": "Q"}]}', "dup.json"),
        ('{"other": 1}', "dict.json"),
        ('"just a string"', "str.json"),
        ('[{"id": "a"}]', "nokey.json"),
        ('[1, 2]', "ints.json"),
        ('{"sequences": [\n  {"id": "a", "sequence": "Q"},\n]}', "trailing.json"),
        ('{"sequences": [\n  {"id": "a" "sequence": "Q"}\n]}', "comma.json"),
        ('[{"id": "a\\x"}]', "escape.json"),
        ('[{"id": "\\ud83d\\ude00\\u00e9", "sequence": "Q"}]', "uescape.json"),
        ('[{"id": "a\tb", "sequence": "Q"}]', "control.json"),
        ("[NaN, Infinity]", "nan.json"),
        ("", "empty.json"),
        ("  \n  ", "ws.json"),
        ('[1] x', "extra.json"),
        ("{'a': 1}", "single.json"),
        ('[{"id": "a", "sequence": "Q"}', "open.json"),
        ("x", "book.xlsx"),
    ]
    return cases


def fuzz_read_input_cases(n: int = 400, seed: int = 11) -> list[tuple[str, str | None]]:
    """Valid JSON/CSV/FASTA texts with random characters deleted, inserted or replaced, to
    compare error messages (JSON line/column, csv errors) as well as successful parses."""
    rng = random.Random(seed)
    base = read_input_cases()[:30]
    noise = '",{}[]:\n\r\t \\>x1e-.'
    out = []
    for _ in range(n):
        text, filename = rng.choice(base)
        chars = list(text)
        for _ in range(rng.randint(1, 3)):
            pos = rng.randrange(len(chars) + 1)
            op = rng.random()
            if op < 0.4 and pos < len(chars):
                del chars[pos]
            elif op < 0.7:
                chars.insert(pos, rng.choice(noise))
            elif pos < len(chars):
                chars[pos] = rng.choice(noise)
        out.append(("".join(chars), filename))
    return out


def read_input_diffs(cases: list[tuple[str, str | None]]) -> list[str]:
    answers = run_node([{"op": "readInput", "text": t, "filename": f} for t, f in cases])
    diffs = []
    for (text, filename), js in zip(cases, answers):
        py = py_read_input(text, filename)
        if py != js:
            diffs.append(f"[readInput {filename!r} {text[:30]!r}] python={py!r} js={js!r}")
    return diffs


# --------------------------------------------------------------------------
# Number helpers
# --------------------------------------------------------------------------

def helper_values(rng: random.Random) -> tuple[list, list, list]:
    """Values for round(), repr(float) and sum(): random, exact ties and awkward ones."""
    floats = [rng.uniform(-20, 20) for _ in range(20000)]
    floats += [rng.uniform(-1, 1) * 10 ** rng.randint(-30, 30) for _ in range(5000)]
    floats += [k / 8 for k in range(-200, 200)] + [k / 1000 + 0.0005 for k in range(-3000, 3000)]
    floats += [k / 200 for k in range(-3000, 3000)] + [0.0, -0.0, 1e16, 1e15, 1e-4, 1e-5,
                                                      2.675, 1.005, 0.125, 0.375, 2.5, 3.5]
    rounds = [[x, n] for x in floats for n in (2, 3)] + [[x, 0] for x in floats[:2000]]
    sums = [[rng.choice([1.8, -4.5, -3.5, 2.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8,
                         -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]) for _ in range(rng.randint(1, 130))]
            for _ in range(3000)]
    return rounds, floats, sums


def helper_diffs() -> tuple[int, list[str]]:
    rounds, floats, sums = helper_values(random.Random(7))
    strs = ["*", "'", '"', "'\"", "\\", "﻿", "\xa0", "\x7f", "\x00", "​", "é",
            "\U0001F600", "\U000E0001", "　", "한", "\ud800", " ", " ", "\x85"]
    answers = run_node([
        {"op": "round", "values": rounds},
        {"op": "repr", "values": floats},
        {"op": "sum", "values": sums},
        {"op": "reprStr", "values": strs},
    ])
    expected = [
        [repr(round(x, n)) for x, n in rounds],
        [repr(x) for x in floats],
        [repr(sum(v)) for v in sums],
        [repr(s) for s in strs],
    ]
    names = ["round", "float repr", "sum", "str repr"]
    inputs = [rounds, floats, sums, strs]
    diffs = []
    for name, inp, exp, ans in zip(names, inputs, expected, answers):
        for i, (e, got) in enumerate(zip(exp, ans["values"])):
            if e != got:
                diffs.append(f"[{name}] input={inp[i]!r} python={e} js={got}")
    return sum(len(x) for x in inputs), diffs


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def standard_cases() -> list[tuple[str, list[tuple[str, str]], float]]:
    corpus, edges = corpus_entries(), edge_cases()
    cases = [("corpus pH 7.4", corpus, 7.4), ("edge cases pH 7.4", edges, 7.4)]
    for ph in (5.0, 9.0):
        cases.append((f"corpus pH {ph}", corpus, ph))
        cases.append((f"edge cases pH {ph}", edges, ph))
    # each corpus record alone too, so a summary of one sequence is covered
    cases += [(f"single {sid}", [(sid, seq)], 7.4) for sid, seq in corpus]
    return cases


def fasta_cases(paths: list[Path]) -> list[tuple[str, list[tuple[str, str]], float]]:
    cases = []
    for path in paths:
        entries = fasta_entries(path)
        for i in range(0, len(entries), CHUNK):
            cases.append((f"{path.name}[{i}:{i + CHUNK}]", entries[i:i + CHUNK], 7.4))
    return cases


def main(argv: list[str]) -> int:
    if node_executable() is None:
        print("node 가 PATH 에 없다. Node.js 를 설치한 뒤 다시 실행하라.", file=sys.stderr)
        return 2

    total_diffs: list[str] = []

    n_helpers, diffs = helper_diffs()
    print(f"helpers (round/repr/sum): {n_helpers} values, {len(diffs)} differences")
    total_diffs += diffs

    ri_cases = read_input_cases() + fuzz_read_input_cases()
    diffs = read_input_diffs(ri_cases)
    print(f"readInput: {len(ri_cases)} inputs, {len(diffs)} differences")
    total_diffs += diffs

    n, diffs = analyze_cases(standard_cases())
    print(f"synthetic corpus + edge cases (pH 7.4/5.0/9.0): {n} sequences, {len(diffs)} differences")
    total_diffs += diffs

    paths = [Path(a) for a in argv]
    if paths:
        cases = fasta_cases(paths)
        n_fasta, diffs = 0, []
        for i in range(0, len(cases), 8):  # a few chunks per node run keeps memory small
            k, d = analyze_cases(cases[i:i + 8])
            n_fasta += k
            diffs += d
        print(f"FASTA files ({len(paths)}): {n_fasta} sequences, {len(diffs)} differences")
        total_diffs += diffs

    for line in total_diffs[:MAX_SHOWN]:
        print("  " + line)
    if len(total_diffs) > MAX_SHOWN:
        print(f"  ... and {len(total_diffs) - MAX_SHOWN} more")
    print(f"total: {len(total_diffs)} differences")
    return 1 if total_diffs else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
