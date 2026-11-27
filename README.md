# ab-seq-selfcheck

A screening tool that annotates **CDR positions** and flags **developability
liabilities** in antibody variable domain (VH/VL) sequences. Built for the
candidate-triage step: a 30-second read on whether a sequence is worth carrying
forward.

No dependencies — Python 3.9+ standard library only. It makes no network calls,
so proprietary sequences never leave the machine.

**English** · [한국어](README.ko.md)

```bash
python3 -m abcheck.cli candidates.fasta --html
```

> The tool's own output (reports, error messages) is in Korean — it was built for
> Korean research teams. Step-by-step usage docs are in Korean too:
> [user manual PDF](docs/manual/abcheck_manual.pdf) (22 pages, real outputs and
> screenshots) and [TUTORIAL.md](TUTORIAL.md).

## Validation

Measured on **7,000 sequences** from the Observed Antibody Space (OAS) public
repertoire database, across 7 datasets spanning naive B, memory B, plasma B, plasmablast and
unsorted B cells in human and mouse. The numbers are reproducible →
[`scripts/validate_on_oas.py`](scripts/validate_on_oas.py)

| Metric | Result |
|---|---|
| Chain (H/L) assignment | **0 incorrect** / 6,984 |
| CDR annotation success | 98.7% |
| **CDR boundaries** vs ANARCI ground truth, germline family and IMGT length fixed (7 datasets, 16,498 comparisons) | **0.35% deviation** |
| FR3 motif → conserved Cys distance = 29 (measured before the cross-check existed) | 3,429 / 3,443 — all 14 exceptions were motifs matched in the wrong place |

OAS ships ANARCI annotations (IMGT numbering) alongside each sequence, which makes
CDR boundaries checkable against ground truth. Because this tool approximates
**Kabat** boundaries, string equality is the wrong metric — the comparison asks
whether the *offset* between the two schemes is consistent. Within a germline family
and IMGT CDR length it is constant →
[`scripts/compare_with_anarci.py`](scripts/compare_with_anarci.py)

**Running on real data surfaced seven bugs** that the synthetic test corpus could not
have caught, because the synthetic data was built from the same assumptions as the
code. One fix was **rejected once and reversed after comparing with ground truth** —
a frequency distribution tells you what is common, not what is correct. The full
account is in [VALIDATION.md](VALIDATION.md).

## Design principles

- **Never guess.** When an anchor residue is missing, the tool reports
  *not confident* with a reason instead of inventing a position. This is why chain
  assignment has zero incorrect calls across the 6,984 sequences that passed input
  validation.
- **No overall verdict.** The tool does not produce a developability score.
  Weighting one `NG` in a CDR against two `DG` in framework has no agreed answer, and
  once a number exists, people read the number instead of the evidence.
- **Always show the evidence.** Every flag carries an explanation of why it matters.
- **It is a screening heuristic, not a numbering engine.** It does not replace
  ANARCI/IMGT, and every output says so.

### Limitations

- CDR boundaries are approximate. Insertions and deletions can shift them.
- Liability flags look at sequence motifs only. Without structure (solvent
  accessibility), a buried residue that poses no real problem is still flagged.
- pI, GRAVY and net charge are sequence-based estimates and can differ from
  measured values.

### Re-check with ANARCI/IMGT when

1. a sequence comes back `confident=False`;
2. a CDR length falls outside the typical range;
3. the sequence may be a non-standard scaffold (nanobody/VHH, scFv, bispecific linker);
4. **you are about to drop a candidate based on this output.**

## Quick start

Run from the repository root (on Windows, use `py` instead of `python3`). To use it from
any folder, run `python3 -m pip install .` once and then call `abcheck`.

Input can be a pasted sequence, FASTA, a CSV/TSV saved from Excel (`id`/`sequence`
columns), or JSON; UTF-8, UTF-8 with BOM, UTF-16 and cp949 files are all read.

**1. One sequence, printed as a table**

```bash
python3 -m abcheck.cli "YITTIETSSFHYGRVRHEIEQCALLKGGLSRAHSIWVRQALKLPYITVTKDNREFVYFQEVGQYDKATLSHQSVSESLQHSLQSSLPVVRVVYYCKLAPSGSYFGRVWGQGLGSSGGF"
```

**2. Many sequences, single-file HTML report** (`sample_output/report.html`)

```bash
python3 -m abcheck.cli data/synthetic/sequences.json --html
```

FASTA files work the same way:

```bash
python3 -m abcheck.cli my_sequences.fasta --html --out-dir results
```

**3. JSON for scripts or spreadsheets**

```bash
python3 -m abcheck.cli data/synthetic/sequences.json --json > result.json
```

| Option | Meaning | Default |
|---|---|---|
| `--html` | write an HTML report | off (terminal table) |
| `--json` | print JSON to stdout | off |
| `--out-dir` / `--out-name` | where to write the HTML report | `sample_output` / `report.html` |
| `--ph` | pH for net charge | 7.4 |

Exit codes: `0` processed (sequences flagged *not confident* still count as
processed), `1` some input failed validation, `2` input could not be read.

## Browser app

`web/abcheck.html` is the same analysis as one self-contained page. Open the file in a
browser and paste sequences or drop files; the page contains no network call (a test
checks this), so sequences stay on the machine.

- **Same engine.** `web/abcheck-engine.js` is a JavaScript port of `abcheck/`.
  `scripts/js_parity.py` runs both on the synthetic corpus and on edge cases and compares
  every output field; `tests/test_js_parity.py` fails on any difference (it needs `node`
  and is skipped without it).
- **Around the engine** (`web/app-lib.js`, tested through `node` in `tests/test_web_app.py`):
  rows pasted from a spreadsheet, several files at once, triage rule sets kept in the
  browser (the default set reproduces the engine's triage; a set can be handed to a
  colleague as one line of text), VH/VL partner chains found by name and shown as one
  antibody with each chain's own result, CDR comparison against a reference candidate
  (a dot where a residue is the same; no alignment is made), groups of candidates whose
  CDR3 have the same length and differ in a few residues, a table of how the motifs
  spread over the panel (candidates per motif and region; a number opens those
  candidates), a substitution trial in the candidate drawer (`N56Q` is applied and run
  through the same engine and rule set, shown beside the parent, and can be added to the
  list; places are 1-based places in the sequence as entered, not Kabat numbers), a note
  per candidate, a pI × GRAVY scatter, a one-page sheet for picked candidates, and
  exports (table for a spreadsheet, mail text, FASTA, CSV, JSON; picks and notes travel
  with the mail text, the sheet and the tables). With a reference set, every residue in
  which another candidate of that chain differs is listed over the whole domain, piece
  by piece (FR1 · CDR1 · … · FR4, cut at the engine's CDR borders; a piece of another
  length is named and left out, no alignment is made), and that line goes out with the
  sheet, the mail text and the tables. The pH can be changed for the list on screen
  (only net charge depends on it). One line per detected motif can be copied for a
  spreadsheet. The JSON result carries the review (picks, notes, references under a
  key `review`, which also keeps the text of an input that was rejected, so the file can
  be opened again; without marks and rejected inputs the file equals the CLI's `--json`
  output) and gives it back when it is opened, so a review can be handed to a colleague
  or another PC.
- **Shell.** The layout values are shared with two sibling tools and measured by one
  script (`series-check.mjs`, kept outside this repository); the only addresses in the
  page are two links to those tools and the link to `VALIDATION.md`, each followed only
  when the viewer presses it.
- **Build.** `python3 web/build.py` inlines `web/app.html`, `web/app-lib.js`, the engine
  and the bundled fonts into `web/abcheck.html`. `python3 web/build.py --check` fails when
  the committed page is stale; the test suite runs the same check.
- [CHANGELOG.md](CHANGELOG.md) lists each version of the page with the steps a task took
  before and after.

## What it reports

### CDR annotation

Conserved-residue anchors approximate Kabat boundaries. Chain type is inferred from
the terminal motif (`WGxG` → heavy, `FGxG` → light). Coordinates are 1-based and
inclusive.

| Chain | Search order |
|---|---|
| VH | `Cys22` → **CDR-H1** (from Cys22+9) → `Trp36` (followed by `[VIA][RKQ]Q`) → **CDR-H2** (from Trp36+14) → FR3 motif `[RKQ].[TSAV][IMLFV][STAN]` → `Cys92` (preceded by `Y/F`) → **CDR-H3** (from Cys92+3) → `WGxG` |
| VL | `Cys23` → **CDR-L1** (from Cys23+1) → `Trp35` (followed by `[YF][QL]Q`) → **CDR-L2** (from Trp35+15, 7 residues) → `Cys88` (preceded by `Y/F`) → **CDR-L3** (from Cys88+1) → `FGxG` |

If the FR3 motif is missing (typically broken by somatic hypermutation), the
boundary is back-calculated 29 residues from the conserved Cys, and the result is
marked *not confident*. The same distance cross-checks every motif match: in real VH
it was 29 in 3,429 of 3,443 cases, and all 14 exceptions were motifs matched in the
wrong place, which the annotator now re-anchors.

A CDR length outside the typical range keeps its position but drops confidence.
The range answers "can this annotation be trusted?", not "is this antibody unusual?".

| CDR | Typical length |
|---|---|
| CDR-H1 | 5-7 |
| CDR-H2 | 16-20 |
| CDR-H3 | 3-30 |
| CDR-L1 | 7-17 |
| CDR-L2 | 7 |
| CDR-L3 | 5-15 |

### Liability motifs

Risk depends on **where** a motif sits: inside a CDR it is escalated one level, in
framework it is de-escalated — the same `DG` is *high* in CDR-H3 and *medium* in
framework.

| Liability | Motif | Base risk |
|---|---|---|
| N-glycosylation | `N[^P][ST]` | high |
| Deamidation | `NG` | high |
| Deamidation | `NS` `NT` `NN` `NH` | medium |
| Isomerization | `DG` | high |
| Isomerization | `DS` `DT` `DD` | medium |
| Fragmentation | `DP` | medium |
| Oxidation | `M` `W` | low |
| Unpaired Cys | odd Cys count in the domain | domain-level flag |

### Physicochemical properties

pI (EMBOSS pKa, bisection), GRAVY (Kyte-Doolittle) and net charge at pH 7.4, for the
whole domain and per CDR. Per-CDR pI treats the fragment as a free peptide — use it
for relative comparison between regions only.

### Triage

User-defined rules mark sequences as *needs review*. Defaults: a high-risk,
top-of-class motif (`N[^P][ST]`, `NG`, `DG`, `DP`) inside a CDR, or a high-risk
unpaired Cys. Results are *needs review* / *no rule matched* / *cannot evaluate*
(input or annotation failed, or no rule matched but the CDR boundaries are not confident —
never silently treated as clean). Rules live in
`abcheck/triage.py` (`DEFAULT_RULES`).

### Input validation

Rejected with a reason: non-amino-acid characters, ambiguous codes (`X B Z J U O`),
embedded FASTA headers, lengths outside 70-200 residues (split scFv/Fab into VH and
VL), and DNA/RNA pasted by mistake. Whitespace and lowercase are normalized.

## Synthetic data

All sequence data in this repository is **synthetic** — no therapeutic antibody or
germline database records. The framework is weighted random residues with the
conserved anchors fixed in place; CDRs are inserted at known positions, so every
sequence carries ground-truth CDR coordinates. See
[data/synthetic/README.md](data/synthetic/README.md).

```bash
python3 scripts/make_synthetic_sequences.py
```

## Tests and CI

```bash
python3 -m pytest -q
```

`.github/workflows/ci.yml` runs on pull requests into `main` and pushes to `main`:

| Job | Checks |
|---|---|
| `pytest (Python 3.9–3.13)` | the full test suite, including that the committed sample report and user manual match current output |
| `Works on windows-latest` / `macos-latest` | the full test suite on real Windows and macOS, then `pip install .` and `abcheck` from another folder |
| `No runtime dependencies` | `pyproject.toml` still declares no runtime dependencies |
| `Synthetic corpus is reproducible` | the committed corpus equals the generator output |
| `CLI smoke test` | the single-sequence, HTML and JSON examples run on the synthetic corpus; invalid input exits with exactly 1 |

To make failing checks block merges, enable branch protection on `main` and require
the ten checks above (`pytest (Python 3.9)` … `pytest (Python 3.13)`,
`Works on windows-latest`, `Works on macos-latest`, `No runtime dependencies`,
`Synthetic corpus is reproducible`, `CLI smoke test`).
The exact `gh api` command is in [README.ko.md](README.ko.md#ci-와-차단).

A local pre-push hook runs the tests and the corpus check before every push:

```bash
git config core.hooksPath .githooks
```

## Repository layout

```
abcheck/          validate · annotate · liabilities · properties · triage · report · render · cli
scripts/          synthetic corpus generator · OAS → FASTA converter · two validation scripts
data/synthetic/   24 synthetic sequences with ground-truth CDR coordinates
sample_output/    example HTML report
docs/manual/      user manual PDF (Korean) and the script that builds it from live output
tests/            test suite
web/              browser app: app.html + app-lib.js + abcheck-engine.js → build.py → abcheck.html
CHANGELOG.md      browser app versions: what changed, why, steps before → after
VALIDATION.md     measurement method, results, known limitations
DECISIONS.md      design decisions and their rationale
README.ko.md      Korean README
TUTORIAL.md       step-by-step guide (Korean)
```

## Out of scope

Structure prediction, docking, affinity prediction. No external database or API calls.

## License

MIT — see [LICENSE](LICENSE).
