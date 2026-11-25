# ab-seq-selfcheck

A screening tool that annotates **CDR positions** and flags **developability
liabilities** in antibody variable domain (VH/VL) sequences. Built for the
candidate-triage step: a 30-second read on whether a sequence is worth carrying
forward.

No dependencies — Python 3.10+ standard library only. It makes no network calls,
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
repertoire database, across 7 datasets spanning naive B, memory B, plasma B and
plasmablast cells in human and mouse. The numbers are reproducible →
[`scripts/validate_on_oas.py`](scripts/validate_on_oas.py)

| Metric | Result |
|---|---|
| Chain (H/L) assignment | **0 incorrect** / 6,984 |
| CDR annotation success | 98.3% |
| **CDR boundaries** vs ANARCI ground truth (7 datasets, 16,427 comparisons) | **0.35% deviation** |
| Internal consistency (FR3 ↔ conserved Cys distance) | 100% |

OAS ships ANARCI annotations (IMGT numbering) alongside each sequence, which makes
CDR boundaries checkable against ground truth. Because this tool approximates
**Kabat** boundaries, string equality is the wrong metric — the comparison asks
whether the *offset* between the two schemes is consistent. Within a germline family
and IMGT CDR length it is constant →
[`scripts/compare_with_anarci.py`](scripts/compare_with_anarci.py)

**Running on real data surfaced six bugs** that the synthetic test corpus could not
have caught, because the synthetic data was built from the same assumptions as the
code. One fix was **rejected once and reversed after comparing with ground truth** —
a frequency distribution tells you what is common, not what is correct. The full
account is in [VALIDATION.md](VALIDATION.md).

## Design principles

- **Never guess.** When an anchor residue is missing, the tool reports
  *not confident* with a reason instead of inventing a position. This is why chain
  assignment has zero incorrect calls across 7,000 sequences.
- **No overall verdict.** The tool does not produce a developability score.
  Weighting one `NG` in a CDR against two `DG` in framework has no agreed answer, and
  once a number exists, people read the number instead of the evidence.
- **Always show the evidence.** Every flag carries an explanation of why it matters.
- **It is a screening heuristic, not a numbering engine.** It does not replace
  ANARCI/IMGT, and every output says so.

### Re-check with ANARCI/IMGT when

1. a sequence comes back `confident=False`;
2. a CDR length falls outside the typical range;
3. the sequence may be a non-standard scaffold (nanobody/VHH, scFv, bispecific linker);
4. **you are about to drop a candidate based on this output.**

## Quick start

Run from the repository root.

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
boundary is back-calculated 29 residues from the conserved Cys — a distance that was
constant across 3,443 real VH sequences — and the result is marked *not confident*.

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
(annotation failed — never silently treated as clean). Rules live in
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

`.github/workflows/ci.yml` runs on pull requests and pushes to `main`:

| Job | Checks |
|---|---|
| `pytest (Python 3.10–3.13)` | the full test suite, including that the committed sample report and user manual match current output |
| `No runtime dependencies` | `pyproject.toml` still declares no runtime dependencies |
| `Synthetic corpus is reproducible` | the committed corpus equals the generator output |
| `CLI smoke test` | the three Quick start examples run; invalid input exits with 1 |

To make failing checks block merges, enable branch protection on `main` and require
the seven checks above (`pytest (Python 3.10)` … `pytest (Python 3.13)`,
`No runtime dependencies`, `Synthetic corpus is reproducible`, `CLI smoke test`).
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
VALIDATION.md     measurement method, results, known limitations
DECISIONS.md      design decisions and their rationale
README.ko.md      Korean README
TUTORIAL.md       step-by-step guide (Korean)
```

## Out of scope

Structure prediction, docking, affinity prediction. No external database or API calls.

## License

MIT — see [LICENSE](LICENSE).
