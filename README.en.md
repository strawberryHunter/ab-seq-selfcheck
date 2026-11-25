# ab-seq-selfcheck

A screening tool that annotates **CDR positions** and flags **developability
liabilities** in antibody variable domain (VH/VL) sequences. Built for the
candidate-triage step: a 30-second read on whether a sequence is worth carrying
forward.

No dependencies — Python 3.10+ standard library only. It makes no network
calls, so proprietary sequences never leave the machine.

*The full documentation is in Korean. This page is a summary.*

```bash
python3 -m abcheck.cli candidates.fasta --html
```

## What it does

1. **CDR annotation** — locates CDR1/2/3 using conserved-residue anchors
   (Cys22/Cys92, Trp36, the FR3 motif, the terminal `WGxG`/`FGxG`), approximating
   Kabat boundaries. Chain type (H/L) is inferred from the terminal motif.
2. **Liability motifs** — N-glycosylation `N[^P][ST]`, deamidation `NG` > `NS/NT/NN/NH`,
   isomerization `DG` > `DS/DT/DD`, fragmentation `DP`, oxidation `M`/`W`, unpaired Cys.
   **Risk is adjusted by region**: a motif inside a CDR is escalated one level, in
   framework it is de-escalated — the same `DG` is *high* in CDR-H3 and *medium* in framework.
3. **Physicochemical properties** — pI (EMBOSS pKa, bisection), GRAVY (Kyte-Doolittle),
   net charge at pH 7.4; computed for the whole domain and per CDR.
4. **Triage** — user-defined rules mark sequences as *needs review*. The tool applies
   the rules; it does not invent them.
5. **Output** — terminal table, single-file HTML report, or JSON.

## Validation

Measured on **7,000 sequences** from the Observed Antibody Space (OAS) public
repertoire database, across 7 datasets spanning naive B, memory B, plasma B and
plasmablast cells in human and mouse.

| Metric | Result |
|---|---|
| Chain (H/L) assignment | **0 incorrect** / 6,984 |
| CDR annotation success | 98.3% |
| **CDR boundaries** vs ANARCI ground truth (16,427 comparisons) | **0.35% deviation** |
| Internal consistency (FR3↔conserved Cys distance) | 100% |

Both measurements are reproducible:

```bash
python3 scripts/oas_to_fasta.py unit.csv.gz candidates.fasta --count 500
python3 scripts/validate_on_oas.py candidates.fasta
python3 scripts/compare_with_anarci.py unit.csv.gz --count 500
```

OAS ships ANARCI annotations (IMGT numbering) alongside each sequence, which makes
CDR boundaries verifiable against ground truth. Since this tool approximates **Kabat**
boundaries, string equality is not the right metric — the comparison measures whether
the *offset* between the two schemes is consistent. Within a germline family and IMGT
CDR length it is constant, which is what "no boundary error" means here.

## Design principles

- **Never guess.** When an anchor is missing, the tool reports `not confident`
  with a reason instead of inventing a position. This is why chain assignment has
  zero incorrect calls across 7,000 sequences.
- **No overall verdict.** The tool does not produce a developability score.
  Weighting `NG` in a CDR against two `DG` in framework has no agreed answer, and
  once a number exists people read the number instead of the evidence.
- **Always show the evidence.** Every flag carries an explanation of why it matters.
- **It is a screening heuristic, not a numbering engine.** It does not replace
  ANARCI/IMGT, and says so in every output.

## What the real data taught me

Six bugs were found by running on real repertoire data — none of them would have
surfaced from the synthetic test corpus, because the synthetic data was built from
the same assumptions as the code.

- Conserved Trp inside CDR-H1 (`SYWMQ`) was taken as the FR2 anchor, truncating CDR-H1
- An extra Cys in FR3 was taken as the Cys92/Cys88 anchor — one case passed as
  *confident* with a plausible-looking but wrong CDR-H3
- Germline variants of the FR3 motif (`QVTIS`, `RITIN`) fell outside the regex
- Somatic hypermutation broke the FR3 motif entirely (mouse plasmablast)
- A decoy FR3 motif inside CDR-H2 shifted the boundary
- Two Trp residues in IGKV1 light chains (`...QGISS**W**LA**W**YQQ...`) shifted CDR-L1/L2 by 3

The last one had been **rejected once**: judging from frequency distributions alone,
the deviating 6% looked like normal germline variation. Comparing against ANARCI ground
truth showed that 6% *was* the error. Distribution tells you what is common, not what is
correct. The reasoning — including the reversal — is kept in
[VALIDATION.md](VALIDATION.md).

## Repository

```
abcheck/        validate · annotate · liabilities · properties · report · render · cli
scripts/        synthetic corpus generator · OAS converter · two validation scripts
data/synthetic/ 24 synthetic sequences with ground-truth CDR coordinates
tests/          204 tests
```

All sequence data in this repository is **synthetic**. No therapeutic antibody or
germline database records are included; the synthetic framework is weighted random
residues with conserved anchors fixed in place.

- [README.md](README.md) — full documentation (Korean)
- [TUTORIAL.md](TUTORIAL.md) — step-by-step guide (Korean)
- [VALIDATION.md](VALIDATION.md) — measurement method, results, known limitations (Korean)
- [DECISIONS.md](DECISIONS.md) — design decisions and their rationale (Korean)

## Out of scope

Structure prediction, docking, affinity prediction. No external database or API calls.

## License

MIT — see [LICENSE](LICENSE).
