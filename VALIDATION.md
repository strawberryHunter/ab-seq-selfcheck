# Validation record

These are the results of measuring the annotator on public repertoire data. **With synthetic
data, the generation rules are the ground truth, so it is only correct within "the assumptions
I made."** Real sequences contain patterns those assumptions never anticipated, and that is
exactly how every bug fixed in this repository came to light.

To reproduce:

```bash
# 1) Download a paired data unit (.csv.gz) from OAS — https://opig.stats.ox.ac.uk/webapps/oas/
python3 scripts/oas_to_fasta.py downloaded.csv.gz candidates.fasta --count 500
python3 scripts/validate_on_oas.py candidates.fasta
```

---

## Data

7 paired data units from OAS (Observed Antibody Space, CC-BY 4.0), chosen to differ in
cell type and degree of somatic hypermutation (SHM). 500 pairs each (VH+VL = 1,000 sequences),
**7,000 sequences** in total.

| Dataset | Species | Cells | Notes |
|---|---|---|---|
| Jaffe 2022 | Human | naive B | Almost no SHM |
| Alsoiussi 2020 | Mouse | plasmablast | Heavy SHM |
| King 2020 | Human | memory B | Tonsil, pediatric |
| Ota 2024 | Human | memory B | Peanut allergy |
| Phad 2022 | Human | plasma B | |
| Setliff 2019 | Human | unsorted B | HIV — long CDR-H3 |
| Sokal 2021 | Human | memory B | COVID-19, includes truncated reads |

The real sequences are not committed to the repository. Only the conversion and measurement
scripts are.

## Results

| Dataset | Input validation | CDR annotation succeeded | Confident | Chain identification |
|---|---|---|---|---|
| Jaffe 2022 (human naive) | 100% | 99.9% | 99.8% | 999/1000 |
| Alsoiussi 2020 (mouse) | 100% | **100%** | 99.2% | **1000/1000** |
| King 2020 | 100% | 99.1% | 94.7% | 994/1000 |
| Ota 2024 | 100% | 99.0% | 95.0% | 991/1000 |
| Phad 2022 | 100% | 99.0% | 93.1% | 990/1000 |
| Setliff 2019 (HIV) | 100% | 99.9% | 97.7% | 999/1000 |
| Sokal 2021 | 98.4% | 91.2% | 89.4% | 975/984 |
| **Total** | **99.8%** | **98.3%** | **95.6%** | **6948/6984** |

### CDR boundaries — compared with ANARCI ground truth

OAS data ships with **CDRs assigned by ANARCI using IMGT numbering**. We use these as ground
truth.

Caveat: this tool approximates **Kabat** boundaries, so the strings cannot be identical, and
should not be. We therefore check not the match rate but **whether the offset is constant**.
A constant offset is just a difference between numbering schemes; the cases that vary are the
real errors.

```
python3 scripts/compare_with_anarci.py downloaded.csv.gz --count 500
```

Human naive B, 500 pairs (2,996 comparisons):

| CDR | Compared | Modal offset | Agreement |
|---|---|---|---|
| CDR-H1 | 500 | (+5, +2) | 100.0% |
| CDR-H2 | 500 | (−1, +8) | 100.0% |
| **CDR-H3** | 500 | (+2, +0) | **100.0%** |
| CDR-L1 | 499 | (−3, +2) | 100.0% |
| CDR-L2 | 498 | (+0, +4) | 99.2% |
| **CDR-L3** | 499 | (+0, +0) | **100.0%** |
| **Total** | **2,996** | | **99.87%** |

**Holding germline family and IMGT CDR length fixed, the mismatch rate is 0.00% (0/2,996).**

All 7 datasets (400 pairs each, 16,427 comparisons):

| Dataset | Mismatch with family and length fixed |
|---|---|
| Jaffe 2022 (human naive) | **0.00%** (0/2,396) |
| Alsoiussi 2020 (mouse) | **0.00%** (0/2,385) |
| Setliff 2019 (HIV) | 0.25% (6/2,387) |
| Sokal 2021 | 0.28% (6/2,169) |
| King 2020 | 0.46% (11/2,369) |
| Ota 2024 | 0.63% (15/2,363) |
| Phad 2022 | 0.85% (20/2,358) |
| **Total** | **0.35%** (58/16,427) |

The CDR1/CDR2 offsets differ between families because each family places its insertions
differently, which changes the Kabat↔IMGT relationship. Within the same family and the same
length, the offset is perfectly constant.

CDR-H3 and CDR-L3 have a single offset regardless of family or length. In other words, the
boundaries are systematically consistent on the loops that matter most for binding.

### Chain identification — 0 misclassifications

In OAS paired data it is known which chain is the heavy chain, so this is **validation against
ground truth**.

```
체인 판별 정확   6948 / 6984  (99.5%)   오판 0,  판별 실패 36
```

(Gloss: chain identified correctly 6948 / 6984 (99.5%), misclassified 0, identification failed 36.)

**It never gave a wrong answer.** All 36 cases it cannot handle are reported as
"identification failed." They are truncated reads missing the terminal motif (`WGxG`/`FGxG`).

### Internal consistency — FR3↔conserved Cys distance

The distance from the FR3 start motif to the conserved Cys should be structurally constant.
Measured:

```
29 residues   3443 / 3443   100.0%
```

The annotator **uses this distance as a cross-check.** If the motif matches in the wrong place,
the distance comes out wrong, and the annotator then re-anchors from the Cys. Before this check
was added the figure was 3429/3443, and all 14 mismatches were cases where CDR-H2 came out
as 15 or 29 residues.

### Reasons for low confidence (311 of 7,000)

| Count | Reason | Nature |
|---|---|---|
| 145 | FR3 motif not found; back-calculated from Cys | Hypermutation. Position recovered but flagged conservatively |
| 39 | Conserved Cys23 not found | N-terminally truncated read |
| 36 | Terminal motif not found | C-terminally truncated read |
| 25 | Conserved Cys22 not found | N-terminally truncated read |
| 12 | End of CDR-H2 cannot be determined | Both anchors failed |
| 12 | CDR-H2 length out of range (15 or 21) | Atypical germline |
| 5 | FR3 motif inconsistent with distance → corrected | Caught by the cross-check |
| 6 | CDR-H1 length out of range (9–17) | Conserved Trp mutated |

**Cell type decides the outcome.** Human naive B (almost no SHM) reaches 99.8%; heavily
hypermutated memory/plasma B cells reach 93–95%. The FR3 motif assumes germline conservation,
and hypermutation breaks that assumption.

---

## Bugs found and fixed using real data

None of these would have been caught with synthetic data alone. Each is pinned as a regression
case.

### 1. A Trp inside CDR-H1 mistaken for the anchor

```
...CKASGYTFT SYWMQ W VKQ...
              ↑       ↑real Trp36
```

Using "the first Trp after CDR1" picked up a Trp inside CDR-H1 and truncated CDR-H1 to 2
residues. Trp is common in the CDR-H1 of real antibodies. Fixed by skipping the Kabat minimum
length (5) before searching.
→ The synthetic data was built without W in the CDRs, so this case **never occurred** there.
Regression: `SYN-VH-009`, `SYN-VL-006`

### 2. An extra Cys in FR3 mistaken for the anchor — a silent error

```
...LTIS C LQSEDFATYY C QQYYSYPPT FGGGTK...
       ↑fake          ↑real Cys88
```

CDR-L3 started 11 residues too early. **One of these (MOU-259-VH) had a CDR-H3 length of 25,
within the normal range, so it slipped past the length check and passed as `confident`** — a
silently wrong answer. The residue preceding the conserved Cys is measured as Y 83.2% + F 16.1%
= 99.3%, so the search was changed to use the `[YF]C` context.
Regression: `SYN-VL-007`

### 3. Germline variants of the FR3 motif not detected

The IGHV2 family's `QVTIS` (Q in the first position) and `RITIN` (N in the last position) fell
outside the regex. 26 of 1,000 VH sequences failed CDR-H2/H3 annotation. 4 candidate patterns
were compared on real data, and the motif was widened to `[RKQ].[TSAV][IMLFV][STAN]` — failures
went from 26 to 8, **with no change to the 1,000 previously successful annotations.**
Regression: `SYN-VH-011`, `SYN-VH-012`

### 4. Motif broken by hypermutation

`EATLT`, `RLNIS`, `KAILT` — substitutions in position 1 or 3. Rather than widening the regex
further, the annotator back-calculates from the fact that the FR3↔conserved Cys distance was
29 in all 992 measured cases. The position is recovered, but because the evidence is weaker it
is flagged as `low confidence`.
Regression: `SYN-VH-013`

### 5. A spurious FR3 motif inside CDR-H2

The motif matched first inside CDR-H2, producing a CDR-H2 of 15 or 29 residues. The distance-29
cross-check corrects this.
Regression: `SYN-VH-014`

### 6. A Trp inside CDR1 mistaken for the anchor — revealed by the ANARCI comparison

```
IGKV1   ...SCRAS QGISSW LA W YQQKPG...
                      ↑        ↑real Trp35
```

In germlines with two Trp residues, one at the end of CDR-L1 and one in FR2, the annotator took
the first one, shifting every subsequent boundary by 3 residues. Fixed by preferring the FR2
context after the conserved Trp (`W[VIA][RKQ]Q` / `W[YF][QL]Q`).

Across roughly 12,000 comparisons on 5 datasets, **CDR boundary mismatches fell from 3.83% to
0.31%**, and on the 2,996 human naive B comparisons the rate became **0.00% with family and
length fixed**.

> **I wrongly rejected this fix once.** Looking only at the distribution, without ground truth,
> I concluded that "the 6% with `LAW` is normal germline, so this pattern must not be used." When
> I checked against ANARCI ground truth, that 6% turned out to be exactly my error. The lesson I
> keep from this: a distribution alone cannot tell "common" apart from "correct."

Regression: `test_trp_anchor_uses_fr2_context_for_light_chain`

---

## Known limitations — not fixed on purpose

These were measured and **deliberately left unfixed**. Fixing them would break something bigger.

### Sequences with a mutated conserved Trp

```
KING-133-VL   ...RASQTVNINLA F YQQKPGQ...      (normally W YQQ)
SOKA-351-VH   CDR-H1 comes out as 9 residues
```

When Trp35/Trp36 is substituted by another residue, the anchor cannot be found and the annotator
picks up the wrong Trp.

**I measured whether the FR2 context (the 3 residues after `W`) could recover it:**

| Chain | Candidate pattern | Hit rate |
|---|---|---|
| VH | `[VI][RK]Q` | 85.6% |
| VL | `[YF][QL]Q` | 83.1% |

**I first rejected this, then reversed the decision after comparing with ANARCI ground truth.**
Looking only at the distribution, I assumed "the 6–8% that don't match must be normal germline,"
but when matched against ground truth, a substantial share of them were my errors. Adding the
context-first rule brought CDR boundary mismatches down from 3.83% to 0.31% (item 6 above).

However, sequences where **the conserved Trp itself is substituted** by another residue still
cannot be recovered, because the FR2 context pattern itself contains the Trp. Such sequences
(about 1 in 3,426) are still rejected as `low confidence`.

**The current behavior is correct — it does not guess; it rejects as `low confidence`.**

### Truncated reads

In data that mixes in reads truncated at the 5'/3' end, such as Sokal 2021, the anchors cannot
be found (annotation success on that dataset: 91.2%). Even if a sequence passes input length
validation (70-200), if part of the domain is missing the anchor is simply absent. Rejection is
the correct behavior here too.

### Nanobodies (VHH) — runs, but outside the validated scope

I took sequences from 5 real nanobody structures in the PDB and ran them. **Annotations are
produced** — the conserved Cys, the Trp and the terminal `WGxG` are all found, and the values
look plausible: CDR-H1 5 residues, CDR-H2 17 residues, CDR-H3 15-20 residues. For sequences
with an N-terminal expression tag (`MA…`), the anchors shift together, so the annotation stays
self-consistent.

However, **the FR2 context warning always appears.**

```
참고: 보존 Trp36 을 FR2 문맥(W[VIA][RKQ]Q)으로 찾지 못해 38번의 Trp 을 그대로 썼다.
```

(Gloss: "Note: conserved Trp36 could not be found via the FR2 context (W[VIA][RKQ]Q), so the Trp
at position 38 was used as is.")

This is due to the VHH hallmark substitution (V37F) — the `W V R Q` of a typical VH is `W F R Q`
in a VHH. This position is what makes a VHH a VHH, so "fixing" it by adding `F` to the context
pattern would be wrong. The warning is, if anything, an accurate signal.

**I could not validate this.** None of the 91 OAS studies include camelid data, so there was no
way to compare against ANARCI ground truth. Accordingly, no accuracy claim is made for VHH. If
you use it, confirm the boundaries with ANARCI first.

### Why a failed FR2 context lookup for the conserved Trp does not lower confidence

A back-calculated FR3 motif lowers the result to `low confidence`, but a failed Trp FR2 context
lookup only leaves a warning. That looked inconsistent, so I measured it.

```
Trp FR2 context failures among 6,881 conventional antibodies: 503 (7.31%)
```

Yet on the same data, the CDR boundary mismatch against ANARCI is 0.35%. In other words,
**even without the context, it picks the correct Trp in most cases.** Downgrading all 7.31% to
`low confidence` would add 500 false alarms and dull the warning. The FR3 back-calculation, by
contrast, rests on one more level of indirect evidence, so it is downgraded.

### No structural information

The liability flags look only at sequence motifs. Sites that are buried and pose no real problem
are flagged too. This is outside the design scope and is stated in the README and in every output.

---

## Summary

| Item | Value |
|---|---|
| Sequences measured | 7,000 (7 public datasets) |
| **CDR boundaries (vs. ANARCI ground truth)** | **0.35% mismatch** (16,427 comparisons, 7 datasets) |
| Chain identification **misclassifications** | **0** |
| CDR annotation succeeded | 98.3% |
| Internal consistency (FR3↔Cys) | 100% |
| Bugs found and fixed using real data | 6 |
| Fixes rejected once, then reversed after ground-truth comparison | 1 |

The design goal of this tool was **to state the reason rather than guess when it cannot do
something.** 0 chain misclassifications across 7,000 sequences, and a Korean-language reason
attached to every failure, confirm that goal.
