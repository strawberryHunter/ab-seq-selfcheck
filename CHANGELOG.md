# Changelog — browser app (web/abcheck.html)

One line per published version: what changed, why, and the steps a task takes before → after.
Steps are counted as clicks or drops from an open page; typing a name counts as one step.

| Version | Changed | Why | Task: steps before → after |
|---|---|---|---|
| v4 (2026-11-26) | Series design system (green tokens, bundled Pretendard, cards + drawer); hero card; rule sets saved in the browser; several files at once; pick → one-page sheet; mail text copy; last analysis, filters and theme restored on reload | The app has to shorten real triage work, not only look tidy | Add candidates from N files: N×3 → 2 · Apply the team's criteria: re-read every sequence → 0 (saved rule set applies itself) · One-page meeting sheet: build a table by hand → 3 |
| v5 (2026-11-26) | Rows copied from Excel (name + sequence, with or without a header), a single column of sequences, CSV text and JSON text can be pasted directly; a save that cannot happen says so and copies the content instead of claiming "saved" | Round 1: pasting Excel rows failed with a misleading "ambiguous code B" error although the guide said it worked | Add candidates from an Excel sheet: 6 (save as CSV UTF-8 → new analysis → open file → pick) → 3 (new analysis → paste → analyse) |
| v6 (2026-11-26) | Partner chains found by name (mAb01_VH / mAb01_VL): "검토 필요 넣기" brings the partner chain along, the sheet lists each antibody's chains together, the drawer links to the partner; "FASTA 복사" copies the sequences of the list on screen | Round 2: a meeting sheet built from flagged chains lacked the other chain of each antibody, and rechecking 판정 불가 sequences in ANARCI meant copying them one by one | One-page sheet with both chains: 3 + 3 per missing partner → 3 · Sequences to recheck: 2 per sequence → 2 in total |
| v7 (2026-11-26) | CDR comparison shows 200 rows per chain with "더 보기"; deleting a rule set can be undone from the toast | Round 3 (1,200 candidates at once: analyse + first screen 0.14 s, re-triage 0.05 s): the comparison page drew 1,100 rows (67,000 px), and a deleted rule set was gone for good | Compare a large batch: scroll 67,000 px → search or page · Recover a deleted rule set: rebuild it (1 + toggles) → 1 |
| v8 (2026-11-26) | Home list has a 서열 / 항체 switch (shown when partner chains exist, remembered): in 항체 an antibody is one card with its VH and VL, each chain keeping its own classification, CDR3 and numbers — nothing is merged; a chain that passes the filter brings its partner, and copy / CSV / drawer steps follow the cards on screen; one bookmark picks both chains; chain cards show a "짝 VL" link; the mail text lists an antibody's chains on adjacent lines | Round 4 (12 antibodies pasted from a sheet, heavy chains first): the two chains of an antibody sat far apart in the list and in the mail text, and a flagged VH said nothing about its VL | See both chains of an antibody: 2 per antibody (open drawer → partner link) → 1 (항체, stays on) · Pick an antibody for the sheet: 2 → 1 · Mail text grouped by antibody: reorder lines by hand → 0 |

## Next candidates (not built)

- Direct print of the one-page sheet (the artifact viewer blocks the print dialog; today: save HTML, open, print).
- Rule results in the JSON export (today the JSON is the engine's default-rule output; CSV carries the active rule set).
- Scatter labels for more than 15 candidates.
- Antibody view as the default when most chains have a partner (today: the switch is off until chosen once).
