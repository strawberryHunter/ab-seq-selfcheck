# Changelog — browser app (web/abcheck.html)

One line per published version: what changed, why, and the steps a task takes before → after.
Steps are counted as clicks or drops from an open page; typing a name counts as one step.

| Version | Changed | Why | Task: steps before → after |
|---|---|---|---|
| v4 (2026-11-26) | Series design system (green tokens, bundled Pretendard, cards + drawer); hero card; rule sets saved in the browser; several files at once; pick → one-page sheet; mail text copy; last analysis, filters and theme restored on reload | The app has to shorten real triage work, not only look tidy | Add candidates from N files: N×3 → 2 · Apply the team's criteria: re-read every sequence → 0 (saved rule set applies itself) · One-page meeting sheet: build a table by hand → 3 |
| v5 (2026-11-26) | Rows copied from Excel (name + sequence, with or without a header), a single column of sequences, CSV text and JSON text can be pasted directly; a save that cannot happen says so and copies the content instead of claiming "saved" | Round 1: pasting Excel rows failed with a misleading "ambiguous code B" error although the guide said it worked | Add candidates from an Excel sheet: 6 (save as CSV UTF-8 → new analysis → open file → pick) → 3 (new analysis → paste → analyse) |
