"""Check that the README has not drifted from the code and data.

Researchers copy and paste straight from the README -- that is the tool's basic usage --
so tests guard against examples going stale and no longer running.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

from abcheck.annotate import TYPICAL_LENGTHS
from abcheck.liabilities import MOTIF_SPECS
from abcheck.validate import MAX_DOMAIN_LENGTH, MIN_DOMAIN_LENGTH

REPO_ROOT = Path(__file__).resolve().parents[1]
README_EN = REPO_ROOT / "README.md"
README_KO = REPO_ROOT / "README.ko.md"


def _readme() -> str:
    """The Korean README (README.ko.md) -- the full reference for researchers."""
    return README_KO.read_text(encoding="utf-8")


def _readme_en() -> str:
    """The English README (README.md) -- the repository landing page."""
    return README_EN.read_text(encoding="utf-8")


BOTH_READMES = pytest.mark.parametrize("read", [_readme, _readme_en], ids=["ko", "en"])


def _bash_blocks(text: str) -> list[str]:
    return [b.strip() for b in re.findall(r"```bash\n(.*?)```", text, re.DOTALL)]


# --------------------------------------------------------------------------
# Acceptance criterion: 3 run examples
# --------------------------------------------------------------------------

@BOTH_READMES
def test_readme_has_three_runnable_cli_examples(read):
    blocks = _bash_blocks(read())
    cli_examples = [b for b in blocks if "abcheck.cli" in b]
    assert len(cli_examples) >= 3, f"CLI 실행 예시가 {len(cli_examples)}개뿐이다"


def test_readme_examples_section_covers_three_modes():
    text = _readme()
    assert "## 실행 예시 3개" in text
    section = text.split("## 실행 예시 3개")[1].split("\n---")[0]
    # Terminal / HTML / JSON
    assert "--html" in section
    assert "--json" in section
    cli_blocks = [b for b in _bash_blocks(section) if "abcheck.cli" in b]
    assert len(cli_blocks) >= 3


@BOTH_READMES
def test_pasted_sequence_example_actually_runs(records, read):
    """The pasted sequence in example 1 is in the current corpus and actually gets annotated."""
    text = read()
    matches = re.findall(r'abcheck\.cli "([ACDEFGHIKLMNPQRSTVWY]+)"', text)
    assert matches, "README 에 붙여넣기용 서열 예시가 없다"

    corpus_sequences = {r["sequence"] for r in records}
    for seq in matches:
        assert seq in corpus_sequences, (
            "README 예시 서열이 코퍼스에 없다 — 코퍼스를 다시 만든 뒤 README 를 "
            "갱신하지 않은 것으로 보인다."
        )

    result = subprocess.run(
        [sys.executable, "-m", "abcheck.cli", matches[0]],
        cwd=REPO_ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr
    assert "CDR-H1" in result.stdout


def test_corpus_html_example_runs(tmp_path):
    """Example 2 actually produces HTML."""
    result = subprocess.run(
        [sys.executable, "-m", "abcheck.cli", "data/synthetic/sequences.json",
         "--html", "--out-dir", str(tmp_path)],
        cwd=REPO_ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "report.html").exists()


# --------------------------------------------------------------------------
# Do the docs and code state the same values?
# --------------------------------------------------------------------------

@pytest.mark.parametrize("read, marker", [
    (_readme, "| CDR | 통상 길이 |"),
    (_readme_en, "| CDR | Typical length |"),
], ids=["ko", "en"])
def test_documented_cdr_lengths_match_code(read, marker):
    """Does the typical CDR length table state the same values as TYPICAL_LENGTHS in the code?

    Strings like CDR-H1 also appear as column headers in the summary table, so the
    search is narrowed to the 'typical length' table.
    """
    text = read()
    assert marker in text, "README 에 통상 CDR 길이 표가 없다"
    table = text.split(marker)[1].split("\n\n")[0]

    for name, (lo, hi) in TYPICAL_LENGTHS.items():
        expected = str(lo) if lo == hi else f"{lo}-{hi}"
        row = re.search(rf"\|\s*{re.escape(name)}\s*\|\s*([^|]+?)\s*\|", table)
        assert row, f"통상 길이 표에 {name} 행이 없다"
        assert row.group(1) == expected, (
            f"{name}: README 는 {row.group(1)}, 코드는 {expected}"
        )


def test_documented_length_limits_match_code():
    assert f"{MIN_DOMAIN_LENGTH}-{MAX_DOMAIN_LENGTH} 잔기" in _readme()
    assert f"{MIN_DOMAIN_LENGTH}-{MAX_DOMAIN_LENGTH} residues" in _readme_en()


@BOTH_READMES
def test_documented_motifs_match_code(read):
    """The README's liability table includes every motif in the code."""
    text = read()
    for spec in MOTIF_SPECS:
        if spec.kind == "oxidation":
            continue                      # The table lists `M` `W` together
        assert f"`{spec.pattern}`" in text or all(
            f"`{alt}`" in text for alt in _expand(spec.pattern)
        ), f"README 에 {spec.pattern} 가 없다"


def _expand(pattern: str) -> list[str]:
    """Expand a pattern like `N[STNH]` into `NS` `NT` ... (to match the README notation)."""
    m = re.fullmatch(r"([A-Z])\[([A-Z]+)\]", pattern)
    if not m:
        return [pattern]
    first, options = m.groups()
    return [first + ch for ch in options]


# --------------------------------------------------------------------------
# Is the limitations notice in the README?
# --------------------------------------------------------------------------

def test_readme_states_it_is_not_a_numbering_engine():
    text = _readme()
    assert "numbering 엔진이 아니다" in text
    assert "ANARCI" in text


def test_readme_lists_when_to_recheck_with_anarci():
    text = _readme()
    assert "재확인해야 하는 경우" in text
    assert "confident=False" in text
    assert "후보를 탈락시키려는" in text


def test_readme_states_no_overall_verdict():
    text = _readme()
    assert "종합 판정을 내리지 않는다" in text


def test_readme_states_synthetic_data_policy():
    text = _readme()
    assert "전부 합성" in text
    assert "germline" in text


def test_readme_states_out_of_scope():
    text = _readme()
    assert "구조 예측" in text
    assert "도킹" in text
    assert "결합력 예측" in text


# --------------------------------------------------------------------------
# TUTORIAL.md
# --------------------------------------------------------------------------

TUTORIAL = REPO_ROOT / "TUTORIAL.md"


@BOTH_READMES
def test_tutorial_exists_and_is_linked_from_readme(read):
    assert TUTORIAL.exists()
    assert "TUTORIAL.md" in read(), "README 에서 튜토리얼로 가는 링크가 없다"


def test_tutorial_example_sequence_is_in_corpus(records):
    """The tutorial's paste example has not drifted from the current corpus."""
    text = TUTORIAL.read_text(encoding="utf-8")
    matches = re.findall(r'abcheck\.cli "([ACDEFGHIKLMNPQRSTVWY]+)"', text)
    assert matches, "TUTORIAL 에 붙여넣기용 서열 예시가 없다"

    corpus = {r["sequence"] for r in records}
    for seq in matches:
        assert seq in corpus, (
            "TUTORIAL 예시 서열이 코퍼스에 없다 — 코퍼스를 다시 만든 뒤 "
            "TUTORIAL 을 갱신하지 않은 것으로 보인다."
        )


def test_tutorial_quoted_output_matches_real_output(records):
    """Does the CDR table pasted into the tutorial match the actual output?

    Keeps the numbers in the docs from going stale. Changing the output format or annotation rules trips this.
    """
    from abcheck.annotate import annotate

    text = TUTORIAL.read_text(encoding="utf-8")
    seq = re.findall(r'abcheck\.cli "([ACDEFGHIKLMNPQRSTVWY]+)"', text)[0]
    ann = annotate(seq)

    for name, cdr in ann.cdrs.items():
        row = re.search(
            rf"{re.escape(name)}\s*\|\s*(\d+)-(\d+)\s*\|\s*(\d+)\s*\|\s*(\w+)", text
        )
        assert row, f"TUTORIAL 에 {name} 행이 없다"
        start, end, length, seq_text = row.groups()
        assert (int(start), int(end), int(length), seq_text) == (
            cdr.start, cdr.end, cdr.length, cdr.seq
        ), f"{name}: 문서 {row.groups()} vs 실제 {(cdr.start, cdr.end, cdr.length, cdr.seq)}"


def test_tutorial_covers_the_common_errors():
    text = TUTORIAL.read_text(encoding="utf-8")
    for topic in ["보존 Cys22 를 찾지 못했다", "너무 길다", "DNA/RNA", "모호한 아미노산"]:
        assert topic in text, f"TUTORIAL 에 '{topic}' 안내가 없다"


def test_tutorial_states_no_network_and_no_scoring():
    text = TUTORIAL.read_text(encoding="utf-8")
    assert "외부 API 를 호출하지 않는다" in text
    assert "점수를 매기지 않는" in text


# --------------------------------------------------------------------------
# VALIDATION.md
# --------------------------------------------------------------------------

VALIDATION = REPO_ROOT / "VALIDATION.md"


@BOTH_READMES
def test_validation_doc_exists_and_is_linked(read):
    assert VALIDATION.exists()
    assert "VALIDATION.md" in read(), "README 에서 검증 기록으로 가는 링크가 없다"


def test_validation_doc_states_how_to_reproduce():
    """Numbers without a way to reproduce them are not validation."""
    text = VALIDATION.read_text(encoding="utf-8")
    assert "scripts/oas_to_fasta.py" in text
    assert "scripts/validate_on_oas.py" in text
    assert "opig.stats.ox.ac.uk" in text


def test_validation_doc_records_known_limitations():
    """What we decided not to fix, and why, must stay documented."""
    text = VALIDATION.read_text(encoding="utf-8")
    assert "## Known limitations" in text
    assert "rejected" in text
    assert "reversed" in text, "the record of a fix rejected and then reversed is missing"
    assert "the conserved Trp itself is substituted" in text, "the remaining limitation is missing"


def test_validation_script_runs_on_the_synthetic_fasta():
    """Does the validation script actually run? Smoke test with a synthetic FASTA."""
    result = subprocess.run(
        [sys.executable, "scripts/validate_on_oas.py", "data/synthetic/sequences.fasta"],
        cwd=REPO_ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    assert result.returncode == 0, result.stderr
    assert "CDR 주석 성공" in result.stdout
    assert "FR3↔보존Cys" in result.stdout


# --------------------------------------------------------------------------
# What it needs as a portfolio piece
# --------------------------------------------------------------------------

def test_license_exists():
    license_file = REPO_ROOT / "LICENSE"
    assert license_file.exists(), "LICENSE 가 없으면 남이 쓸 수 없는 코드다"
    assert "MIT License" in license_file.read_text(encoding="utf-8")


def test_readmes_link_to_each_other():
    assert "README.ko.md" in _readme_en(), "no link from the English README to the Korean one"
    assert "(README.md)" in _readme(), "no link from the Korean README to the English one"


def test_english_readme_carries_the_same_headline_numbers():
    """The headline figures must not drift between the two READMEs and VALIDATION.md."""
    validation = VALIDATION.read_text(encoding="utf-8")
    for token in ["7,000", "0.35%", "16,498", "98.7%"]:
        assert token in _readme_en(), f"README.md is missing {token}"
        assert token in _readme(), f"README.ko.md is missing {token}"
        assert token in validation, f"VALIDATION.md is missing {token}"


def test_english_readme_states_the_limits():
    """The English README also states it is not a numbering engine and gives no verdict."""
    english = _readme_en()
    assert "not a numbering engine" in english
    assert "No overall verdict" in english
    assert "confident=False" in english
    assert "about to drop a candidate" in english
    assert "**synthetic**" in english
    for item in ["Structure prediction", "docking", "affinity prediction"]:
        assert item in english


def test_readme_leads_with_validation():
    """Reviewers look for 30 seconds -- validation figures must come before usage."""
    text = _readme()
    assert text.index("## 검증") < text.index("## 실행 예시 3개")
    english = _readme_en()
    assert english.index("## Validation") < english.index("## Quick start")


def test_korean_summary_example_matches_real_output():
    """The summary-table rows quoted in README.ko.md match the CLI's actual summary rows.

    The example once lacked the triage column and showed totals that matched neither the
    quoted rows nor the real corpus.
    """
    text = _readme().split("### 0. CDR별 liability 집계")[1].split("### 0-2.")[0]
    quoted = [
        [c.strip().strip("*") for c in line.strip("|").split("|")]
        for line in text.splitlines()
        if line.startswith("| SYN-") or line.startswith("| **합계**")
    ]
    assert len(quoted) >= 4, "the Korean summary example has no rows"

    result = subprocess.run(
        [sys.executable, "-m", "abcheck.cli", "data/synthetic/sequences.json"],
        cwd=REPO_ROOT, capture_output=True, text=True, encoding="utf-8",
    )
    actual = {}
    for line in result.stdout.splitlines():
        cells = [c.strip() for c in line.split("|")]
        if len(cells) > 3 and (cells[0].startswith("SYN-") or cells[0] == "합계"):
            actual.setdefault(cells[0], cells)   # first occurrence = summary table
    for row in quoted:
        key = "합계" if row[0] == "합계" else row[0]
        assert key in actual, f"{key} is not in the real summary output"
        assert row == actual[key], f"README.ko.md {row} vs real {actual[key]}"
