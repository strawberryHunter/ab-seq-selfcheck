"""README 가 코드·데이터와 어긋나지 않는지 확인한다.

연구원이 README 를 그대로 복사해 붙여넣는 것이 이 도구의 기본 사용법이므로,
예시가 낡아서 안 돌아가는 상황을 테스트로 막는다.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from abcheck.annotate import TYPICAL_LENGTHS
from abcheck.liabilities import MOTIF_SPECS
from abcheck.validate import MAX_DOMAIN_LENGTH, MIN_DOMAIN_LENGTH

REPO_ROOT = Path(__file__).resolve().parents[1]
README = REPO_ROOT / "README.md"


def _readme() -> str:
    return README.read_text(encoding="utf-8")


def _bash_blocks(text: str) -> list[str]:
    return [b.strip() for b in re.findall(r"```bash\n(.*?)```", text, re.DOTALL)]


# --------------------------------------------------------------------------
# 완료 기준: 실행 예시 3개
# --------------------------------------------------------------------------

def test_readme_has_three_runnable_cli_examples():
    blocks = _bash_blocks(_readme())
    cli_examples = [b for b in blocks if "abcheck.cli" in b]
    assert len(cli_examples) >= 3, f"CLI 실행 예시가 {len(cli_examples)}개뿐이다"


def test_readme_examples_section_covers_three_modes():
    text = _readme()
    assert "## 실행 예시 3개" in text
    section = text.split("## 실행 예시 3개")[1].split("\n---")[0]
    # 터미널 / HTML / JSON 세 가지
    assert "--html" in section
    assert "--json" in section
    cli_blocks = [b for b in _bash_blocks(section) if "abcheck.cli" in b]
    assert len(cli_blocks) >= 3


def test_pasted_sequence_example_actually_runs(records):
    """예시 1의 붙여넣기 서열이 현재 코퍼스에 있고, 실제로 주석된다."""
    text = _readme()
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
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "CDR-H1" in result.stdout


def test_corpus_html_example_runs(tmp_path):
    """예시 2가 실제로 HTML 을 만든다."""
    result = subprocess.run(
        [sys.executable, "-m", "abcheck.cli", "data/synthetic/sequences.json",
         "--html", "--out-dir", str(tmp_path)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "report.html").exists()


# --------------------------------------------------------------------------
# 문서와 코드가 같은 값을 말하는지
# --------------------------------------------------------------------------

def test_documented_cdr_lengths_match_code():
    text = _readme()
    for name, (lo, hi) in TYPICAL_LENGTHS.items():
        expected = str(lo) if lo == hi else f"{lo}-{hi}"
        row = re.search(rf"\|\s*{re.escape(name)}\s*\|\s*([^|]+?)\s*\|", text)
        assert row, f"README 에 {name} 길이 행이 없다"
        assert row.group(1) == expected, (
            f"{name}: README 는 {row.group(1)}, 코드는 {expected}"
        )


def test_documented_length_limits_match_code():
    text = _readme()
    assert f"{MIN_DOMAIN_LENGTH}-{MAX_DOMAIN_LENGTH} 잔기" in text


def test_documented_motifs_match_code():
    """README 의 liability 표가 코드의 모티프를 모두 담고 있다."""
    text = _readme()
    for spec in MOTIF_SPECS:
        if spec.kind == "oxidation":
            continue                      # 표에 `M` `W` 로 묶어 적는다
        assert f"`{spec.pattern}`" in text or all(
            f"`{alt}`" in text for alt in _expand(spec.pattern)
        ), f"README 에 {spec.pattern} 가 없다"


def _expand(pattern: str) -> list[str]:
    """`N[STNH]` 같은 패턴을 `NS` `NT` ... 로 펼친다 (README 표기와 맞추기 위함)."""
    m = re.fullmatch(r"([A-Z])\[([A-Z]+)\]", pattern)
    if not m:
        return [pattern]
    first, options = m.groups()
    return [first + ch for ch in options]


# --------------------------------------------------------------------------
# 한계 고지가 README 에 있는지
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
