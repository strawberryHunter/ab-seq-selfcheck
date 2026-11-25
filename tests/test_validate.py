"""입력 검증 테스트 — 거부 사유가 한글로 나오는지까지 본다."""

from __future__ import annotations

import pytest

from abcheck.validate import (
    MAX_DOMAIN_LENGTH,
    MIN_DOMAIN_LENGTH,
    SequenceError,
    ensure_valid,
    validate,
)


def _has_korean(text: str) -> bool:
    return any("가" <= ch <= "힣" for ch in text)


def test_valid_sequence_passes(records):
    for rec in records:
        result = validate(rec["sequence"])
        assert result.ok, f"{rec['id']}: {result.errors}"
        assert result.cleaned == rec["sequence"]


def test_lowercase_and_whitespace_are_normalized(records):
    seq = records[0]["sequence"]
    wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)).lower()

    result = validate(wrapped)
    assert result.ok
    assert result.cleaned == seq
    assert any("공백" in w for w in result.warnings)


def test_empty_input_rejected():
    for bad in ["", "   ", "\n\n"]:
        result = validate(bad)
        assert not result.ok
        assert _has_korean(result.errors[0])


def test_fasta_header_rejected(records):
    seq = records[0]["sequence"]
    result = validate(f">SYN-VH-001 설명\n{seq}")

    assert not result.ok
    assert "FASTA" in result.errors[0]
    assert _has_korean(result.errors[0])


def test_non_amino_acid_characters_rejected():
    result = validate("S" * 50 + "123" + "S" * 50)
    assert not result.ok
    joined = " ".join(result.errors)
    assert "아미노산이 아닌 문자" in joined


def test_dna_like_input_rejected():
    """DNA 서열은 전부 유효한 아미노산 문자지만, 핵산 오입력으로 걸러야 한다."""
    result = validate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 3)
    assert not result.ok
    joined = " ".join(result.errors)
    assert "DNA/RNA" in joined
    assert "번역한 단백질 서열" in joined


def test_homopolymer_is_not_mistaken_for_dna():
    """poly-A 는 핵산 판정에서 제외한다 (ACGT 중 3종 이상일 때만 핵산으로 본다)."""
    result = validate("A" * 120)
    assert result.ok, result.errors


def test_ambiguous_codes_rejected_with_explanation():
    result = validate("S" * 50 + "X" + "S" * 50)
    assert not result.ok
    joined = " ".join(result.errors)
    assert "모호한 아미노산 코드" in joined
    assert "불특정 잔기" in joined
    assert "추측하지 않는다" in joined


def test_each_ambiguous_code_is_named():
    for code, label in [("B", "Asp 또는 Asn"), ("Z", "Glu 또는 Gln"), ("J", "Leu 또는 Ile")]:
        result = validate("S" * 50 + code + "S" * 50)
        assert not result.ok
        assert label in " ".join(result.errors)


def test_too_short_rejected():
    result = validate("ACDEFGHIKL")
    assert not result.ok
    joined = " ".join(result.errors)
    assert "너무 짧다" in joined
    assert str(MIN_DOMAIN_LENGTH) in joined
    assert "CDR 조각" in joined


def test_too_long_rejected_with_scfv_hint():
    result = validate("S" * (MAX_DOMAIN_LENGTH + 1))
    assert not result.ok
    joined = " ".join(result.errors)
    assert "너무 길다" in joined
    assert "scFv" in joined


def test_boundary_lengths_accepted():
    assert validate("S" * MIN_DOMAIN_LENGTH).ok
    assert validate("S" * MAX_DOMAIN_LENGTH).ok


def test_ensure_valid_raises_korean_message():
    with pytest.raises(SequenceError) as exc:
        ensure_valid("ACDEF")
    assert _has_korean(str(exc.value))


def test_ensure_valid_returns_cleaned(records):
    seq = records[0]["sequence"]
    assert ensure_valid(seq.lower() + "\n") == seq


def test_error_messages_only_mention_options_that_exist():
    """오류 메시지가 없는 옵션을 안내하면 안 된다.

    한때 FASTA 헤더 오류가 존재하지 않는 '--fasta' 옵션을 안내했다. 매뉴얼을 만들며
    실제 메시지를 옮기다 발견했다. CLI 의 실제 옵션과 대조한다.
    """
    import re

    from abcheck.cli import _build_parser

    known = {a for action in _build_parser()._actions for a in action.option_strings}
    vh = "Q" * 100
    inputs = [">x\n" + vh, "S" * 10, "S" * 300, "S" * 50 + "X" + "S" * 50,
              "S" * 50 + "1" + "S" * 50, "ACGT" * 30]
    for raw in inputs:
        for message in validate(raw).errors:
            for opt in re.findall(r"--[a-z][a-z-]*", message):
                assert opt in known, f"없는 옵션 {opt} 를 안내한다: {message}"
