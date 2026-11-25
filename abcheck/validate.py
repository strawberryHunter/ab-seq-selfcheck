"""Input validation — rejects non-amino-acid characters, FASTA headers and out-of-range lengths.

All error messages are in Korean. Researchers who don't code must be able to read them
and fix the input right away, so each one packs "what is wrong + where + how to fix it"
into a single sentence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: The 20 standard amino acids
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Ambiguity codes. Kept separate so they get their own guidance message.
AMBIGUOUS_AA: dict[str, str] = {
    "X": "불특정 잔기",
    "B": "Asp 또는 Asn",
    "Z": "Glu 또는 Gln",
    "J": "Leu 또는 Ile",
    "U": "셀레노시스테인",
    "O": "피롤리신",
}

#: Letters used to detect nucleotide input by mistake. All are also valid amino acid codes.
NUCLEOTIDE_LETTERS = frozenset("ACGTUN")

#: A sequence at least this long, made only of the letters above, and containing 3+
#: of ACGT is treated as mistaken nucleotide input. Without the 3-letter rule,
#: single-residue sequences like poly-A would be false positives.
_NUCLEOTIDE_MIN_LENGTH = 30
_NUCLEOTIDE_MIN_DISTINCT = 3

#: Length range for a single variable domain. VH is usually 115-130 residues, VL 105-115.
MIN_DOMAIN_LENGTH = 70
MAX_DOMAIN_LENGTH = 200

#: Whitespace is stripped silently, since most users paste line-wrapped sequences as-is.
_STRIPPABLE = frozenset(" \t\r\n\v\f")


class SequenceError(ValueError):
    """Raised when an input sequence is unusable. The message is in Korean."""


@dataclass
class ValidationResult:
    ok: bool
    cleaned: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def validate(raw: str) -> ValidationResult:
    """Validate a sequence and return the cleaned result. Never raises.

    Passes when:
      - no FASTA header (a line starting with '>') is mixed in
      - every non-whitespace character is one of the 20 standard amino acids (case-insensitive)
      - the length is within the variable-domain range (70-200 residues)
    """
    errors: list[str] = []
    warnings: list[str] = []

    if raw is None or not str(raw).strip():
        return ValidationResult(False, "", ["서열이 비어 있다. 아미노산 한 글자 코드로 붙여넣어라."])

    text = str(raw)

    # --- FASTA header ---
    header_lines = [ln.strip() for ln in text.splitlines() if ln.strip().startswith(">")]
    if header_lines:
        sample = header_lines[0][:40]
        errors.append(
            f"FASTA 헤더가 섞여 있다 (예: '{sample}'). "
            "'>' 로 시작하는 줄을 지우고 아미노산 서열만 남겨라. "
            "여러 서열을 FASTA 로 넣으려면 파일로 저장한 뒤 파일 경로를 넣으면 된다 "
            "(예: python3 -m abcheck.cli 후보.fasta)."
        )
        return ValidationResult(False, "", errors)

    # --- strip whitespace, uppercase ---
    cleaned = "".join(ch for ch in text if ch not in _STRIPPABLE).upper()
    if len(cleaned) != len(text):
        warnings.append("공백과 줄바꿈은 제거하고 처리했다.")

    if not cleaned:
        return ValidationResult(False, "", ["서열이 비어 있다. 아미노산 한 글자 코드로 붙여넣어라."])

    # --- character check ---
    bad: list[tuple[int, str]] = [
        (i + 1, ch) for i, ch in enumerate(cleaned) if ch not in STANDARD_AA
    ]
    if bad:
        ambiguous = sorted({ch for _, ch in bad if ch in AMBIGUOUS_AA})
        others = sorted({ch for _, ch in bad if ch not in AMBIGUOUS_AA})

        if ambiguous:
            desc = ", ".join(f"{ch}({AMBIGUOUS_AA[ch]})" for ch in ambiguous)
            errors.append(
                f"모호한 아미노산 코드가 있다: {desc}. "
                "이 도구는 모호 코드를 추측하지 않는다. 확정된 잔기로 바꿔라."
            )
        if others:
            spots = ", ".join(f"{pos}번째 '{ch}'" for pos, ch in bad[:5] if ch in others)
            errors.append(
                f"아미노산이 아닌 문자가 있다: {', '.join(repr(c) for c in others)} ({spots}). "
                "DNA 서열이나 숫자·구두점이 섞이지 않았는지 확인하라."
            )
        return ValidationResult(False, cleaned, errors, warnings)

    # --- nucleotide sequence entered by mistake ---
    # A/C/G/T/U/N are all valid one-letter amino acid codes too, so the character
    # check alone can't catch DNA. A long sequence made only of these six letters is
    # practically never a protein, so we treat it as mistaken nucleotide input.
    if (
        len(cleaned) >= _NUCLEOTIDE_MIN_LENGTH
        and set(cleaned) <= NUCLEOTIDE_LETTERS
        and len(set(cleaned) & frozenset("ACGT")) >= _NUCLEOTIDE_MIN_DISTINCT
    ):
        errors.append(
            f"DNA/RNA 서열로 보인다 (A/C/G/T/U/N 으로만 이뤄진 {len(cleaned)} 글자). "
            "이 도구는 아미노산 서열만 받는다. 번역한 단백질 서열을 넣어라."
        )
        return ValidationResult(False, cleaned, errors, warnings)

    # --- length check ---
    n = len(cleaned)
    if n < MIN_DOMAIN_LENGTH:
        errors.append(
            f"서열이 너무 짧다 ({n} 잔기). 가변도메인 하나는 보통 {MIN_DOMAIN_LENGTH} 잔기 이상이다. "
            "CDR 조각만 넣은 것이 아닌지 확인하라 — 이 도구는 도메인 전체를 받는다."
        )
    elif n > MAX_DOMAIN_LENGTH:
        errors.append(
            f"서열이 너무 길다 ({n} 잔기). 이 도구는 VH 또는 VL 도메인 하나만 받는다 "
            f"(허용 범위 {MIN_DOMAIN_LENGTH}-{MAX_DOMAIN_LENGTH} 잔기). "
            "scFv·Fab·전장 항체라면 VH 와 VL 을 각각 나눠서 넣어라."
        )

    if errors:
        return ValidationResult(False, cleaned, errors, warnings)

    return ValidationResult(True, cleaned, errors, warnings)


def ensure_valid(raw: str) -> str:
    """Return the cleaned sequence if it passes validation; raise SequenceError otherwise."""
    result = validate(raw)
    if not result.ok:
        raise SequenceError(" / ".join(result.errors))
    return result.cleaned
