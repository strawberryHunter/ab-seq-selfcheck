"""입력 검증 — 아미노산 아닌 문자, FASTA 헤더, 길이 범위 이탈을 거부한다.

오류 메시지는 전부 한글이다. 코딩을 모르는 연구원이 읽고 바로 고칠 수 있어야 하므로,
"무엇이 문제인지 + 어디가 문제인지 + 어떻게 하면 되는지"를 한 문장에 담는다.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: 표준 아미노산 20종
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: 모호 코드. 따로 안내 메시지를 내기 위해 분리해 둔다.
AMBIGUOUS_AA: dict[str, str] = {
    "X": "불특정 잔기",
    "B": "Asp 또는 Asn",
    "Z": "Glu 또는 Gln",
    "J": "Leu 또는 Ile",
    "U": "셀레노시스테인",
    "O": "피롤리신",
}

#: 핵산 오입력 판정에 쓰는 문자 집합. 전부 유효한 아미노산 코드이기도 하다.
NUCLEOTIDE_LETTERS = frozenset("ACGTUN")

#: 이 길이 이상이면서 위 집합만으로 이뤄졌고, ACGT 중 3종 이상이 나오는 서열을
#: 핵산 오입력으로 본다. 3종 조건이 없으면 poly-A 같은 단일 잔기 서열이 오탐된다.
_NUCLEOTIDE_MIN_LENGTH = 30
_NUCLEOTIDE_MIN_DISTINCT = 3

#: 가변도메인 하나의 길이 범위. VH 는 보통 115-130, VL 은 105-115 잔기다.
MIN_DOMAIN_LENGTH = 70
MAX_DOMAIN_LENGTH = 200

#: 공백류는 조용히 제거한다. 줄바꿈된 서열을 그대로 붙여넣는 경우가 대부분이기 때문이다.
_STRIPPABLE = frozenset(" \t\r\n\v\f")


class SequenceError(ValueError):
    """입력 서열을 쓸 수 없을 때 올린다. 메시지는 한글이다."""


@dataclass
class ValidationResult:
    ok: bool
    cleaned: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def validate(raw: str) -> ValidationResult:
    """서열을 검증하고 정리한 결과를 돌려준다. 예외를 올리지 않는다.

    통과 조건:
      - FASTA 헤더('>' 로 시작하는 줄)가 섞여 있지 않다
      - 공백류를 뺀 모든 문자가 표준 아미노산 20종이다 (대소문자 무관)
      - 길이가 가변도메인 범위(70-200 잔기) 안이다
    """
    errors: list[str] = []
    warnings: list[str] = []

    if raw is None or not str(raw).strip():
        return ValidationResult(False, "", ["서열이 비어 있다. 아미노산 한 글자 코드로 붙여넣어라."])

    text = str(raw)

    # --- FASTA 헤더 ---
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

    # --- 공백류 제거, 대문자화 ---
    cleaned = "".join(ch for ch in text if ch not in _STRIPPABLE).upper()
    if len(cleaned) != len(text):
        warnings.append("공백과 줄바꿈은 제거하고 처리했다.")

    if not cleaned:
        return ValidationResult(False, "", ["서열이 비어 있다. 아미노산 한 글자 코드로 붙여넣어라."])

    # --- 문자 검증 ---
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

    # --- 핵산 서열 오입력 ---
    # A/C/G/T/U/N 은 전부 유효한 아미노산 한 글자 코드이기도 하므로 문자 검증만으로는
    # DNA 를 걸러낼 수 없다. 이 여섯 글자만으로 이뤄진 긴 서열은 단백질일 가능성이
    # 사실상 없으므로 핵산 오입력으로 본다.
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

    # --- 길이 검증 ---
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
    """검증을 통과한 정리된 서열을 돌려준다. 실패하면 SequenceError 를 올린다."""
    result = validate(raw)
    if not result.ok:
        raise SequenceError(" / ".join(result.errors))
    return result.cleaned
