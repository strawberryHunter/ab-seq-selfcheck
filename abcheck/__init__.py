"""abcheck — a self-check tool for antibody variable-domain sequences.

This tool is a **screening heuristic**, not a numbering engine like ANARCI.
It approximates Kabat boundaries by regex-matching conserved anchor residues, and
when an anchor is broken it reports failure instead of guessing positions.
See the README and LIMITATION_NOTICE for the full list of limitations.
"""

from __future__ import annotations

__version__ = "0.1.0"

#: Limitation notice that must appear verbatim in every output (terminal/HTML/JSON).
LIMITATION_NOTICE: tuple[str, ...] = (
    "이 도구는 보존 잔기 앵커 기반 스크리닝 휴리스틱이다. ANARCI/IMGT 같은 numbering 엔진이 아니다.",
    "CDR 경계는 Kabat 을 근사한 것이며, 번호를 매기지 않는다. 삽입/결실이 있는 서열에서는 경계가 밀릴 수 있다.",
    "liability 플래그는 서열 모티프만 본다. 구조적 노출도(용매 접근성)를 보지 않으므로 실제 위험도는 다를 수 있다.",
    "pI/GRAVY/net charge 는 서열 기반 추정치다. 실측값과 다를 수 있다.",
    "다음 경우에는 ANARCI/IMGT 로 반드시 재확인할 것: "
    "confident=False 인 서열, CDR 길이가 통상 범위를 벗어난 서열, "
    "비(非)항체 또는 비정형 스캐폴드(나노바디·scFv·이중특이체 링커 포함)로 의심되는 서열, "
    "그리고 이 결과를 근거로 후보를 탈락시키려는 모든 경우.",
)

#: Why no overall verdict is given — also included in the output.
NO_VERDICT_NOTICE: str = (
    "이 도구는 '개발 가능성 높음/낮음' 같은 종합 판정을 내리지 않는다. "
    "플래그와 그 근거만 제시한다. 판단은 사람이 한다."
)
