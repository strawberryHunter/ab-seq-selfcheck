# ab-seq-selfcheck

항체 가변도메인(VH 또는 VL) 아미노산 서열을 붙여넣으면 **CDR 위치**와
**개발 가능성(developability) 점검 결과**를 보여주는 도구다. 후보 선별 단계에서
"이 서열 일단 괜찮은지" 빠르게 보는 용도.

설치할 것이 없다. 파이썬 3.10 이상이면 바로 돌아간다 (표준 라이브러리만 쓴다).

---

## 먼저 읽을 것 — 이 도구의 한계

**이 도구는 스크리닝 휴리스틱이지 numbering 엔진이 아니다.**

보존 잔기(Cys, Trp)와 모티프를 정규식으로 찾아 Kabat 경계를 **근사**할 뿐,
번호를 매기지 않는다. ANARCI 나 IMGT 를 대체하지 않는다.

- CDR 경계는 근사치다. 삽입/결실이 있는 서열에서는 경계가 밀릴 수 있다.
- liability 플래그는 **서열 모티프만** 본다. 구조적 용매 접근성을 보지 않으므로,
  묻혀 있어 실제로는 문제되지 않는 잔기도 그대로 플래그된다.
- pI / GRAVY / net charge 는 서열 기반 추정치다. 실측값과 다를 수 있다.
- **종합 판정을 내리지 않는다.** "개발 가능성 높음/낮음" 같은 결론은 이 도구가
  하지 않는다. 플래그와 그 근거만 보여주고, 판단은 사람이 한다.

### ANARCI / IMGT 로 반드시 재확인해야 하는 경우

1. 결과가 `confident=False` 로 나온 서열
2. CDR 길이가 통상 범위를 벗어난다고 표시된 서열
3. 비(非)항체이거나 비정형 스캐폴드로 의심되는 서열 — 나노바디(VHH), scFv,
   이중특이체 링커가 섞인 서열 등
4. **이 결과를 근거로 후보를 탈락시키려는 모든 경우**

앵커를 못 찾으면 이 도구는 위치를 **추측하지 않는다.** `confident=False` 로 표시하고
왜 실패했는지 한글로 알려준다. 결과가 비어 있다면 그건 버그가 아니라 설계다.

---

## 실행 예시 3개

아래 명령은 이 저장소 폴더에서 그대로 복사해 붙여넣으면 된다.

### 1. 서열 하나를 붙여넣어 터미널에서 보기

따옴표 안에 서열을 넣는다. 대소문자와 줄바꿈은 알아서 정리한다.

```bash
python3 -m abcheck.cli "YITTIETSSFHYGRVRHEIEQCALLKGGLSRAHSIWALKLPYITVTHQSKDNREFVYFQEVGQYDKATLSVSESLQHSLQSSLPVVRVVKLCLGAPSGSYFGRVWGQGSSGGFQY"
```

CDR 위치, liability 플래그, 물성이 표로 나온다. 서열 아래 `^` 표시가 CDR 구간이다.

### 2. 여러 서열을 HTML 리포트로 만들기

```bash
python3 -m abcheck.cli data/synthetic/sequences.json --html
```

`sample_output/report.html` 이 만들어진다. 브라우저로 열면 된다.
FASTA 파일도 같은 방식으로 넣을 수 있다.

```bash
python3 -m abcheck.cli my_sequences.fasta --html --out-dir results
```

### 3. 결과를 JSON 으로 받기

다른 스크립트나 엑셀로 넘길 때 쓴다.

```bash
python3 -m abcheck.cli data/synthetic/sequences.json --json > result.json
```

---

## 무엇을 보여주는가

### 0. CDR별 liability 집계 (리포트 맨 앞)

여러 후보를 나란히 비교할 수 있게, 서열 × 구간 × 위험도 건수를 교차 집계한 표가
리포트 맨 앞에 나온다. **중쇄와 경쇄는 표를 나눠서** 낸다 — 그래야 열 이름을
`CDR-H1` / `CDR-L1` 로 정확히 쓸 수 있다.

**중쇄 (VH)**

| 서열 | 주석 | CDR-H1 | CDR-H2 | CDR-H3 | framework | 도메인 플래그 |
|---|---|---|---|---|---|---|
| SYN-VH-003 | 신뢰 가능 | - | 높음 1 | - | 낮음 2 | - |
| SYN-VH-005 | 신뢰 가능 | - | - | - | 낮음 2 | 비결합 Cys 가능성 (높음) |
| SYN-VH-902 | 신뢰 불가 | 주석 실패 | 주석 실패 | 주석 실패 | 낮음 2 | 비결합 Cys 가능성 (중간) |
| **합계** | | **높음 1** | **높음 1** | **높음 1, 중간 2** | **중간 1, 낮음 21** | **2건** |

**경쇄 (VL)** 도 같은 모양으로 `CDR-L1` / `CDR-L2` / `CDR-L3` 열을 쓴다.

**체인 판별 실패 · 입력 오류** 는 따로 모은다. CDR 구간을 가를 수 없으므로
CDR 열 자체를 두지 않고 검출된 모티프 건수만 보여준다.

| 서열 | 상태 | 검출 모티프 | 도메인 플래그 |
|---|---|---|---|
| SYN-VH-901 | 체인 판별 실패 | 낮음 1 | - |

**건수 집계일 뿐 종합 판정이 아니다.** 점수를 매기거나 순위를 내지 않고, 입력 순서도
바꾸지 않는다. 건수가 적다고 안전한 것도, 많다고 탈락시킬 근거도 되지 않는다.

체인은 판별됐지만 주석이 실패한 서열(`SYN-VH-902`)은 해당 체인 표에 남되 CDR 칸에
`주석 실패` 로 표시하고, 검출된 모티프는 `framework` 칸에 모은다. 실제로는 CDR 안일
수 있으니 과소평가하면 안 된다.

HTML 리포트에서는 서열 ID 를 누르면 해당 서열의 상세 섹션으로 이동한다.
JSON 에서는 `summary` 키로 같은 집계를 받을 수 있다.

### 1. CDR 주석

보존 잔기 앵커로 Kabat 경계를 근사한다. 체인(H/L)은 말단 모티프로 자동 판별한다
(`WGxG` → 중쇄, `FGxG` → 경쇄).

| 체인 | 찾는 순서 |
|---|---|
| VH | `Cys22` → **CDR-H1**(Cys22+9부터) → `Trp36` → **CDR-H2**(Trp36+14부터) → FR3 모티프 `[RK].[TSAV][IMLFV][STA]` → `Cys92` → **CDR-H3**(Cys92+3부터) → `WGxG` |
| VL | `Cys23` → **CDR-L1**(Cys23+1부터) → `Trp35` → **CDR-L2**(Trp35+15부터, 7잔기) → `Cys88` → **CDR-L3**(Cys88+1부터) → `FGxG` |

좌표는 모두 **1-based, 양끝 포함**이다.

CDR 길이가 통상 범위를 벗어나면 위치는 그대로 보여주되 `confident=False` 로 떨어뜨린다.

| CDR | 통상 길이 |
|---|---|
| CDR-H1 | 5-7 |
| CDR-H2 | 16-19 |
| CDR-H3 | 3-25 |
| CDR-L1 | 7-17 |
| CDR-L2 | 7 |
| CDR-L3 | 5-13 |

### 2. liability 모티프

**CDR 안에 있는 것과 framework 에 있는 것의 위험도를 구분한다.** 같은 모티프라도
CDR 안이면 결합에 직접 영향을 주므로 더 위험하다.

| 항목 | 모티프 | 기본 위험도 |
|---|---|---|
| N-결합 당화 | `N[^P][ST]` | 높음 |
| deamidation | `NG` | 높음 |
| deamidation | `NS` `NT` `NN` `NH` | 중간 |
| isomerization | `DG` | 높음 |
| isomerization | `DS` `DT` `DD` | 중간 |
| fragmentation | `DP` | 중간 |
| 산화 | `M` `W` | 낮음 |
| 비결합 Cys | 도메인 내 Cys 개수가 홀수 | 도메인 단위 플래그 |

**등급 조정 규칙**: CDR 안이면 한 단계 **올리고**, framework 면 한 단계 **내린다**
(높음/중간/낮음 세 단계에서 양 끝은 그대로).

→ CDR 안의 `DG` = **높음**, framework 의 `DG` = **중간**
→ CDR 안의 `M` = **중간**, framework 의 `M` = **낮음**

CDR 주석이 실패한 서열에서는 구간을 구분할 수 없으므로 등급 조정을 하지 않고
`구간 미확정` 으로 표시한다. 과소평가를 피하기 위한 선택이다.

> 참고: 보존 Trp 앵커(Trp36, `WGxG` 의 Trp)도 산화 모티프라 그대로 플래그된다.
> 실제로는 대개 묻혀 있어 위험이 낮지만, 검출 자체를 숨기지는 않는다.

### 3. 물성

도메인 전체와 **CDR별로 각각** 계산한다.

- **pI** — EMBOSS pKa 값, bisection 으로 탐색
- **GRAVY** — Kyte-Doolittle hydropathy 평균 (양수=소수성, 음수=친수성)
- **net charge** — 기본 pH 7.4 (`--ph` 로 변경 가능)

CDR 조각의 pI 는 그 조각을 독립된 펩타이드로 보고 계산한 값이다. 실제 CDR 은
도메인 안에 박혀 있어 자유 말단이 없으므로, **구간 간 상대 비교에만** 쓸 것.

### 4. 입력 검증

다음은 거부하고 한글로 이유를 알려준다.

- 아미노산이 아닌 문자 (숫자, 구두점 등)
- 모호 코드 (`X` `B` `Z` `J` `U` `O`) — 추측하지 않는다
- FASTA 헤더가 섞인 입력 — 파일로 넣으면 된다
- 길이 범위 이탈 (70-200 잔기 밖). scFv·Fab 은 VH/VL 로 나눠서 넣어야 한다
- DNA/RNA 서열 오입력 (`ACGT` 는 전부 유효한 아미노산 코드라 따로 걸러낸다)

공백과 줄바꿈은 자동으로 제거한다. 소문자도 받는다.

---

## 합성 데이터

`data/synthetic/` 에 테스트용 서열 16개가 있다. **전부 합성이다** — 실제 치료용
항체나 germline 데이터베이스 레코드를 저장소에 넣지 않았다. framework 는 V 도메인의
잔기 조성을 흉내낸 가중 난수에 보존 잔기만 고정한 스캐폴드다.

CDR 을 삽입하는 시점에 위치를 알 수 있으므로 **정답 CDR 좌표를 함께 저장**한다.
테스트는 이 정답을 전부 복원하는지로 주석 로직을 검증한다.

다시 만들려면:

```bash
python3 scripts/make_synthetic_sequences.py
```

자세한 내용은 [`data/synthetic/README.md`](data/synthetic/README.md) 참고.

---

## 테스트

```bash
python3 -m pytest -q
```

핵심 테스트 두 개:

- `tests/test_annotate.py::test_ground_truth_cdr_positions` — 합성 서열의 정답 CDR
  위치를 전부 복원하는지
- `tests/test_annotate.py::test_broken_anchors_are_not_confident` — 앵커가 깨진
  서열에서 `confident=False` 가 나오고 위치를 추측하지 않는지

---

## 구조

```
abcheck/
  validate.py     입력 검증 (한글 오류 메시지)
  annotate.py     CDR 주석 — 앵커 탐색, confident 판정
  liabilities.py  liability 모티프 검출, CDR/framework 위험도 구분
  properties.py   pI / GRAVY / net charge
  report.py       위 결과를 한 덩어리로 합친다 (한계 고지 포함)
  render.py       터미널 표 / HTML / JSON
  cli.py          명령줄 진입점
scripts/
  make_synthetic_sequences.py   합성 코퍼스 생성기
data/synthetic/   합성 서열 + 정답 CDR 좌표
sample_output/    HTML 리포트 예시
tests/
DECISIONS.md      설계 판단 기록
```

## 범위 밖

구조 예측, 도킹, 결합력 예측은 하지 않는다. 외부 DB·API 를 호출하지 않는다.
