# ab-seq-selfcheck

항체 가변도메인(VH/VL) 서열의 **CDR 위치**와 **developability liability** 를 보여주는
스크리닝 도구. 후보 선별 단계에서 "이 서열 일단 괜찮은지" 30초 만에 보는 용도다.

의존성 없음 — 파이썬 3.10+ 표준 라이브러리만. 외부 API 를 호출하지 않으므로 사내
서열도 밖으로 나가지 않는다.

*[English summary](README.en.md)*

```bash
python3 -m abcheck.cli 후보.fasta --html
```

## 검증

공개 레퍼토리(OAS) **7,000 서열**로 측정했다. 수치는 재현 가능하다 →
[`scripts/validate_on_oas.py`](scripts/validate_on_oas.py)

| 항목 | 결과 |
|---|---|
| 체인(H/L) 판별 | **오판 0건** / 6,984 |
| CDR 주석 성공 | 98.3% |
| **CDR 경계** (ANARCI 정답 대비, 7개 데이터셋 16,427건) | **불일치 0.35%** |
| 내부 일관성 (FR3↔보존Cys 거리) | 100% |

CDR 경계는 OAS 에 들어 있는 **ANARCI 주석을 정답으로** 비교했다. Kabat 과 IMGT 는
경계 정의가 다르므로 일치율이 아니라 오프셋의 일관성을 본다 →
[`scripts/compare_with_anarci.py`](scripts/compare_with_anarci.py)

**실제 데이터로 버그 6개를 찾아 고쳤다.** 합성 데이터만으로는 하나도 못 잡았을 것들이다.
그중 하나는 **한 번 잘못 기각했다가 정답 비교 후 번복**했다 — 분포만으로는 "흔한 것"과
"맞는 것"을 구분할 수 없다는 걸 배웠다. 과정은 [VALIDATION.md](VALIDATION.md) 에 그대로 남겼다.

## 설계 원칙

- **못 하면 추측하지 않는다.** 앵커를 못 찾으면 위치를 지어내지 않고 `신뢰 불가` 와
  한글 사유를 낸다. 7,000건에서 체인 판별 오판이 0건인 이유다.
- **종합 판정을 내리지 않는다.** "개발 가능성 높음/낮음" 같은 점수를 만들지 않는다.
  가중치에 근거가 없고, 숫자가 생기면 사람은 근거 대신 숫자만 본다.
- **근거를 같이 낸다.** 모든 플래그에 왜 문제인지가 붙는다.

---

> 처음 쓰신다면 **[TUTORIAL.md](TUTORIAL.md)** — 설치부터 출력 읽는 법까지 단계별로.

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
python3 -m abcheck.cli "YITTIETSSFHYGRVRHEIEQCALLKGGLSRAHSIWVRQALKLPYITVTKDNREFVYFQEVGQYDKATLSHQSVSESLQHSLQSSLPVVRVVYYCKLAPSGSYFGRVWGQGLGSSGGF"
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

### 0-2. 선별 규칙 — '검토 필요' 표시

후보가 수십 개일 때 사람이 먼저 볼 것을 골라낸다. 기본 규칙은 두 개다.

- CDR 안에 위험도 **'높음'**이면서 각 분류의 **최고 위험 모티프**(`N[^P][ST]` 당화, `NG`, `DG`, `DP`)가 1건 이상 있으면 **검토 필요**
- 비결합 Cys 등 도메인 플래그가 위험도 **'높음'**이면 **검토 필요**

결과는 세 가지로 나온다.

| 분류 | 뜻 |
|---|---|
| `검토 필요` | 위 규칙에 걸렸다. **항체가 나쁘다는 뜻이 아니다** |
| `기준 해당 없음` | 이 규칙에 걸리지 않았다. **안전하다는 뜻이 아니다** |
| `판정 불가` | CDR 주석이 실패해 규칙을 적용할 수 없다. 사람이 직접 봐야 한다 |

`검토 필요`에는 **어느 규칙에 왜 걸렸는지**가 같이 나온다.

```
선별 분류: 검토 필요
  → CDR-H1 32번에 N-결합 당화 서열 (NLS) — 위험도 높음
```

주석이 실패한 서열을 `기준 해당 없음`으로 넘기지 않고 `판정 불가`로 따로 두는 것이
중요하다. 못 본 것을 괜찮다고 적으면 그게 가장 위험하다.

> **왜 최고 위험 모티프만 보나** — CDR 안이면 위험도를 한 단계 올리므로, 기본 '중간'인
> 약한 모티프(`NS/NT/NN/NH`, `DS/DT/DD`)까지 '높음'이 된다. 실제 OAS 레퍼토리 16건으로
> 돌려 보니 CDR 안 '높음' 25건 중 13건이 `NS/NT/NN/NH` 였고 16건 중 12건(75%)이 검토
> 대상이 됐다. 그만큼 걸리면 "먼저 볼 것을 고른다"는 목적을 못 한다. 최고 위험 모티프로
> 한정하니 6건(37.5%)으로 줄었다. 약한 모티프도 표와 상세에는 '높음'으로 그대로 나오므로
> 정보가 사라지지는 않는다.

**규칙을 바꾸려면** `abcheck/triage.py` 의 `DEFAULT_RULES` 를 고친다. 기준이 코드에
명시적으로 남으므로 나중에 "왜 이 후보를 먼저 봤는지" 추적할 수 있다.

> 이 도구는 점수를 매기지 않는다. "CDR 의 NG 하나"와 "framework 의 DG 둘" 중 무엇이
> 더 나쁜지에는 합의된 답이 없고, 투여 경로·보관 조건·타깃에 따라 달라진다. 임의의
> 가중치로 숫자를 만들면 그 숫자는 근거가 아니라 추측이며, 사람은 근거 대신 숫자만
> 보게 된다. 그래서 규칙을 눈에 보이게 적어 두고 걸린 규칙을 그대로 보여준다.

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
| CDR-H2 | 16-20 |
| CDR-H3 | 3-30 |
| CDR-L1 | 7-17 |
| CDR-L2 | 7 |
| CDR-L3 | 5-15 |

이 범위는 "이 항체가 이상한가"가 아니라 **"주석을 믿을 수 있는가"**의 기준이다. OAS 실제 서열 2,000건의 실측 분포에 맞춰 상한을 잡았다. 하한은 넓히지 않았다 — 앵커를 잘못 집으면 CDR 이 짧아지므로 하한이 그 신호다.

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

## CI 와 차단

`.github/workflows/ci.yml` 이 `main` 으로 가는 PR 과 `main` 에 들어온 커밋을 검사한다.

| job | 검사 |
|---|---|
| `pytest (Python 3.10~3.13)` | 전체 테스트. 커밋된 `sample_output/report.html` 이 현재 출력과 같은지(날짜 제외)도 여기서 본다 |
| `런타임 의존성 없음` | `pyproject.toml` 에 런타임 의존성이 조용히 추가됐는지 |
| `합성 코퍼스 재현성` | 저장된 코퍼스가 생성기 출력과 같은지 |
| `CLI 동작` | README 실행 예시 3개가 실제로 돌아가는지, 잘못된 입력이 종료코드 1 인지 |

출력 양식을 바꿨다면 예시 리포트를 다시 만들어 함께 커밋해야 한다. 안 그러면 테스트가 막는다.

```bash
python3 -m abcheck.cli data/synthetic/sequences.json --html
```

### 실패한 코드를 막으려면 (저장소 설정 — 한 번만)

CI 는 실패를 **알려줄** 뿐이고, 실제로 **막는** 것은 GitHub 의 브랜치 보호다.
저장소 관리자가 한 번 켜야 한다.

`Settings → Branches → Add branch ruleset` (또는 `Add rule`) 에서:

1. 대상 브랜치: `main`
2. **Require status checks to pass** 체크
3. 아래 7개를 필수 검사로 추가

   ```
   pytest (Python 3.10)
   pytest (Python 3.11)
   pytest (Python 3.12)
   pytest (Python 3.13)
   런타임 의존성 없음
   합성 코퍼스 재현성
   CLI 동작
   ```
4. **Require branches to be up to date before merging** 체크
5. 관리자에게도 적용하려면 **Do not allow bypassing the above settings** 체크

`gh` CLI 로 하려면:

```bash
gh api -X PUT repos/strawberryHunter/ab-seq-selfcheck/branches/main/protection --input - <<'JSON'
{
  "required_status_checks": {
    "strict": true,
    "checks": [
      {"context": "pytest (Python 3.10)"},
      {"context": "pytest (Python 3.11)"},
      {"context": "pytest (Python 3.12)"},
      {"context": "pytest (Python 3.13)"},
      {"context": "런타임 의존성 없음"},
      {"context": "합성 코퍼스 재현성"},
      {"context": "CLI 동작"}
    ]
  },
  "enforce_admins": true,
  "required_pull_request_reviews": null,
  "restrictions": null,
  "allow_force_pushes": false,
  "allow_deletions": false
}
JSON
```

`checks` 는 객체의 배열이라 `gh api -f key=value` 로는 못 만든다. 반드시 `--input` 으로
JSON 을 통째로 넘겨야 한다.

켜졌는지 확인:

```bash
gh api repos/strawberryHunter/ab-seq-selfcheck/branches/main --jq .protected
```

> 비공개 저장소에서 브랜치 보호를 쓰려면 GitHub Pro/Team 이상이 필요하다.
> 무료 플랜이면 저장소를 공개로 바꾸거나, 아래 pre-push 훅으로 대신한다.

### push 전에 로컬에서 막기 (선택)

브랜치 보호와 별개로, 각자 클론에서 한 번 설정하면 테스트가 깨진 상태로 push 하는 것을
로컬에서 먼저 막는다.

```bash
git config core.hooksPath .githooks
```

`.githooks/pre-push` 가 `pytest` 와 코퍼스 재현성 검사를 돌리고, 실패하면 push 를 중단한다.
의도적으로 건너뛰려면 `git push --no-verify`.

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
  oas_to_fasta.py               OAS 데이터 → FASTA 변환
  validate_on_oas.py            공개 데이터로 정확도 측정
  compare_with_anarci.py        ANARCI 정답과 CDR 경계 비교
data/synthetic/   합성 서열 + 정답 CDR 좌표
sample_output/    HTML 리포트 예시
tests/
README.en.md      English summary
VALIDATION.md     검증 기록 (공개 데이터 7,000 서열)
DECISIONS.md      설계 판단 기록
```

## 범위 밖

구조 예측, 도킹, 결합력 예측은 하지 않는다. 외부 DB·API 를 호출하지 않는다.

## 라이선스

MIT — [LICENSE](LICENSE) 참고.
