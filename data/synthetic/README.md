# 합성 서열 코퍼스

`scripts/make_synthetic_sequences.py` 가 생성한다. 다시 만들려면:

```bash
python3 scripts/make_synthetic_sequences.py
```

기존 파일과 내용이 같은지만 확인하려면:

```bash
python3 scripts/make_synthetic_sequences.py --check
```

생성기는 시드가 고정돼 있어 같은 입력에서 항상 같은 서열이 나온다.

## 데이터 출처

**전부 합성이다.** 실제 치료용 항체 서열이나 특정 germline 데이터베이스 레코드를 복사해
오지 않았다. framework 부분은 V 도메인에서 흔히 보이는 잔기 조성을 흉내낸 가중 난수에
보존 잔기(Cys, Trp, WGxG/FGxG)만 정해진 자리에 고정해 만든 "germline 유사" 스캐폴드이며,
어떤 실제 서열과도 대응되지 않는다. 주석 로직의 회귀 테스트용 입력이고, 생물학적 해석의
대상이 아니다.

## 만드는 방식

"스캐폴드에 CDR 끼워 넣기" 방식이다. CDR 을 삽입하는 시점에 그 위치를 알 수 있으므로
정답 CDR 좌표를 함께 저장한다. 모든 VH 케이스는 같은 framework 스캐폴드를, 모든 VL
케이스는 같은 VL 스캐폴드를 공유한다. CDR 길이가 달라지면 그 뒤 앵커 위치가 전부
밀리므로(`SYN-VL-002` 는 CDR-L1 이 짧아 Cys2 가 85번에 온다), 주석기가 고정 번호가 아니라
앵커 상대 거리로 찾는지 검증할 수 있다.

### VH 좌표 설계 (1-based 선형 위치, Kabat 경계 모사)

| 구간 | 위치 |
|---|---|
| FR1 | 1..21 |
| 보존 Cys (앵커) | 22 |
| FR2a | 23..30 |
| CDR-H1 | 31.. (= Cys22 + 9) |
| 보존 Trp (앵커) | CDR-H1 직후 (통상 36) |
| FR2b | Trp+1 .. Trp+13 |
| CDR-H2 | Trp+14 .. |
| FR3 모티프 `[RK].[TSAV][IMLFV][STA]` (앵커) | CDR-H2 직후 (통상 66) |
| FR3 나머지 | 모티프+5 부터 21 잔기 |
| 보존 Cys (앵커) | 통상 92 |
| FR3 꼬리 | Cys+1 .. Cys+2 |
| CDR-H3 | Cys+3 .. (= Cys92 + 3) |
| `WGxG` (앵커) | CDR-H3 직후 |
| FR4 | 그 뒤 7 잔기 |

### VL 좌표 설계

| 구간 | 위치 |
|---|---|
| FR1 | 1..22 |
| 보존 Cys (앵커) | 23 |
| CDR-L1 | 24.. (= Cys23 + 1) |
| 보존 Trp (앵커) | CDR-L1 직후 (통상 35) |
| FR2 | Trp+1 .. Trp+14 |
| CDR-L2 (7 잔기 고정) | Trp+15 .. |
| FR3 | CDR-L2 직후 31 잔기 |
| 보존 Cys (앵커) | 통상 88 |
| CDR-L3 | Cys+1 .. |
| `FGxG` (앵커) | CDR-L3 직후 |
| FR4 | 그 뒤 7 잔기 |

## 생성기의 자기 검증

조립한 서열이 아래를 모두 만족할 때까지 CDR 내용을 재추첨한다. 통제된 입력만 코퍼스에
들어가도록 하는 장치다.

1. 말단 `WGxG`(또는 `FGxG`)가 서열의 **마지막** 매치이고, 그 뒤에 반대쪽 모티프가 없다
   — 체인 판별이 모호해지지 않게.
2. 15..30(VL 은 16..30) 창 안의 Cys 가 앵커 하나뿐이다.
3. CDR1 시작 이후 **첫** Trp 가 의도한 Trp 앵커다.
4. CDR-H2 시작 +9 이후 **첫** FR3 모티프 매치가 의도한 FR3 앵커다 (CDR-H2 안에서
   우연히 먼저 걸리지 않는다).
5. FR3 모티프(VL 은 CDR-L2) 이후 **첫** Cys 가 의도한 두 번째 Cys 앵커다.
6. 각 CDR 안의 liability 모티프가 **주입한 조각 안에 완전히 포함된 것만** 존재한다
   — 그래서 "이 케이스의 CDR liability 는 정확히 이것뿐"을 테스트에서 단정할 수 있다.

framework filler 알파벳에서는 Cys/Trp/Met/Asn/Asp 를 뺐다. 보존 앵커를 유일하게
유지하고, framework 에 우연한 liability 가 생기지 않게 하기 위한 것이다. framework
liability 가 필요한 케이스(`SYN-VH-006`)는 명시적으로 주입한다.

주의: 보존 Trp 앵커(Trp36, `WGxG` 의 Trp)는 성질상 산화 모티프에 해당하므로,
liability 검출기는 이들을 **framework 내 산화 부위**로 보고한다. 정상적인 동작이고,
실제 항체에서도 이 Trp 들은 대개 묻혀 있어 위험도가 낮다.

## 케이스 목록

| ID | 체인 | expected.confident | 설명 |
|---|---|---|---|
| SYN-VH-001 | H | True | 정상 VH, 주입 liability 없음 |
| SYN-VH-002 | H | True | 정상 VH, CDR-H3 18 잔기 |
| SYN-VL-001 | L | True | 정상 VL, kappa 유사 (L1 11 / L3 9) |
| SYN-VL-002 | L | True | 정상 VL, CDR-L1 8 잔기 — 이후 앵커 전체가 밀린다 |
| SYN-VH-003 | H | True | CDR-H2 에 `NG` (deamidation 최고) |
| SYN-VH-004 | H | True | CDR-H1 에 `NLS` (N-glycosylation) |
| SYN-VL-003 | L | True | CDR-L3 에 `NG` |
| SYN-VH-005 | H | True | 도메인 내 Cys 3개 — 비결합 Cys |
| SYN-VH-006 | H | True | framework 에만 `DG` + `M` — CDR 건과 위험도 비교용 |
| SYN-VH-007 | H | True | CDR-H3 에 `DG` (isomerization 최고) |
| SYN-VH-008 | H | True | CDR-H3 에 `M` 과 `W` (산화) |
| SYN-VL-004 | L | True | CDR-L1 에 `DP` (fragmentation) |
| SYN-VL-005 | L | True | CDR 에 Lys/Arg 집중 — pI / net charge 검증용 |
| SYN-VH-901 | H | **False** | 말단 `WGxG` 의 Trp→Ser — 체인 판별부터 실패해야 한다 |
| SYN-VH-902 | H | **False** | Cys22→Ser — CDR-H1 앵커를 못 찾아야 한다 |
| SYN-VH-903 | H | **False** | CDR-H1 9 잔기 — 앵커는 찾히지만 통상 범위(5-7) 밖 |

`SYN-VH-901` / `SYN-VH-902` 는 `expected.cdrs` 가 `null` 이다. 주석기는 위치를
추측하지 말고 `confident=False` 와 한글 사유를 내야 한다. `SYN-VH-903` 은 앵커가
멀쩡하므로 정답 CDR 좌표가 있지만, 길이 이탈 때문에 `confident=False` 가 되어야 한다.

## 파일

- `sequences.json` — 정답 CDR 좌표와 기대값이 붙은 정식 코퍼스
- `sequences.fasta` — 눈으로 보기 쉬운 사본 (정답 좌표는 없다)

### `sequences.json` 레코드 구조

```
id          서열 ID
chain       "H" | "L"
description 한글 설명
sequence    아미노산 서열
length      길이
notes       의도적으로 주입한 liability 목록
expected
  cdrs      {"CDR-H1": {"start", "end", "seq"}, ...} 또는 null(앵커 깨짐)
            start/end 는 1-based, 양끝 포함
  chain     기대 체인 또는 null
  confident 주석기가 confident=True 를 내야 하는지
  anchors_broken  앵커가 일부러 깨진 케이스인지
  length_outlier  CDR 길이가 통상 범위를 벗어나는 케이스인지
anchors_intended  의도한 앵커 위치 (1-based)
```
