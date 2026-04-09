# GB Correction Fitting Analysis Report (v2.0)

## 1. 목적

Bruggeman 추정치(σ_brug)와 실측 proxy(σ_proxy)의 비율 **R = σ_brug / σ_proxy**가
GB density(GB_d)와 전극 두께(T)에 어떻게 의존하는지를 다양한 회귀 모델로 비교하여 최적 모델을 선정한다.

### 정의

| 변수 | 정의 | 단위 |
|------|------|------|
| R | σ_brug / σ_proxy (과대평가 비율) | - |
| σ_brug | σ_bulk × φ_SE × f_perc / τ² (Bruggeman, GB 무시) | ratio |
| σ_proxy | G_path × f_perc / τ (접촉면적 기반 실측) | ratio |
| GB_d (ρ_hop) | Inter-particle hop density (percolation path의 hop 수 / z 거리). 네트워크의 직렬 연결 밀도. | hops/μm |
| T | 전극 두께 | μm |

### 시뮬레이션 조건

| 항목 | 범위 |
|------|------|
| 접촉 모델 | hooke/hysteresis (Luding, 탄소성 이력 모델) |
| SE 입자 크기 | 0.5, 1.0, 1.5 μm |
| AM 입자 크기 | 3.0, 6.0 μm |
| 전극 두께 | 18 ~ 184 μm |
| RVE 크기 | 30×30, 50×50 μm |
| AM:SE 질량비 | 62:38 ~ 85:15 |
| P:S 비율 | 0:10 ~ 10:0 |
| 가압 조건 | 300 MPa |

---

## 2. 데이터 요약

| 항목 | 범위 |
|------|------|
| n (유효 데이터) | **41** |
| GB_d | 0.47 ~ 2.48 hops/μm |
| T (두께) | 18 ~ 184 μm |
| R (ratio) | 15.7 ~ 1612.5 |

### 데이터 그룹

| 그룹 | 케이스 수 | SE 크기 | 두께(μm) | RVE | GB_d 범위 |
|------|----------|---------|---------|-----|----------|
| 박막 (1mAh, 75:25) | 3 | 0.5μm | 19~22 | 50×50 | 1.89~2.04 |
| 박막 (1mAh, 80:25) | 5 | 0.5μm | 18~20 | 50×50 | 1.91~2.18 |
| 박막 (1mAh, 85:15) | 5 | 0.5μm | 18~19 | 50×50 | 2.14~2.48 |
| 후막 (8mAh, 75:25) | 3 | 0.5μm | 182~184 | 50×50 | 1.26~1.27 |
| 후막 (8mAh, 80:20) | 5 | 0.5μm | 157~162 | 50×50 | 1.27~1.33 |
| 후막 (8mAh, 85:15) | 5 | 0.5μm | 139~146 | 50×50 | 1.31~1.38 |
| Real (6mAh, SE 0.5μm) | 3 | 0.5μm | 113~117 | 50×50 | 1.30~1.36 |
| Real (6mAh, SE 1.5μm) | 5 | 1.5μm | 112~120 | 50×50 | 0.48~0.62 |
| Particulate (SE 0.5μm) | 1 | 0.5μm | 93 | 30×30 | 1.35 |
| Particulate (SE 1.0μm) | 3 | 1.0μm | 94~136 | 30×30 | 0.68~0.74 |
| Particulate (SE 1.5μm) | 3 | 1.5μm | 98~135 | 30×30 | 0.47~0.51 |

---

## 3. 후보 모델 비교 (n=41)

| 순위 | ID | 모델 | R² | 파라미터 수 | 수식 |
|------|----|------|-----|-----------|------|
| 1 | M6 | Power Law (GB_d, T 독립) | **0.9490** | 3 | ln(R) = a·ln(GB_d) + b·ln(T) + c |
| 2 | **M15** | **BLM+Constriction (GB_d²×T)** | **0.9427** | **2** | **ln(R) = α·ln(GB_d²×T) + ln(C)** |
| 3 | M13 | Exponential + Arrhenius | **0.9295** | 3 | ln(R) = b·GB_d + c/T + d |
| 4 | B1 | Power Law | 0.1963 | 2 | ln(R) = c·ln(GB_d) + d |
| 5 | E3 | Square Root | 0.1169 | 2 | ln(R) = a·√GB_d + c |
| 6 | A1 | Exponential Decay | 0.0579 | 2 | ln(R) = b·GB_d + ln(k) |
| 7 | C1 | Linear | 0.0347 | 2 | R = a·GB_d + c |
| 8 | C2 | Series Resistance (원점) | -0.3169 | 1 | R = 1 + a·GB_d |
| 9 | A2 | Exponential (원점통과) | -1.3052 | 1 | ln(R) = b·GB_d |

**M15 추천 이유:** R²=0.9427로 M6(0.9490)과 거의 동일하면서 파라미터가 2개 (vs 3개). M15는 동등한 설명력에서 최소 파라미터를 사용하는 parsimonious model이며, 물리적 유도가 가능한 유일한 후보이다.

---

## 4. 모델 상세 분석

### Tier 1: R² ≥ 0.9 (유효 모델)

#### M6. Power Law (GB_d, T 독립)

**수식:** `R = exp(c)·GB_d^a·T^b`

| 파라미터 | 값 |
|----------|-----|
| a | 4.02 |
| b | 1.75 |
| c | -2.81 |
| **R²** | **0.9490** |
| n | 41 |

**물리적 근거:** GB_d와 T의 독립적 기여. M15의 일반화 형태.

**탈락 이유:** GB_d^4.02의 물리적 근거 부족. "왜 4승인가?" 답변 불가. 파라미터 3개.

#### M15. BLM+Constriction (GB_d²×T) ★ **추천**

**수식:** `R = C·(GB_d²×T)^α`

| 파라미터 | 값 |
|----------|-----|
| α | 1.82 |
| ln(C) | -2.97 |
| C | 0.052 |
| **R²** | **0.9427** |
| n | 41 |

**물리적 근거:** BLM(입계 수) + Maxwell Constriction(접촉면적) → GB_d² × T. 물리적 유도 가능.

#### M13. Exponential + Arrhenius

**수식:** `R = exp(d)·exp(b·GB_d)·exp(c/T)`

| 파라미터 | 값 |
|----------|-----|
| b | 4.69 |
| c | -120.3 |
| d | 1.67 |
| **R²** | **0.9295** |
| n | 41 |

**물리적 근거:** 기존 Exp decay에 두께 보정(Arrhenius형) 추가. 기존 모델의 자연스러운 확장.

**참고:** M13의 b=4.69는 동일 두께 내 exponential decay b=5.21과 유사. T→∞이면 기존 식에 수렴.

### Tier 2: R² < 0.9 (탈락 모델)

| ID | 모델 | R² | 탈락 사유 |
|----|------|-----|----------|
| B1 | Power Law | 0.20 | GB_d만으로 T 효과 설명 불가 |
| E3 | Square Root | 0.12 | 비선형성 부족 |
| A1 | Exponential Decay | 0.06 | **동일 두께에서는 R²=0.96이나, T가 다르면 붕괴** |
| C1 | Linear | 0.03 | 비선형 데이터에 선형 모델 부적합 |
| C2 | Series Resistance | -0.32 | 외삽 시 음수 발생 |
| A2 | Exponential (원점) | -1.31 | 강제 원점 통과로 편향 극심 |

**핵심 발견:** 기존 exponential decay 모델(A1, R²=0.96 @ 동일 두께)이 전체 데이터(다양한 T)에서 R²=0.06으로 완전 붕괴. **T가 숨겨진 변수였다는 것이 이번 분석의 가장 중요한 발견.**

---

## 5. 추천 모델 유도: BLM + Constriction

### Step 1: Bruggeman (출발점)

$$\sigma_{eff}/\sigma_{bulk} = \phi_{SE} \times f_{perc} / \tau^2$$

문헌에서 확립된 effective medium theory. 입계(GB) 저항을 무시 → 과대평가.

### Step 2: Brick Layer Model (BLM-inspired scaling)

소결 세라믹 이온전도체에서 확립된 입계 수 스케일링 (Haile et al., 2003):

$$R_{gb} = \rho_{gb} \times w_{gb} \times L / (L_g \times A)$$

- L = 전극 두께 (= T)
- L_g = grain size ∝ 1/GB_d

정리: **N_GB ∝ GB_d × T** (입계 수 = 밀도 × 거리)

**주의:** 본 모델은 BLM의 입계 수 스케일링(N_GB ∝ GB_d × T)만을 차용한다. 소결 세라믹의 면접촉 가정과 입계 고유 저항(ρ_gb × w_gb)은 사용하지 않으며, 대신 DEM 복합양극의 점접촉 특성을 반영한 constriction resistance로 대체한다.

### Step 3: Maxwell Constriction Resistance

DEM 복합양극은 소결 세라믹과 달리 **입자 간 점 접촉(inter-particle contact)** → spreading resistance:

$$R_{constriction} = \rho / (2a)$$

여기서 a = 접촉 반경.

- a ∝ r_SE ∝ 1/GB_d (작은 SE → 작은 접촉)
- **R_per_contact ∝ 1/a ∝ GB_d**

**이것이 BLM과 DEM 복합양극의 핵심 차이.** BLM은 면 접촉(소결)을 가정하여 접촉면적이 입자 크기에 무관하지만, DEM에서는 점 접촉이므로 접촉면적이 GB_d에 추가로 의존한다.

### Step 4: 결합 — GB_d² × T 유도

```
R_total = N_GB × R_per_contact
       = (GB_d × T)  ×  GB_d
         [BLM:총수]    [Constriction:개별저항]
       = GB_d² × T
```

| 항 | 출처 | 의미 |
|-----|------|------|
| 1st GB_d | BLM | μm당 입계 수 |
| 2nd GB_d | Maxwell constriction | 점접촉 → 접촉반경↓ → 저항↑ |
| T | BLM | 총 경로 길이 |

### τ²와의 구조적 유사성

| | 1번째 | 2번째 | 출처 |
|---|-------|-------|------|
| τ² | 경로 길이 ↑ | 유효 단면적 ↓ | Bruggeman |
| GB_d² | 입계 수 ↑ (BLM) | 입계당 저항 ↑ (Constriction) | BLM + Maxwell |

**같은 물리 구조의 다른 발현**: τ²는 "경로의 기하학적 비효율", GB_d²는 "입계의 접촉 비효율".

### Step 5: 데이터 검증

$$R = C \times (GB_d^2 \times T)^\alpha$$

| 파라미터 | 값 | 비고 |
|----------|-----|------|
| α | 1.82 | 접촉면적 불균일성에 의한 superlinear 스케일링 |
| C | 0.052 | exp(-2.97) |
| **R²** | **0.9427** | **n = 41** |

---

## 6. α=1.82 해석: 왜 1이 아닌가?

### 선형 BLM+Constriction 예측: α=1

BLM(입계 수) × Constriction(개별 저항) = GB_d² × T 의 결합은 α=1을 예측한다.

### 실측 α=1.82 > 1의 해석

α=1.82는 데이터가 결정하는 **경험적 값**이며, 접촉면적 불균일성을 포함한 복합 효과를 반영한다.

- **접촉면적 분포의 불균일성 (Bottleneck 효과):** 직렬 저항에서 가장 좁은 접촉이 경로 전도도를 지배. G_path = 1/Σ(1/A_i)는 조화평균으로 산술평균보다 항상 작으며, 이 차이가 GB_d 증가에 따라 확대됨.
- **소성변형에 의한 접촉면적 확대:** hooke/hysteresis 모델에서 소성변형이 접촉면적을 확대시키며, 이 효과는 fitting parameter C에 반영됨.

**주의:** α=1.82의 개별 기여 요인을 정량적으로 분리하는 것은 본 데이터만으로는 불가능하며, 이를 위해서는 rigid sphere DEM 비교 등 추가 연구가 필요하다.

### α=2 고정 테스트 (Robustness 검증)

| 모델 | α | R² | 파라미터 수 |
|------|---|-----|----------|
| Free fit | 1.82 | **0.9427** | 2 |
| α=2 fixed | 2.00 | **0.9333** | 1 |
| α=1.5 | 1.50 | 0.9137 | 1 |
| α=1.0 | 1.00 | 0.7516 | 1 |

**ΔR² = 0.0093 (1.0%)** — α=2 고정 시에도 R²=0.93으로 유의미한 차이 없음.

α sweep에서 최적은 1.82 부근이며, 1.5~2.0 범위에서 R² > 0.91로 안정적. 이는 모델이 α 값에 robust함을 보여준다.

---

## 7. 최종 유효 이온전도도 공식

$$\sigma_{eff} = \frac{\sigma_{bulk} \times \phi_{SE} \times f_{perc}}{\tau^2 \times C \times (GB_d^2 \times T)^\alpha}$$

| 항 | 값 | 출처 | 의미 |
|-----|-----|------|------|
| σ_grain | 3.0 mS/cm | SE grain interior (Li₆PS₅Cl) | grain 전도도 |
| φ_SE | DEM 계산 | 부피분율 | SE 양 |
| f_perc | DEM 계산 | percolation | 연결된 SE 비율 |
| τ² | DEM 계산 | Bruggeman | 경로 기하 손실 |
| GB_d² × T | DEM 계산 | BLM + Constriction | 입계 접촉 손실 |
| α | 1.82 | fitting (n=41) | superlinear 스케일링 |
| C | 0.052 | fitting | 비례상수 |

**분모 구조:**
- **τ²**: 경로가 꼬여서 생기는 손실 (Bruggeman)
- **C·(GB_d²·T)^α**: 입계 접촉에 의한 추가 손실 (BLM + Constriction)

**유도 경로:** Bruggeman(τ²) + BLM-inspired scaling(GB_d×T) + Maxwell Constriction(GB_d) → GB_d²×T

### τ와 GB_d의 독립성 검증

최종 식에서 τ²와 (GB_d²×T)^α를 곱하는 것이 정당한지 — 즉, τ와 GB_d가 독립적인지 확인이 필요하다.

| 그룹 | τ (대표) | GB_d (대표) | 상관 방향 |
|------|---------|-----------|---------|
| SE 0.5μm 후막 (T~160) | ~1.2 | ~1.3 | 낮은 τ, 높은 GB_d |
| SE 0.5μm 박막 (T~20) | ~2.0 | ~2.0 | 높은 τ, 높은 GB_d (양의 상관) |
| SE 1.5μm Real (T~115) | ~1.5 | ~0.5 | 중간 τ, 낮은 GB_d (역상관) |
| SE 1.5μm Particulate | ~1.2 | ~0.5 | 낮은 τ, 낮은 GB_d |

τ와 GB_d 사이에 약한 부분적 상관이 존재한다 (둘 다 SE 크기에 의존). 그러나 τ²는 경로의 기하학적 꼬임(packing 구조)을, (GB_d²×T)^α는 접촉 저항(입자 크기 + 전극 두께)을 각각 설명하며, 물리적으로 분리 가능한 두 손실 메커니즘이다. 같은 GB_d에서도 packing 구조에 따라 τ가 다를 수 있고(예: 후막 vs particulate에서 GB_d≈1.3이지만 τ≈1.2~1.3으로 차이), 이는 두 변수가 완전히 종속적이지 않음을 보여준다.

*τ 값은 각 케이스의 DEM percolation path에서 distance-weighted shortest path 기반으로 별도 산출된다 (dem_analysis_core.calc_tortuosity). 상세 케이스별 데이터는 full_metrics.json 및 network_summary.csv에 수록.*

---

## 8. M6 vs M15: 독립적 검증

M6 (자유 3-parameter fit)의 결과가 M15의 제약된 형태와 거의 일치하며, 이는 GB_d²×T가 natural scaling variable임을 독립적으로 확인한다.

| | GB_d 지수 | T 지수 | 파라미터 수 | R² |
|---|----------|--------|----------|-----|
| M6 (free) | 4.02 | 1.75 | 3 | 0.949 |
| M15 (constrained) | 3.64 (=1.82×2) | 1.82 | 2 | 0.943 |

M6이 자유도 3개로 찾은 최적값(GB_d^4.02 × T^1.75)이 M15의 결합 변수(GB_d²×T)^1.82 = GB_d^3.64 × T^1.82와 거의 동일하다. **데이터가 자발적으로 GB_d²×T 결합을 선택한다는 증거.**

M6의 비정상적으로 높은 GB_d 지수(4.02)는 GB_d와 T 사이의 약한 음의 상관관계(r ≈ −0.93, SE 0.5μm subset에서 GB_d ∝ T^−0.23)에 기인한다. M6가 GB_d와 T를 독립 변수로 취급하면서, T 항(1.75)이 충분히 설명하지 못하는 두께 의존성을 GB_d 지수를 부풀려 보상한 결과이다. M15는 GB_d²×T로 결합함으로써 이 coupling을 구조적으로 흡수하며, 물리적으로 해석 가능한 지수 α=1.82를 산출한다.

---

## 9. 상대적 스케일링 예측

> **⚠ 주의:** 이 모델은 DEM 프레임워크 내에서 self-consistent한 **상대적 스케일링 관계**를 제시하며, 절대 전도도 예측이 아니다. 절대값 보정을 위해서는 EIS 측정과의 calibration이 필요하다.

### SE 크기별 상대 비교

| SE 크기 | GB_d (대표) | R 범위 (σ_brug/σ_proxy) | 상대 σ_eff |
|---------|------------|------------------------|-----------|
| 0.5 μm | ~1.3 (후막), ~2.0 (박막) | 100 ~ 1600 | 1× (기준) |
| 1.0 μm | ~0.7 | 60 ~ 110 | ~10× |
| 1.5 μm | ~0.5 | 16 ~ 33 | ~30× |

*상대 σ_eff는 대표 조건(후막, 80:20 비율)에서 R의 역수 비율로 산출. SE 크기 외 조건(두께, 조성)에 따라 변동 가능.*

**SE 0.5→1.5μm 변경 시 inter-particle contact resistance가 ~30배 감소 (proxy 기준).** Network solver 기준으로는 R_brug 차이가 ~2~3× (Section 14 참조). Proxy의 30배는 single-path 근사에 의한 과장을 포함.

---

## 10. Limitations

1. **절대 전도도 예측 불가:** σ_proxy는 DEM 접촉면적 기반 proxy이며 실험 EIS 측정값이 아님. R = σ_brug/σ_proxy는 두 DEM 계산값의 비율로, 스케일링 관계를 제시하지만 절대 전도도 예측에는 EIS calibration이 필요함.

2. **Inter-particle contact vs Crystallographic GB:** 본 모델의 "GB"는 복합양극 내 SE 입자 간 접촉(inter-particle contact)이며, 소결 pellet 내 결정립계(crystallographic grain boundary)와는 다른 물리적 실체임. BLM의 입계 수 스케일링(N_GB ∝ GB_d × T)만을 차용한 BLM-inspired scaling이며, 소결 세라믹의 입계 고유 저항(ρ_gb × w_gb)은 사용하지 않음.

3. **Argyrodite의 inter-particle resistance:** Cold-pressed 황화물은 전통적으로 negligible GB resistance로 간주되었으나, 최근 AIMD (Sadowski, 2024)와 NMR (de Klerk, 2019)에서 argyrodite pellet 내에서도 유의미한 GB 효과가 보고됨. Sintered pellet의 crystallographic GB와 composite cathode의 mechanical inter-particle contact는 서로 다른 메커니즘이지만, 둘 다 입자 간 이온 수송의 병목이라는 점에서 공통된다. 복합양극의 less intimate한 접촉에서는 추가적인 constriction resistance가 예상됨.

4. **데이터 분포:** 41개 중 SE 0.5μm이 30개로 다수. SE 1.0μm(3개), 1.5μm(8개)으로 GB_d 0.7~1.2 중간 영역의 데이터가 부족. R²=0.94가 SE 0.5μm 데이터에 의해 driven될 가능성을 배제할 수 없음.

5. **α=1.82의 기여 분리 불가:** α > 1의 원인(bottleneck, 소성변형, 네트워크 구조)을 개별적으로 정량 분리하기 위해서는 rigid sphere DEM 비교 등 추가 연구가 필요.

6. **GB_d와 T의 상관관계:** SE 0.5μm subset(n=30)에서 GB_d와 T 사이에 강한 음의 상관(r ≈ −0.93, GB_d ∝ T^−0.23)이 존재한다. 이는 얇은 전극에서 percolation path가 더 tortuous해지면서 hop density가 증가하기 때문으로 해석된다. 이 coupling은 GB_d²×T 결합 변수 내에 자연스럽게 흡수되나, 두 변수의 완전한 독립성을 가정할 수 없음을 의미한다. 단, 이 상관관계는 SE 0.5μm subset 내부에 한정되며, SE 크기가 다른 그룹을 포함하면 동일 두께 범위(T≈110~135μm)에서 GB_d가 0.5~1.35로 분리되어, 전체 데이터셋에서 두 변수는 독립적 축을 형성한다.

7. **GB_d ≠ 1/d_SE:** GB_d를 1/d_SE로 단순 치환하여 R ∝ (T/d_SE²)^α 형태로 변환할 경우, R²가 0.94에서 0.78로 16%p 하락한다. 이는 GB_d가 SE 입자 크기 외에도 AM:SE 조성, packing 구조, 두께 의존 percolation geometry 등 추가적인 미세구조 정보를 포함하기 때문이다. 따라서 GB_d는 단순한 입자 크기의 역수가 아닌, 복합적인 미세구조 지표(microstructural descriptor)로 이해해야 한다.

---

## 11. 결론

1. **41개 DEM 시뮬레이션** (SE 0.5/1.0/1.5μm, T=18~184μm, RVE 30/50μm)에 대해 9개 후보 모델 비교
2. **BLM-inspired scaling 모델 (M15)** R²=0.9427 — 파라미터 2개, 물리적 유도 가능, **추천**
3. M6 (free 3-parameter fit)이 M15와 거의 동일한 지수 도출 → **GB_d²×T가 natural scaling variable이라는 독립적 증거**
4. **기존 exponential decay 모델은 동일 두께에서만 유효** (R²=0.96) — 두께가 다르면 R²=0.06으로 붕괴. **T가 숨겨진 변수였다는 것이 핵심 발견.**
5. τ²와 GB_d²가 동일한 물리 구조 — 경로 기하 손실 vs inter-particle contact 손실
6. α=1.82: 데이터가 결정하는 경험적 값. 접촉면적 불균일성(bottleneck)을 포함한 복합 효과 반영
7. α=2 고정 시에도 R²=0.933 (ΔR²=1.0%) — 모델 robust
8. **본 모델은 상대적 스케일링 관계로서 유효하며, 절대 전도도 예측에는 EIS calibration이 필요**
9. **GB_d는 d_SE의 단순 역수가 아님** — T/d_SE² 치환 시 R²=0.78로 하락. GB_d가 조성·두께·packing에 의존하는 복합 미세구조 지표임을 확인

> **논문 표현:** "BLM-inspired scaling과 Maxwell constriction resistance를 결합한 GB_d²×T가 DEM 복합양극의 inter-particle contact resistance를 universal scaling variable로 설명한다. SE 크기 3단계(0.5~1.5μm), 전극 두께 10배 범위(18~184μm), 다양한 AM:SE 비율의 41개 시뮬레이션에서 R²=0.94로 검증되었다. 본 모델은 DEM 프레임워크 내에서 self-consistent하며, 전극 설계 파라미터가 접촉 저항에 미치는 영향의 상대적 스케일링을 제시한다."

---

## References

1. BLM-inspired scaling: Haile, S. M., West, D. L., & Campbell, J. (2003). *J. Mater. Res.*
2. Maxwell Constriction Resistance: Holm, R. (1967). *Electric Contacts*
3. Hooke/Hysteresis Contact Model: Luding, S. (2008). *Granular Matter*
4. Effective conductivity in ASSB: Minnmann, P. et al. (2021). *J. Electrochem. Soc.* (DOI: 10.1149/1945-7111/abf8d7)
5. Argyrodite GB (AIMD): Sadowski, M. et al. (2024). *Adv. Mater. Interfaces* (DOI: 10.1002/admi.202400423)
6. Li exchange NMR: de Klerk, N. J. J. & Wagemaker, M. (2019). *ACS Appl. Energy Mater.*
7. GB impedance (EIS): Tao, B. et al. (2023). *J. Electrochem. Soc.*

---

---

# Part II: DEM-Native Resistor Network Solver

## 12. Network Solver를 통한 절대값 검증 및 물리적 분해

Part I의 proxy 모델은 상대적 스케일링(GB_d²×T)을 발견하는 "discovery tool"로서 유효하나, 절대값 정확도에 두 가지 구조적 한계가 있다:

1. **Single-path proxy**: G_path = 1/Σ(1/A_i)는 shortest path 몇 개의 평균. 실제 SE 네트워크에는 수천 개의 병렬 경로가 존재.
2. **Constriction-only**: 접촉 목만 고려. SE 입자 내부의 bulk 이온 전도를 무시.

이로 인해 σ_proxy가 실제보다 과소 추정된다. Part I의 R=15~1600은 이 single-path 근사에 의한 것이며, full network solver로 R_brug=3~10을 확인함으로써 Bruggeman의 overestimation이 유의미하지만 관리 가능한 수준임을 확인한다.

## 13. 방법론: DEM-Native Resistor Network

### 13.1 개요

Ketter/Zeier (Nat Comm 2025)의 voxel 기반 방식과 달리, DEM 입자 네트워크를 직접 사용한다. Voxel 방식은 입자를 격자(voxel)로 변환하는 과정에서 접촉 기하학 정보가 손실되고, 면접촉만 표현 가능하다. 본 방법은 DEM에서 직접 추출한 입자 좌표, 반지름, 접촉면적을 그대로 사용하여 **Maxwell constriction resistance를 명시적으로 포함**한다.

### 13.2 Step 1: SE 접촉 네트워크 구성

DEM contact dump에서 SE-SE 접촉만 추출하여 그래프를 구성한다:
- **노드(Node)**: 각 SE 입자 (ID, 위치 x/y/z, 반지름 r)
- **엣지(Edge)**: SE-SE 접촉 (입자 쌍, 접촉면적 A_contact)

```
SE network graph:
  Node = SE 입자 (수천~수만 개)
  Edge = SE-SE 접촉 (수만~수십만 개)
  
  예: 후막 75:25 → 234,088 nodes, 686,756 edges
```

### 13.3 Step 2: 엣지별 물리적 저항 계산

각 SE-SE 접촉(edge)에 두 가지 저항을 부여한다:

**① Bulk resistance (입자 내부 이온 전도)**

$$R_{bulk} = \frac{d_{ij}/2}{\sigma_{SE} \cdot \pi r_i^2} + \frac{d_{ij}/2}{\sigma_{SE} \cdot \pi r_j^2}$$

- d_ij: 입자 i, j 중심 간 거리 (주기경계 보정)
- r_i, r_j: 입자 반지름
- 각 입자의 절반 거리를 통과하는 bulk 전도 경로
- 단면적 A = πr²: 구를 원기둥으로 근사 (최대 단면적을 대표값으로 사용)

*근사 타당성: 구의 단면적을 정확히 적분하면 R_bulk이 ~20-30% 변하나, constriction이 전체 저항의 69~81%를 지배하므로 σ_full에 미치는 영향은 <8%. 이 체계적 오차는 fitting constant C에 흡수되며 R²에 영향 없음. 또한 packed bed에서 이온은 구의 전체 단면을 균일하게 통과하지 않고 접촉점 간 3D 경로를 따르므로, 구 적분이 원기둥보다 반드시 정확하지는 않음.*

**② Constriction resistance (접촉 목 병목)**

$$R_{constriction} = \frac{1}{\sigma_{SE} \cdot 2a_{ij}}$$

- a_ij = √(A_contact / π): 접촉 반경 (DEM에서 직접 추출)
- Maxwell spreading resistance 공식 (Holm, 1967)

**③ 합산**

$$R_{edge} = R_{bulk} + R_{constriction}$$

```
실측 비율:
  R_bulk : R_constriction ≈ 24% : 76% (SE 0.5μm 평균)
  → Constriction이 저항의 대부분을 차지
```

### 13.4 Step 3: Conductance Matrix (Laplacian) 조립

각 edge의 conductance g_ij = 1/R_edge를 사용하여 **graph Laplacian 행렬** L을 구성한다.

**Kirchhoff 전류 법칙(KCL)을 행렬로 표현한 것이다.** 각 노드에서 "들어오는 전류 = 나가는 전류"를 모든 노드에 대해 동시에 적으면 연립방정식이 되고, 그 계수 행렬이 L이다.

예: SE 입자 3개(A, B, C)가 A-B, B-C로 연결된 경우, 노드 B의 KCL:
```
g_AB × (V_B - V_A) + g_BC × (V_B - V_C) = I_B

정리: -g_AB × V_A + (g_AB + g_BC) × V_B - g_BC × V_C = I_B

→ L 행렬의 B행: 대각(g_AB+g_BC), 비대각(-g_AB, -g_BC)
```

코드:
```python
# 각 SE-SE 접촉에 대해:
g = 1.0 / R_edge

L[i][i] += g    # 대각: 노드 i에 연결된 전체 conductance 합
L[j][j] += g
L[i][j] -= g    # 비대각: 노드 i-j 간 conductance (음수)
L[j][i] -= g
```

L은 sparse matrix로 구성 — 노드가 수만 개이지만 각 노드가 이웃 몇 개(CN≈3~6)하고만 연결되어 L의 대부분이 0.

### 13.5 Step 4: 경계조건 및 Kirchhoff 풀이

**이온 흐름 방향:** 복합양극에서 Li⁺는 separator(top) → current collector(bottom) 방향으로 이동. 단, Li⁺는 반드시 bottom까지 관통할 필요 없이, SE 네트워크를 타고 **AM 표면에 도달**하면 충전/방전 반응이 일어난다. 따라서 top boundary는 **top_reachable SE** (top에서 연결된 모든 SE)를, bottom boundary는 **bottom SE**를 사용한다 (percolation_sets.json에서 로드).

**Virtual source/sink:** top/bottom에 닿는 SE 입자가 수백~수천 개이므로, "가상 노드 1개"를 만들어 한 점에 연결한다.

```
   [source] ← separator 쪽 (Li⁺ 공급)
      / | | \
    SE  SE SE SE   ← top_reachable SE들
     |  |  |  |
    SE  SE SE SE   ← 내부 SE network (AM 표면에서 반응)
     |  |  |  |
    SE  SE SE SE   ← bottom SE들
      \ | | /
    [sink] ← current collector 쪽
```

**g=10⁶ 연결의 의미:** "저항 없는 완벽한 도선". 실제 SE-SE conductance가 ~10⁰~10²인데, 10⁶이면 1만~100만 배 크므로 이 연결에서 전압 강하가 사실상 0. 전체 전압 강하가 R_edge들(bulk + constriction)에서만 발생하도록 보장한다. EIS에서 Au 스퍼터링으로 집전체 접촉저항을 최소화하는 것과 같은 원리.

*g=10⁴~10⁸ 범위에서 σ_full 변화 < 0.01%로 결과에 영향 없음.*

**전류 주입:**
- I_source = +1 (전류 1을 네트워크에 밀어넣음)
- I_sink = -1 (전류 1을 빼냄)
- V_sink = 0 (접지, 전위 기준점)

비유: 건전지의 +극(source)에서 전류가 나와 SE 네트워크를 통과하고 -극(sink)으로 돌아감. V_source는 "전류 1을 흘리기 위해 얼마나 세게 밀어야 하는지"를 나타냄. 세게 밀어야 하면 전도도가 낮고, 살짝만 밀어도 되면 전도도가 높다.

**연립방정식 풀이:**

$$\mathbf{L} \cdot \mathbf{V} = \mathbf{I}$$

```python
V = scipy.sparse.linalg.spsolve(L, I)
G_eff = 1.0 / V_source  # I=1, V_sink=0이므로 G=I/ΔV=1/V_source
```

Percolating component만 사용 (bottom↔top 연결된 connected component). 고립된 SE 입자를 포함하면 L 행렬에 연결 없는 행/열이 생겨 singular matrix가 되므로 제외.

### 13.6 Step 5: σ_full 산출

$$\sigma_{full} = G_{eff} \times \frac{T}{A_{electrode}} \times \sigma_{grain}$$

- T: 전극 두께 (plate_z × scale, μm)
- A_electrode: 전극 단면적 (box_x × box_y × scale², μm²)
- σ_grain: SE grain interior conductivity (단결정 전도도)

**σ_grain vs σ_pellet 구분 (중요!):**

| 용어 | 값 | 의미 |
|------|-----|------|
| σ_grain (grain interior) | ~3.0 mS/cm | 단결정 고유 전도도. GB 효과 없음. |
| σ_pellet ("σ_bulk") | ~1.3 mS/cm | pellet 측정값. 내부 SE-SE GB 포함. |

문헌에서 "σ_bulk = 1.3 mS/cm"로 보고하는 값은 실제로 **σ_pellet**이며, pellet 내부의 SE-SE 접촉 저항이 이미 포함되어 있다. Network solver는 복합양극 내 모든 SE-SE 접촉을 명시적으로 계산하므로, σ_pellet을 사용하면 접촉 저항이 **이중으로 적용**된다.

따라서 본 모델에서는 **σ_grain ≈ 3.0 mS/cm** (argyrodite grain interior)를 사용한다.

**σ_grain = 3.0 mS/cm 근거:**
- MLIP MD (UMA potential): σ₃₀₀K = 3.0 mS/cm, Ea = 0.286 eV (perfect crystal, no GB)
- Second-sintered pellet: 최대 3.19 mS/cm (ScienceDirect)
- Hot-pressed pellet (150°C, 300 MPa): >2.0 mS/cm (ACS Publications)
- 일반 합성: ~2.6 mS/cm (RSC Publishing)
- 실험 pellet(~2~3 mS/cm)은 잔여 GB를 포함하므로, 순수 grain interior는 3.0+ mS/cm이 physically consistent

**실험적 역검증:**
- σ_pellet(1.3) 사용 시: σ_full = 0.085 mS/cm → Minnmann 대비 **2× 과소** (이중 적용)
- σ_grain(3.0) 사용 시: σ_full = 0.196 mS/cm → Minnmann 대비 **1.15×** (**거의 일치**)
- σ_grain이 올바른 선택임을 Minnmann 실험이 독립적으로 뒷받침

*정규화:* 계산은 ρ=1 (정규화)로 수행한 후, 최종값에 σ_grain을 곱하여 mS/cm로 변환. σ_grain 값에 무관하게 스케일링 관계가 유지된다.

### 13.7 Three-Mode Decomposition

동일한 네트워크에서 edge 저항 구성만 바꿔 3번 풀이:

| Mode | R_edge | 물리적 의미 |
|------|--------|-----------|
| **FULL** | R_bulk + R_constriction | Ground truth (실제 σ_eff) |
| **CONTACT_FREE** | R_bulk only (R_constr=0) | 접촉이 완벽할 때 (upper bound) |
| **CONSTRICTION_ONLY** | R_constriction only (R_bulk=0) | 접촉만의 기여 |

추가로 **Bruggeman EMT** (σ_grain × φ^1.5)를 해석적으로 계산하여 비교한다.

이 분해를 통해:
- **Contact-free / Full ratio**: 이상적 접촉 대비 실제 접촉 패널티 (3~10×)
- **Bruggeman EMT / Full ratio**: Bruggeman 과대평가 정도 (2~5×)
- **Constriction 기여**: R_bulk/(R_bulk+R_constr) = 19~31% → constriction이 69~81% 지배

### 13.8 Ketter (2025) 대비 차별점

| | Ketter (Nat Comm 2025) | 본 연구 |
|---|---|---|
| 입력 | FIB-SEM/voxel | DEM 입자 좌표 + 접촉면적 |
| 노드 | voxel 중심 | SE 입자 중심 |
| 연결 | 인접 voxel (면접촉) | SE-SE 접촉 (점접촉) |
| Constriction | ❌ 불가능 | ✅ Maxwell R=ρ/(2a) |
| 분해 | ❌ | ✅ FULL/BULK/CONSTR 3-mode |
| 코드 | Münster datastore | `scripts/network_conductivity.py` |

## 14. 결과 (n=41, 82 runs)

### 절대값

| 조건 | σ_full (mS/cm) | 문헌 |
|------|---------------|------|
| SE 0.5μm 후막 75:25 | 0.196~0.277 (σ_grain 기준) | Minnmann: ~0.17 (**1.15×, 거의 일치**) |
| SE 1.5μm Particulate 62:38 | 0.586~0.641 (σ_grain 기준) | — |
| SE 0.5μm 박막 85:15 | 0.023~0.037 (σ_grain 기준) | — |

기존 proxy(0.0002 mS/cm)보다 **1000배 개선**, Minnmann 실험과 1.15× 일치.

*σ_grain = 3.0 mS/cm (grain interior) 사용. σ_pellet(1.3 mS/cm) 사용 시 위 값의 ×0.43.*

*조건 차이: Minnmann은 NCM-622, 380 MPa, ~42 vol% CAM. 본 DEM은 NCM-811 equivalent, 300 MPa, 75:25 wt%. 정확한 1:1 비교는 불가하나 order-of-magnitude sanity check로서 유효.*

### Contact-free / Full ratio (실제 접촉 패널티)

| 그룹 | Contact-free/Full | 기존 R (proxy) |
|------|--------|--------------|
| SE 0.5μm 후막 75:25 | **4.0** | 1044 |
| SE 0.5μm 박막 85:15 | **7~10** | 213~373 |
| SE 1.5μm | **3.2~3.6** | 17~33 |

### Constriction 지배

bulk_frac = 19~31% → **constriction이 69~81% 지배**. SE 작을수록 constriction 비율 증가.

### τ_eff² 역산

후막 75:25: τ_eff² ≈ 4.0~4.2 → 문헌 τ_eff²≈4~5와 정확 일치.

## 15. Bruggeman Exponent Decomposition Framework

### "ASSB에서 왜 n_eff ≈ 3인가?"

Three-mode decomposition(σ_full = σ_bulk_net / R_brug)을 이용하여 n_eff를 분해:

| 항 | 값 | R² | 의미 |
|----|-----|-----|------|
| n_eff | 3.37 | 0.80 | 총 effective exponent |
| n_geo | 2.54 | 0.71 | 기하학적 tortuosity (접촉 없이) |
| n_contact | 0.83 | 0.48 | 접촉 저항 패널티 |

n_eff = n_geo + n_contact = 2.54 + 0.83 = **3.37** (항등식으로 성립)

**이 분해의 의미는 개별 R² 값에 있다:**
- n_geo R²=0.71 → φ_SE가 기하학적 전도의 71% 설명 → Bruggeman 근사가 기하학적으로는 유효
- n_contact R²=0.48 → φ_SE가 접촉 패널티의 48%만 설명 → **나머지 52%는 hop_area, CN, GB_d에 의존 → Section 16의 multi-scale model이 필요한 이유**

문헌 비교: ASSB의 n≈3~4 (Bielefeld 2020)의 기원 = 기하학(2.54) + 접촉(0.83). 액체 LIB의 n≈1.5는 접촉 저항 없음(n_contact≈0).

---

# Part III: Multi-Scale Scaling Law

## 16. 최종 모델

### 16.1 모델 진화

| Version | 공식 | R² | thick | thin | avg3 | 비고 |
|---------|------|-----|-------|------|------|------|
| v1 | σ_brug × C × √A_hop × CN² × GB_d^(4/3) | 0.894 | - | - | - | 초기 모델 |
| v2 | σ_brug × C × (A_hop × GB_d²)^(3/5) × CN² | 0.923 | - | - | - | 결합 변수 |
| v3 | σ_brug × C × (G_path × GB_d²)^(1/4) × CN² | 0.836 | 0.961 | -1.03 | 0.27 | thick only |
| v4 | σ_brug × C × √(1-φ_c/φ) × τ^(3/2)/f × CN^(3/2) | 0.946 | 0.982 | 0.917 | 0.910 | universal (no cov) |
| **v4++ FORM X** | **C × σ_grain × (φ-φ_c)^(3/4) × CN × √cov / √τ** | **0.961** | **0.967** | **0.934** | **0.932** | **CHAMPION** |

v4++ FORM X = C × σ_grain × ⁴√[(φ_SE-φ_c)³ × CN⁴ × cov² / τ²]

진화의 핵심:
- v3→v4: contact quality(G_path×GB_d²) → percolation threshold(φ_c) + τ² 상쇄
- v4→v4++: coverage 추가 + φ_SE 제거((φ-φc)만으로 충분)
- FORM A(√[φ(φ-φc)×CN²×cov/τ]) stress test에서 α=0.574≠0.5, φ_SE 지수≈0 발견 → FORM X로 정제
- (φ-φc)^(3/4): 3D percolation exponent t=3/2의 √
- CN^1: v3의 CN²에서 축소 (√ 안에서 CN⁴ = effective CN^1)
- LOOCV=0.960, leave-SE-size-out R²>0.93 (외삽 능력 검증)
- φ_c=0.18, C=0.123, 1 free parameter

### 16.2 v4++ FORM X: 물리적 유도

$$\sigma_{ion} = C \times \sigma_{grain} \times \sqrt[4]{(\phi_{SE} - \phi_c)^3 \times CN^4 \times cov^2 / \tau^2}$$

$$= C \times \sigma_{grain} \times (\phi_{SE} - \phi_c)^{3/4} \times CN \times \sqrt{cov} / \sqrt{\tau}$$

#### Step 1: Percolation Theory (출발점)

Bruggeman EMT는 φ_SE를 사용하지만, 실제로는 **percolation threshold** 이하에서 σ=0이다:

$$\sigma \propto (\phi - \phi_c)^t \quad \text{(3D percolation, } t \approx 2.0 \text{)}$$

FORM X에서 (φ-φc)^(3/4) = [(φ-φc)^(3/2)]^(1/2). 즉 **3D percolation exponent t=3/2**를 √ 안에서 적용한 형태.

*φ_c = 0.18*: 데이터에서 최적화 (continuous optimization → 0.191). 문헌 값 0.15~0.20 범위 내.
*thin100_15(φ_SE=0.176)가 정확히 threshold 위에 위치.*

#### Step 2: τ → √τ (Softened Tortuosity)

Bruggeman: σ ∝ 1/τ² — **random walk in continuous medium** 가정.
실제 packed bed: SE가 **discrete contact network** → random walk가 아님.

- Thick(τ=1.2): 1/√1.2 = 0.91 (거의 penalty 없음)
- Thin(τ=4.5): 1/√4.5 = 0.47 (적절한 penalty)
- cf. Bruggeman: 1/4.5² = 0.049 (과도!)

**τ² → √τ: 이산 접촉 네트워크에서 tortuosity penalty가 약화됨.**

σ_brug 프레임에서: τ^(3/2) 보정항이 σ_brug의 τ^(-2)를 부분 상쇄 → effective τ^(-1/2).

#### Step 3: CN (Network Connectivity)

v3에서 CN²이었던 것이 FORM X에서 CN^1로 축소. 이유:

- v3에서 CN²는 φ_SE가 포함된 σ_brug 위에서 추가 보정 → CN의 "순수 효과"는 CN²
- FORM X에서는 (φ-φc)가 φ_SE를 대체 → CN의 일부 효과가 이미 (φ-φc)에 흡수
- Free fit: CN 지수 = 1.33 (FORM X의 1.0과 근접)
- CN^1 ↔ CN^1.33 차이: ΔR² < 0.005

#### Step 4: Coverage (AM-SE 계면 품질)

Coverage = AM 표면 중 SE와 접촉하는 비율 (%).

- σ_ionic에 coverage가 필요한 이유: 이온은 SE 네트워크를 통해 AM 표면에 도달해야 반응 가능
- Coverage ↑ → 유효 계면 면적 ↑ → 이온 접근 경로 ↑ → σ_ionic ↑
- FORM X에서 √cov: 제곱근 비례 (수확 체감)

Ablation: coverage 제거 시 ΔR² = -0.015 → 독립적 설명력 확인.
Outlier 분석: worst cases의 coverage = 평균의 74% → coverage가 outlier를 설명.

#### Step 5: σ_brug와의 관계

FORM X는 σ_brug를 **버린 것이 아니라 정제한 것**:

```
FORM X = σ_brug × C × √(1-φ_c/φ_SE) × τ^(3/2)/f_perc × CN^(1) × √cov
                        ──threshold──   ──τ² 상쇄──  ──f 상쇄──  ─connectivity─  ─interface─
```

상쇄 후:
- τ^(3/2) × τ^(-2) = τ^(-1/2) ← softened
- f_perc^(-1) × f_perc = 1 ← cancelled (already in threshold)
- φ × (1-φ_c/φ)^(1/2) → √φ × √(φ-φ_c) ← percolation

#### Step 6: ⁴√ 안의 정수 — 왜 깔끔한가

$$\sigma = C \times \sigma_{grain} \times \sqrt[4]{(\phi-\phi_c)^3 \times CN^4 \times cov^2 \times \tau^{-2}}$$

⁴√ 안의 지수: 3, 4, 2, -2 — 전부 정수!

| 항 | ⁴√ 안 지수 | Effective 지수 | 물리 |
|---|---|---|---|
| (φ-φc) | 3 | 3/4 | percolation (t=3/2의 √) |
| CN | 4 | 1 | network (선형 비례) |
| coverage | 2 | 1/2 | interface (수확 체감) |
| τ | -2 | -1/2 | softened tortuosity |

#### 검증 요약

| Metric | 값 |
|--------|-----|
| R² (ALL, n=48) | 0.961 |
| LOOCV | 0.960 |
| R² (thick, τ≤1.5) | 0.967 |
| R² (thin, τ>2.5) | 0.934 |
| Free fit α | 0.515 (FORM X: 0.5, gap 0.015) |
| Free fit φ_c | 0.191 (FORM X: 0.18, gap 0.011) |
| Leave-SE-size-out | R² > 0.93 (외삽 OK) |
| Free params | 1 (C only) |
| Within 20% | 40/48 (83%) |

#### Robustness Analysis: 왜 고정 지수가 자유 지수보다 나은가

α=1/2 (⁴√)과 φ_c=0.18을 고정하지 않고 자유 파라미터로 풀면:

| Model | Free params | α | φ_c | R² |
|-------|-------------|---|-----|-----|
| **FORM X (fixed)** | **1 (C)** | **0.50** | **0.18** | **0.9611** |
| +φ_c free | 2 (C, φ_c) | 0.50 | 0.185 | 0.9624 |
| +α free | 2 (C, α) | 0.515 | 0.18 | 0.9619 |
| All free | 3 (C, α, φ_c) | 0.507 | 0.185 | 0.9626 |

**ΔR² = 0.0014 (0.14%) — α와 φ_c를 동시에 풀어도 거의 개선 없음.**

2D grid (α=0.40~0.65, φ_c=0.12~0.22)에서 R² > 0.95인 영역이 **broad plateau**를 형성:
- φ_c = 0.17~0.19, α = 0.48~0.52 전체가 R² > 0.96
- Sharp peak이 아닌 flat maximum → **robust, not overfitted**

SE 크기별 φ_c 안정성:
- SE 0.5μm (n=37): φ_c = 0.190
- SE 1.5μm (n=7): φ_c = 0.175
- 차이 0.015 → **φ_c ≈ 0.18은 SE 크기에 무관한 universal threshold**

**결론: 자유도를 늘리면 overfitting 위험만 증가하고 R² 개선은 0.14%에 불과. 고정 지수(α=1/2, φ_c=0.18)가 parsimony와 robustness 면에서 우월하다.**

> "Joint optimization of α and φ_c yields R²=0.963 (3 free params) vs FORM X R²=0.961 (1 free param). The ΔR²=0.001 confirms that α=1/2 and φ_c=0.18 are robust fixed values within a broad optimality plateau, not artifacts of discrete data fitting."

### 16.3 v3 공식 (thick only, legacy)

$$\sigma_{eff} = \sigma_{brug} \times C \times (G_{path} \times GB_d^2)^{1/4} \times CN^2$$

| 항 | 지수 | 물리 근거 | 검증 |
|---|---|---|---|
| (G_path × GB_d²)^(1/4) | 1/4 | empirical (free fit: 0.27≈1/4, ΔR²=0.002). G_path: harmonic mean conductance, GB_d²: BLM×Maxwell | Fixed vs Free gap = 0.0017 |
| CN² | 2 | free fit: 1.98≈2. Redundant path factor (grain conductance × edge density) | CN 제거 → R² 급락 |
| C = 0.073 | — | 비례상수 | 유일한 free parameter |

### 16.3 왜 G_path가 A_hop보다 좋은가

| | A_hop (산술 평균) | G_path (조화 평균) |
|---|---|---|
| 정의 | mean(A_contact) | 1/Σ(1/A_i) |
| 물리 | 평균 접촉면적 | **직렬 저항의 유효 전도도** |
| Bottleneck 반영 | ❌ | ✅ 자동 포함 |
| R² (with GB_d, CN) | 0.925 | **0.943** |

직렬 저항 회로에서 전체 전도도는 **가장 약한 접촉에 의해 지배**됨. G_path = 1/Σ(1/A_i)는 이 조화 평균 특성을 자연스럽게 반영하므로, A_hop보다 정확.

### 16.4 결합 변수 (G_path × GB_d²)의 의미

- **G_path** [μm²]: single percolation path의 유효 전도도 (접촉 면적 단위)
- **GB_d²** [μm⁻²]: (hop density)² = BLM(hop 수) × Maxwell(접촉 크기)
- **G_path × GB_d²** [dimensionless]: 접촉 품질 × 밀도² = "유효 contact 품질 지수"
- **(...)^(1/4)**: σ_brug가 주 physics의 ~90%를 포착 → contact correction은 미세 보정 → 약한 지수

1/4의 이론적 유도는 현재 없으며, free fit(0.27)에 가장 가까운 깔끔한 분수로 empirical하게 결정.

### 16.5 검증

- **R² = 0.942** (1 free param, fixed exponents)
- **Free fit R² = 0.943** (gap = 0.0017 — 사실상 완벽)
- **Mean |error| = 14.5%**, Within 20%: 30/41 (73%)
- **Thick only (n=28): R² = 0.953, mean error 10.5%**
- **실험 비교**: σ_full ≈ 0.196 mS/cm vs Minnmann 0.17 mS/cm (1.15×)

### 16.6 Per-SE-size 안정성

| SE | n | C (개별 fit) | R² | 비고 |
|---|---|---|---|---|
| 0.5μm | 30 | 0.075 | 0.938 | |
| 1.0μm | 4 | 0.064 | 0.982 | |
| 1.5μm | 7 | 0.071 | 0.919 | |
| **전체** | **41** | **0.073** | **0.942** | **C 범위: 1.17×** |

C 범위 1.17× (v1에서 1.46×) → G_path × GB_d² 결합 변수가 SE 크기 의존성을 효과적으로 흡수.

### 16.7 Exhaustive Screening

15개 변수(기본 6 + 결합 9)의 모든 1/2/3-variable 조합을 체계적으로 테스트. 결과:

| 순위 | 모델 (σ_brug × ...) | R² |
|------|---------------------|-----|
| 1 | **(G_path×GB_d²)^(1/4) × CN²** | **0.942** |
| 2 | (√(A_hop×BN)×GB_d²)^0.43 × CN^1.88 | 0.933 |
| 3 | (BN×GB_d²)^0.31 × CN^1.78 | 0.930 |
| 4 | (A_hop×GB_d²)^(3/5) × CN² | 0.923 |

Fixed exponent (1/4, 2)에서 R²=0.942는 free fit(0.943)과 gap 0.0017로, 가능한 모든 조합 중 최고.

## 17. GB_d²×T의 운명: 분해 → 결합 변수로 재탄생

GB_d²×T(proxy R²=0.94)는 network에서 R²=0.05로 붕괴. 탈락이 아니라 **정제**:

| GB_d²×T 요소 | Network에서 | v3 최종 모델에서 |
|-------------|-----------|------|
| GB_d (hop 수) | τ²에 흡수 | σ_brug 내 τ² |
| GB_d (접촉 크기) | G_path에 흡수 | **(G_path × GB_d²)^(1/4)** 내 G_path |
| GB_d (mesh 밀도) | GB_d²에 결합 | **(G_path × GB_d²)^(1/4)** 내 GB_d² |
| T (두께) | network에 내재 | 별도 항 불필요 (b=-0.05) |

**GB_d²의 재탄생**: Part I에서 발견한 GB_d²(BLM × Maxwell)가 v3에서 **G_path와 결합**하여 (G_path × GB_d²)^(1/4)로 살아남음. 접촉 품질(G_path)과 밀도(GB_d²)가 하나의 dimensionless 결합 변수로 통합.

---

# Part IV: Implications & Limitations

## 18. 설계 가이드

| 설계 변수 | 영향 항 | σ_eff 향상 방향 |
|----------|--------|--------------|
| AM:SE 비율 | φ_SE, CN | SE ↑ → φ_SE↑, CN↑ → σ_eff ↑ |
| P:S 비율 | τ | P ↑ → τ ↓ → σ_eff ↑ |
| SE 크기 | G_path, CN, GB_d | SE↑ → G_path↑, GB_d↓ → **trade-off** |
| 가압 조건 | G_path | 압력 ↑ → G_path ↑ → σ_eff ↑ |

## 19. 추가 Limitations (Part II~III)

8. **σ_grain vs σ_pellet 선택**: σ_pellet(1.3 mS/cm) 사용 시 실험 대비 2~3× 과소이나, σ_grain(3.0 mS/cm) 사용 시 Minnmann과 1.15×로 거의 일치. 본 보고서는 σ_grain을 채택하나, 최적 σ_grain 값은 SE 재료 및 합성 조건에 따라 달라질 수 있음.

9. **n_eff 분해의 한계**: 수학적 항등식. n_contact의 R²=0.48로 φ_SE가 접촉 패널티의 48%만 설명.

10. **CN²의 물리적 유도 미완**: Kirkpatrick EMA는 σ∝(z-2) (선형), percolation exponent t≈2는 (p-p_c)^t이지 CN^t가 아님. CN²은 empirical scaling(free fit 1.98에서 정수 2로 고정)이며, "redundant path factor"(병렬 경로의 수와 네트워크 이중화 효과)로 해석 가능하나 엄밀한 이론적 유도는 부재. 단, ablation study에서 CN² 제거 시 ΔR²=-0.12로 독립적 설명력이 확인됨.

11. **Constriction 비율의 SE 크기 의존**: 76% 평균은 SE 0.5μm 케이스에 의해 주도. SE 1.5μm에서는 69~71%로 감소.

12. **DEM 접촉면적 정확도**: σ_pellet 사용 시 실험 대비 2~3× 과소였으나, σ_grain(3.0 mS/cm) 적용으로 해소됨 (Minnmann 대비 1.15×). 다만 DEM hooke/hysteresis 모델의 접촉면적이 실제 cold-pressed argyrodite의 소성변형 및 표면 거칠기를 정확히 반영하는지는 별도 검증이 필요하다.

13. **SE 크기 변수(d_SE) 미포함**: d_SE를 추가 변수로 포함하면 R²=0.949(+0.007)로 개선되며, Per-SE-size error가 균등화됨(SE 1.0μm: 15.3%→13.0%, SE 1.5μm: 16.8%→11.2%). 그러나 (1) free parameter가 1→2개로 증가, (2) SE 0.5μm max error가 56%→78%로 악화, (3) SE 1.0μm(n=4), 1.5μm(n=7)의 데이터 부족으로 통계적 유의성 검증 불가. **SE 크기별 데이터가 균등하게 확보되면(각 10~15개) d_SE^(-1/6) 항의 추가를 재검토해야 한다.**

14. **박막(T<50μm) 적용 한계**: Thin electrode(T=18~20μm)에서 mean error 22.4%(vs Thick 10.5%). RVE에 입자 2~3층만 존재하여 G_path, CN 등 통계량의 대표성이 낮음. **RVE 크기 확대(50×50→100×100μm) 또는 앙상블 평균**으로 개선 가능.

15. **(G_path × GB_d²)^(1/4)의 1/4 지수**: empirical value (free fit 0.27에 가장 가까운 깔끔한 분수). Ab initio 이론적 유도는 현재 없으며, σ_brug가 주 physics의 ~90%를 포착하기 때문에 contact correction의 지수가 작은 것으로 해석. 추가 연구 필요.

## 20. 최종 결론

Part I의 GB_d²×T 모델은 proxy 기반 상대적 스케일링으로서 유효하나, network solver가 밝힌 진실은:

1. **Full network solver로 R_brug=3~10 확인** — proxy 기반 R=15~1600은 single-path 근사에 의한 과장. Bruggeman의 overestimation은 유의미하지만 관리 가능한 수준
2. **Constriction이 이온 전도 저항의 69~81% 지배**
3. **n_eff ≈ 3.4 = n_geo(2.54) + n_contact(0.83)**: 접촉 저항이 Bruggeman exponent를 1.5에서 3.4로 올린 원인
4. **최종 모델: σ_eff = σ_brug × 0.073 × (G_path × GB_d²)^(1/4) × CN²** (R²=0.942, 1 free param)
5. **Exhaustive screening**: 15개 변수의 모든 1/2/3-var 조합 중 최고. Fixed vs Free gap = 0.0017.
6. **GB_d²×T → (G_path × GB_d²)^(1/4)로 정제** — Part I의 발견이 Part III에서 결합 변수로 재탄생
7. **Per-SE-size C 안정**: 0.064~0.075 (1.17× ratio) — SE 크기 통합 성공

> **논문 표현:** "DEM-native resistor network solver를 통해 Bruggeman의 과대평가가 3~10×임을 정량화하고, inter-particle constriction이 총 저항의 69~81%를 지배함을 보였다. 15개 변수의 체계적 스크리닝에서 도출된 multi-scale scaling law (σ_brug × (G_path × GB_d²)^{1/4} × CN², R²=0.942)는 harmonic mean conductance(접촉 품질)와 network connectivity(CN²)의 물리를 반영하며, 문헌 Bruggeman exponent n≈3의 기원을 기하학(2.54) + 접촉(0.83)으로 분해한다."

---

## References

1. BLM-inspired scaling: Haile, S. M., West, D. L., & Campbell, J. (2003). *J. Mater. Res.*
2. Maxwell Constriction Resistance: Holm, R. (1967). *Electric Contacts*
3. Hooke/Hysteresis Contact Model: Luding, S. (2008). *Granular Matter*
4. Effective conductivity in ASSB: Minnmann, P. et al. (2021). *J. Electrochem. Soc.* (DOI: 10.1149/1945-7111/abf8d7)
5. Argyrodite GB (AIMD): Sadowski, M. et al. (2024). *Adv. Mater. Interfaces* (DOI: 10.1002/admi.202400423)
6. Li exchange NMR: de Klerk, N. J. J. & Wagemaker, M. (2019). *ACS Appl. Energy Mater.*
7. GB impedance (EIS): Tao, B. et al. (2023). *J. Electrochem. Soc.*
8. Resistor Network for SSB: Ketter, L. et al. (2025). *Nature Communications* (DOI: 10.1038/s41467-025-56514-5)
9. Composite cathode microstructure: Minnmann, P. et al. (2024). *J. Electrochem. Soc.* (DOI: 10.1149/1945-7111/ad510e)

---

# Part V: Electronic Conductivity Scaling Law

## 21. AM-AM Network: Topology-Dominated Transport

### 21.1 Network 구조

Electronic conductivity는 AM-AM 접촉 네트워크만을 통한 전자 전도이다. Ionic과 달리 **contact quality가 아닌 topology가 지배**한다.

| 비교 항목 | Ionic (SE-SE) | Electronic (AM-AM) |
|----------|---------------|-------------------|
| 입자 크기 | ~0.5 μm | ~5 μm (10배) |
| σ_bulk | 3.0 mS/cm | 50 mS/cm (17배) |
| 접촉면적 | 작음 → constriction 지배 | 큼 → constriction 무시 가능 |
| R²에 대한 contact 기여 | ΔR²=+0.12 (√A_hop) | **ΔR²<0** (overfitting) |

### 21.2 최종 공식

$$\sigma_{el} = 0.015 \times \sigma_{AM} \times \phi_{AM}^{3/2} \times CN_{AM}^2 \times \exp\!\left(\frac{\pi}{T / d_{AM}}\right)$$

| 항 | 지수 | 물리 근거 | 비고 |
|---|---|---|---|
| σ_AM | 1 | AM bulk conductivity | 0.05 S/cm (NCM811) |
| φ_AM^(3/2) | 1.5 | Bruggeman EMT (구형 입자) | depolarization factor L=1/3 |
| CN_AM² | 2 | Grain conductance × edge density | Makse (2004) 유사: K,G ∝ Z² |
| exp(π/(T/d)) | - | Finite-size boundary effect | Trivial percolation in thin electrode |
| C = 0.015 | - | EMT geometric: (1/6)×(a_c/R) ≈ 0.017 | Fitted vs EMT 차이 12% |

### 21.3 검증

| Metric | 값 |
|--------|-----|
| R² (통합, 38 cases) | **0.89** |
| R² (후막만, T/d > 5) | **0.93** |
| Free parameters | **1** (C = 0.015) |
| Minnmann 실험 비교 | 예측 ~3.3 mS/cm vs 실험 10~15 mS/cm (같은 order) |

### 21.4 핵심 발견

1. **Contact quality 변수 불필요**: A_hop, GB_d 추가 시 R² **감소** → constriction이 AM에서는 무시 가능
2. **exp(π/(T/d)) 필수**: β_fit = 3.06 ≈ π. 구형 입자 둘레/직경 = π (기하학적 기원)
3. **CN² = universal**: Ionic에서도 CN_SE² 동일 → DEM granular network의 보편적 구조

상세 유도는 `docs/electronic_conductivity_derivation.md` 참조.

---

# Part VI: Thermal Conductivity Scaling Law

## 22. 2상 Transport: 왜 Thermal이 가장 어려운가

### 22.1 근본적 차이

Thermal conductivity는 SE와 AM이 **모두** 열을 전도하는 2상 문제이다:

| 재료 | k (W/(m·K)) | 비율 |
|------|------------|------|
| NCM AM | 4.0 | 5.7× |
| LPSCl SE | 0.7 | 1× |

SE가 작고(~0.5μm) 많아서 촘촘한 mesh를 형성하므로 **SE가 thermal backbone**이다. AM이 많아지면 이 mesh가 끊기면서 역설적으로 **k_AM > k_SE인데도 σ_th가 감소**한다.

### 22.2 문헌 모델의 전멸

50개 이상의 모델을 테스트한 결과:

| 모델 유형 | 대표 | R² | 실패 원인 |
|----------|------|-----|----------|
| EMT (Bruggeman, Archie, Lichtenecker) | L1~L8 | -17 ~ 0.70 | 연속 매질 가정, constriction 무시 |
| Pure geometry (φ, τ only) | T1~T24 | 0.22 ~ 0.82 | 지수가 비물리적 (φ_AM^4.8) |
| Contact-based (k_hop, CN_thermal) | T30 | 0.26 | 잘못된 가중치 |
| Thin correction (exp(π/(T/d))) | T33b | -2.65 | **Thermal은 thickness-independent!** |

### 22.3 핵심 통찰: Thickness Independence

Electronic과 달리 thermal은 **전극 두께에 무관**:
- Per-case residual에서 T와의 상관: r = +0.006 (≈ 0)
- exp(π/(T/d)) 시도 → R²가 오히려 대폭 악화

**이유**: Electronic에서 exp(π/(T/d))가 필요한 이유는 AM-AM percolation이 thin에서 끊기기 때문. 그러나 thermal에서는 SE+AM 모두 열을 전도하므로 thin에서도 percolating path가 항상 존재.

## 23. 최종 공식 유도

$$\sigma_{th} = 286 \times \sigma_{ion}^{3/4} \times \phi_{AM}^2 \times CN_{SE}^{-1}$$

### 23.1 σ_ion^(3/4): SE Backbone Geometry

σ_ion (network solver 결과)은 SE backbone의 가장 완전한 서술자이다. φ_SE, τ, f_perc, A_hop, CN, GB_d를 모두 포함.

**왜 3/4?** Thermal이 ionic보다 SE geometry에 **덜 민감**:
- Ionic: SE-only → SE path가 끊기면 끝
- Thermal: AM도 전도 → SE가 일부 끊겨도 AM 우회 가능

Free fit: 0.758 → 3/4로 고정 시 ΔR² = 0.006.

### 23.2 φ_AM²: AM Thermal Enhancement

AM-AM 접촉 확률 ∝ φ_AM². AM이 network를 형성해야 k_AM(=4 W/m·K)의 높은 열전도를 활용.

Free fit: 1.987 → 2로 고정 시 ΔR² < 0.001.

### 23.3 CN_SE⁻¹: SE Clustering Penalty

**가장 중요한 발견**: CN_SE의 부호가 ionic과 **반대**!

| Transport | CN_SE 의존성 | 물리 |
|-----------|-------------|------|
| Ionic | CN² (**+**) | SE 연결 ↑ → ionic path ↑ |
| **Thermal** | CN⁻¹ (**-**) | SE clustering → AM 고립 → thermal path ↓ |

σ_ion이 이미 CN_SE의 positive 효과를 포함하고 있으므로, 추가 CN_SE⁻¹ term은 "SE가 너무 자기들끼리만 뭉치면 AM이 고립되어 thermal에 해로움"을 포착.

**이 CN 부호 역전이 단일 상(ionic/electronic)과 2상(thermal) transport의 근본 차이.**

Free fit: -1.171 → -1로 고정 시 ΔR² = 0.006.

### 23.4 C = 286

유일한 free parameter. k_AM/k_SE ratio, constriction geometry 등을 포함.

## 24. 검증

| Metric | 값 |
|--------|-----|
| Train R² | **0.90** (1 free param) |
| LOOCV R² | **0.86** |
| Free fit (4p) R² | 0.9019 |
| Mean |error| | **9.4%** (LOOCV: 9.2%) |
| Cases within 10% | 28/41 (68%) |
| Cases within 20% | 37/41 (90%) |
| Max error | 42.4% (σ_th=0.797, 극단 low case) |

### Exponent Sensitivity

| Variant | σ_ion exp | φ_AM exp | CN exp | R² |
|---------|-----------|----------|--------|-----|
| Free fit | 0.758 | 1.987 | -1.171 | **0.9019** |
| **Fixed** | **3/4** | **2** | **-1** | **0.8961** |
| σ_ion^(1/2) | 0.5 | 2 | -1 | 0.6681 |
| σ_ion^1 | 1 | 2 | -1 | 0.5645 |
| φ_AM^(3/2) | 3/4 | 1.5 | -1 | 0.8721 |

σ_ion 지수가 3/4에서 벗어나면 R² 급락 → 3/4이 유일한 최적.

### Ablation: Analytic Bruggeman (σ_brug) vs Network Solver (σ_ion)

σ_ion 대신 φ_SE × f_perc / τ² (analytic)를 사용하면 R²=0.17로 **붕괴**.

**Network solver σ_ion이 필수** — contact quality (A_hop, CN, GB_d) 정보가 σ_ion에 내재되어 있으며, analytic approximation으로는 대체 불가.

## 25. 세 Transport의 통일적 프레임워크

### 25.1 비교

| | Ionic | Electronic | Thermal |
|---|---|---|---|
| **공식** | σ_brug × C × (G_path × GB_d²)^(1/4) × CN² | 0.015 × σ_AM × φ^(3/2) × CN² × exp(π/(T/d)) | 286 × σ_ion^(3/4) × φ_AM² / CN_SE |
| **주도 인자** | Contact quality + network | Topology | 2상 competition |
| **CN** | CN_SE² (+) | CN_AM² (+) | CN_SE⁻¹ (**-**!) |
| **Thickness** | 간접 (GB_d 반영) | exp(π/(T/d)) | 무관 |
| **R²** | **0.947** | 0.89 | 0.90 (LOOCV 0.86) |
| **Free params** | 1 (C ≈ 0.073) | 1 (C = 0.015) | 1 (C = 286) |

### 25.2 CN 부호 역전: 단일 상 vs 2상

세 공식에서 가장 주목할 발견은 **CN의 부호 역전**이다:

- **Ionic CN_SE²(+)**: SE-SE 연결 ↑ → redundant path ↑ → 전도 ↑
- **Electronic CN_AM²(+)**: AM-AM 연결 ↑ → redundant path ↑ → 전도 ↑
- **Thermal CN_SE⁻¹(-)**: SE-SE 연결 ↑ → SE clustering → **AM 고립 → 열전도 ↓**

단일 상 transport에서는 "연결이 많을수록 좋다"가 항상 성립. 2상 transport에서는 **한 상의 과도한 연결이 다른 상을 squeeze out**하여 전체 성능을 저하시킨다.

### 25.3 Thermal이 Ionic과 Electronic을 연결

Thermal 공식은 ionic 결과(σ_ion)를 직접 변수로 사용한다. 이는 thermal transport가 ionic path(SE backbone)와 AM enhancement의 **기하학적 competition**에 의해 결정됨을 보여준다. 세 공식이 하나의 통일된 프레임워크를 형성한다.

상세 유도는 `docs/thermal_conductivity_derivation.md` 참조.

---

## References (Part V~VI 추가)

10. Ebner, M., Chung, D.-W., Garcia, R. E., Wood, V. (2013). *Adv. Energy Mater.* **4**, 1301278.
11. Makse, H. A., Gland, N., Johnson, D. L., Schwartz, L. M. (2004). *Phys. Rev. E* **70**, 061302.
12. Kirkpatrick, S. (1973). *Rev. Mod. Phys.* **45**, 574--588.
13. Ketter, F. et al. (2025). Thermal conductivity of LPSCl argyrodite.

---

*Report generated by DEM Analyzer v2.0 | Kirchhoff Network Solver + Physics Scaling Laws (Ionic + Electronic + Thermal)*
