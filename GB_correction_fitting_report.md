# GB Correction Fitting Analysis Report

## 1. 목적

Bruggeman 추정치(σ_brug)와 실측 proxy(σ_proxy)의 비율 **R = σ_brug / σ_proxy**가
GB density(GB_d)와 전극 두께(T)에 어떻게 의존하는지를 다양한 회귀 모델로 비교하여 최적 모델을 선정한다.

### 정의

| 변수 | 정의 | 단위 |
|------|------|------|
| R | σ_brug / σ_proxy (과대평가 비율) | - |
| σ_brug | σ_bulk × φ_SE × f_perc / τ² (Bruggeman, GB 무시) | ratio |
| σ_proxy | G_path × f_perc / τ (접촉면적 기반 실측) | ratio |
| GB_d | 입계 밀도 (percolation path의 hop 수 / z 거리) | hops/μm |
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

**M15 추천 이유:** R²=0.9427로 M6(0.9490)과 거의 동일하면서 파라미터가 2개 (vs 3개). 물리적 유도가 가능한 유일한 모델.

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
| σ_bulk | 1.3 mS/cm | SE 물성 (Li₆PS₅Cl, cold-pressed) | 벌크 전도도 |
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

---

## 8. M6 vs M15: 독립적 검증

M6 (자유 3-parameter fit)의 결과가 M15의 제약된 형태와 거의 일치하며, 이는 GB_d²×T가 natural scaling variable임을 독립적으로 확인한다.

| | GB_d 지수 | T 지수 | 파라미터 수 | R² |
|---|----------|--------|----------|-----|
| M6 (free) | 4.02 | 1.75 | 3 | 0.949 |
| M15 (constrained) | 3.64 (=1.82×2) | 1.82 | 2 | 0.943 |

M6이 자유도 3개로 찾은 최적값(GB_d^4.02 × T^1.75)이 M15의 결합 변수(GB_d²×T)^1.82 = GB_d^3.64 × T^1.82와 거의 동일하다. **데이터가 자발적으로 GB_d²×T 결합을 선택한다는 증거.**

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

**SE 0.5→1.5μm 변경 시 inter-particle contact resistance가 ~30배 감소.** 원인: GB_d² 감소에 의한 접촉 저항 저감.

---

## 10. Limitations

1. **절대 전도도 예측 불가:** σ_proxy는 DEM 접촉면적 기반 proxy이며 실험 EIS 측정값이 아님. R = σ_brug/σ_proxy는 두 DEM 계산값의 비율로, 스케일링 관계를 제시하지만 절대 전도도 예측에는 EIS calibration이 필요함.

2. **Inter-particle contact vs Crystallographic GB:** 본 모델의 "GB"는 복합양극 내 SE 입자 간 접촉(inter-particle contact)이며, 소결 pellet 내 결정립계(crystallographic grain boundary)와는 다른 물리적 실체임. BLM의 입계 수 스케일링(N_GB ∝ GB_d × T)만을 차용한 BLM-inspired scaling이며, 소결 세라믹의 입계 고유 저항(ρ_gb × w_gb)은 사용하지 않음.

3. **Argyrodite의 inter-particle resistance:** Cold-pressed 황화물은 전통적으로 negligible GB resistance로 간주되었으나, 최근 AIMD (Sadowski, 2024)와 NMR (de Klerk, 2019)에서 argyrodite pellet 내에서도 유의미한 GB 효과가 보고됨. Sintered pellet의 crystallographic GB와 composite cathode의 mechanical inter-particle contact는 서로 다른 메커니즘이지만, 둘 다 입자 간 이온 수송의 병목이라는 점에서 공통된다. 복합양극의 less intimate한 접촉에서는 추가적인 constriction resistance가 예상됨.

4. **데이터 분포:** 41개 중 SE 0.5μm이 30개로 다수. SE 1.0μm(3개), 1.5μm(8개)으로 GB_d 0.7~1.2 중간 영역의 데이터가 부족. R²=0.94가 SE 0.5μm 데이터에 의해 driven될 가능성을 배제할 수 없음.

5. **α=1.82의 기여 분리 불가:** α > 1의 원인(bottleneck, 소성변형, 네트워크 구조)을 개별적으로 정량 분리하기 위해서는 rigid sphere DEM 비교 등 추가 연구가 필요.

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

*Report generated by DEM Analyzer | BLM-inspired Inter-particle Contact Resistance Scaling Model*
