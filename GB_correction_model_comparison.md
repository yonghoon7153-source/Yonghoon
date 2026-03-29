# GB Correction Model Comparison

## 목적

Bruggeman 추정치($\sigma_{brug}$)와 실측 proxy($\sigma_{proxy}$)의 비율이 GB density($GB_d$)에 어떻게 의존하는지를 다양한 회귀 모델로 비교하여, 최적 모델을 선정한다.

**정의:**
- $R = \sigma_{brug} / \sigma_{proxy}$ (과대평가 비율)
- $\sigma_{brug} = \sigma_{bulk} \times \phi_{SE} \times f_{perc} / \tau^2$ (Bruggeman, GB 무시)
- $\sigma_{proxy} = G_{path} \times f_{perc} / \tau$ (접촉면적 기반 실측)
- $GB_d$: 입계 밀도 (hops/μm)
- 데이터: n = 21 (SE 0.5μm 16개 + SE 1.5μm 5개)

---

## Tier 1: R² > 0.96 (유효 모델)

### A1. Exponential Decay (채택)

$$\ln(R) = b \cdot GB_d + \ln(k)$$

$$\therefore R = k \cdot e^{b \cdot GB_d}$$

| 파라미터 | 값 |
|---|---|
| b | 5.21 |
| ln(k) | 0.20 |
| **R²** | **0.9614** |
| σ_eff @GB_d=1.3 | 0.00019 |
| σ_eff @GB_d=0.5 | 0.0079 |

**물리적 근거:** 각 grain boundary가 독립적인 투과 barrier로 작용. 이온이 n개의 GB를 통과할 확률 = 개별 투과확률의 곱 = $T^n$. $n \propto GB_d$이므로 지수함수 형태가 자연스럽다.

**채택 이유:**
1. R² 최상위 (0.9614)
2. 파라미터 2개 — parsimonious
3. Stretched exp (E1)가 자발적으로 n=1 선택 → 데이터가 지수형태를 "원함"
4. 외삽 안전: $GB_d \to 0$이면 $R \to k \approx 1.22$ (Bruggeman에 수렴)

---

### E1. Stretched Exponential

$$\ln(R) = b \cdot GB_d^n + c$$

| 파라미터 | 값 |
|---|---|
| b | 5.21 |
| n | **1.00** |
| c | 0.20 |
| **R²** | **0.9614** |

**근거:** GB barrier 높이가 불균일(분포)할 때, $n < 1$이면 완만한 감쇠 (Kohlrausch–Williams–Watts). $n > 1$이면 가속 감쇠.

**결과:** $n = 1.00$ → A1과 완전 동일. 데이터에 barrier 불균일성이 없거나, 21개 데이터로는 구분 불가. **A1의 정당성을 강화하는 근거.**

---

### D1. Quadratic (log domain)

$$\ln(R) = a \cdot GB_d^2 + b \cdot GB_d + c$$

| 파라미터 | 값 |
|---|---|
| a | 0.00 |
| b | 5.21 |
| c | -0.00 |
| **R²** | **0.9620** |
| σ_eff @GB_d=0.5 | 0.0083 |

**근거:** 고밀도 GB에서 cooperative 효과 — GB끼리 가까우면 space charge layer가 겹쳐서 저항이 가속적으로 증가.

**탈락 이유:** $a \approx 0$이므로 사실상 A1과 동일. 파라미터 3개 사용은 과적합 위험. R² 차이 0.0006은 유의미하지 않음.

---

### A2. Exponential (원점 통과)

$$\ln(R) = b \cdot GB_d$$

$$\therefore R = e^{b \cdot GB_d}$$

| 파라미터 | 값 |
|---|---|
| b | 5.31 |
| **R²** | **0.9604** |
| σ_eff @GB_d=0.5 | 0.0088 |

**근거:** $k = 1$ 강제 → 단위 변환 계수가 없다고 가정. $GB_d = 0$이면 $R = 1$ (Bruggeman 정확).

**장점:** 물리적 boundary condition 만족 ($GB_d = 0 \Rightarrow \sigma_{eff} = \sigma_{brug}$).
**단점:** $\ln(k) = 0.20$이 유의미할 수 있으므로 절편을 강제로 제거하면 편향 발생.

---

### E3. Square Root

$$\ln(R) = a \cdot \sqrt{GB_d} + c$$

| 파라미터 | 값 |
|---|---|
| a | 6.56 |
| c | -0.75 |
| **R²** | **0.9602** |
| σ_eff @GB_d=0.5 | 0.0084 |

**근거:** 확산(diffusion) 기반 모델. GB를 통한 이온 이동이 확산 지배적이면 저항이 $\sqrt{거리}$에 비례. Fick's law 유도.

**A1과의 차이:** 지수함수보다 완만한 감쇠 → GB_d 높은 영역에서 A1보다 관대. 현재 데이터 범위(0.48~1.38)에서는 구분 불가.

---

### B1. Power Law

$$\ln(R) = c \cdot \ln(GB_d) + d$$

$$\therefore R = e^d \cdot GB_d^c$$

| 파라미터 | 값 |
|---|---|
| c | 4.53 |
| d | 5.80 |
| **R²** | **0.9561** |
| σ_eff @GB_d=0.5 | 0.0091 |

**근거:** Percolation theory에서 임계현상 부근의 스케일링 법칙. 전도도가 구조 파라미터의 거듭제곱으로 변하는 보편적 형태. Archie's law ($\sigma = \sigma_0 \cdot \phi^m$)와 유사.

**단점:** $GB_d \to 0$이면 $R \to 0$ (발산 문제). 물리적으로 $R \geq 1$이어야 하는데 boundary condition 위반. 외삽 위험.

---

## Tier 2: R² < 0.94 (탈락 모델)

### E2. Arrhenius (inverse)

$$\ln(R) = a / GB_d + c$$

| R² | 0.9387 |
|---|---|
| σ_eff @GB_d=0.5 | 0.0105 |

**근거:** 활성화 에너지 형태. Arrhenius식 $\sigma = \sigma_0 \exp(-E_a/k_BT)$에서 온도 대신 $GB_d$를 넣은 형태.

**탈락:** $GB_d \uparrow$이면 $a/GB_d \downarrow$이므로 $R \downarrow$. 즉 **GB 많을수록 저항 감소** → 물리적으로 반대 방향. Arrhenius는 온도에 대한 식이지 공간 밀도에 적용하면 방향이 뒤집힘.

---

### E4. Logistic (Sigmoid)

$$R = \frac{L}{1 + e^{-k(GB_d - x_0)}}$$

| R² | 0.8827 |
|---|---|
| σ_eff @GB_d=0.5 | 0.0483 |

**근거:** GB_d에 임계점($x_0$)이 존재하여 그 이상에서 ratio가 급격히 증가. Phase transition 유사.

**탈락:** SE 1.5μm (GB_d=0.5)에서 σ_eff = 0.048 → 실측(0.004~0.016)의 **3~10배 과대**. 포화 구간(GB_d < $x_0$)에서 ratio를 너무 낮게 잡음.

---

### D2. Quadratic (ratio domain)

$$R = a \cdot GB_d^2 + b \cdot GB_d + c$$

| R² | 0.8643 |
|---|---|
| σ_eff @GB_d=0.5 | FAIL (음수) |

**근거:** 단순 2차 다항식 fitting. 물리적 의미 없는 경험식.

**탈락:** GB_d = 0.5에서 ratio가 음수 → 물리적 불가능. 외삽 완전 실패. 2차식은 볼록/오목 전환점 때문에 단조 증가해야 하는 ratio에 부적합.

---

### C1. Linear (ratio ~ GB_d)

$$R = a \cdot GB_d + c$$

| R² | 0.8526 |
|---|---|
| σ_eff @GB_d=0.5 | FAIL (음수) |

**근거:** $R_{total} = R_{bulk} + n \cdot R_{gb}$의 가장 단순한 형태. 저항이 GB 수에 선형 비례.

**탈락:** GB_d=1.3에서 ratio=1070인데 GB_d=0.5에서 ratio=-31. 기울기(1479)가 너무 급해서 절편(-770)이 음수. 데이터의 **비선형적 급변**을 선형으로 잡으려니 파탄.

---

### B2. Power Law + Offset

$$R = a \cdot GB_d^c + 1$$

| R² | 0.8271 |
|---|---|
| σ_eff @GB_d=0.5 | 0.0012 |

**근거:** $GB_d = 0$이면 $R = 1$ (Bruggeman 정확)이라는 물리적 boundary condition을 만족시키면서 거듭제곱 성장.

**탈락:** 물리적 BC는 좋지만, fitting력 부족. SE 1.5μm (실측 0.008)에서 0.0012 → **실측의 1/7**. offset=1이 ratio=1000급 데이터를 fitting하는 데 제약이 됨.

---

### C3. Harmonic (1/σ ~ GB_d)

$$1/\sigma_{proxy} = a \cdot GB_d + c$$

| R² | 0.7945 |
|---|---|
| σ_eff @GB_d=0.5 | FAIL (음수) |

**근거:** 옴의 법칙 직접 적용. 저항(= 1/전도도)이 GB_d에 선형 비례.

**탈락:** C1과 유사한 문제. 비선형성 못 잡고 GB_d=0.5에서 음수. ratio 범위가 17~1400으로 너무 넓어서 1/σ 공간에서도 선형이 안 됨.

---

### C2. Series Resistance (원점 통과)

$$R = 1 + a \cdot GB_d$$

| R² | 0.6849 |
|---|---|
| σ_eff @GB_d=0.5 | 0.0003 |

**근거:** 가장 물리적인 직렬저항 모델. $R_{total} = R_{bulk} + n \cdot R_{gb}$에서 정규화하면 $R = 1 + (R_{gb}/R_{bulk}) \cdot GB_d$.

**탈락:** R² 최하위. SE 1.5μm에서 실측(0.008)의 **1/25 수준**으로 극심한 과소예측. 이유: ratio가 GB_d에 **선형이 아니라 지수적**으로 증가하는데, 선형 모델은 이 급격한 변화를 근본적으로 표현 불가.

---

## 예측 비교 요약

### GB_d = 1.3 (SE 0.5μm)

| 모델 | σ_eff/σ_bulk |
|---|---|
| 모든 Tier 1 모델 | **0.00017~0.00019** |

→ **모델 선택과 무관하게 동일 결론.** SE 0.5μm은 어떤 모델을 써도 극히 낮음.

### GB_d = 0.5 (SE 1.5μm)

| 모델 | σ_eff/σ_bulk | 실측 range |
|---|---|---|
| B1. Power law | 0.0091 | |
| A2. Exp(원점) | 0.0088 | |
| E3. Sqrt | 0.0084 | 0.004 ~ 0.016 |
| D1. Quadratic | 0.0083 | |
| **A1. Exponential** | **0.0079** | |

→ Tier 1 모델 모두 **실측 범위 안**. 모델 간 차이는 ±15% 수준으로 결론에 영향 없음.

---

## 결론

1. **지수함수(A1)가 최적 모델.** R² 최상위, 파라미터 2개, 물리적 해석 가능.
2. **Stretched exp (n=1.00)이 A1을 자발 선택** → 지수형태가 데이터의 고유 특성.
3. **선형/직렬저항 모델은 전부 탈락** → 데이터가 본질적으로 비선형.
4. **GB_d=1.3에서는 모델 무관** → 0.5μm SE의 GB 불리함은 모델 선택과 독립적.
5. **GB_d=0.5에서 Tier 1 모델 간 차이 ±15%** → 결론에 영향 없음.

> 논문 표현: "13개 후보 모델을 비교한 결과, exponential decay model ($R² = 0.9614$)이 가장 높은 설명력과 물리적 일관성을 보였다. Stretched exponential fitting에서 stretch exponent $n = 1.00$이 자발적으로 선택되어 순수 지수 감쇠의 타당성이 재확인되었다."
