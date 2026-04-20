# v29 FINAL Formula — 각 항의 물리적 의미와 함수 형태

> 구현 single source of truth: `scripts/generate_comparison_plots.py`의
> `_formx_v29_predict()` (line 2740) 및 `_fit_at()` (line 1112).
> 이 문서는 코드 재검증 완료(2026-04-20).

---

## 전체 공식

$$
\sigma_{\mathrm{ionic}} = e^{C_{\mathrm{corr}}}\,\cdot\,C_{\mathrm{blend}}(\tau)\,\cdot\,\sigma_{\mathrm{grain}}\,\cdot\,(\phi_{SE}-\phi_c)^{1/2}\,\cdot\,\mathrm{CN}^{3/2}\,\cdot\,\mathrm{cov}^{2/5}\,\cdot\,f_{\mathrm{perc}}^{3}
$$

여기서 $C_{\mathrm{corr}}$는 3-term residual correction (P:S sigmoid + τ-bump × p + gb sigmoid).

**고정 상수** (v12-clean v3, `_formx_v29_params`):

| 상수 | 값 | 비고 |
|---|---|---|
| $\sigma_{\mathrm{grain}}$ | 3.0 mS/cm | Li₆PS₅Cl grain interior |
| $\phi_c$ | 0.20 | Percolation threshold (data-native) |
| $\alpha$ (φ 지수) | 1/2 | Mean-field / quasi-2D percolation |
| $\beta$ (CN 지수) | 3/2 | Kirkpatrick |
| $\gamma$ (cov 지수) | 2/5 | AM-SE interface |
| $\delta$ ($f_{\mathrm{perc}}$ 지수) | 3 | v12 free fit 2.90 → 3 반올림 |

---

## 1. $f_{\mathrm{perc}}^{3}$ — Percolation fraction cube

### 정의

$$
f_{\mathrm{perc}} = \frac{\text{percolating cluster에 속한 SE 개수}}{\text{전체 SE 개수}}
$$

코드: `f_perc = percolation_pct / 100`, `max(f_perc, 0.01)` floor 적용.

### 물리 의미

**percolation 네트워크의 fragmentation penalty**.

- $f_{\mathrm{perc}} = 1$: SE가 완전히 연결된 하나의 cluster — penalty 없음
- $f_{\mathrm{perc}} = 0.5$: 절반만 percolating — $(0.5)^3 = 0.125$로 세제곱 급감
- $f_{\mathrm{perc}} \to 0$: SE fragmented, 이온 경로 차단 — σ → 0 ✓

### monomodal에서 문제 없는 이유

$f_{\mathrm{perc}}$는 **fine fraction이 아니라 percolation fraction**. monomodal AM + SE 시스템에서도 SE가 percolate하면 $f_{\mathrm{perc}} \approx 1$. monomodal ↔ bimodal 구분은 이 항이 아니라 residual correction의 $p_{\mathrm{frac}}$(P:S ratio) 항이 담당.

### 왜 3제곱?

**경험적**. v12 데이터-네이티브 free fit에서 δ=2.90 추정, 단순 분수로 반올림하여 3. "v12 자유도 7개 재fit 시 3이 1/2·1/2의 noise 범위 내 최적" (`ionic_scaling_law_experiments.md:135`).

물리적 해석 시도:
- 3D cluster size distribution: finite cluster의 평균 크기 $\sim f_{\mathrm{perc}}^{\nu d}$ with $\nu d \approx 3$ (3D lattice percolation)
- 연속 Monte-Carlo가 이 지수를 $2.5{-}3.0$으로 확인하긴 함 — 확정적 derivation은 아님

---

## 2. $C_{\mathrm{blend}}(\tau)$ — Tortuosity blend (regime-aware)

### 구조

**두 개의 서로 다른 τ-함수를 smooth blend**:

$$
\ln C_{\mathrm{blend}}(\tau) = (1-w_{\mathrm{bl}}(\tau))\cdot \ln C_{v5}(\tau) + w_{\mathrm{bl}}(\tau)\cdot \ln C_{p3}(\tau)
$$

#### $C_{v5}(\tau)$ — thick↔thin asymptote (내부 sigmoid)

$$
\ln C_{v5}(\tau) = \ln C_t + (\ln C_n - \ln C_t)\cdot\sigma_{k_\tau}(\tau - \tau_{0})
$$

- $\tau_0 = 2.1$, $k_\tau = 5.0$ (고정)
- $C_t$, $C_n$: thick/thin asymptote, LOOCV fit (≈ 0.035 / 0.020)

#### $C_{p3}(\tau)$ — 극단 thin의 poly3 in $\ln\tau$

$$
\ln C_{p3}(\tau) = a_0 + a_1\ln\tau + a_2(\ln\tau)^2 + a_3(\ln\tau)^3
$$

계수 (기본): $(-3.80,\ +2.38,\ -5.58,\ +2.81)$.

#### Blending weight

$$
w_{\mathrm{bl}}(\tau) = \sigma_{k_{\mathrm{bl}}}(\tau - \tau_{c,\mathrm{bl}}) = \frac{1}{1+e^{-k_{\mathrm{bl}}(\tau - \tau_{c,\mathrm{bl}})}}
$$

- $k_{\mathrm{bl}}, \tau_{c,\mathrm{bl}}$: Nelder-Mead로 LOOCV 최대화
- 실측 최적 ≈ $k_{\mathrm{bl}}\sim 5$, $\tau_{c,\mathrm{bl}}\sim 2.0$

### 물리 의미

- **moderate τ** ($\tau < \tau_{c,\mathrm{bl}}$): $C_{v5}$가 thick↔thin 전이 담당
- **extreme thin** ($\tau > \tau_{c,\mathrm{bl}}$): poly3가 $\ln\tau$에 대한 큰 곡률 커버
- 두 함수 어느 쪽도 전 영역에서 단독으로 fit 불가 → 경계 sigmoid로 부드럽게 이음

---

## 3. Residual correction 3항

최종 $C_{\mathrm{corr}}$는 아래 3항의 **합** (linear regression on residual):

$$
C_{\mathrm{corr}} = \beta_{pf}(w_{pf}-\bar{w}_{pf}) + \beta_{\mathrm{lin}}(p\cdot w_{\mathrm{win}} - \overline{p\cdot w_{\mathrm{win}}}) + \beta_{gb}(w_{gb}-\bar{w}_{gb})
$$

각 항은 "평균값 빼기"로 centering — 기존 prefactor에 영향 없이 순수 residual만 fit.

### 3.1 $\beta_{pf}\cdot w_{pf}$ — P:S ratio sigmoid

$$
w_{pf}(p) = \sigma_{k_{pf}}(p - p_c) = \frac{1}{1+e^{-k_{pf}(p-p_c)}}
$$

- $p = p_{\mathrm{frac}} = P/(P+S)$, P:S 질량비에서 유도
- 기본값: $k_{pf}=50$, $p_c = 0.598$, $\beta_{pf} \approx -0.10$
- 해석: P가 많으면 ($p > p_c$) 모델이 과대예측 → $\exp(\beta_{pf}) = 0.905$만큼 깎음
- "P 입자 많음 = particulate 접촉 empty contact 많음" 보정

### 3.2 $\beta_{\mathrm{lin}}\cdot p\cdot w_{\mathrm{win}}$ — Mixed-regime τ-bump

$$
w_{\mathrm{win}}(\tau) = \exp\left[-\frac{1}{2}\left(\frac{\tau-\tau_{c,\mathrm{win}}}{\sigma_{\tau,\mathrm{win}}}\right)^2\right]
$$

- Gaussian bump, center $\tau_{c,\mathrm{win}} \approx 2.0$, width $\sigma_{\tau,\mathrm{win}} \approx 0.15$
- $p_{\mathrm{frac}}$와 곱해서 "τ≈2 근처 + P-heavy"일 때만 활성
- $\beta_{\mathrm{lin}} \approx -0.49$
- 해석: thin regime 전이 구간에서 particulate 케이스가 특히 과대예측

### 3.3 $\beta_{gb}\cdot w_{gb}$ — gb_density sigmoid (v29 신규)

$$
w_{gb}(\rho_{gb}) = \sigma_{k_{gb}}(\ln\rho_{gb} - \overline{\ln\rho_{gb}})
$$

- log-space sigmoid, center = $\ln\rho_{gb}$의 **median** (데이터 의존)
- $k_{gb} = 4.0$ (고정)
- $\beta_{gb} \approx +0.043$ (작은 값, bounded)
- 해석: GB density가 높으면 constriction resistance 증가 — Bruggeman이 놓치는 $R_{\mathrm{constr}}$ 보정

---

## 요약 표

| 항 | 정의 | 형태 | 파라미터 |
|---|---|---|---|
| $(\phi_{SE}-\phi_c)^{1/2}$ | SE 부피분율 excess | 멱함수 | $\phi_c=0.20$ 고정, 지수 1/2 |
| $\mathrm{CN}^{3/2}$ | SE-SE 접촉 수 | 멱함수 | 지수 3/2 (Kirkpatrick) |
| $\mathrm{cov}^{2/5}$ | AM-SE 표면 coverage | 멱함수 | 지수 2/5 (0.25→0.40 enhanced) |
| $f_{\mathrm{perc}}^{3}$ | SE percolation fraction | 멱함수 | 지수 3 (v12 free fit 2.90 반올림) |
| $C_{\mathrm{blend}}(\tau)$ | v5 sigmoid ⊕ poly3 blend | 이중 sigmoid blend | $k_{\mathrm{bl}}, \tau_{c,\mathrm{bl}}$ + v5 asymptote 2개 + poly3 계수 4개 |
| $\beta_{pf}\cdot w_{pf}$ | P:S ratio 보정 | sigmoid | $k_{pf}, p_c, \beta_{pf}$ |
| $\beta_{\mathrm{lin}}\cdot p\cdot w_{\mathrm{win}}$ | thin 전이 + P-heavy | Gaussian × linear | $\tau_{c,\mathrm{win}}, \sigma_{\tau,\mathrm{win}}, \beta_{\mathrm{lin}}$ |
| $\beta_{gb}\cdot w_{gb}$ | GB density 보정 | log-sigmoid | $k_{gb}=4, \beta_{gb}$ |

**총 free parameter** (LOOCV 최적화 대상): $k_{\mathrm{bl}}, \tau_{c,\mathrm{bl}}, k_{pf}, p_c, \tau_{c,\mathrm{win}}, \sigma_{\tau,\mathrm{win}}$ 6개 + residual $\beta$ 3개 + v5 asymptote ($C_t, C_n$) 2개 + poly3 계수 4개. 나머지($\phi_c$, 4개 멱지수, $k_{gb}$)는 **고정**.

**성능** (n=57): R² = 0.9813, LOOCV = 0.9791, ±20% = 51/57, median |err| = 6.2%.

---

## 물리적 해석 요점

1. **Base Kirkpatrick scaling** (네 멱함수 곱): $(\phi-\phi_c)^{1/2}\,\mathrm{CN}^{3/2}\,\mathrm{cov}^{2/5}\,f_{\mathrm{perc}}^3$
   — percolation theory + 접촉 역학 + SE 네트워크 연결성 반영
2. **$C_{\mathrm{blend}}(\tau)$**: 단일 $\tau$-멱함수로는 moderate/thin regime 동시 커버 불가 → 두 함수 blend. τ² (Bruggeman)는 thin 과벌점, τ^(1/2)는 thick 과소벌점 → regime-aware 필수
3. **3-term residual**: base formula의 잔차에서 **통계적으로 유의한 신호 3개만** 남김 (v30 phase-split은 noise σ=0.0018 아래로 기각)

---

## 확인 방법

```python
# 구현 확인
grep -n "_formx_v29_predict\|_fit_at\|PHI_C\|K_PF\|B_PF\|B_LIN\|B_GB" \
  scripts/generate_comparison_plots.py

# 파라미터 기본값
grep -A 5 "_formx_v29_params" scripts/generate_comparison_plots.py | head -30
```

---

## 변경 이력

- **2026-04-20**: 전면 개정. 원 문서의 오류 수정:
  - `f_p^3` 라벨을 "fine-particle volume fraction"(틀림) → "percolation fraction"(맞음)으로 정정
  - Furnas-McGeary 해석 제거 (이 항의 물리 아님 — `docs/Packing_Regime_Analysis.md`의 porosity 해석과 혼동)
  - $C_{\mathrm{blend}}(\tau)$를 단일 sigmoid 보정(틀림) → v5 sigmoid ⊕ poly3 blend(맞음)로 재기술
  - $G(\tau,p)$ 라벨을 "bilinear sigmoid coupling"(틀림) → "Gaussian τ-bump × linear $p$"(맞음)로 정정
  - $C_{gb}$ 항 추가 (v29의 신규 항, 원 문서 누락)
  - $\phi_c$ 값을 0.185 → 0.20으로 수정 (v9 → v12-clean v3 production)
  - 기본 상수값을 코드와 일치시킴 ($\beta_{pf}=-0.10$, $p_c=0.598$, $k_{pf}=50$ 등)
