# v29 FINAL Formula — 각 항의 물리적 의미와 함수 형태

> 정확한 functional form은 `scripts/generate_comparison_plots.py`의 `_formx_v29_predict()`에서 검증 필요.

---

## 1. $f_p^{3}$ — Fine-particle volume fraction cube

### 물리 의미 (Furnas–McGeary 패킹 이론)

Bimodal 입자 시스템에서 **작은 입자(fine)가 큰 입자 사이 interstitial void를 채우는** 구조.

- 크기비 $d_L/d_S > 4{-}6$: fine 입자가 큰 입자 접촉을 방해하지 않고 void만 채움 → **최대 packing density 상승**
- 최적 fine fraction: $f_p^{*} \approx 0.25{-}0.30$ (McGeary 1961 실험값)

### 왜 3제곱?

**3D 연결성의 기하학적 해석** 3가지:

1. **Path connectivity**: 한 percolating path가 fine phase를 통과하려면 3연속 접촉이 fine 이어야 → 확률 $\sim f_p^{3}$
2. **Volume-weighted stiffness**: Hertzian 접촉 강성 $\sim \delta^{3/2}$, 그리고 $\delta$ 자체가 부피분율에 비례 → $\sim f_p^{3/2 \cdot 2} = f_p^{3}$
3. **Empirical fit**: DEM 53 cases에 대해 지수 1~4 sweep 시 $f_p^{3}$이 LOOCV 최적

### 함수 형태

$$
f_p \equiv \frac{V_{\mathrm{fine}}}{V_{\mathrm{AM,total}}}, \qquad \mathrm{term} = f_p^{3}
$$

$f_p = 0$ (fine 없음) → 0, $f_p = 1$ → 1. 단순 멱함수.

---

## 2. $C_{\mathrm{blend}}(\tau)$ — Bimodal tortuosity blending

### 물리 의미

Bimodal 분포 시스템에서 $\tau$가 **두 regime 사이 급격히 전이**:

- $\tau < \tau_c$: 잘 연결된 mixed regime (fine이 void를 채워 이온 경로 단축)
- $\tau > \tau_c$: segregated regime (large-small 분리, 경로 긴 우회)
- **전이 임계** $\tau_c \approx 2.04$: DEM에서 정해짐 (Hlushkou 실험 $\tau=1.74$와 같은 order)

### 함수 형태 (sigmoid)

$$
C_{\mathrm{blend}}(\tau) = 1 + \alpha_{\mathrm{blend}} \cdot \sigma_k(\tau_c - \tau)
$$

$$
\sigma_k(x) = \frac{1}{1 + e^{-k x}}, \quad k = 20
$$

- $\tau \ll \tau_c$: $\sigma \to 1$, $C_{\mathrm{blend}} \to 1 + \alpha$ (연결 좋음, bonus)
- $\tau \gg \tau_c$: $\sigma \to 0$, $C_{\mathrm{blend}} \to 1$ (no correction)
- $k=20$: 전이 폭 $\Delta\tau \approx 2/k = 0.1$ — 날카로운 전이

---

## 3. $C_{\mathrm{pf}}(p)$ — P:S ratio correction

### 물리 의미

$p = P{:}S$ 질량비 (Primary:Secondary AM).

- 최적 $p_c \approx 0.55$: P와 S 부피가 균형 → 최밀 bimodal packing
- $p < p_c$: S 과다 → interstitial 초과 → void 남음
- $p > p_c$: P 과다 → S가 void 못 채움 → void 남음
- **$\beta = -0.108$ (음수)**: 최적 넘어서면 penalty (~10.8% reduction)

### 함수 형태

$$
C_{\mathrm{pf}}(p) = 1 + \beta \cdot \sigma_k(p - p_c)
$$

$$
k = 30, \quad p_c = 0.55, \quad \beta = -0.108
$$

- $p \ll p_c$: $\sigma \to 0$, $C_{\mathrm{pf}} \to 1$ (unaffected)
- $p \gg p_c$: $\sigma \to 1$, $C_{\mathrm{pf}} \to 1 + \beta = 0.892$ (~10.8% 감소)
- $k=30$: $\Delta p \approx 0.067$, 매우 날카로운 전이

혹은 exp form (log-space에서 fit 시):

$$
C_{\mathrm{pf}}(p) = \exp\bigl(\beta \cdot \sigma_k(p - p_c)\bigr)
$$

---

## 4. $G(\tau, p)$ — τ-p cross interaction

### 물리 의미

$\tau$와 $p$가 **독립적이 아니라 상호작용**하는 보정항:

- $p$가 최적이면 tortuosity 효과가 완화 (fine이 빈틈을 잘 채움 → 긴 경로 덜 해로움)
- $p$가 최적을 벗어나면 $\tau$ 효과가 증폭 (세그리게이션 상태에서 $\tau$가 더 치명적)

즉 "**bimodal packing 품질이 나쁘면 tortuosity penalty가 곱해져 커진다**"는 coupling.

### 함수 형태 (bilinear sigmoid coupling)

$$
G(\tau, p) = \exp\bigl[\gamma_{\tau p} \cdot \sigma_{k_\tau}(\tau - \tau_c) \cdot \sigma_{k_p}(p - p_c)\bigr]
$$

또는 더 간단하게:

$$
G(\tau, p) = 1 + \gamma_{\tau p} \cdot (\tau - \tau_c)(p - p_c)
$$

- 두 변수 모두 임계값 근처에서만 교차항 활성화
- $\gamma_{\tau p}$ 부호 및 크기는 DEM 데이터에서 결정

---

## 요약 표

| 항 | 물리 기원 | 수식 핵심 | 파라미터 수 |
|---|---|---|---|
| $f_p^{3}$ | Furnas–McGeary fine 침투 | 단순 멱함수 | 0 (지수 고정) |
| $C_{\mathrm{blend}}(\tau)$ | Bimodal τ 전이 | $1 + \alpha\,\sigma_k(\tau_c{-}\tau)$ | 2 ($k, \tau_c$) + 1 ($\alpha$) |
| $C_{\mathrm{pf}}(p)$ | P:S packing 최적점 | $1 + \beta\,\sigma_k(p{-}p_c)$ | 3 ($k, p_c, \beta$) |
| $G(\tau,p)$ | 교차 coupling | $\exp[\gamma\,\sigma_\tau\,\sigma_p]$ | 1 ($\gamma$) |

**총 자유 파라미터 ≈ 7개** (물리적으로 해석 가능한 sigmoid 모양 변수).

확인하려면:

```bash
grep -n "C_blend\|C_pf\|sigmoid\|tau_c\|pc\s*=" scripts/generate_comparison_plots.py
```

에서 실제 구현 확인 가능.
