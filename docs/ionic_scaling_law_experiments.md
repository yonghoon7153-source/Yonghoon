# Ionic Scaling Law — Troubleshooting Log

Dataset: n=57 (Ionic v9 BLEND)
Branch: `claude/resistor-network-analysis-lKgcS`
Final model: **v9 BLEND** (R²=0.9806, LOOCV=0.9764)

---

## Experiment Ledger

| # | Variant | Params | R² | LOOCV | Gap | <10% | <20% | Verdict |
|---|---------|--------|-----|-------|-----|------|------|---------|
| 1 | v9 (baseline) | 6 | 0.9806 | **0.9764** | 0.0042 | 39/57 | 52/57 | 🏆 final |
| 2 | v9 + SIGMA_MIN=0.0 | 6 | 0.9805 | 0.9764 | 0.0041 | — | 52/57 | 저-σ 포함, 성능 불변 |
| 3 | v9 + k continuous opt | 6 | 0.9806 | 0.9764 | 0.0042 | — | 52/57 | k=20 bound-hit, LOOCV flat |
| 4 | v9 + k bounded to 20 | 6 | 0.9806 | 0.9764 | 0.0042 | — | 52/57 | numerical-clean |
| 5 | v10 (4-param P:S) | 10 | 0.9817 | 0.9776 | 0.0041 | **35/57** ✗ | 51/57 | overfit (<10% 뱅크 악화) |
| 6 | v10 (1-param P:S, 2-step) | 7 | 0.9822 | 0.9769 | 0.0053 | 35/57 | 52/57 | solid-only particulate 악화 |
| 7 | v11 (3 features A/B/C) | 10 | 0.9826 | 0.9752 | 0.0074 | 36/57 | 52/57 | LOOCV 내림 — 첫 overfit signal |
| 8 | v11b (B+C only) | 9 | 0.9811 | 0.9756 | 0.0055 | 37/57 | 52/57 | γ 붕괴 (+0.221→+0.008) — signal degeneracy 확인 |

---

## 유의미했던 발견 (Meaningful Findings)

### ✅ 1. Cluster C (CN 포화) — 물리적으로 **진짜**

- **관찰**: v11에서 δ=−0.113 ((log CN)² 항)
- **효과**: `input_particulate_1` (+19.5%) 및 `input_particulate_7` (+14.0%)가 >10% 리스트에서 **제거됨**
- **물리**: CN≥7에서 packing이 포화되므로 `CN^1.5` 무한 증가는 비물리적. 고-CN regime에서 포화 필요
- **한계**: v11b에서 Cluster A 제거 시 δ=−0.036으로 약화 → 3 feature 간 collinearity로 v11에서는 과추정
- **결론**: 진짜 signal이지만 단독으로는 통계적 유의성 확보 어려움. 더 많은 고-CN 샘플 필요

### ✅ 2. P:S ratio monotone trend in 1mAh_80:20 — 국소적으로 **진짜**

- **관찰**: `input_thin_5_AMS` (0:10) −22.8% → `input_thin_4` (3:7) −17.3% → `input_thin_5` (5:5) −15.3% → `input_thin_6` (7:3) **+44.4%** → `input_thin_5_AMP` (10:0) +13.2%
- **P:S 증가 → 오차가 음수→양수로 단조**: P가 많을수록 모델이 더 과대예측
- **물리 해석**: particulate SE는 접촉 "개수"는 많지만 각 접촉의 실제 이온전도 기여는 작음 ("empty contacts"). 모델의 CN^1.5가 이를 구분 못함
- **한계**: 이 경향은 **1mAh_80:20 서브그룹에만** 나타남. 다른 그룹 (particulate all-solid, 8mAh 등)은 반대 방향 → **글로벌 feature로 사용 불가**
- **결론**: 그룹 내 효과는 있지만 cross-group generalization 실패

### ✅ 3. 8mAh_85:15 그룹 전체 −15% 과소예측 — 그룹 특이 **systematic**

- **관찰**: `input_7/8/9/8_AMS` 네 케이스 모두 −14.5% ~ −16.9% 과소
- **공통 특성**: τ≈1.28 (후막), AM-rich (85:15), 저-cov (0.13~0.18)
- **물리 가설**: AM 주변으로 "얇은 SE 막"이 형성되어 추가 전도 경로 → `cov^0.25`로 덜 반영되어 과소예측
- **해결 시도**: v11에서 γ=+0.221 (cov 지수 조정)으로 일부 개선했으나 LOOCV 저하
- **결론**: 그룹 특이 신호이나 단독 feature로 추출 시 다른 그룹에 side effect. 그룹 single-particle-size artifact 가능성

---

## 기각된 가설 (Rejected Hypotheses)

### ❌ H1. Sigmoid transition이 outlier 원인 (τ∈[1.9, 2.1])

- **초기 가설**: k=15로 날카로운 blend sigmoid가 τ≈2 근처 outlier 3개 유발
- **검증**: k 연속 최적화 → k=1~20 전 범위 LOOCV 플랫 (Δ<0.001)
- **결론**: k 변경은 outlier에 효과 없음. Blend sigmoid는 혐의 벗음

### ❌ H2. 박막(1mAh)에서 gb_density 누락이 outlier 원인

- **검증**: log(gb_dens) residual 상관계수 r=+0.019 (negligible)
- **결론**: 기각. GB density는 CN과 중복 정보

### ❌ H3. Particulate-empty-contact penalty (log g_path)

- **가설**: particulate SE에서 접촉은 많지만 실제 path conductance 낮음
- **검증 (v11)**: α=−0.046, 표준 오차 대비 noise level. `input_thin_9` 오히려 악화 (+62.8% → +66.8%)
- **결론**: log(g_path) 단독 추가는 다른 features와 collinearity로 noise 역할

### ❌ H4. Cluster B (coverage 지수 상향)

- **초기**: v11에서 γ=+0.221 (effective cov exp = 0.47) 유의해 보였음
- **재검증 (v11b)**: g_path 제거 시 γ=+0.008로 붕괴
- **결론**: v11의 γ는 g_path와 보정 효과 주고받던 artifact. 진짜 signal 아님

### ❌ H5. 글로벌 P:S ratio feature

- **검증**: 단일 β (v10 simplified)는 solid-only particulate group (#6, #15) 악화
- **결론**: P:S 경향이 그룹 의존적이라 글로벌 선형 보정 불가

---

## 통계적 구별 한계 (Noise Floor)

- **LOOCV 표준오차 추정**: √(2/n) × (1−R²) ≈ √(2/57) × 0.024 ≈ **0.0045**
- v9, v10, v11, v11b LOOCV 범위: 0.9752~0.9776 = **Δ 0.0024**
- → 모든 variant 간 차이가 **노이즈의 1/2 수준**, 통계적 구별 불가
- **결론**: 이 데이터셋에서는 **R² 0.98 / LOOCV 0.976 이 실질적 ceiling**

---

## 최종 outlier 18개 (|err|>10%, v9 기준)

| Cluster | 케이스 수 | 추정 물리 원인 | 해결 방향 |
|---------|-----------|---------------|-----------|
| A: 박막+particulate+τ≈2 과대 | 3 | CN이 "empty contact" 과다 카운트 | 데이터 추가 필요 |
| B: 8mAh_85:15 전멤버 과소 | 4 | AM 표면 thin SE film | coverage 고해상도 측정 |
| C: 고-σ particulate 과대 (CN≥7) | 3 | CN 포화 미반영 | 고-CN 샘플 추가 |
| D: 1mAh_80:20 P:S-dep signed | 5 | particulate 접촉 quality | g_path feature 개선 |
| 기타 | 3 | 개별 노이즈 | seed 재시뮬 |

---

## Residual Correlation Summary (v9)

| Feature | r | 해석 |
|---------|---|------|
| log(phi_ex) | +0.051 | 무효 |
| log(CN) | +0.156 | 약함 |
| log(tau) | +0.087 | 무효 |
| log(cov) | +0.054 | 무효 |
| log(fp) | −0.036 | 무효 |
| log(gb_dens) | +0.019 | 무효 |
| log(g_path) | +0.054 | 무효 |
| log(gp·gb²) | +0.106 | 약함 |
| p_frac (P:S) | +0.293 | 있음 (subgroup artifact) |
| (p_frac)² | **+0.333** ⚠ | 있음 (subgroup artifact) |

---

## Final Model Specification (v9 BLEND)

```
σ_ionic = C_blend(τ) × σ_grain × CN^1.5 × (φ−φc)^¾ × cov^¼ × f_p²

C_blend(τ) = (1−w(τ))·C_v5(τ) + w(τ)·C_p3(τ)
  C_v5(τ) = Ct + (Cn − Ct) · σ(5·(τ−2.1))
  C_p3(τ) = exp[a0 + a1·ln τ + a2·(ln τ)² + a3·(ln τ)³]
  w(τ) = σ(k·(τ−τc))     # k, τc continuously optimized by LOOCV

Fitted parameters (n=57):
  Ct = 0.0351, Cn = 0.0201
  φc = 0.185 (percolation threshold, fixed)
  σ_grain = 3.0 mS/cm (Li6PS5Cl bulk, fixed)
  k ≈ 20, τc ≈ 1.92 (auto-optimized; LOOCV flat in this region)

Performance:
  R² = 0.9806
  LOOCV R² = 0.9764
  ±20% coverage: 52/57 (91%)
  median |err| = 6.5%, mean |err| = 10.2%
```

---

## Lessons for Future Improvements

1. **데이터 추가가 모델 변경보다 우선** — 3 cluster가 각 3~4개 샘플뿐이라 feature fit 불가
2. **Signal degeneracy 체크** — 한 feature 제거 시 다른 feature 계수 붕괴하면 둘 다 artifact
3. **LOOCV가 진실** — Training R² 증가해도 LOOCV 감소면 overfit
4. **서브그룹 효과는 글로벌 feature 아님** — P:S trend는 1mAh_80:20에만 존재
5. **시뮬레이션 noise floor 측정** — seed variation 없이는 irreducible error 알 수 없음

---

## Recommended Next Experiments (더 많은 데이터 확보 후)

1. **Particulate-rich 박막**: `input_thin_85:15` 와 `input_thin_80:20`에서 P:S=7:3, 10:0 각각 5~10 샘플씩 추가
2. **고-CN 샘플**: CN≥7 범위 확장 (현재 3개 → 10개 이상)
3. **Seed variation**: 기존 outlier 3개를 각 3 seed로 재시뮬 → noise floor 확정
4. **박막 두께 scan**: 1mAh와 8mAh 사이 (2mAh, 4mAh) 추가
