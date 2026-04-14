# Scaling Law Formula Catalog: 전체 시도 수식 목록

**총 시도 수식**: ~200+ (50개 스크립트에서 추출)

---

## A. Ionic Conductivity — 시도한 모든 수식

### A.1 Bruggeman 기반 (v2.0~v2.0c)

```
[BASE] σ_brug = σ_grain × φ_SE × f_perc / τ²                          R²≈0.30

[v2.0-1]  σ_brug × C × (G_path × GB_d²)^(1/4) × CN²                  R²=0.940 ★ v3 champion
[v2.0-2]  σ_brug × C × GB_d^α                                         R²≈0.84
[v2.0-3]  σ_brug × C × G_path^α                                       R²≈0.80
[v2.0-4]  σ_brug × C × CN^α                                           R²≈0.75
[v2.0-5]  σ_brug × C × (GB_d × T)^α                                   R²≈0.82

[v2b-1]   σ_brug × (G×d²)^¼ × CN² × φ_AM^c                           R²≈0.934
[v2b-2]   C × φ_SE^a × φ_AM^b × (G×d²)^c                              R²≈0.92
[v2b-3]   σ_brug × (G×d²)^0.189 × CN^2.149 × φ_AM^(-0.278)           R²≈0.94 (free)

[v2c-1]   σ_grain × φ_SE × f_perc / τ^3 × (G×d²)^a × CN^b            R²≈0.88
[v2c-2]   σ_grain × φ_SE × f_perc / τ^1.5 × CN^b                      R²≈0.90
[v2c-3]   C × φ_SE^a × φ_AM^b × CN^c × exp(-d×τ)                      R²≈0.89
[v2c-4]   C × (φ_SE×f_perc)^a / τ^b × CN^c                            R²≈0.92
```

### A.2 φ_SE × φ_AM 결합 (v2.0d~v2.0f)

```
[v2d-1]   φ_SE^7.0 × φ_AM^5.5 × CN^(-0.7)                            R²=0.95 (overfitting!)
[v2d-2]   σ_brug^1.01 × φ_AM^2.05 × CN^(-0.68)                       R²≈0.94
[v2d-3]   φ_SE^a × φ_AM^b (a=2~8, b=0~7.5)                           R²≈0.92 max

[v2e-1]   σ_soft(τ^n) × φ_AM^a × CN^b                                 R²≈0.93
[v2e-2]   (φ×f_perc)^5.0 × φ_AM^2.0                                   R²≈0.93
[v2e-3]   φ^1.5 × CN^2 × φ_AM^b                                       R²≈0.91

[v2f-1]   (φ/τ)^a × φ_AM^b × CN^c                                     R²≈0.93 balanced
[v2f-2]   (φ×φ_AM)^a × (CN/τ)^b                                       R²≈0.91
[v2f-3]   σ_brug × (CN/τ)^k × φ_AM^m                                  R²≈0.92
```

### A.3 ⚡ Percolation Threshold 도입 (v2.0g~v2.0i)

```
[v2g-1]   σ_grain × (φ-0.12)^1.5 × CN                                 R²≈0.92
[v2g-2]   σ_grain × (φ-0.15)^1.5 × φ_AM × CN                          R²≈0.93
[v2g-3]   σ_grain × (φ-0.18)^1.5 × φ_AM × CN                          R²≈0.94 ★★
[v2g-4]   σ_grain × (φ-0.20)^1.5 × φ_AM × CN                          R²≈0.93

[v2h-1]   σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.5                  R²=0.943
[v2h-2]   σ_brug × (1-φc/φ)^0.5 × τ^b × CN^d                         R²≈0.94
[v2h-3]   σ_brug × √(1-φc/φ) × τ^1.5 / f × CN^1.5 × (G×d²)^c       R²≈0.944

[v2i-1]   σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.5                  R²=0.944 LOOCV=0.943
[v2i-2]   (위 + φ_AM^0.25)                                              R²≈0.945
[v2i-3]   (위 + G_path^0.125)                                           R²≈0.944
```

### A.4 Coverage + Fine-tuning (v2.0j~v2.0k)

```
[v2j-1]   σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.25 × cov^0.375    R²=0.945 LOOCV=0.953
[v2j-2]   σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.5 × cov^0.25      R²≈0.944
[v2j-3]   σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.5 × cov^0.5       R²≈0.943

[v2k-1]   √perc × τ^(3/2) × CN^(5/4) × ⁴√cov                         R²=0.960
[v2k-2]   √perc × τ^(3/2) × CN^(5/4) × cov^(3/8) × G_path^0.125     R²≈0.961
[v2k-3]   √perc × τ^(3/2) × CN^(3/2) × cov^(1/4)                     R²≈0.955
```

### A.5 Elegant Forms (v2.0l)

```
[FORM A]  √[φ(φ-φc) × CN² × cov / τ]                                  R²=0.934
[FORM B]  √[φ(φ-φc)/τ] × CN × ⁴√cov                                   R²=0.935
[FORM C]  √[φ(φ-φc)/τ] × CN^(3/2) × cov^(3/8)                        R²=0.940
[FORM D]  √[φ(φ-φc)/τ] × CN^(3/2)                                     R²=0.925
```

### A.6 FORM A Stress Test → FORM X (v2.0m~v2.0n)

```
[v2m] FORM A tests:
  α=0.3~0.7 sweep → optimal α=0.486 (≈0.5 ✓)
  φ_c=0.14~0.22 → optimal φ_c=0.176 (≈0.18 ✓)
  (φ-φc)² → R²=0.930 (worse)
  (φ-φc)^1.5 → R²=0.936 ★ (better than FORM A!)

[v2n] HEAD-TO-HEAD:
  v3:     σ_brug × (G×d²)^¼ × CN²                                     R²=0.926
  FORM A: √[φ(φ-φc) × CN² × cov / τ]                                  R²=0.934
  FORM Z: (φ-φc) × CN × √cov / √τ                                     R²=0.939
  FORM W: (φ-φc)^¾ × CN^(3/2) × cov^(3/8) / √τ                       R²=0.938
  FORM V: (φ-φc) × CN^(5/4) × √cov / √τ                               R²=0.939
  FORM X: (φ-φc)^¾ × CN × √cov / √τ                                   R²=0.942 ★★★ CHAMPION
  FORM Y: ⁴√[(φ-φc)³ × CN⁴ × cov² / τ²]                              R²=0.942 (=FORM X)
```

### A.7 φ_c Optimization (screening_ionic_perfect.py)

```
  φ_c=0.150: R²=0.916     φ_c=0.170: R²=0.935
  φ_c=0.175: R²=0.939     φ_c=0.180: R²=0.943
  φ_c=0.185: R²=0.944 ★★★ OPTIMAL
  φ_c=0.190: R²=0.944     φ_c=0.200: R²=0.918
```

### A.8 ★ FINAL CHAMPION

```
σ_ionic = 0.1275 × σ_grain × (φ_SE - 0.185)^(3/4) × CN × √cov / √τ

         = 0.1275 × 3.0 × ⁴√[(φ_SE - 0.185)³ × CN⁴ × cov² / τ²]

R² = 0.944, LOOCV = 0.940, φ_c = 0.185
```

---

## B. Electronic Conductivity — 시도한 모든 수식

### B.1 초기 탐색 + exp(π/(T/d))

```
[el-1]    σ_AM × φ^1.5 × CN² × exp(π/(T/d))                           R²=0.900 (unified)
[el-2]    σ_AM × φ^3.45 × CN^1.72 (free fit)                          R²≈0.91
[el-3]    σ_AM × φ^1.5 × CN² × exp(0.7π/(T/d))                        R²≈0.91
[el-4]    σ_AM × φ^2.5 × CN² × sigmoid(T/d)                           R²≈0.90
[el-5]    σ_AM × φ^1.5 × CN² × tanh(T/d)                              R²≈0.88
```

### B.2 2-Regime 발견 (T/d=10 분할)

```
THICK (T/d ≥ 10):
[tk-1]    φ^4 × CN^1.5                                                  R²=0.963
[tk-2]    φ^4 × CN^1.5 × cov                                           R²=0.963
[tk-3]    φ^4 × CN^1.5 × cov × √τ                                     R²=0.969 ★
[tk-4]    φ^4 × CN^2 × cov                                             R²=0.952
[tk-5]    φ^3.5 × CN^1.5 × cov                                         R²=0.949
[tk-6]    φ^4.5 × CN^1.5 × cov                                         R²=0.964
[tk-7]    φ^4 × CN^1.5 × cov × √P_contact                             R²=0.963
[tk-8]    φ^4 × CN^1.5 × cov × √δ                                     R²=0.936
[tk-9]    φ^4 × CN^1.5 × cov × √a_c                                   R²=0.952
[tk-10]   (φ-0.15)^3 × CN^1.5 × cov                                   R²=0.963
[tk-11]   (φ-0.20)^2 × CN^2 × cov                                     R²=0.956

THIN (T/d < 10) — 초기:
[tn-1]    CN^1.25 × por^3 × cov^1.25 × φ_SE^1.25 / (T/d)^0.5        R²=0.918
[tn-2]    CN × por^2.5 × cov^1.5 / √(T/d)                             R²=0.899
[tn-3]    P_contact^3 × CN / √(T/d)                                    R²=0.934
[tn-4]    CN × por^2 × φ_SE^1.5 / √(T/d)                              R²=0.894
[tn-5]    CN × por^2 × cov / √(T/d)                                    R²=0.884
```

### B.3 Universal (exp) 시도

```
[uni-1]   φ² × CN² × √cov × exp(π/(T/d))                              R²=0.913
[uni-2]   φ^1.5 × CN² × exp(π/(T/d))                                  R²=0.900
[uni-3]   φ^1.5 × CN² × (1 + 9.64/(T/d))                              R²=0.892
[uni-4]   φ^1.5 × CN² × (T/d)^(-0.59)                                 R²=0.865
[uni-5]   φ^1.5 × CN² × exp(β/(T/d)), β=3.096≈π                      R²=0.900 ★ β≈π!
```

### B.4 접촉 역학 도입

```
[cm-1]    δ^0.5 × CN × (T/d)^(-0.5)                                    R²=0.891
[cm-2]    (δ²/A)^0.75 × CN / √(T/d)                                    R²=0.903
[cm-3]    (δ²/A)^0.5 × CN × √hop / √(T/d)                             R²=0.919
[cm-4]    CN × (δ²/A)^(5/8) × √hop / √(T/d)                           R²=0.922
[cm-5]    CN^0.875 × (δ²/A)^0.75 × hop^0.25 / (T/d)^0.5              R²=0.936
[cm-6]    P^3 × CN / √(T/d)                                            R²=0.934 (P 불안정)
[cm-7]    √(CN×A) × por / √(T/d)                                       R²=0.681
```

### B.5 δ/R Hertz + τ proxy

```
[hz-1]    φ^4 × (δ/R) × CN × √τ                                       R²=0.694
[hz-2]    φ^4 × (δ/R)^(3/4) × √τ                                      R²=0.788
[hz-3]    φ^4 × √τ × δ² / A                                           R²=0.864
[hz-4]    φ^3.75 × √τ × δ² / A^0.75                                   R²=0.894
[hz-5]    φ^3.5 × (T/d)^(-0.25) × δ^2.25 / A                          R²=0.904
```

### B.6 T/d + 접촉역학 짬뽕

```
[mx-1]    φ^3.5 × (T/d)^(-0.25) × δ^2.25 / A                          R²=0.904
[mx-2]    φ^3 × (T/d)^(-0.25) × (δ²/A) × hop^0.5 × CN^0.5           R²=0.910
[mx-3]    φ^2 × (T/d)^(-0.25) × (δ²/A) × hop^0.5 × CN^0.5           R²=0.905
[mx-4]    φ^4 × exp(π/ξ) × δ² / A                                     R²=0.710
```

### B.7 Thin 개형 최적화 — Spearman ρ 기반

```
Intra-group Spearman ρ (개형 정확도):
  por×cov:    ρ=+0.806 CONSISTENT [+1.00, +0.82, +0.60] ★★
  A:          ρ=+0.727 CONSISTENT
  δ:          ρ=+0.679 CONSISTENT
  hop:        ρ=+0.645 CONSISTENT
  CN:         ρ=-0.512 NOT consistent (75:25에서 반전!)
  τ:          ρ=-0.135 mixed
  cov:        ρ=+0.033 거의 무관

[shape-1]   por × cov / √(T/d)                                         ρ=0.81, R²=0.39
[shape-2]   por × cov / (φ_SE × √(T/d))                                ρ=0.81, R²=-0.001
[shape-3]   √(por×cov) / (φ_SE × √(T/d))                              ρ=0.81, R²=0.269
[shape-4]   ⁴√(por×cov) / (φ_SE × √(T/d))                             ρ=0.81, R²=0.364
[shape-5]   por^0.6 × cov^0.3 / (√φ_SE × √(T/d))                      ρ≈0.80, R²=0.138
```

### B.8 Thin 접촉역학 조합 (Spearman 스크리닝)

```
[sp-1]    d^0.3 × CN^0.2 / φ_SE                                        ρ=0.917, R²=0.771
[sp-2]    hop^0.3 × CN^0.3 / φ_SE                                      ρ=0.934, R²=0.723
[sp-3]    d^0.3 × τ^0.5 / φ_SE                                         ρ=0.946, R²=0.658
[sp-4]    hop^0.5 × τ^0.7 / √φ_SE                                      ρ=0.963, R²=0.593
[sp-5]    d^0.5 × CN^0.3 / √φ_SE                                       ρ=0.848, R²=0.817 ★
```

### B.9 Thin por×cov vs CN 갈등

```
Simpson's Paradox 발견:
  intra-group: cov↑ → σ↑ (ρ=+0.81)
  inter-group: cov↑ → σ↓ (75:25 cov높은데 σ낮음)

[pc-1]    por × cov² / √(T/d)                                          R²=-8.599 (level 반전!)
[pc-2]    por × cov / (φ_SE × √(T/d))                                  level 맞지만 개형 급함
[pc-3]    por^1.0 × cov^0.2 / (√φ_SE × √(T/d))                        ρ=0.922, range ±139%
[pc-4]    por^0.5 × cov^0.2 / (√φ_SE × √(T/d))                        range ±더 좁음
[pc-5]    por^0.7 × cov^0.2 / (√φ_SE × √(T/d))                        중간
```

### B.10 Thin hop × CN × δ 최종 조합

```
[final-1]  hop^0.3 × CN^0.3 / (φ_SE × √(T/d))                         R²=0.556
[final-2]  hop^0.3 × CN^0.3 / (φ_SE^1.25 × √(T/d))                    R²=0.349
[final-3]  hop^0.3 × CN^0.3 / (φ_SE^1.5 × √(T/d))                     R²=-0.021
[final-4]  hop^0.2 × CN^0.2 / (φ_SE × √(T/d))                         R²=0.584, 12/12 within 20%
[final-5]  hop^0.35 × CN^0.4 / (φ_SE^0.85 × √(T/d))                   R²=0.612
[final-6]  hop^0.25 × CN^0.4 × δ^0.2 × f_p^0.15 / (φ_SE^0.85 × √ξ)  R²=0.811, 18/20 ★★
```

### B.11 ★ THICK FINAL

```
σ_el(thick) = 0.79 × σ_AM × φ⁴ × CN^(3/2) × cov × √τ

R² = 0.969, LOOCV = 0.967, 22/24 within 20%
```

### B.12 ★ THIN (진행 중)

```
σ_el(thin) = C × σ_AM × hop^0.25 × CN^0.4 × δ^0.2 × f_p^0.15 / (φ_SE^0.85 × √ξ)

R² = 0.811 (전체 thin 20개), 18/20 within 20%
```

---

## C. Thermal Conductivity

```
σ_thermal = 286 × σ_ionic^(3/4) × φ_AM² / CN_SE
R² ≈ 0.88
```

---

## D. 실패한 접근들 (교훈)

| 시도 | 왜 실패 | 교훈 |
|------|---------|------|
| φ_SE^7 × φ_AM^5.5 | 과적합, 물리 없음 | 지수 > 3은 위험 |
| τ^4, τ^5 보정 | thin/thick 동시 만족 불가 | τ² 자체가 문제 |
| CN^(-0.7) | thin에서만 유효 | 편상관 ≠ 인과관계 |
| exp(π/(T/d)) × 접촉역학 | exp가 thin에서 폭발 | exp는 그룹 간 안정성 파괴 |
| CN × A (총 접촉면적) | 스파이크 유발 | A 변동 너무 큼 |
| CN × √Gc (path cond) | 스파이크 | mean이 network 거동 못 잡음 |
| φ^(-3) × CN × por² | CN 방향 뒤집음 | φ의 음수 지수가 CN과 싸움 |
| cn_std (CN 표준편차) | inter-group R² 악화 | intra-group ρ ≠ inter-group R² |
| (φ-φc) for electronic | thick에서 효과 없음 | AM percolation ≠ SE percolation |

---

*Total formulas cataloged: ~120 unique expressions*
*From 50 screening scripts (15 ionic + 35 electronic)*
