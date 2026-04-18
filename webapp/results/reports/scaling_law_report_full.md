# DEM-Based Ionic Conductivity Scaling Law: Complete Report

**Project**: ASSB Composite Cathode DEM Analysis
**Generated**: 2026-04-14
**Scripts Analyzed**: 50 screening scripts (15 ionic + 35 electronic)
**Final Formula**: FORM X — σ = C × σ_grain × (φ_SE - φc)^¾ × CN × √cov / √τ

---

## Part 1: Network Solver Methodology

### 1.1 Kirchhoff Resistor Network

DEM 시뮬레이션의 입자 접촉 정보를 기반으로 Kirchhoff 저항 네트워크를 구성.
각 SE-SE (또는 AM-AM) 접촉이 하나의 저항기(edge)가 됨.

```
R_total = R_bulk + R_constriction  (직렬 연결)
```

**R_bulk** (벌크 저항):
- R_bulk = d / (σ × π × r²)
- d: 입자 중심 간 거리 (hop distance)
- r: 입자 반경
- 입자 내부를 통과하는 전류의 저항

**R_constriction** (Maxwell 확산 저항, Holm 1967):
- R_constr = 1 / (2σa)
- a = √(A_contact / π): 접촉 반경
- A_contact: 접촉 면적 (DEM에서 직접 계산)
- **접촉이 작을수록 저항 급증** → bottleneck

### 1.2 Three Decomposition Runs

| Run | R_bulk | R_constriction | 의미 | 용도 |
|-----|--------|----------------|------|------|
| FULL | ✓ | ✓ | σ_full (ground truth) | 실제 전도도 |
| CONTACT_FREE | ✓ | ✗ (R_constr=0) | σ_cf (upper bound) | 이상적 접촉 |
| CONSTRICTION_ONLY | ✗ (R_bulk=0) | ✓ | σ_constr | 확산 저항 한계 |

**핵심 결과**:
- Bulk resistance fraction: ~20-30%
- **Constriction이 70-80% 지배** → 접촉 면적이 전도도의 핵심
- R_brug/R_full = 3~10× (Bruggeman 과대추정)

### 1.3 Three Transport Modes

| Mode | Network | σ_bulk | 물리 |
|------|---------|--------|------|
| Ionic | SE-SE contacts | 3.0 mS/cm (LPSCl grain) | 이온 전도 |
| Electronic | AM-AM contacts | 50 mS/cm (NCM) | 전자 전도 |
| Thermal | ALL contacts | k_AM=4e-2, k_SE=0.7e-2 W/(cm·K) | 열전도 |

### 1.4 Solver Implementation

1. NetworkX 그래프 구성 (atoms → nodes, contacts → edges)
2. Bottom/Top 경계 입자 식별
3. Kirchhoff 행렬 (Laplacian) 구성
4. scipy.sparse.linalg.spsolve로 전위 분포 계산
5. Ohm's law로 σ_eff 산출: σ_eff = G_eff × L / A

---

## Part 2: Ionic Scaling Law Evolution (15 Scripts)

### 2.1 출발점: Bruggeman EMT

```
σ_brug = σ_grain × φ_SE × f_perc / τ²
```
- R² ≈ 0.30 (Network solver 대비)
- **문제**: 접촉 저항 완전 무시, τ² 과도한 벌점

### 2.2 v2.0 — Exhaustive Screening (screening_v2.py)

**목표**: 모든 1/2/3변수 조합을 σ_brug 기반으로 탐색

**테스트한 변수들**:
- GB density (GB_d), Path conductance (G_path)
- SE-SE CN, φ_AM, τ, f_perc
- Contact area, Bottleneck area

**결과**:
```
Champion: σ_brug × C × (G_path × GB_d²)^(1/4) × CN²
R² = 0.940 (전체), R²_thin ≈ 0.45
```

**발견**: GB density와 path conductance가 Bruggeman 보정에 효과적
**문제**: thin electrode (τ > 2.5)에서 R² < 0.50

### 2.3 v2.0b — φ_AM Addition (screening_v2b.py)

**목표**: thin100 케이스 커버를 위해 φ_AM 추가

**테스트**:
- Strategy A: σ_brug × (G×d²)^¼ × CN² × φ_AM^c
- Strategy B: C × φ_SE^a × φ_AM^b × contact^c
- Strategy D: free fit

**결과**:
```
Best (fixed): R² ≈ 0.934, thin R² ≈ 0.50
Free fit: σ_brug × (G×d²)^0.189 × CN^2.149 × φ_AM^-0.278
```

**발견**: φ_AM이 thin regime에서 필수적
**문제**: τ > 2.5에서 fixed exponent R² < 0.5

### 2.4 v2.0c — τ 극복 (screening_v2c.py)

**목표**: 비표준 τ 보정 (τ³, τ⁴, exp(-dτ), CN/τ^n)

**테스트한 형태들**:
- σ_grain × φ_SE × f_perc / τ^n (n=2~5)
- C × φ_SE^a × φ_AM^b × CN^c × exp(-d×τ)
- C × (φ_SE×f_perc)^a / τ^b × CN^c

**결과**:
```
Best universal: R² ≈ 0.92 (free 5-variable)
Best thin specialist: R² ≈ 0.65 (여전히 낮음)
```

**핵심 발견**: Softened Bruggeman (τ^1.5 또는 τ^1)이 thin에 유리하지만 thick을 악화
**결론**: 단일 τ 지수로는 모든 regime 커버 불가능

### 2.5 v2.0d — φ_SE × φ_AM Deep Dive (screening_v2d.py)

**목표**: φ_SE와 φ_AM의 결합 효과 + CN 편상관 분석

**핵심 테스트**:
- 2변수: φ_SE^a × φ_AM^b
- 3변수 free: φ_SE^7.0 × φ_AM^5.5 × CN^-0.7 (R²=0.95!)

**발견**:
- **CN의 음의 편상관**: φ_SE 통제 후 CN이 σ와 음의 관계
- 이유: CN↑ → 접촉 많지만 각 접촉이 약해짐 (constriction 효과)
- φ_SE^7은 비현실적 → 물리적 해석 필요

### 2.6 v2.0e — Physics-Driven (screening_v2e.py)

**목표**: 물리적 분해 (softened Bruggeman + φ_AM + percolation weighting)

**테스트**:
- σ_soft(τ^n) × φ_AM^a × CN^b
- (φ×f_perc)^a × φ_AM^b
- Physical: φ^1.5 × CN² × φ_AM^b

**결과**:
```
Best: (φ×f_perc)^5.0 × φ_AM^2.0 → R² ≈ 0.93
```

**발견**: Percolation-weighted φ_SE가 별도 항보다 효과적

### 2.7 v2.0f — Balanced Optimization (screening_v2f.py)

**목표**: √(R²_thick × R²_thin) 균형 메트릭 도입

**접근**: 6가지 형태 동시 비교
- (φ/τ)^a × φ_AM^b × CN^c
- σ_brug × (CN/τ)^k × φ_AM^m
- Percolation (φ-φc)^t

**결과**: Best balanced ≈ 0.93
**핵심**: 균형 메트릭이 essential → thick만 좋은 식은 가치 없음

### 2.8 v2.0g — Percolation Threshold 도입 (screening_v2g.py)

**⚡ BREAKTHROUGH MOMENT ⚡**

**목표**: Bruggeman 구조에 percolation threshold 추가

**핵심 전환**:
```
기존: σ ∝ φ_SE / τ²
새로운: σ ∝ (φ_SE - φ_c)^n / τ^b
```

**테스트**: φ_c ∈ [0.12, 0.15, 0.18, 0.20]

**결과**:
```
Best: σ_grain × (φ-0.18)^1.5 × φ_AM × CN → R² = 0.94+
```

**핵심 발견**: **(φ-φc) framework가 plain φ_SE보다 근본적으로 우수**
→ Bruggeman의 φ를 (φ-φc)로 대체하는 패러다임 전환

### 2.9 v2.0h — σ_brug Hybrid (screening_v2h.py)

**목표**: σ_brug 구조 안에서 percolation 통합

**Framework**:
```
σ = σ_brug × C × R_perc × R_comp × CN² × R_contact
R_perc = √(1 - φc/φ)
R_comp = τ^1.5 / f_perc
```

**결과**: R² = 0.943 전체, thick 0.88, mid 0.90, thin 0.65
**발견**: τ^1.5 보정이 σ_brug의 τ^-2를 부분 상쇄 → 실효 τ^-0.5

### 2.10 v2.0i — 미세조정 (screening_v2i.py)

**목표**: 0.25 단위 지수 최적화

**테스트**: φ_c, p, τ_b, f_c, CN_exp, φ_AM 전부 0.25 step

**결과**:
```
Best: σ_brug × √(1-0.18/φ) × τ^1.5 / f × CN^1.5
R² = 0.944, LOOCV = 0.943
```

**결론**: 지수 미세조정으로 ΔR² < 0.03 → 구조가 중요, 소수점이 아님

### 2.11 v2.0j — Coverage 발견 (screening_v2j.py)

**목표**: Coverage 변수 추가 + outlier 분석

**Coverage 테스트**: 지수 0~0.5 sweep

**결과**:
```
+cov^0.375: R² = 0.945 → LOOCV = 0.953!
```

**Ablation**: 모든 항 제거 시 R² 유의미 하락 → 모든 항 필수
**Per-case**: 86% within 30%, 98% within 50%

### 2.12 v2.0k — 극한 Tuning (screening_v2k.py)

**목표**: 0.125 단위 + 2nd 추가변수

**결과**:
```
Clean: √perc × τ^(3/2) × CN^(5/4) × ⁴√cov
R² = 0.960 전체, mean |err| = 7.6%
```

**발견**: 2nd 변수 추가 시 ΔR² < 0.005 → 한계 도달

### 2.13 v2.0l — 깔끔한 지수 (screening_v2l.py)

**목표**: n/2, n/4 지수만 사용하는 elegant form 비교

**후보 15개 비교**:
| Form | 식 | R² |
|------|-----|-----|
| FORM A | √[φ(φ-φc) × CN² × cov / τ] | 0.934 |
| FORM B | √[φ(φ-φc)/τ] × CN × ⁴√cov | 0.935 |
| FORM C | √[φ(φ-φc)/τ] × CN^(3/2) × cov^(3/8) | 0.940 |
| FORM D | √[φ(φ-φc)/τ] × CN^(3/2) | 0.925 |

### 2.14 v2.0m — FORM A Stress Test (screening_v2m.py)

**목표**: FORM A의 robustness 검증

**7가지 테스트**:
1. α=0.5 최적성 → ✓ (gap=0.014)
2. φ_c=0.18 최적성 → ✓ (gap=0.004)
3. f_perc 추가 → marginal
4. Free fit → distance 0.058 from FORM A
5. Leave-SE-size-out → R²=0.91~0.94 ✓
6. Outlier robustness → C 변화 1.2% ✓
7. **Percolation term → (φ-φc)^1.5 > φ×(φ-φc)!**

**⚡ FORM A BROKEN ⚡**: (φ-φc)^1.5가 φ×(φ-φc)보다 우수!

### 2.15 v2.0n — FORM X 탄생 (screening_v2n.py)

**⚡ FINAL CHAMPION ⚡**

**Head-to-head 비교**:
| Formula | R² | LOOCV | avg3 |
|---------|-----|-------|------|
| v3 | 0.926 | 0.910 | 0.917 |
| FORM A | 0.934 | 0.932 | 0.929 |
| FORM Z | 0.939 | 0.937 | 0.938 |
| **FORM X** | **0.942** | **0.940** | **0.945** |

**Free fit 검증**:
```
Free: (φ-φc)^0.789 × CN^0.968 × cov^0.516 × τ^-0.495
FORM X: (φ-φc)^0.750 × CN^1.000 × cov^0.500 × τ^-0.500
Distance: 0.039 (excellent match!)
```

### 2.16 screening_ionic_perfect.py — 최종 검증

**φ_c Fine Scan (0.005 단위)**:
| φ_c | R² |
|-----|-----|
| 0.175 | 0.939 |
| 0.180 | 0.943 |
| **0.185** | **0.944** |
| 0.190 | 0.944 |
| 0.195 | 0.937 |

**Path metrics 추가**: Gc, BN → 효과 없음 확인
**결론**: FORM X가 최종 확정

---

## Part 3: FORM X Final Specification

### 3.1 공식

```
σ_ionic = 0.1275 × σ_grain × (φ_SE - 0.185)^(3/4) × CN × √cov / √τ

= 0.1275 × σ_grain × ⁴√[(φ_SE - 0.185)³ × CN⁴ × cov² / τ²]
```

### 3.2 매개변수

| Parameter | Value | Description |
|-----------|-------|-------------|
| σ_grain | 3.0 mS/cm | LPSCl grain interior conductivity |
| φ_c | 0.185 | Percolation threshold (optimized) |
| C | 0.1275 | Single free parameter |

### 3.3 성능

| Metric | Value |
|--------|-------|
| R² (log-space) | 0.944 |
| LOOCV | 0.940 |
| Mean |error| | ~15% |
| Within 20% | 41/49 (84%) |
| n (training) | 49 cases |

### 3.4 각 항의 물리

| 항 | 유래 | 물리 |
|----|------|------|
| (φ_SE-φc)^¾ | 3D percolation theory (t≈2, softened to ¾) | SE가 threshold 이상일 때만 이온 전도 |
| CN | SE-SE coordination number | 접촉 수 = 병렬 이온 경로 수 |
| √cov | AM-SE surface coverage | AM-SE 계면 품질 → 접촉 안정성 |
| 1/√τ | Softened tortuosity | Bruggeman의 1/τ² → 1/√τ (과벌점 완화) |

### 3.5 Ablation Study

| 제거 항 | R² | ΔR² |
|---------|-----|------|
| Full (baseline) | 0.944 | — |
| Remove cov | ~0.91 | -0.03 |
| Remove τ | ~0.93 | -0.01 |
| Remove CN | ~0.85 | -0.09 |
| Remove (φ-φc) | ~0.38 | **-0.56** |

---

## Part 4: Network Solver — Electronic & Thermal

### 4.1 Electronic (AM-AM Network)

Network Solver가 AM-AM 접촉으로 Kirchhoff 네트워크를 구성하여 전자 전도도 계산.
σ_AM = 50 mS/cm (NCM, discharged state)

**2-Regime Structure**:
- T/d ≥ 10 (Thick): topology 지배 → R²=0.97
- T/d < 10 (Thin): contact mechanics + finite-size → R²=0.81

(Electronic scaling law 상세는 별도 보고서에서 다룸)

### 4.2 Thermal (ALL Contacts)

모든 접촉(AM-AM, AM-SE, SE-SE)을 포함한 열전도 네트워크.
- k_AM = 4.0×10⁻² W/(cm·K)
- k_SE = 0.7×10⁻² W/(cm·K)
- AM-SE 접촉: harmonic mean

```
σ_thermal ≈ 286 × σ_ionic^(3/4) × φ_AM² / CN_SE
R² ≈ 0.88
```

---

## Part 5: Evolution Summary

### 5.1 Ionic — 15 Scripts, 4 Phases

**Phase 1 (v2.0~v2.0c): σ_brug 기반 탐색**
- σ_brug × (G_path × GB_d²)^¼ × CN² → R²=0.94
- τ² 과벌점 문제 발견
- Thin (τ>2.5) R² < 0.50

**Phase 2 (v2.0d~v2.0f): 변수 결합 탐색**
- φ_SE × φ_AM 조합
- CN 음의 편상관 발견
- 균형 메트릭 도입

**Phase 3 (v2.0g~v2.0k): Percolation Breakthrough**
- (φ-φc) 도입 → R² quantum jump
- Coverage 발견
- 0.125 단위 극한 tuning

**Phase 4 (v2.0l~v2.0n): FORM X 확립**
- FORM A stress test → broken
- (φ-φc)^¾ > φ×(φ-φc) 발견
- FORM X = final champion

### 5.2 Key Milestones

| Milestone | Script | R² | Description |
|-----------|--------|-----|-------------|
| σ_brug champion | v2.0 | 0.940 | GB density 보정 |
| τ 문제 인식 | v2.0c | 0.920 | τ² 과벌점 |
| CN 편상관 | v2.0d | 0.950 | CN^-0.7 (overfitting) |
| **(φ-φc) 도입** | **v2.0g** | **0.938** | **패러다임 전환** |
| Coverage 발견 | v2.0j | 0.945 | +cov^0.375 |
| FORM A broken | v2.0m | 0.934 | (φ-φc)^1.5 발견 |
| **FORM X 탄생** | **v2.0n** | **0.942** | **최종 champion** |
| φ_c 최적화 | ionic_perfect | 0.944 | φ_c=0.185 |

---

*Report generated from analysis of 50 screening scripts*
*DEM Analyzer v2.0*
