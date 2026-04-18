# DEM Simulation Methodology & References Database

## 1. Contact Model & Material Properties

### 1.1 SE Effective Young's Modulus Calibration

**Method**: SE-only pellet 실험 porosity(~10% at 300 MPa)에 DEM을 매칭하여 E_eff 결정

```
실험: LPSCl pellet, 300 MPa cold-pressing → porosity ≈ 10%
DEM:  E_SE 변경하면서 선형회귀
      → E_eff = 1.35 GPa (문헌 E_bulk = 24 GPa)
      → 소성변형 + 소결 효과를 effective modulus로 반영
```

**DEM Parameters:**
| Parameter | AM (NCM) | SE (LPSCl) | 비고 |
|-----------|----------|------------|------|
| E (Young's) | 140 GPa | 1.35 GPa | SE는 calibrated |
| E_literature | 140 GPa | 22-24 GPa | |
| Contact model | hooke/hysteresis (Luding) | 같음 | 탄소성 이력 |
| Pressure | 300 MPa | | 일축 가압 |
| Scale factor | 1000 | | 계산 효율 |

### 1.2 왜 E_eff가 E_bulk보다 낮은가

- LPSCl은 300 MPa에서 소성변형 + 부분 소결 발생
- DEM의 hooke/hysteresis 모델은 제한적 소성만 표현
- E_eff를 낮춰서 동일 하중에서 더 큰 변형 → 실제 소성변형 효과 재현
- **실험 porosity 매칭으로 검증된 접근**

---

## 2. Packing Theory References

### 2.1 Binary Particle Packing

**Furnas, C.C.** (1931). "Grading Aggregates: I—Mathematical Relations for Beds of Broken Solids of Maximum Density." *Industrial & Engineering Chemistry*, 23(9), 1052-1058.
- 최적 소입자 분율 x_S ≈ 0.27 (binary packing)
- ε_min = ε₁ × ε₂ (이론 최소 porosity)
- x_S > optimal → matrix disruption → porosity 증가

**McGeary, R.K.** (1961). "Mechanical Packing of Spherical Particles." *Journal of the American Ceramic Society*, 44(10), 513-522.
- 크기비 7:1에서 최고 packing density (~86%)
- 7:1 이상 plateau (10:1, 15:1 유사)
- 소입자 과잉 시 density 감소 실험 확인

**Westman, A.E.R. & Hugill, H.R.** (1930). "The Packing of Particles." *Journal of the American Ceramic Society*, 13(10), 767-779.
- "Matrix disruption" 개념 최초 제안
- x_S > optimal: 소입자가 큰 입자 골격 파괴

**Stovall, T., de Larrard, F., & Buil, M.** (1986). "Linear Packing Density Model of Grain Mixtures." *Powder Technology*, 48(1), 1-12.
- Compressible Packing Model (CPM)
- Binary packing 수학적 모델

**Cumberland, D.J. & Crawford, R.J.** (1987). *The Packing of Particles*. Elsevier.
- 교과서, binary/ternary packing 이론 종합

**German, R.M.** (2014). *Sintering: From Empirical Observations to Scientific Principles*. Elsevier.
- Chapter 3: Bimodal packing theory
- Critical size ratio, interstitial vs replacement

### 2.2 Packing Regime Transition

**Interstitial regime** (d_L/d_S > 5~7):
- 소입자가 큰 입자 사이 void에 들어감
- 과잉 소입자 → void 넘침 → porosity ↑
- tetrahedral void 입구 ≈ 0.155 × d_L

**Replacement regime** (d_L/d_S < 5):
- 소입자가 void에 못 들어감
- 두 종류 입자가 함께 packing
- 더 많은 입자 → 더 높은 packing diversity → porosity ↓

**전환점**: d_L/d_S ≈ 5~7 (void 입구 크기 = 소입자 크기)

---

## 3. DEM Calibration References

**Coetzee, C.J.** (2017). "Review: Calibration of the Discrete Element Method." *Powder Technology*, 310, 104-142.
- DEM 파라미터 캘리브레이션 방법론 총 정리
- E_eff를 실험 매칭으로 결정하는 표준 절차
- "effective E는 소성변형을 포함한 bulk behavior를 재현"

**Marigo, M. & Stitt, E.H.** (2015). "Discrete Element Method (DEM) for Industrial Applications." *KONA Powder and Particle Journal*, 32, 236-252.
- E_eff calibration 방법론
- 실험 매칭 기반 검증

---

## 4. SSB Composite Cathode References

### 4.1 Ionic Transport

**Bielefeld, A., Weber, D.A., & Janek, J.** (2019). "Modeling Effective Ionic Conductivity and Binder Influence in Composite Cathodes for All-Solid-State Batteries." *ACS Applied Materials & Interfaces*, 11(34), 30765-30778.
- Bruggeman vs percolation theory for SSB
- φ_c (percolation threshold) 개념
- AM:SE ratio 최적화

**Bielefeld, A., Weber, D.A., & Janek, J.** (2022). "Influence of Lithium Ion Kinetics, Particle Morphology and Voids on the Electrochemical Performance." *Journal of The Electrochemical Society*, 169, 020539.
- DEM + ionic transport
- E_eff calibration 유사 접근

**Minnmann, P., Quillman, L., Burkhardt, S., Richter, F.H., & Janek, J.** (2021). "Quantifying the Impact of Charge Transport Bottlenecks in Composite Cathodes of All-Solid-State Batteries." *Journal of The Electrochemical Society*, 168(4), 040537.
- Electronic percolation in SSB
- σ_ionic ≈ 0.17 mS/cm (실험 reference)
- AM:SE ratio 효과

### 4.2 SE Pellet Properties

**Doux, J.M., Nguyen, H., Tan, D.H.S., et al.** (2020). "Stack Pressure Considerations for Room-Temperature All-Solid-State Lithium Metal Batteries." *Advanced Energy Materials*, 10(1), 1903253.
- LPSCl pellet 압축 실험
- 250~400 MPa에서 porosity 5~15%
- 300 MPa에서 ~10% → **우리 calibration reference**

**Kraft, M.A., Culver, S.P., Caldwell, M., et al.** (2017). "Influence of Lattice Polarizability on the Ionic Conductivity in the Lithium Superionic Argyrodites." *Journal of the American Chemical Society*, 139(31), 10909-10918.
- Li₆PS₅Cl σ_grain = 3.0 mS/cm
- σ_pellet < σ_grain (grain boundary 효과)

**Janek, J. & Zeier, W.G.** (2016). "A solid future for battery development." *Nature Energy*, 1, 16141.
- SE pellet 제조 조건 overview

### 4.3 DEM for SSB

**Shi, T., Tu, Q., Tian, Y., et al.** (2020). "High Active Material Loading in All-Solid-State Battery Electrode via Particle Size Optimization." *Advanced Energy Materials*, 10(1), 1902881.
- DEM으로 SSB cathode packing
- 크기비 효과, packing 최적화
- E_eff calibration 유사 접근

---

## 5. Contact Mechanics References

**Holm, R.** (1967). *Electric Contacts: Theory and Application*. 4th ed., Springer.
- Maxwell constriction resistance: R = 1/(2σa)
- 접촉 반경 a에 반비례

**Thornton, C. & Ning, Z.** (1998). "A theoretical model for the stick/bounce behaviour of adhesive, elastic-plastic spheres." *Powder Technology*, 99(2), 154-162.
- 탄성-소성 접촉 이론
- yield pressure 기반 permanent deformation

**Thornton, C. & Antony, S.J.** (2000). "Quasi-static deformation of particulate media." *Philosophical Transactions of the Royal Society A*, 356, 2763-2782.
- 소성변형 packing

**Luding, S.** (2008). "Cohesive, frictional powders: contact models for tension." *Granular Matter*, 10, 235-246.
- hooke/hysteresis 모델 원논문
- 이력 루프로 소성변형 에너지 반영

---

## 6. Percolation Theory

**Kirkpatrick, S.** (1973). "Percolation and Conduction." *Reviews of Modern Physics*, 45(4), 574-588.
- σ ∝ (p - p_c)^t
- 3D percolation exponent t ≈ 2.0

**Stauffer, D. & Aharony, A.** (1994). *Introduction to Percolation Theory*. 2nd ed., Taylor & Francis.
- 교과서, φ_c 개념, universality

---

## 7. Key Equations Summary

### Network Solver
```
R_total = R_bulk + R_constriction
R_bulk = d / (σ × π × r²)
R_constriction = 1 / (2σa)         [Holm 1967]
a = √(A_contact / π)
σ_eff = G_eff × L / A              [Kirchhoff]
```

### Ionic FORM X
```
σ_ionic = 0.1275 × σ_grain × (φ_SE - 0.185)^(3/4) × CN × √cov / √τ
σ_grain = 3.0 mS/cm [Kraft 2017]
φ_c = 0.185 [optimized, consistent with Bielefeld 2019]
R² = 0.944
```

### Binary Packing
```
x_S,optimal = ε₁ / (1 + ε₁ - ε₂) ≈ 0.27    [Furnas 1931]
d_L/d_S optimal ≈ 7:1                         [McGeary 1961]
void opening = 0.155 × d_L (tetrahedral)
ε_min = ε₁ × ε₂ ≈ 0.13                       [Furnas 1931]
```

### DEM Calibration
```
E_eff(SE) = 1.35 GPa
Calibrated to: SE pellet porosity ≈ 10% at 300 MPa
Reference: Doux et al. (2020), Coetzee (2017)
```

---

## 8. DEM Calibration Philosophy & Limitations

### 8.1 캘리브레이션 접근법

**Coetzee, C.J.** (2017). "Review: Calibration of the Discrete Element Method." *Powder Technology*, 310, 104-142.
- DEM 파라미터 캘리브레이션 방법론 종합 리뷰
- **핵심**: "DEM 파라미터는 실험 매칭을 통해 결정하며, effective parameter는 bulk behavior를 재현하는 것이 목적"
- E_eff를 실험 porosity에 매칭하는 접근이 표준 절차
- Micro-level 파라미터 → Macro-level 응답 매칭 (angle of repose, porosity, bulk density)
- **우리 접근의 근거**: E_eff=1.35 GPa는 SE pellet porosity ~10% (실험)에 캘리브레이션

**Coetzee, C.J.** (2016). "Calibration of the Discrete Element Method and the Effect of Particle Shape." *Powder Technology*, 297, 50-70.
- 입자 형상이 캘리브레이션에 미치는 영향
- 구형 가정의 한계와 보정 방법

### 8.2 Reduced Young's Modulus 접근

**Lommen, S., Schott, D., & Lodewijks, G.** (2014). "DEM speedup: Stiffness effects on behavior of bulk material." *Particuology*, 12, 107-112.
- E 감소 → timestep 증가 → 계산 효율 향상
- E를 1/100로 줄여도 bulk flow는 유사하게 유지
- **단, 접촉 면적은 영향 받음** → contact-sensitive 분석 시 주의

**Yan, Z., Wilkinson, S.K., Stitt, E.H., & Marigo, M.** (2015). "Discrete Element Modelling of Fluid Bed Granulation." *Journal of Computational Physics*, 302, 245-264.
- E_eff 사용 시 computational efficiency 대 accuracy trade-off
- Bulk behavior (porosity, flow) 보존, contact-level behavior (면적, 응력) 변화

### 8.3 SSB DEM 캘리브레이션 — 최신 연구

**Alabdali, M., Zanotto, F.M., et al.** (2022). "Contact model for DEM simulation of compaction and sintering of all-solid-state battery electrodes." *MethodsX*, 9, 101861.
- ASSB 전극용 DEM 접촉 모델 (탄성+소성+점탄성)
- Maxwell 점탄성 모델로 소결 효과 시뮬레이션
- **핵심 한계 인정**: "spherical overlap approximation에 의한 접촉 면적은 실제 소성변형보다 과소평가"
- 소성변형 시 2차 접촉 발생 → 구형 overlap 모델이 감지 못함

**Parameter sensitivity analysis and calibration of a DEM model for optimizing ASSB cathode microstructures** (2025). *Electrochimica Acta*.
- DEM 파라미터 민감도 분석 (friction이 가장 영향 큼)
- **높은 AM loading에서 DEM 정확도 저하** → 우리 82:18 케이스 관련
- 낮은 CAM loading에서 높은 정확도, 높은 loading에서 discrepancy 증가
- 추가 개선 필요성 언급

**Shi, T., Tu, Q., Tian, Y., et al.** (2020). "High Active Material Loading in All-Solid-State Battery Electrode via Particle Size Optimization." *Advanced Energy Materials*, 10(1), 1902881.
- DEM으로 SSB cathode packing 시뮬레이션
- AM loading은 크기비에 의해 결정됨
- **Percolation 문제**: 높은 cathode loading에서 separator 근처만 활성화
- 크기비 크게 → cathode utilization 향상

### 8.4 DEM 접촉 면적 한계 — 물리적 설명

```
DEM 접촉:     a = √(R* × δ)        → 구-구 Hertz 접촉, 원형 patch
실제 접촉:     소성유동 + 소결       → AM 곡면에 conform, 넓은 비정형 면적

결과:
  - Porosity:  E_eff로 캘리브레이션 → 실험과 일치 ✓
  - Coverage:  구-구 접촉 면적 공식 → 실제보다 과소평가 △
  - CN:        접촉 유무(topology) → E_eff 무관, 정확 ✓
  - σ_ionic:   Network solver = f(접촉 면적) → 상대적 경향 정확, 절대값 과소 △
```

**캘리브레이션 가능 vs 불가능:**

| 지표 | 캘리브레이션 | 이유 |
|------|:-:|------|
| Porosity | ✓ | 부피 지표 — E_eff로 전체 변형량 매칭 |
| CN, Percolation | ✓ | 위상 지표 — packing geometry에 의존 |
| Tortuosity | ✓ | 경로 지표 — 위상에 의존 |
| Coverage | △ | 접촉 면적 — 구-구 공식 한계 |
| σ_ionic (절대값) | △ | Coverage 의존 → 상대적 경향만 신뢰 |
| Contact Pressure | △ | 접촉 면적 과소 → 압력 과대 |

**논문 표현:**
> "The DEM model with calibrated effective Young's modulus (E_eff = 1.35 GPa) reproduces experimental porosity (~10% for SE-only pellet at 300 MPa, Doux et al. 2020). While the spherical overlap approximation inherently underestimates contact areas compared to actual plastic deformation and partial sintering of sulfide SE (Alabdali et al. 2022), the relative trends in ionic conductivity, coverage, and coordination number across different compositions remain physically valid. This calibration approach follows the established DEM methodology where effective parameters capture bulk-level response rather than individual contact mechanics (Coetzee 2017)."

---

*Database compiled for DEM-Based ASSB Composite Cathode Analysis*
*Last updated: 2026-04-14*
