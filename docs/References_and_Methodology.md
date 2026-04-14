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

*Database compiled for DEM-Based ASSB Composite Cathode Analysis*
*Last updated: 2026-04-14*
