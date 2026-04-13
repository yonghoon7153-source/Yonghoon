# Adhesion v5 — Crystalline Slab xy-shift 방법론 상세 보고서

## 1. 배경 및 동기

### 왜 adhesion energy를 계산하는가?
ASSB에서 SE/NCM 계면 접착은 전지 성능과 수명을 직접 좌우한다.
충방전 시 NCM 부피 변화 → 계면 응력 → 약한 접착이면 박리 → 용량 열화.
SE 조성(Cl/Br, vacancy)이 계면 접착에 미치는 영향 정량화가 목표.

### 이전 방법 한계
| Version | 방법 | 문제 |
|---------|------|------|
| v1 | Random SE + melt | Seed 의존성 큼 |
| v2 | CIF + 3000K melt | Vacancy 파괴! Li5.4 불가 |
| v3 | 800K melt, 5x1x1 | Mismatch ±20% → PBC artifact |
| v4 | 1500K melt | Framework 파괴 + mismatch |

핵심 딜레마: 3000K melt → mismatch 해소 but vacancy 파괴.
v5 해결: supercell repeat로 mismatch 해소 → melt 불필요 → vacancy 보존!

## 2. Cell 설계

### Lattice Matching
- comp1/2B: SE 2x2=20.11A, NCM 7x7=20.15A → strain +0.2%
- comp3/4/5: SE 2x2=14.23A, NCM 5x5=14.39A → strain +1.1%

### SE 두께 매칭
- comp1/2B: cubic a=10A → 2x2x1=10A(너무 얇!) → 2x2x3=30A
- comp3/4/5: rhombo a3=29A → 2x2x1=29A (적당)
- 30A vs 29A ≈ 동일

### SE Orthogonalization
- Rhombo a3=[-3.558,-2.054,29.047] → z축에서 35도 기울어짐
- a3에서 a1,a2 성분 제거 → a3=[0,0,29.047] → z축 정렬
- _slab.xyz 파일로 저장

### 최종 Cell
| Family | SE repeat | SE at | SE thick | NCM | NCM at | A (A2) | strain |
|--------|-----------|-------|----------|-----|--------|--------|--------|
| Li6 | (2,2,3) | 624 | 30A | 7x7x1 | 196 | 351.5 | +0.2% |
| Li5.4 | (2,2,1) | 248 | 29A | 5x5x1 | 100 | 179.3 | +1.1% |

## 3. xy-shift Sampling

### z-shift 실패
z-shift=0.4 → SE가 z=0.6에서 잘림 → dangling bond → Wad=0.006~9.252 (비물리적)

### xy-shift 원리
```python
frac_se[:,0] = (frac_se[:,0] + dx) % 1.0
frac_se[:,1] = (frac_se[:,1] + dy) % 1.0
# z 안 건드림! SE slab 연속 유지!
```
PBC로 xy 주기적 → shift해도 slab 무결성 유지 → NCM-SE registry만 변경

비유: z-shift = 빵 자르기, xy-shift = 빵 밀기

## 4. Wad 계산

### E_int
1. NCM + SE stacking (gap=2.5A)
2. cell_z = SE_max + 30A (UMA vacuum sensitivity → 30A 최적!)
3. LBFGS relax (fmax=0.01, steps=200)
4. E_int = interface.get_potential_energy()

MQA 미사용 이유: 800K에서 Li 58/248 interdiffusion, 500K에서도 z_boundary shift.

### E_sep
1. cell_z += 30A (필수! 안 하면 PBC wrap)
2. SE atoms z += 30A
3. new_calc() (필수! 같은 calc → shape mismatch)
4. E_sep = single point (relax 안 함! strain 해소 → Wad 음수)

### Wad
Wad = (E_sep - E_int) / A * 16.0218 (eV/A2 → J/m2)

## 5. 실패한 대안들

| Method | Wad | Problem |
|--------|-----|---------|
| Isolated slab + relax | 10 | relax releases strain |
| Isolated slab, no relax | 75 | UMA 60A vacuum → E=+248eV |
| Sep 30A, cell fixed | 40 | PBC wrap |
| Independent slab cells | 17 | pre-MQA structure |
| MQA 800K | 62 | Li interdiffusion 58/248 atoms |
| MQA 500K + element sep | 51 | z_boundary shift, Li ambiguous |
| Sep + relax | -1.9 | strain release → negative Wad |
| Vacuum 60A | 24.5 | UMA out of training distribution |
| **Sep 30A + cell+30A** | **1.25** | **CORRECT** |

## 6. 최종 결과

| Comp | Family | Wad (J/m2) | std | Expt (aJ) |
|------|--------|-----------|-----|-----------|
| comp1 | Li6 | **1.433** | 0.288 | 194 |
| comp2B | Li6 | **1.244** | 0.356 | 180 |
| comp3 | Li5.4 | **2.361** | 0.41 | 316 |
| comp4 | Li5.4 | **2.202** | 0.33 | 298 |
| comp5 | Li5.4 | **2.037** | 0.44 | 249 |

Trends:
- Li6: comp1(1.43) > comp2B(1.24) → Br↑ Wad↓
- Li5.4: comp3(2.36) > comp4(2.20) > comp5(2.04) → Br↑ Wad↓
- Cross: Li5.4(2.2) > Li6(1.3) → vacancy enhances adhesion
- Expt order: comp3>comp4>comp5>comp1>comp2 = PERFECT MATCH!

### Vacancy 효과 반증
3000K melt comp3: Wad≈1.0 = comp1(1.1) → vacancy 파괴 → 차이 소멸!
→ vacancy가 adhesion 향상의 핵심 원인 증명!

## 7. 논문 서술

Methods: "Interfacial adhesion energies were computed using crystalline SE/NCM
slab models with UMA. Five configurations per composition via random
xy-translation of the SE slab. Wad = (E_sep - E_int)/A."

Results: "Wad decreases with Br in both families. Li5.4 >> Li6 due to
vacancy-mediated interfacial anchoring. Control: 3000K melt destroys
vacancy → Wad drops to Li6 level, confirming vacancy is essential."

## 8. 기술적 주의사항
1. UMA vacuum = 30A (60A에서 비정상)
2. Separation 시 cell도 확장 (cell_z += 30A)
3. Separation 후 relax 안 함 (single point only)
4. 새 calculator 매번 생성 (new_calc())
5. xy-shift만 사용 (z-shift → slab 파괴)
6. SE 두께 매칭 (Li6: 2x2x3, Li5.4: 2x2x1)
7. a3 orthogonalization 필수 (rhombo cell)
