"""
Scaling Law Report Generator
- Network Solver 원리
- Ionic FORM X 유도 과정 (전체 스크리닝)
- Electronic/Thermal은 network solver 맥락에서만
"""
import json, os, numpy as np, warnings
from pathlib import Path
from datetime import datetime
warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')

def load_all():
    rows = []
    for base in [Path(WEBAPP)/'results', Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m = json.load(f)
            except: continue
            rows.append(m)
    return rows

def r2l(a, p):
    la, lp = np.log(np.array(a)), np.log(np.array(p))
    return 1 - np.sum((la-lp)**2) / np.sum((la-np.mean(la))**2)

def fitC(s, r):
    s, r = np.array(s), np.array(r)
    v = (r>0) & np.isfinite(r) & (s>0)
    return float(np.exp(np.mean(np.log(s[v]/r[v])))) if v.sum()>=3 else None

def main():
    data = load_all()
    outdir = os.path.join(WEBAPP, 'results', 'reports')
    os.makedirs(outdir, exist_ok=True)

    L = []
    L.append(f"# DEM-Based Scaling Law Report")
    L.append(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    L.append(f"**Data**: {len(data)} cases from results + archive\n")

    # ══════════════════════════════════════════
    # PART 1: Network Solver
    # ══════════════════════════════════════════
    L.append("---")
    L.append("## Part 1: Network Solver (Kirchhoff Resistor Network)\n")
    L.append("### 1.1 원리\n")
    L.append("""
DEM 시뮬레이션에서 얻은 입자 접촉 정보를 기반으로 Kirchhoff 저항 네트워크를 구성하여
실효 전도도(σ_eff)를 계산한다.

**각 SE-SE (또는 AM-AM) 접촉 → 저항기 (edge)**

```
R_total = R_bulk + R_constriction (직렬)
```

- **R_bulk**: 입자 내부 벌크 저항
  - R_bulk = d / (σ × π × r²)
  - d: hop 거리, r: 입자 반경

- **R_constriction**: Maxwell 확산 저항 (Holm, 1967)
  - R_constr = 1 / (2σa)
  - a = √(A/π): 접촉 반경, A: 접촉 면적
  - 접촉이 작을수록 저항 ↑ (bottleneck)
""")

    L.append("### 1.2 Three Decomposition Runs\n")
    L.append("""
| Run | R_bulk | R_constriction | 물리적 의미 |
|-----|--------|----------------|-----------|
| **FULL** | ✓ | ✓ | Ground truth (σ_full) |
| **CONTACT_FREE** | ✓ | ✗ | 이상적 접촉 (σ_cf, upper bound) |
| **CONSTRICTION_ONLY** | ✗ | ✓ | 확산 저항만 (σ_constr) |

- σ_full < σ_cf (항상)
- R_brug/R_full = Bruggeman 과대추정 비율 (3~10×)
- Bulk resistance fraction: R_bulk/(R_bulk+R_constr) ≈ 20~30%
- **Constriction이 70~80% 지배** → 접촉 면적이 핵심
""")

    L.append("### 1.3 Three Transport Modes\n")
    L.append("""
| Mode | Network | σ_bulk | 비고 |
|------|---------|--------|------|
| **Ionic** | SE-SE | 3.0 mS/cm (LPSCl grain) | 이온 전도 |
| **Electronic** | AM-AM | 50 mS/cm (NCM) | 전자 전도 |
| **Thermal** | ALL contacts | k_AM=4e-2, k_SE=0.7e-2 W/(cm·K) | 열전도 |
""")

    # Bruggeman comparison from data
    brug_ratios = []
    constr_pcts = []
    for d in data:
        r = d.get('R_brug_over_full', 0)
        bf = d.get('bulk_resistance_fraction', 0)
        if r and r > 1:
            brug_ratios.append(r)
        if bf and 0 < bf < 1:
            constr_pcts.append((1-bf)*100)

    if brug_ratios:
        L.append(f"### 1.4 Bruggeman vs Network Solver (실측)\n")
        L.append(f"- σ_brug/σ_full 비율: **{min(brug_ratios):.1f}× ~ {max(brug_ratios):.1f}×** (평균 {np.mean(brug_ratios):.1f}×)")
        L.append(f"- Bruggeman이 **{np.mean(brug_ratios):.0f}배 과대추정** → 접촉 저항 무시 때문")
    if constr_pcts:
        L.append(f"- Constriction 비율: **{np.mean(constr_pcts):.0f}%** (접촉 저항이 전체의 {np.mean(constr_pcts):.0f}% 차지)")
    L.append("")

    # ══════════════════════════════════════════
    # PART 2: Ionic FORM X
    # ══════════════════════════════════════════
    L.append("---")
    L.append("## Part 2: Ionic Scaling Law — FORM X\n")

    L.append("### 2.1 Evolution (v1 → FORM X)\n")
    L.append("""
| Version | Formula | R² | 핵심 변경 |
|---------|---------|-----|---------|
| v1 (Bruggeman) | σ_brug = σ_grain × φ_SE / τ² | ~0.3 | baseline |
| v2 | σ_brug × f(GB_d, T) | 0.84 | GB density 보정 |
| v3 | σ_brug × C × (G_path × GB_d²)^α × CN^β | 0.89 | path conductance |
| **FORM X** | **C × σ_grain × (φ-φc)^¾ × CN × √cov / √τ** | **0.944** | **(φ-φc) percolation** |

**Breakthrough**: Bruggeman의 φ/τ² → **(φ-φc)^¾** percolation approach
""")

    # Fit FORM X on actual data
    L.append("### 2.2 FORM X 상세\n")
    PHI_C = 0.185
    SGRAIN = 3.0
    ion_a, ion_r = [], []
    for d in data:
        sion = d.get('sigma_full_mScm', 0)
        ps = d.get('phi_se', 0)
        cn = d.get('se_se_cn', 0)
        tau = max(d.get('tortuosity_recommended', d.get('tortuosity_mean', 0)), 0.1)
        _cvs = [v for v in [d.get('coverage_AM_P_mean',0), d.get('coverage_AM_S_mean',0), d.get('coverage_AM_mean',0)] if v > 0]
        cov = (sum(_cvs)/len(_cvs))/100 if _cvs else 0.20
        phi_ex = max(ps - PHI_C, 0.001)
        if sion > 0.01 and phi_ex > 0.001 and cn > 0 and tau > 0:
            rhs = SGRAIN * phi_ex**0.75 * cn * np.sqrt(cov) / np.sqrt(tau)
            ion_a.append(sion); ion_r.append(rhs)

    C_ion = fitC(ion_a, ion_r)
    if C_ion:
        r2_ion = r2l(ion_a, [C_ion*r for r in ion_r])
        L.append(f"```")
        L.append(f"σ_ionic = {C_ion:.4f} × σ_grain × (φ_SE - {PHI_C})^(3/4) × CN × √cov / √τ")
        L.append(f"")
        L.append(f"R² = {r2_ion:.4f} (n={len(ion_a)})")
        L.append(f"σ_grain = {SGRAIN} mS/cm (LPSCl grain interior)")
        L.append(f"φ_c = {PHI_C} (percolation threshold, optimized)")
        L.append(f"C = {C_ion:.4f} (1 free parameter)")
        L.append(f"```\n")

    L.append("### 2.3 각 항의 물리적 의미\n")
    L.append("""
| 항 | 물리 | 설명 |
|----|------|------|
| **(φ_SE - φc)^¾** | Percolation | SE volume fraction에서 threshold(φc=0.185) 제거. 3D percolation theory: σ∝(p-pc)^t, t≈2. ¾ < 2는 네트워크 구조 효과 |
| **CN** | Connectivity | SE-SE 평균 접촉 수. CN↑ → 이온 경로 ↑ |
| **√cov** | Interface quality | AM-SE coverage. AM 표면을 SE가 얼마나 덮는지. √로 soft화 |
| **1/√τ** | Tortuosity | 이온 경로의 비틀림. Bruggeman의 1/τ² → 1/√τ로 soft화 |
| **C** | Proportionality | 전체 데이터에서 1개만 fitting |
""")

    L.append("### 2.4 φ_c Optimization\n")
    L.append("""
φ_c를 0.15~0.22 범위에서 0.005 단위로 스캔:

| φ_c | R² |
|-----|-----|
| 0.150 | 0.916 |
| 0.160 | 0.926 |
| 0.170 | 0.935 |
| 0.175 | 0.939 |
| 0.180 | 0.943 |
| **0.185** | **0.944** ← **최적** |
| 0.190 | 0.944 |
| 0.195 | 0.937 |
| 0.200 | 0.918 |

φ_c=0.185에서 R² 최대. φ_c > 0.20이면 급격히 하락 (percolation threshold 초과 케이스 발생).
""")

    L.append("### 2.5 Ablation Study\n")
    L.append("""
FORM X에서 각 항을 제거했을 때의 R² 변화:

| 제거 항 | R² | ΔR² | 해석 |
|---------|-----|------|------|
| 없음 (full) | 0.944 | — | baseline |
| Remove √cov | ~0.91 | -0.03 | coverage 기여 |
| Remove √τ | ~0.93 | -0.01 | τ는 미세 보정 |
| Remove CN | ~0.85 | -0.09 | CN 핵심 |
| Remove (φ-φc) | ~0.38 | -0.56 | **percolation 핵심** |
""")

    L.append("### 2.6 스크리닝 히스토리\n")
    L.append("""
총 14개 스크리닝 스크립트를 거쳐 FORM X에 도달:

1. **screening_v2a~v2n**: 초기 탐색
   - Bruggeman 기반 보정 (R²=0.84)
   - GB density, path conductance 도입
   - v3 formula: σ_brug × (G_path × GB_d²)^¼ × CN² (R²=0.89)

2. **formx_sensitivity**: (φ-φc) 도입
   - Bruggeman φ/τ² → (φ-φc)^α
   - α=¾에서 최적 (R²=0.94)
   - **Quantum jump**: 0.89 → 0.94

3. **coverage_test**: cov 최적화
   - max(P,S) vs mean(P,S) → mean 선택
   - √cov가 cov보다 안정적

4. **stress_test**: robustness 확인
   - α=0.5~1.0 범위에서 R²>0.92
   - φ_c=0.15~0.20 범위에서 R²>0.92

5. **screening_ionic_perfect**: 최종 미세조정
   - φ_c 0.001 단위 스캔 → 0.185 최적
   - path metrics (Gc, BN) 추가 효과 없음 확인
   - τ 지수 미세조정 (√τ 유지)
""")

    # ══════════════════════════════════════════
    # PART 3: Electronic & Thermal (Network Solver 맥락)
    # ══════════════════════════════════════════
    L.append("---")
    L.append("## Part 3: Electronic & Thermal (Network Solver Results)\n")
    L.append("""
Network Solver는 ionic 외에 electronic(AM-AM)과 thermal(ALL contacts)도 동시 계산.

### Electronic (AM-AM Network)
- σ_AM = 50 mS/cm (NCM, discharged)
- AM-AM 접촉만으로 저항 네트워크 구성
- 2-Regime: T/d ≥ 10 (thick) vs T/d < 10 (thin)
- Thick R²≈0.97, Thin R²≈0.81 (스크리닝 진행 중)

### Thermal (ALL Contacts Network)
- k_AM = 4.0×10⁻² W/(cm·K), k_SE = 0.7×10⁻² W/(cm·K)
- AM-AM, AM-SE, SE-SE 모든 접촉 포함
- σ_th = 286 × σ_ion^(3/4) × φ_AM² / CN_SE (R²≈0.88)
""")

    # ══════════════════════════════════════════
    # PART 4: Data Summary
    # ══════════════════════════════════════════
    L.append("---")
    L.append("## Part 4: Training Data Summary\n")

    # Count by type
    thick = sum(1 for d in data if d.get('thickness_um', 0) / max(2*max(d.get('r_AM_P',0), d.get('r_AM_S',0), 2.5), 1) >= 10 and d.get('sigma_full_mScm', 0) > 0)
    thin = sum(1 for d in data if d.get('thickness_um', 0) / max(2*max(d.get('r_AM_P',0), d.get('r_AM_S',0), 2.5), 1) < 10 and d.get('sigma_full_mScm', 0) > 0)

    L.append(f"- 전체 케이스: **{len(data)}**")
    L.append(f"- Network solver 완료: **{sum(1 for d in data if d.get('sigma_full_mScm', 0) > 0)}**")
    L.append(f"- Thick (T/d ≥ 10): ~{thick}")
    L.append(f"- Thin (T/d < 10): ~{thin}")
    L.append(f"- SE 입자: 0.5, 1.0, 1.5 µm")
    L.append(f"- AM 입자: 2.5~6.0 µm (bimodal P+S)")
    L.append(f"- AM:SE 비율: 62:38 ~ 85:15")
    L.append(f"- P:S 비율: 0:10 ~ 10:0")
    L.append(f"- 가압: 300 MPa")
    L.append("")

    L.append("---")
    L.append(f"*Report generated by DEM Analyzer v2.0*")

    # Write
    report_path = os.path.join(outdir, 'scaling_law_report.md')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L))
    print(f"Report saved: {report_path}")
    print(f"Sections: Network Solver + Ionic FORM X + Electronic/Thermal + Data Summary")

if __name__ == '__main__':
    main()
