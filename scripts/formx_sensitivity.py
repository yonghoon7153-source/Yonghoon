"""
FORM X sensitivity: α와 φ_c를 유동적으로
=============================================
1. α, φ_c 동시 최적화 (3 free params: C, α, φ_c)
2. φ_c가 SE 크기에 따라 다른가?
3. α가 데이터 범위에 따라 바뀌는가?
4. 2D heatmap: α vs φ_c → R²
"""
import json,os,numpy as np,warnings
from pathlib import Path
from scipy.optimize import minimize
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SG=3.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sn=m.get('sigma_full_mScm',0)
            if not sn or sn<0.001: continue
            ps=m.get('phi_se',0)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cn=m.get('se_se_cn',0);gd=m.get('gb_density_mean',0)
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            if t<=0 or ps<=0 or cn<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'tau':t,'cn':cn,'gd':gd,'T':T,'cov':cov,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p);return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r);return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None

def main():
    rows=load_data();n=len(rows);print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows])
    tau=np.array([r['tau'] for r in rows])
    cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows])
    cov=np.array([r['cov'] for r in rows])
    thick=np.array([r['tau']<=1.5 for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])

    # ═══════════════════════════════════════
    # 1. α, φ_c 동시 최적화
    # ═══════════════════════════════════════
    print("="*80)
    print("1. α, φ_c 동시 최적화 (joint optimization)")
    print("="*80)

    def neg_r2(params):
        alpha, phic = params
        if phic <= 0.05 or phic >= 0.25 or alpha <= 0.1 or alpha >= 1.0:
            return 0
        phi_ex = np.clip(ps - phic, 0.001, None)
        inner = phi_ex**1.5 * cn**2 * cov / tau
        rhs = SG * inner**alpha
        C = fitC(sn, rhs)
        if C is None: return 0
        return -r2l(sn, C*rhs)

    # Grid search first
    best_r2 = -999; best_params = (0.5, 0.18)
    print(f"\n  2D Grid (α × φ_c):")
    print(f"  {'φ_c↓ / α→':>10s}", end='')
    for a in np.arange(0.40, 0.66, 0.05):
        print(f"  {a:5.2f}", end='')
    print()

    for phic in np.arange(0.12, 0.23, 0.01):
        print(f"  {phic:10.2f}", end='')
        for a in np.arange(0.40, 0.66, 0.05):
            phi_ex = np.clip(ps-phic, 0.001, None)
            inner = phi_ex**1.5 * cn**2 * cov / tau
            rhs = SG * inner**a
            C = fitC(sn, rhs)
            if C is None: print(f"  {'--':>5s}", end=''); continue
            r2 = r2l(sn, C*rhs)
            flag = '★' if r2 > best_r2 else ' '
            if r2 > best_r2:
                best_r2 = r2; best_params = (a, phic)
            print(f" {flag}{r2:.3f}", end='')
        print()

    # Scipy optimize
    res = minimize(neg_r2, x0=[0.5, 0.18], method='Nelder-Mead',
                   options={'xatol':0.001,'fatol':0.0001})
    opt_a, opt_phic = res.x
    print(f"\n  Grid best: α={best_params[0]:.2f}, φ_c={best_params[1]:.2f}, R²={best_r2:.4f}")
    print(f"  Scipy opt: α={opt_a:.4f}, φ_c={opt_phic:.4f}, R²={-res.fun:.4f}")
    print(f"  FORM X:    α=0.5000, φ_c=0.1800, R²=", end='')

    phi_ex = np.clip(ps-0.18, 0.001, None)
    rhs_X = SG * (phi_ex**1.5 * cn**2 * cov / tau)**0.5
    C_X = fitC(sn, rhs_X); print(f"{r2l(sn, C_X*rhs_X):.4f}")

    print(f"\n  Gap: α={abs(opt_a-0.5):.4f}, φ_c={abs(opt_phic-0.18):.4f}")

    # ═══════════════════════════════════════
    # 2. φ_c가 SE 크기에 따라 다른가?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ_c per SE size group")
    print("="*80)

    groups = {
        'SE~0.5μm (GB_d>1.0)': gd > 1.0,
        'SE~1.0μm (0.6<GB_d≤1.0)': (gd > 0.6) & (gd <= 1.0),
        'SE~1.5μm (GB_d≤0.6)': gd <= 0.6,
    }

    for gname, gmask in groups.items():
        if gmask.sum() < 5: continue
        sn_g = sn[gmask]; ps_g = ps[gmask]; tau_g = tau[gmask]
        cn_g = cn[gmask]; cov_g = cov[gmask]

        best_phic_g = 0.18; best_r2_g = -999
        for phic in np.arange(0.05, 0.25, 0.005):
            phi_ex_g = np.clip(ps_g - phic, 0.001, None)
            rhs_g = SG * (phi_ex_g**1.5 * cn_g**2 * cov_g / tau_g)**0.5
            C_g = fitC(sn_g, rhs_g)
            if C_g is None: continue
            r2_g = r2l(sn_g, C_g*rhs_g)
            if r2_g > best_r2_g:
                best_r2_g = r2_g; best_phic_g = phic

        print(f"  {gname}: n={gmask.sum()}, best φ_c={best_phic_g:.3f}, R²={best_r2_g:.4f}")

    # ═══════════════════════════════════════
    # 3. FORM X with free α, φ_c → R² gain?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. Free params comparison")
    print("="*80)

    # 1 free param (C only): FORM X
    print(f"  1 free (C): α=0.50, φ_c=0.18  → R²={r2l(sn, C_X*rhs_X):.4f}")

    # 2 free (C + φ_c): α=0.5 fixed
    best_phic2 = 0.18; best_r2_2 = -999
    for phic in np.arange(0.10, 0.25, 0.005):
        phi_ex2 = np.clip(ps-phic, 0.001, None)
        rhs2 = SG * (phi_ex2**1.5 * cn**2 * cov / tau)**0.5
        C2 = fitC(sn, rhs2)
        if C2 is None: continue
        r2_2 = r2l(sn, C2*rhs2)
        if r2_2 > best_r2_2: best_r2_2 = r2_2; best_phic2 = phic
    print(f"  2 free (C+φ_c): α=0.50, φ_c={best_phic2:.3f} → R²={best_r2_2:.4f}")

    # 2 free (C + α): φ_c=0.18 fixed
    best_a3 = 0.5; best_r2_3 = -999
    for a in np.arange(0.35, 0.70, 0.005):
        phi_ex3 = np.clip(ps-0.18, 0.001, None)
        rhs3 = SG * (phi_ex3**1.5 * cn**2 * cov / tau)**a
        C3 = fitC(sn, rhs3)
        if C3 is None: continue
        r2_3 = r2l(sn, C3*rhs3)
        if r2_3 > best_r2_3: best_r2_3 = r2_3; best_a3 = a
    print(f"  2 free (C+α): α={best_a3:.3f}, φ_c=0.18 → R²={best_r2_3:.4f}")

    # 3 free (C + α + φ_c)
    print(f"  3 free (C+α+φ_c): α={opt_a:.3f}, φ_c={opt_phic:.3f} → R²={-res.fun:.4f}")

    print(f"\n  ΔR² for each extra free param:")
    r2_1 = r2l(sn, C_X*rhs_X)
    print(f"    +φ_c free: ΔR²={best_r2_2 - r2_1:+.4f}")
    print(f"    +α free:   ΔR²={best_r2_3 - r2_1:+.4f}")
    print(f"    +both:     ΔR²={-res.fun - r2_1:+.4f}")

    # ═══════════════════════════════════════
    # 4. 결론
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("CONCLUSION")
    print("="*80)
    gain = -res.fun - r2_1
    print(f"""
  FORM X (fixed α=0.5, φ_c=0.18): R²={r2_1:.4f}
  Full optimization (free α, φ_c):  R²={-res.fun:.4f}

  Total gain from freeing α and φ_c: ΔR²={gain:+.4f}

  {'→ FORM X의 고정 지수는 충분히 최적에 가깝다 (ΔR²<0.005)' if gain < 0.005 else '→ 유동적 α/φ_c가 유의미한 개선을 준다'}
  {'→ α=0.5 (⁴√)과 φ_c=0.18은 물리적으로 깔끔하고 R² 손실이 미미' if gain < 0.005 else ''}
    """)

if __name__=='__main__':
    main()
