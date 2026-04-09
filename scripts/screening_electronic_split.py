"""
Electronic: σ=0 제거 + thick/thin 분리 스크리닝
===============================================
1. σ=0 빼고 전체 회귀
2. thick only (T/d > 10)
3. thin only (T/d < 10)
4. 어디서 개형이 안 맞는지 per-case 분석
5. 모든 parameter 총동원
"""
import json,os,numpy as np,warnings
from pathlib import Path
from itertools import combinations
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SAM=50.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sel=m.get('electronic_sigma_full_mScm',0)
            if sel is None: sel = 0
            # σ=0도 포함! (percolation 실패 케이스)
            pa=max(m.get('phi_am',0),0.01);ps=m.get('phi_se',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01);cn=m.get('se_se_cn',0)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            el_act=max(m.get('electronic_active_fraction',0),0.01)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'por':por,'el_act':el_act,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['T']:.1f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p);return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r);return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None

def screen(sel, var_dict, label, threshold=0.85):
    """Power law screening — σ>0 cases only for log-space R²."""
    pos = sel > 0
    n_pos = pos.sum()
    if n_pos < 5:
        print(f"  Not enough positive cases ({n_pos})")
        return []
    log_sel = np.log(sel[pos])
    var_pos = {k: v[pos] for k, v in var_dict.items()}
    n = n_pos
    results = []

    for nvars in [3, 4]:
        for vnames in combinations(var_pos.keys(), nvars):
            X = np.column_stack([var_pos[v] for v in vnames] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sel, rcond=None)
                pred = np.exp(X @ coefs); r2 = r2l(sel[pos], pred)
                if r2 > threshold:
                    terms = ' × '.join(f'{v}^{coefs[i]:.2f}' for i, v in enumerate(vnames))
                    results.append({'name': f'{terms}', 'r2': r2, 'C': np.exp(coefs[-1]), 'nvars': nvars, 'coefs': coefs, 'vnames': vnames})
            except:
                pass

    results.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 10 ({label}, n={n}):")
    for i, r in enumerate(results[:10]):
        print(f"    #{i+1} R²={r['r2']:.4f}  {r['name']}")
    return results

def main():
    rows = load_data(); n = len(rows)
    sel = np.array([r['sel'] for r in rows])
    pa = np.array([r['pa'] for r in rows]); ps = np.array([r['ps'] for r in rows])
    am_cn = np.array([r['am_cn'] for r in rows]); cn = np.array([r['cn'] for r in rows])
    T = np.array([r['T'] for r in rows]); d_am = np.array([r['d_am'] for r in rows])
    tau = np.array([r['tau'] for r in rows]); cov = np.array([r['cov'] for r in rows])
    gd = np.array([r['gd'] for r in rows]); gp = np.array([r['gp'] for r in rows])
    ratio = T / d_am; el_act = np.array([r['el_act'] for r in rows])
    por = np.array([r['por'] for r in rows])
    names = [r['name'] for r in rows]

    thick = ratio >= 10
    thin = ratio < 10

    zero_mask = sel == 0
    print(f"n={n} (σ>0: {(~zero_mask).sum()}, σ=0: {zero_mask.sum()})")
    print(f"thick (T/d≥10): {thick.sum()}")
    print(f"thin (T/d<10): {thin.sum()}")

    # Variable pool
    var_pool = {
        'φ_AM': np.log(pa), 'φ_SE': np.log(ps), 'AM_CN': np.log(am_cn),
        'SE_CN': np.log(cn), 'τ': np.log(tau), 'T': np.log(T),
        'd_AM': np.log(d_am), 'T/d': np.log(ratio), 'cov': np.log(cov),
        'GB_d': np.log(gd), 'el_act': np.log(el_act),
        'por': np.log(np.clip(por, 0.1, None)),
        'φ_AM×φ_SE': np.log(pa * ps),
    }

    # ═══════════════════════════════════════
    # PART 1: v1 per-regime R²
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 1: v1 per-regime")
    print("="*80)
    rhs_v1 = SAM * pa**1.5 * am_cn**2 * np.exp(np.pi / ratio)
    C_v1 = fitC(sel, rhs_v1); pred_v1 = C_v1 * rhs_v1
    print(f"  ALL:   R²={r2l(sel, pred_v1):.4f}")
    print(f"  thick: R²={r2l(sel[thick], pred_v1[thick]):.4f}  (n={thick.sum()})")
    print(f"  thin:  R²={r2l(sel[thin], pred_v1[thin]):.4f}  (n={thin.sum()})")

    # Per-case thin
    print(f"\n  Thin cases (T/d<10):")
    for i in range(n):
        if thin[i]:
            err = abs(sel[i] - pred_v1[i]) / sel[i] * 100
            print(f"    {names[i]:35s} σ_act={sel[i]:.2f} σ_pred={pred_v1[i]:.2f} err={err:.0f}% T/d={ratio[i]:.1f} φ_AM={pa[i]:.3f}")

    # ═══════════════════════════════════════
    # PART 2: ALL data screening (σ>0 only)
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: ALL data (σ>0)")
    print("="*80)
    screen(sel, var_pool, "ALL")

    # ═══════════════════════════════════════
    # PART 3: THICK only
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: THICK only (T/d≥10)")
    print("="*80)
    var_thick = {k: v[thick] for k, v in var_pool.items()}
    screen(sel[thick], var_thick, "thick")

    # ═══════════════════════════════════════
    # PART 4: THIN only
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: THIN only (T/d<10)")
    print("="*80)
    var_thin = {k: v[thin] for k, v in var_pool.items()}
    screen(sel[thin], var_thin, "thin", threshold=0.70)

    # ═══════════════════════════════════════
    # PART 5: Correlation diff thick vs thin
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 5: Correlation — thick vs thin")
    print("="*80)
    log_sel_all = np.log(sel)
    print(f"  {'Variable':>15s} {'corr_ALL':>10s} {'corr_thick':>10s} {'corr_thin':>10s} {'diff':>10s}")
    for vn, vv in sorted(var_pool.items()):
        c_all = np.corrcoef(vv, log_sel_all)[0, 1]
        c_thick = np.corrcoef(vv[thick], np.log(sel[thick]))[0, 1] if thick.sum() > 3 else 0
        c_thin = np.corrcoef(vv[thin], np.log(sel[thin]))[0, 1] if thin.sum() > 3 else 0
        diff = c_thin - c_thick
        flag = '!!!' if abs(diff) > 0.3 else ''
        print(f"  {vn:>15s} {c_all:+10.3f} {c_thick:+10.3f} {c_thin:+10.3f} {diff:+10.3f} {flag}")

if __name__ == '__main__':
    main()
