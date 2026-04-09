"""
Electronic: exp 완전 제거. 처음부터 다시.
=============================================
exp(π/(T/d))에 얽매이지 않는다.
T와 d_AM을 독립 변수로 풀거나, 아예 다른 구조.
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
            if not sel or sel<0.001: continue
            pa=max(m.get('phi_am',0),0.01);ps=m.get('phi_se',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01);cn=m.get('se_se_cn',0)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            el_act=max(m.get('electronic_active_fraction',0),0.01)
            por=m.get('porosity',0)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'el_act':el_act,'por':por,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['T']:.1f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p);return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r);return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn,rhs):
    n=len(sn);la=np.log(sn);lr=np.log(rhs);errs=[]
    for i in range(n):
        m=np.ones(n,bool);m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def main():
    rows=load_data();n=len(rows);print(f"n={n}\n")
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    T=np.array([r['T'] for r in rows]);d_am=np.array([r['d_am'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    ratio=T/d_am;el_act=np.array([r['el_act'] for r in rows])
    log_sel=np.log(sel)

    # ═══════════════════════════════════════
    # PART 1: 완전 자유 — 모든 변수 조합 (power law only)
    # ═══════════════════════════════════════
    print("PART 1: Pure power law — NO exp")
    print("="*80)

    var_pool = {
        'φ_AM': np.log(pa), 'φ_SE': np.log(ps), 'AM_CN': np.log(am_cn),
        'SE_CN': np.log(cn), 'τ': np.log(tau), 'T': np.log(T),
        'd_AM': np.log(d_am), 'T/d': np.log(ratio), 'cov': np.log(cov),
        'GB_d': np.log(gd), 'el_act': np.log(el_act),
        'φ_AM×φ_SE': np.log(pa*ps), '1-φ_SE': np.log(1-ps),
        'φ_AM/φ_SE': np.log(pa/ps), 'AM_CN×φ_AM': np.log(am_cn*pa),
    }

    results = []
    for nvars in [3, 4, 5]:
        for vnames in combinations(var_pool.keys(), nvars):
            X = np.column_stack([var_pool[v] for v in vnames] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sel, rcond=None)
                pred = np.exp(X @ coefs); r2 = r2l(sel, pred)
                if r2 > 0.90:
                    terms = ' × '.join(f'{v}^{coefs[i]:.2f}' for i, v in enumerate(vnames))
                    results.append({
                        'name': f'PL({nvars}v): {terms}',
                        'r2': r2, 'C': np.exp(coefs[-1]),
                        'rhs': np.exp(X @ coefs) / np.exp(coefs[-1]),
                        'nvars': nvars, 'coefs': coefs, 'vnames': vnames,
                    })
            except:
                pass

    # ═══════════════════════════════════════
    # PART 2: Percolation threshold on φ_AM (no exp)
    # (φ_AM - φ_c)^a × AM_CN^b × ...
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: (φ_AM - φ_c)^a NO exp")
    print("="*80)
    for phic in [0, 0.15, 0.20, 0.25, 0.30, 0.35]:
        pa_ex = np.clip(pa - phic, 0.001, None) if phic > 0 else pa
        for a in [1, 1.5, 2, 2.5, 3, 4, 5]:
            for b in [0.5, 1, 1.5, 2]:
                for c in [-1, -0.5, 0, 0.5, 1]:  # T/d exponent (replaces exp)
                    for d in [0, 0.25, 0.5]:  # cov
                        rhs = SAM * pa_ex**a * am_cn**b * ratio**c * cov**d
                        C = fitC(sel, rhs)
                        if C is None or C <= 0: continue
                        pred = C * rhs; r2 = r2l(sel, pred)
                        if r2 > 0.88:
                            results.append({
                                'name': f'PERC_PL: (φ-{phic})^{a}×CN^{b}×(T/d)^{c}×cov^{d}',
                                'r2': r2, 'C': C, 'rhs': rhs,
                            })

    # ═══════════════════════════════════════
    # PART 3: FORM X mirror — ⁴√[...] (no exp)
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: ⁴√[...] like FORM X (no exp)")
    print("="*80)
    for phic in [0, 0.25, 0.30]:
        pa_ex = np.clip(pa - phic, 0.001, None) if phic > 0 else pa
        for a in [3, 4, 5, 6, 8]:
            for b in [3, 4, 5, 6]:
                for c in [0, 1, 2]:  # cov
                    for d in [0, 1, 2, -1, -2]:  # T/d
                        inner = pa_ex**a * am_cn**b * cov**c * ratio**d
                        if np.any(inner <= 0): continue
                        rhs = SAM * inner**0.25
                        C = fitC(sel, rhs)
                        if C is None or C <= 0: continue
                        pred = C * rhs; r2 = r2l(sel, pred)
                        if r2 > 0.88:
                            results.append({
                                'name': f'4RT: ⁴√[(φ-{phic})^{a}×CN^{b}×cov^{c}×(T/d)^{d}]',
                                'r2': r2, 'C': C, 'rhs': rhs,
                            })

    # ═══════════════════════════════════════
    # PART 4: sigmoid / saturating forms
    # 1/(1 + (d/T)^k) instead of exp
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: Sigmoid forms")
    print("="*80)
    for k in [1, 2, 3]:
        sig = 1 / (1 + (d_am / T)**k)
        for a in [2, 2.5, 3]:
            for b in [1.5, 2]:
                rhs = SAM * pa**a * am_cn**b * sig
                C = fitC(sel, rhs)
                if C is None: continue
                pred = C * rhs; r2 = r2l(sel, pred)
                if r2 > 0.88:
                    results.append({'name': f'SIG: φ^{a}×CN^{b}×1/(1+(d/T)^{k})', 'r2': r2, 'C': C, 'rhs': rhs})
                for d in [0.25, 0.5]:
                    rhs2 = rhs * cov**d
                    C2 = fitC(sel, rhs2)
                    if C2 is None: continue
                    pred2 = C2 * rhs2; r2_2 = r2l(sel, pred2)
                    if r2_2 > 0.88:
                        results.append({'name': f'SIG: φ^{a}×CN^{b}×sig({k})×cov^{d}', 'r2': r2_2, 'C': C2, 'rhs': rhs2})

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x: -x['r2'])

    # With exp comparison
    rhs_exp = SAM * pa**2.5 * am_cn**2 * np.exp(np.pi / ratio) * ps**0.5 * cov**0.25
    C_exp = fitC(sel, rhs_exp); r2_exp = r2l(sel, C_exp * rhs_exp)

    print(f"\n{'='*90}")
    print(f"TOP 30 NO-EXP FORMULAS (n={n})")
    print(f"WITH EXP BEST: R²={r2_exp:.4f} for comparison")
    print(f"{'='*90}")

    seen = set(); cnt = 0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']); cnt += 1
        if cnt > 30: break
        cv = loocv_C(sel, r['rhs']) if cnt <= 5 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        beats = "★★" if r['r2'] > r2_exp else "★" if r['r2'] > 0.91 else " "
        print(f"  #{cnt:2d}{beats} R²={r['r2']:.4f}{cv_str}  {r['name']}")

    # Head to head
    print(f"\n{'='*90}")
    print("HEAD TO HEAD: exp vs no-exp vs sigmoid")
    print(f"{'='*90}")
    best_noexp = results[0] if results else None
    print(f"  WITH exp:    R²={r2_exp:.4f}  φ^2.5×CN²×exp(π/(T/d))×√φ_SE×⁴√cov")
    if best_noexp:
        cv_best = loocv_C(sel, best_noexp['rhs'])
        print(f"  NO exp best: R²={best_noexp['r2']:.4f} LOOCV={cv_best:.4f}  {best_noexp['name']}")

if __name__ == '__main__':
    main()
