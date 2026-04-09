"""
Electronic: 비선형 형태 테스트 (exp, sigmoid, threshold)
=====================================================
T/d가 power law가 아닐 수 있다. exp(-k/(T/d))? sigmoid?
φ_AM도 threshold가 있을 수 있다.
"""
import json,os,numpy as np,warnings
from pathlib import Path
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
            am_cn=max(m.get('am_am_cn',0.01),0.01)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,
                'ratio':T/d_am,'por':max(por,0.1),'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['ratio']:.1f}"
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
    rows=load_data()
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    por=np.array([r['por'] for r in rows])

    thick=ratio>=10; thin=ratio<10
    print(f"n={len(rows)}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THIN: T/d 비선형 형태
    # ═══════════════════════════════════════
    sn=sel[thin];cn_n=am_cn[thin];rn=ratio[thin];porn=por[thin];cvn=cov[thin];psn=ps[thin]

    print("="*80)
    print("THIN: T/d 형태 비교 (CN^1.25 × por^3 × cov^1.25 × φ_SE^1.25 × f(T/d))")
    print("="*80)

    base = SAM * cn_n**1.25 * porn**3 * cvn**1.25 * psn**1.25

    forms = [
        # Power law
        ('(T/d)^(-0.25)', rn**(-0.25)),
        ('(T/d)^(-0.5)', rn**(-0.5)),
        ('(T/d)^(-0.75)', rn**(-0.75)),
        ('(T/d)^(-1)', rn**(-1)),
        # Exponential
        ('exp(-1/(T/d))', np.exp(-1/rn)),
        ('exp(-2/(T/d))', np.exp(-2/rn)),
        ('exp(-π/(T/d))', np.exp(-np.pi/rn)),
        ('exp(+π/(T/d))', np.exp(np.pi/rn)),
        ('exp(+1/(T/d))', np.exp(1/rn)),
        ('exp(-0.5/(T/d))', np.exp(-0.5/rn)),
        # Sigmoid
        ('1/(1+(2/(T/d))²)', 1/(1+(2/rn)**2)),
        ('1/(1+(3/(T/d))²)', 1/(1+(3/rn)**2)),
        ('1/(1+(5/(T/d)))', 1/(1+(5/rn))),
        ('(T/d)/(T/d+3)', rn/(rn+3)),
        ('(T/d)/(T/d+5)', rn/(rn+5)),
        ('tanh(T/d)', np.tanh(rn)),
        ('tanh(T/(2d))', np.tanh(rn/2)),
        # Log
        ('ln(T/d)', np.log(rn)),
        ('ln(1+T/d)', np.log(1+rn)),
    ]

    print(f"\n  {'Form':30s} {'R²':>8s} {'LOOCV':>8s}")
    print("  "+"-"*50)
    for label, f_td in sorted(forms, key=lambda x: -r2l(sn, fitC(sn, base*x[1])*(base*x[1]) if fitC(sn,base*x[1]) else np.ones(len(sn)))):
        rhs = base * f_td
        C = fitC(sn, rhs)
        if C is None: continue
        pred = C * rhs; r2 = r2l(sn, pred)
        cv = loocv_C(sn, rhs) if r2 > 0.90 else 0
        cv_str = f"{cv:8.4f}" if cv else "       -"
        flag = '★' if r2 > 0.92 else ' '
        print(f"  {flag}{label:29s} {r2:8.4f} {cv_str}")

    # ═══════════════════════════════════════
    # THIN: 전체 비선형 탐색
    # power + exp/sigmoid 혼합
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN: power × f(T/d) 조합 (best f(T/d) × CN^a × por^b × cov^c × φ_SE^d)")
    print("="*80)

    results = []
    for f_label, f_td in [('exp(π/(T/d))', np.exp(np.pi/rn)),
                           ('(T/d)^(-0.5)', rn**(-0.5)),
                           ('(T/d)/(T/d+5)', rn/(rn+5)),
                           ('tanh(T/(2d))', np.tanh(rn/2)),
                           ('ln(T/d)', np.log(rn))]:
        for a in [1, 1.25, 1.5]:
            for b in [2, 2.5, 3]:
                for c in [0.5, 1, 1.25]:
                    for d in [0.5, 1, 1.25]:
                        rhs = SAM * cn_n**a * porn**b * cvn**c * psn**d * f_td
                        C = fitC(sn, rhs)
                        if C is None or C <= 0: continue
                        pred = C*rhs; r2 = r2l(sn, pred)
                        if r2 > 0.91:
                            results.append({'r2': r2, 'C': C, 'rhs': rhs,
                                'name': f'CN^{a}×por^{b}×cov^{c}×φ_SE^{d}×{f_label}'})

    results.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 15:")
    seen = set()
    for i, r in enumerate(results[:15]):
        if r['name'] in seen: continue
        seen.add(r['name'])
        cv = loocv_C(sn, r['rhs']) if i < 3 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  {r['name']}")

    # ═══════════════════════════════════════
    # THICK: φ_AM 비선형 형태
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THICK: φ_AM 형태 (CN^1.5 × cov × f(φ_AM))")
    print("="*80)

    st=sel[thick];ct=am_cn[thick];cvt=cov[thick];pt=pa[thick]
    base_t = SAM * ct**1.5 * cvt

    forms_phi = [
        ('φ^3', pt**3), ('φ^3.5', pt**3.5), ('φ^4', pt**4), ('φ^4.5', pt**4.5), ('φ^5', pt**5),
        ('(φ-0.25)^2', np.clip(pt-0.25,0.001,None)**2),
        ('(φ-0.25)^3', np.clip(pt-0.25,0.001,None)**3),
        ('(φ-0.30)^2', np.clip(pt-0.30,0.001,None)**2),
        ('(φ-0.30)^3', np.clip(pt-0.30,0.001,None)**3),
        ('exp(5φ)', np.exp(5*pt)),
        ('exp(8φ)', np.exp(8*pt)),
        ('exp(10φ)', np.exp(10*pt)),
    ]

    print(f"\n  {'Form':25s} {'R²':>8s}")
    print("  "+"-"*35)
    for label, f_phi in sorted(forms_phi, key=lambda x: -r2l(st, fitC(st,base_t*x[1])*(base_t*x[1]) if fitC(st,base_t*x[1]) else np.ones(len(st)))):
        rhs = base_t * f_phi
        C = fitC(st, rhs)
        if C is None: continue
        pred = C*rhs; r2 = r2l(st, pred)
        print(f"  {label:25s} {r2:8.4f}")

    # ═══════════════════════════════════════
    # T/d THRESHOLD 탐색
    # thick/thin 임계점이 10이 맞는가?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("T/d THRESHOLD: 어디서 나누는 게 최적?")
    print("="*80)

    for td_cut in np.arange(3, 20, 1):
        tk = ratio >= td_cut; tn = ratio < td_cut
        if tk.sum() < 5 or tn.sum() < 5: continue

        # Thick: φ^4 × CN^1.5 × cov
        rhs_tk = SAM * pa[tk]**4 * am_cn[tk]**1.5 * cov[tk]
        C_tk = fitC(sel[tk], rhs_tk)
        if C_tk is None: continue
        r2_tk = r2l(sel[tk], C_tk * rhs_tk)

        # Thin: CN^1.25 × (T/d)^-0.5 × por^3 × cov^1.25 × φ_SE^1.25
        rhs_tn = SAM * am_cn[tn]**1.25 * ratio[tn]**(-0.5) * por[tn]**3 * cov[tn]**1.25 * ps[tn]**1.25
        C_tn = fitC(sel[tn], rhs_tn)
        if C_tn is None: continue
        r2_tn = r2l(sel[tn], C_tn * rhs_tn)

        avg = (r2_tk + r2_tn) / 2
        flag = '★' if avg > 0.94 else ' '
        print(f"  {flag}T/d={td_cut:5.1f}: thick(n={tk.sum():2d}) R²={r2_tk:.4f} | thin(n={tn.sum():2d}) R²={r2_tn:.4f} | avg={avg:.4f}")

if __name__ == '__main__':
    main()
