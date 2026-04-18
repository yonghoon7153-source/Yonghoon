"""
v2.0l ELEGANT: 깔끔한 지수 후보 비교
"""
import json,os,numpy as np,warnings
from pathlib import Path
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
            ps=m.get('phi_se',0);pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            fp=max(m.get('percolation_pct',0)/100,0.5)
            cn=m.get('se_se_cn',0);gd=m.get('gb_density_mean',0)
            gp=max(m.get('path_conductance_mean',0),1e-6)
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',1))),0.1)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,'gp':gp,'T':T,'cov':cov,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
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
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]);pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    fp=np.array([r['fp'] for r in rows])
    cov=np.array([r['cov'] for r in rows])/100  # fraction
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    sigma_brug=SG*ps*fp/tau**2
    phi_ex=np.clip(ps-0.18,0.001,None)

    print("="*100)
    print("ELEGANT EXPONENT COMPARISON")
    print("All formulas: σ = σ_brug × C × √(1-0.18/φ) × τ^b / f_perc × CN^c × cov^d")
    print("Effective: σ = C × σ_grain × φ^(1/2) × (φ-φ_c)^(1/2) × CN^c × cov^d / τ^(2-b)")
    print("="*100)

    R_perc = np.clip(1-0.18/ps,0.001,None)**0.5

    candidates = [
        # (label, CN_exp, cov_exp, tau_corr_exp, description)
        # Current champions
        ("v4++ best (3/8,5/4)",     1.25,  0.375, 1.5,  "현재 best avg=0.938"),
        ("v4 no cov (3/2)",         1.5,   0,     1.5,  "coverage 없는 v4"),

        # CLEAN: 모든 지수가 n/2 형태
        ("√all: CN¹×√cov",         1.0,   0.5,   1.5,  "σ=C×σ_g×√[φ(φ-φc)×CN²×cov/τ]"),
        ("√all: CN¹×√cov (τ^5/4)", 1.0,   0.5,   1.25, "τ^(3/4) effective"),

        # CLEAN: n/4 형태
        ("¼: CN^(5/4)×⁴√cov",      1.25,  0.25,  1.5,  "coverage를 4제곱근"),
        ("¼: CN×⁴√cov",            1.0,   0.25,  1.5,  "CN도 정수"),
        ("¼: CN^(3/2)×⁴√cov",      1.5,   0.25,  1.5,  "CN^(3/2)"),

        # CLEAN: n/3 형태
        ("⅓: CN^(4/3)×∛cov",       4/3,   1/3,   1.5,  "3제곱근 통일"),
        ("⅓: CN×∛cov",             1.0,   1/3,   1.5,  "CN 정수 + 3제곱근"),
        ("⅓: CN^(4/3)×∛cov (τ^4/3)", 4/3, 1/3,   4/3,  "전부 n/3"),

        # SUPER CLEAN: 정수 + 1/2만
        ("정수: CN×√cov/√τ",        1.0,   0.5,   1.5,  "가장 깔끔"),
        ("정수: CN²×√cov/√τ",       2.0,   0.5,   1.5,  "v3 스타일 CN²"),
        ("정수: CN×cov/√τ",         1.0,   1.0,   1.5,  "coverage 선형"),

        # HALF-INTEGER
        ("반정수: CN^(3/2)×√cov",   1.5,   0.5,   1.5,  "3/2 + 1/2"),
        ("반정수: CN^(3/2)×cov^(3/8)", 1.5, 0.375, 1.5,  "3/2 + 3/8"),

        # τ variations
        ("τ opt: CN^(5/4)×cov^(3/8)×τ^(13/8)", 1.25, 0.375, 1.625, "τ 미세조정"),

        # PHYSICIST's CHOICE: 모든 항이 물리적으로 자연스러운 제곱근
        ("물리: √(φ(φ-φc)/τ) × CN × ⁴√cov", 1.0, 0.25, 1.5, "가장 자연스러운 조합"),
    ]

    print(f"\n  {'Label':42s} {'ALL':>6s} {'LOOCV':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s} {'avg3':>6s} {'C':>7s}")
    print("  "+"-"*110)

    for label, cne, cove, tb, desc in candidates:
        rhs = sigma_brug * R_perc * tau**tb / fp * cn**cne
        if cove > 0:
            rhs = rhs * cov**cove
        C=fitC(sn,rhs)
        if C is None: print(f"  {label:42s} FAILED"); continue
        pred=C*rhs
        r2a=r2l(sn,pred)
        cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3=(r2t+r2m+r2n)/3
        flag='★' if avg3>0.935 else '●' if avg3>0.92 else '○' if avg3>0.90 else ' '
        print(f"  {flag}{label:41s} {r2a:6.4f} {cv:6.4f} {r2t:6.3f} {r2m:6.3f} {r2n:6.3f} {avg3:6.3f} {C:7.4f}")

    # ═══════════════════════════════════════
    # ELEGANT FORMS 정리
    # ═══════════════════════════════════════
    print(f"\n{'='*100}")
    print("ELEGANT FORMS (논문용 표기)")
    print("="*100)

    forms = [
        ("FORM A (single √)", "σ = C × σ_grain × √[ φ_SE(φ_SE-φ_c) × CN² × cov / τ ]",
         1.0, 0.5, 1.5, "모든 항이 √ 안에. 가장 깔끔"),
        ("FORM B (⁴√cov)", "σ = C × σ_grain × √[φ_SE(φ_SE-φ_c)/τ] × CN × ⁴√cov",
         1.0, 0.25, 1.5, "CN 정수, coverage만 4제곱근"),
        ("FORM C (3/2+3/8)", "σ = C × σ_grain × √[φ_SE(φ_SE-φ_c)/τ] × CN^(3/2) × cov^(3/8)",
         1.5, 0.375, 1.5, "현재 best에 가까운 깔끔 버전"),
        ("FORM D (no cov)", "σ = C × σ_grain × √[φ_SE(φ_SE-φ_c)/τ] × CN^(3/2)",
         1.5, 0, 1.5, "coverage 없는 버전"),
    ]

    for fname, formula, cne, cove, tb, note in forms:
        rhs = sigma_brug * R_perc * tau**tb / fp * cn**cne
        if cove > 0:
            rhs = rhs * cov**cove
        C=fitC(sn,rhs)
        if C is None: continue
        pred=C*rhs
        r2a=r2l(sn,pred);cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3=(r2t+r2m+r2n)/3
        err_mean=np.mean(np.abs(sn-pred)/sn*100)
        within20=np.sum(np.abs(sn-pred)/sn<0.2)

        print(f"\n  {fname}: {note}")
        print(f"  {formula}")
        print(f"  C={C:.4f}  R²={r2a:.4f}  LOOCV={cv:.4f}")
        print(f"  thick={r2t:.3f}  mid={r2m:.3f}  thin={r2n:.3f}  avg={avg3:.3f}")
        print(f"  |err|={err_mean:.1f}%  within20%={within20}/{n}")


if __name__ == '__main__':
    main()
