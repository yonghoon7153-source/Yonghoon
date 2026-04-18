#!/usr/bin/env python3
"""
Electronic Conductivity: Network Solver → G_Holm → Scaling Law

DEM 시뮬레이션에서 전자 전도도를 예측하는 전체 파이프라인.
Kirchhoff resistor network → Holm constriction → Percolation scaling law

Reference:
  - Holm, R. (1967). Electric Contacts: Theory and Application.
  - Kirkpatrick, S. (1973). Percolation and Conduction. Rev. Mod. Phys.
"""
import numpy as np
import json
import os
from pathlib import Path


# ═══════════════════════════════════════════════════════════════
# 1. 물리 상수 및 재료 물성
# ═══════════════════════════════════════════════════════════════

SIGMA_AM = 50.0    # AM (NCM) 전자 전도도 [mS/cm]
PHI_C = 0.10       # AM electronic percolation threshold
                   # (ionic φ_c=0.185보다 낮음 — AM이 크기 때문에 적은 양으로도 percolation)


# ═══════════════════════════════════════════════════════════════
# 2. Holm Constriction Resistance
# ═══════════════════════════════════════════════════════════════
#
# 두 구형 입자가 접촉할 때, 전류는 접촉 면적 A를 통해 흐름.
# 접촉 반지름 a = √(A/π) 일 때:
#
#   R_constriction = 1 / (2σ × a)     [Holm 1967]
#   G_constriction = 2σ × a = 2σ × √(A/π)  ∝  √A
#
# 즉, 접촉 면적이 4배 커지면 전도도는 2배 증가 (√A 관계)


def holm_conductance(area, sigma=SIGMA_AM):
    """단일 접촉의 Holm constriction conductance.

    Args:
        area: 접촉 면적 A [µm²]
        sigma: 재료 전도도 [mS/cm]

    Returns:
        G = 2σ√(A/π) [mS·µm]
    """
    a = np.sqrt(area / np.pi)  # 접촉 반지름
    return 2 * sigma * a


# ═══════════════════════════════════════════════════════════════
# 3. G_Holm: Holm Conductance-Weighted Coordination Number
# ═══════════════════════════════════════════════════════════════
#
# 전통적 CN (coordination number):
#   CN = (접촉 수) / (입자 수)  → 개수만 셈, 품질 무시
#
# G_Holm (이 연구에서 제안):
#   G_Holm = Σ(n_type × √A_type) / n_AM
#          = 입자당 평균 Holm 전도도
#
# 왜 G_Holm이 중요한가:
#   - Bimodal cathode에서 P:S 비율 변화 시:
#     CN은 변하지만 (대-소 입자 혼합으로 접촉 수 감소)
#     G_Holm은 거의 일정 (접촉 면적 증가가 보상)
#   - 이것이 bimodal cathode에서 σ_el이 flat한 이유


def compute_g_holm(metrics):
    """full_metrics.json에서 G_Holm 계산.

    Bimodal: PP, PS, SS 접촉 유형별 면적을 가중 합산
    Monomodal: 단일 유형의 접촉 면적 사용

    Args:
        metrics: full_metrics.json 데이터 (dict)

    Returns:
        g_holm: 입자당 평균 Holm 전도도
    """
    n_AM_P = metrics.get('n_AM_P', 0)
    n_AM_S = metrics.get('n_AM_S', 0)
    n_AM = max(n_AM_P + n_AM_S, 1)

    # 접촉 유형별 개수와 평균 면적
    n_pp = metrics.get('area_AM_P_AM_P_n', 0) or 0    # P-P 접촉 수
    n_ps = metrics.get('area_AM_P_AM_S_n', 0) or 0    # P-S 접촉 수
    n_ss = metrics.get('area_AM_S_AM_S_n', 0) or 0    # S-S 접촉 수

    a_pp = metrics.get('area_AM_P_AM_P_mean', 0) or 0  # P-P 평균 면적 [µm²]
    a_ps = metrics.get('area_AM_P_AM_S_mean', 0) or 0  # P-S 평균 면적 [µm²]
    a_ss = metrics.get('area_AM_S_AM_S_mean', 0) or 0  # S-S 평균 면적 [µm²]

    # G_Holm = Σ(n_type × √A_type) / n_AM
    # √A ∝ Holm conductance (R = 1/2σa, a = √(A/π))
    total_holm = (n_pp * np.sqrt(max(a_pp, 0.001)) +
                  n_ps * np.sqrt(max(a_ps, 0.001)) +
                  n_ss * np.sqrt(max(a_ss, 0.001)))

    if total_holm > 0:
        return total_holm / n_AM
    else:
        # Fallback: 일반 CN 사용
        return metrics.get('am_am_cn', 1.0)


# ═══════════════════════════════════════════════════════════════
# 4. Electronic Scaling Law
# ═══════════════════════════════════════════════════════════════
#
# σ_el = C × σ_AM × CN^1.5 × G_Holm^0.25 × (φ_AM - φ_c)² / por^0.35
#
# 각 항의 물리적 의미:
#
#   CN^1.5:
#     - Kirkpatrick percolation theory에서 σ ∝ (p-pc)^t, t≈2 (3D)
#     - CN은 percolation 확률 p의 proxy
#     - 지수 1.5 (< 2.0)은 Holm conductance가 일부 흡수
#
#   G_Holm^0.25:
#     - 접촉 품질 보정 (Holm constriction conductance)
#     - 작은 지수(0.25)는 G_Holm이 CN과 상관되어 있기 때문
#     - Bimodal에서 CN 변화를 보상하여 flat 경향 재현
#
#   (φ_AM - φ_c)²:
#     - Percolation excess (임계점 위 여유)
#     - φ_c = 0.10: AM electronic percolation threshold
#     - 지수 2: 3D percolation universal exponent
#
#   1/por^0.35:
#     - Porosity에 의한 geometric correction
#     - Void가 전도 경로를 방해하는 효과


def predict_electronic_thick(cn, g_holm, phi_am, porosity, C=None):
    """Electronic conductivity scaling law (thick regime, T/d ≥ 8).

    σ_el = C × σ_AM × CN^1.5 × G_Holm^0.25 × (φ_AM - φ_c)² / por^0.35

    Args:
        cn: AM-AM coordination number
        g_holm: Holm conductance-weighted CN
        phi_am: AM volume fraction
        porosity: porosity [%]
        C: fitting constant (None → use default)

    Returns:
        sigma_el: electronic conductivity [mS/cm]
    """
    if C is None:
        C = 7.5  # default from fitting

    phi_excess = max(phi_am - PHI_C, 0.001)

    sigma = (C * SIGMA_AM *
             cn ** 1.5 *
             g_holm ** 0.25 *
             phi_excess ** 2 /
             porosity ** 0.35)

    return sigma


# ═══════════════════════════════════════════════════════════════
# 5. Network Solver vs Scaling Law 비교
# ═══════════════════════════════════════════════════════════════
#
# Network Solver (ground truth):
#   - 각 AM-AM 접촉을 저항기로 모델링
#   - R = R_bulk + R_constriction (Holm)
#   - Kirchhoff 법칙으로 전극 양단 전도도 계산
#   - 정확하지만 계산 비용 높음 (sparse matrix solve)
#
# Scaling Law (이 연구):
#   - Network solver 결과를 3개 핵심 변수로 압축
#   - CN (connectivity) + G_Holm (contact quality) + φ-φc (percolation)
#   - R² = 0.963, 즉석 예측 가능
#
# 핵심 발견:
#   - CN × √A ≈ constant in bimodal (G_Holm ≈ flat)
#   - bimodal에서 CN이 변해도 σ_el이 flat한 이유
#   - G_Holm이 electronic percolation의 진짜 order parameter


def fit_scaling_constant(cases):
    """전체 데이터에서 fitting constant C를 결정.

    log-space에서 OLS fitting:
    log(σ) = log(C) + 1.5×log(CN) + 0.25×log(G_h) + 2×log(φ-φc) - 0.35×log(por)

    Args:
        cases: list of dicts with 'sigma', 'cn', 'g_holm', 'phi_am', 'por'

    Returns:
        C: fitted constant
        r2: R² in log-space
    """
    log_sigma = np.log(np.array([c['sigma'] for c in cases]))
    log_rhs = np.array([
        1.5 * np.log(c['cn']) +
        0.25 * np.log(c['g_holm']) +
        2.0 * np.log(max(c['phi_am'] - PHI_C, 0.001)) -
        0.35 * np.log(c['por'])
        for c in cases
    ])

    # C fitting: log(C) = mean(log(σ) - log(RHS))
    ln_C = np.mean(log_sigma - log_rhs)
    C = np.exp(ln_C)

    # R² calculation
    predicted = ln_C + log_rhs
    ss_res = np.sum((log_sigma - predicted) ** 2)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma)) ** 2)
    r2 = 1 - ss_res / ss_tot

    return C, r2


# ═══════════════════════════════════════════════════════════════
# 6. 실행
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # 데이터 로딩
    base = Path(__file__).resolve().parent.parent / 'webapp'
    cases = []

    for search_dir in [base / 'archive', base / 'results']:
        if not search_dir.is_dir():
            continue
        for mp in sorted(search_dir.rglob('full_metrics.json')):
            try:
                m = json.loads(mp.read_text())
            except:
                continue

            sigma_el = m.get('electronic_sigma_full_mScm', 0)
            if not sigma_el or sigma_el < 0.001:
                continue

            cn = m.get('am_am_cn', 0)
            phi_am = m.get('phi_am', 0)
            por = max(m.get('porosity', 10), 0.1)
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            thickness = m.get('thickness_um', 0)

            if cn <= 0 or phi_am <= 0 or thickness <= 0 or d_am <= 0:
                continue

            ratio = thickness / d_am
            if ratio < 8:  # thick regime only
                continue

            g_holm = compute_g_holm(m)

            cases.append({
                'sigma': sigma_el,
                'cn': cn,
                'g_holm': g_holm,
                'phi_am': phi_am,
                'por': por,
                'path': str(mp.parent),
            })

    # Dedup
    seen = set()
    unique = []
    for c in cases:
        k = f"{c['phi_am']:.4f}_{c['cn']:.2f}"
        if k not in seen:
            seen.add(k)
            unique.append(c)

    print(f"Electronic Thick Scaling Law")
    print(f"{'='*60}")
    print(f"Data: {len(unique)} unique cases (T/d ≥ 8)")

    # Fit C
    C, r2 = fit_scaling_constant(unique)
    print(f"\nσ_el = {C:.2f} × σ_AM × CN^1.5 × G_Holm^0.25 × (φ-φc)² / por^0.35")
    print(f"R² = {r2:.4f}")
    print(f"φ_c = {PHI_C}")

    # Per-case comparison
    print(f"\n{'CN':>6s} {'G_holm':>6s} {'φ_AM':>6s} {'por':>6s} {'pred':>6s} {'actual':>6s} {'err%':>6s}")
    print("-" * 50)

    errors = []
    for c in unique:
        pred = predict_electronic_thick(c['cn'], c['g_holm'], c['phi_am'], c['por'], C)
        err = abs(pred - c['sigma']) / c['sigma'] * 100
        errors.append(err)
        flag = " <--" if err > 20 else ""
        print(f"{c['cn']:6.2f} {c['g_holm']:6.2f} {c['phi_am']:6.3f} {c['por']:6.1f} {pred:6.2f} {c['sigma']:6.2f} {err:5.1f}%{flag}")

    print(f"\n|err| mean = {np.mean(errors):.1f}%")
    print(f"≤20%: {sum(1 for e in errors if e < 20)}/{len(errors)}")
