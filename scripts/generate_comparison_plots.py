#!/usr/bin/env python3
"""
Generate comparison plots from multiple DEM analysis cases.

Usage:
    python generate_comparison_plots.py \
        -i case1/full_metrics.json case2/full_metrics.json \
        -n "post_real_7" "post_real_8" \
        -o ./figures \
        -p porosity am_se_interface se_se_tradeoff se_se_total \
           percolation_tortuosity ionic_active coverage four_panel
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ─── Color palette ────────────────────────────────────────────────────────────
BLUE = "#4472C4"
RED = "#C00000"
GREEN = "#548235"
LIGHT_GREEN = "#A9D18E"
BLACK = "#333333"
GRAY = "#888888"

DPI = 150
FIG_SINGLE = (7, 4.5)
FIG_FOUR = (14, 10)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _apply_style(ax, ylabel, names):
    """Apply common academic style."""
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, color="#CCCCCC", alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis='both', labelsize=9)


def _save(fig, outdir, fname):
    fig.tight_layout(pad=1.5)
    fig.subplots_adjust(right=0.85)  # dual Y-axis 우측 여백 확보
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor='white', pad_inches=0.2)
    plt.close(fig)
    return path


def _get(data, key, default=0.0):
    v = data.get(key, default)
    if v is None:
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def _resolve_am_se(d):
    """Standard: AM-SE → P:S에 맞게 AM_P-SE or AM_S-SE로 분배."""
    am_p = _get(d, "area_AM_P_SE_total")
    am_s = _get(d, "area_AM_S_SE_total")
    if am_p == 0 and am_s == 0:
        total = _get(d, "area_AM_SE_total") or _get(d, "area_AM전체_SE_total")
        ps = d.get("ps_ratio", "")
        if ps in ("P only", "10:0"):
            return total, 0
        else:
            return 0, total
    return am_p, am_s


def _resolve_coverage(d):
    """Standard: coverage_AM → P/S에 맞게 분배."""
    cov_p = _get(d, "coverage_AM_P_mean")
    cov_s = _get(d, "coverage_AM_S_mean")
    std_p = _get(d, "coverage_AM_P_std")
    std_s = _get(d, "coverage_AM_S_std")
    if cov_p == 0 and cov_s == 0:
        cov = _get(d, "coverage_AM_mean")
        std = _get(d, "coverage_AM_std")
        ps = d.get("ps_ratio", "")
        if ps in ("P only", "10:0"):
            return cov, std, 0, 0
        else:
            return 0, 0, cov, std
    return cov_p, std_p, cov_s, std_s


# ─── Individual plot functions ────────────────────────────────────────────────

def plot_porosity(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "porosity") for d in all_data]
    ax.plot(xs, ys, marker="s", markersize=9, color=BLACK, linewidth=1.5,
            markerfacecolor=BLACK, markeredgecolor=BLACK, zorder=3)
    # y range with padding
    if ys:
        ymin, ymax = min(ys), max(ys)
        pad = max((ymax - ymin) * 0.15, 0.5)
        ax.set_ylim(ymin - pad, ymax + pad)
    _apply_style(ax, "Porosity (%)", names)
    ax.set_title("Porosity", fontsize=13, fontweight="bold", pad=10)
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_am_se_interface(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = np.arange(len(names))
    width = 0.5

    resolved = [_resolve_am_se(d) for d in all_data]
    am_p = [r[0] for r in resolved]
    am_s = [r[1] for r in resolved]

    ax.bar(xs, am_p, width, label="AM_P-SE", color=GREEN, zorder=3)
    ax.bar(xs, am_s, width, bottom=am_p, label="AM_S-SE",
           color=LIGHT_GREEN, zorder=3)
    ax.legend(fontsize=10, frameon=False)

    _apply_style(ax, "Interface Area (μm²)", names)
    ax.set_title("AM-SE Interface Area", fontsize=13, fontweight="bold", pad=10)
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_se_se_tradeoff(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)

    xs = np.arange(len(names))
    width = 0.5
    counts = [_get(d, "area_SE_SE_n") for d in all_data]
    means = [_get(d, "area_SE_SE_mean") for d in all_data]

    ax.bar(xs, counts, width, color=BLUE, alpha=0.75, zorder=3, label="Contact Count")
    _apply_style(ax, "SE-SE Contact Count", names)
    ax.set_title("SE-SE Contact: Count vs Mean Area", fontsize=13, fontweight="bold", pad=10)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.legend(loc="upper left", fontsize=9, frameon=False)

    ax2 = ax.twinx()
    ax2.plot(xs, means, marker="s", color=RED, linewidth=2, markersize=8,
             zorder=4, label="Mean Area")
    ax2.set_ylabel("Mean Contact Area (μm²)", fontsize=11, color=RED)
    ax2.tick_params(axis="y", labelcolor=RED, labelsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="upper right", fontsize=9, frameon=False)

    if standalone:
        return _save(fig, "", "")
    return ax


def plot_se_se_total(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "area_SE_SE_total") for d in all_data]
    ax.plot(xs, ys, marker="s", color=GRAY, linewidth=1.5, markersize=8, zorder=3,
            markerfacecolor=GRAY, markeredgecolor=GRAY)

    if ys:
        idx_max = int(np.argmax(ys))
        ax.annotate(f"max: {ys[idx_max]:,.0f}",
                    xy=(idx_max, ys[idx_max]),
                    xytext=(15, 10), textcoords="offset points",
                    fontsize=10, ha="left", color=RED, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1.5))
        ymin, ymax = min(ys), max(ys)
        pad = max((ymax - ymin) * 0.15, 100)
        ax.set_ylim(ymin - pad, ymax + pad)

    _apply_style(ax, "SE-SE Total Contact Area (μm²)", names)
    ax.set_title("SE-SE Total Contact Area", fontsize=13, fontweight="bold", pad=10)
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_percolation_tortuosity(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)

    xs = list(range(len(names)))
    perc = [_get(d, "percolation_pct") for d in all_data]
    tort = [_get(d, "tortuosity_mean") for d in all_data]

    ax.plot(xs, perc, marker="s", color=BLUE, linewidth=2, markersize=8,
            zorder=3, label="Percolation %")
    _apply_style(ax, "Percolation (%)", names)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.set_title("Percolation & Tortuosity", fontsize=13, fontweight="bold", pad=10)
    ax.legend(loc="center left", fontsize=9, frameon=False)

    ax2 = ax.twinx()
    ax2.plot(xs, tort, marker="^", color=RED, linewidth=2, markersize=9,
             zorder=4, label="Tortuosity")
    ax2.set_ylabel("Tortuosity", fontsize=11, color=RED)
    ax2.tick_params(axis="y", labelcolor=RED, labelsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="center right", fontsize=9, frameon=False)

    if standalone:
        return _save(fig, "", "")
    return ax


def plot_ionic_active(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "ionic_active_pct") for d in all_data]

    ax.plot(xs, ys, marker="o", color=GREEN, linewidth=2, markersize=9, zorder=3)
    ax.fill_between(xs, ys, 100, color=RED, alpha=0.12, label="Dead zone")
    ax.fill_between(xs, 0, ys, color=GREEN, alpha=0.12, label="Active zone")
    if ys:
        ymin = max(min(ys) - 3, 0)
        ax.set_ylim(ymin, 101)
    _apply_style(ax, "Ionic Active AM (%)", names)
    ax.set_title("Ionic Active AM Fraction", fontsize=13, fontweight="bold", pad=10)
    ax.legend(fontsize=9, frameon=False, loc="lower left")
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_coverage(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)

    xs = np.arange(len(names))
    width = 0.3

    resolved = [_resolve_coverage(d) for d in all_data]
    mean_p = [r[0] for r in resolved]
    std_p = [r[1] for r in resolved]
    mean_s = [r[2] for r in resolved]
    std_s = [r[3] for r in resolved]

    # Only plot bars where values > 0
    has_p = any(v > 0 for v in mean_p)
    has_s = any(v > 0 for v in mean_s)

    if has_p:
        ax.bar(xs - width / 2, mean_p, width, yerr=std_p, capsize=4,
               color=GREEN, label="AM_P", zorder=3, error_kw=dict(lw=1.2))
    if has_s:
        ax.bar(xs + width / 2, mean_s, width, yerr=std_s, capsize=4,
               color=LIGHT_GREEN, label="AM_S", zorder=3, error_kw=dict(lw=1.2))

    _apply_style(ax, "Coverage (%)", names)
    ax.set_title("AM Coverage (P vs S)", fontsize=13, fontweight="bold", pad=10)
    ax.legend(fontsize=10, frameon=False)
    if standalone:
        return _save(fig, "", "")
    return ax


def _save_csv(outdir, fname, names, all_data, keys):
    """Save plot data as CSV for download."""
    import csv
    path = os.path.join(outdir, fname)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Case'] + [k for k, _ in keys])
        for i, name in enumerate(names):
            row = [name]
            for key, transform in keys:
                val = _get(all_data[i], key)
                if transform:
                    val = transform(all_data[i])
                row.append(val)
            writer.writerow(row)
    return fname


# ─── Four-panel composite ────────────────────────────────────────────────────

def plot_four_panel(all_data, names, outdir):
    fig, axes = plt.subplots(2, 2, figsize=FIG_FOUR)
    plot_porosity(all_data, names, ax=axes[0, 0])
    plot_am_se_interface(all_data, names, ax=axes[0, 1])
    plot_se_se_tradeoff(all_data, names, ax=axes[1, 0])
    plot_percolation_tortuosity(all_data, names, ax=axes[1, 1])

    for i, label in enumerate(['(a)', '(b)', '(c)', '(d)']):
        ax = axes.flat[i]
        ax.text(-0.12, 1.05, label, transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top')

    fig.suptitle("DEM Analysis Comparison", fontsize=16, fontweight="bold", y=1.01)
    return _save(fig, outdir, "four_panel.png")


# ─── Plot dispatch table ─────────────────────────────────────────────────────

PLOT_REGISTRY = {
    "porosity": {
        "func": plot_porosity,
        "file": "porosity.png",
        "title": "Porosity",
        "description": "AM_P 비율 증가에 따른 기공률 변화. 7:3 부근에서 최저 (bimodal packing 효과). V-shape 경향이면 최적 조성 존재.",
        "origin_tip": "Line+Symbol → X: P:S Configuration, Y: Porosity(%).\nSymbol: Square (size 10), Color: Black.\nLine: B-Spline, Width 1.5.\nY축 범위: 자동 ± 1%p 여유.",
    },
    "am_se_interface": {
        "func": plot_am_se_interface,
        "file": "am_se_interface.png",
        "title": "AM-SE Interface Area",
        "description": "AM_P-SE와 AM_S-SE 접촉 면적을 Stacked Bar로 비교.\nAM_P 비율 증가 → AM_P-SE 증가, AM_S-SE 감소.\n전체 AM-SE Total은 S only에서 최대.",
        "origin_tip": "Stacked Bar → X: Configuration, Y: Interface Area (μm²).\nColumn1: AM_P-SE (Dark Green #548235).\nColumn2: AM_S-SE (Light Green #A9D18E).\nLegend 우상단.",
    },
    "se_se_tradeoff": {
        "func": plot_se_se_tradeoff,
        "file": "se_se_tradeoff.png",
        "title": "SE-SE Contact Trade-off",
        "description": "SE-SE 접촉 개수(bar)와 개별 접촉 평균 면적(line)의 trade-off.\n\nAM_P↑ → SE가 넓은 공간에 분산 → 접촉 수↑ + 개별 면적↓\n= 악수 약한 사람 10명 vs 악수 센 사람 3명",
        "origin_tip": "Double-Y Axis → Left: Bar (Contact Count, Blue #4472C4).\nRight: Line+Symbol (Mean Area μm², Red #C00000, Square).\nX: Configuration.",
    },
    "se_se_total": {
        "func": plot_se_se_total,
        "file": "se_se_total.png",
        "title": "SE-SE Total Contact Area",
        "description": "SE-SE 전체 접촉 면적 = N × Mean.\n7:3 부근에서 최대 → N 증가가 Mean 감소를 압도.\n10:0에서 감소 시작 → 과도한 AM_P는 SE 품질 저하.",
        "origin_tip": "Line+Symbol → X: Configuration, Y: Total Area (μm²).\nSymbol: Square (Gray), Line: Solid.\n최대값에 Annotation arrow (Red).",
    },
    "percolation_tortuosity": {
        "func": plot_percolation_tortuosity,
        "file": "percolation_tortuosity.png",
        "title": "Percolation & Tortuosity",
        "description": "이온 전도 경로 형성률(Percolation %)과 경로 꼬임(Tortuosity).\n\nPercolation↑ + Tortuosity↓ = 좋은 이온 경로.\nAM_P↑ → 둘 다 개선 (넓은 SE 공간).",
        "origin_tip": "Double-Y Axis → Left: Percolation % (Blue, Square).\nRight: Tortuosity (Red, Triangle).\n범례: 좌측/우측 분리.",
    },
    "ionic_active": {
        "func": plot_ionic_active,
        "file": "ionic_active.png",
        "title": "Ionic Active AM Fraction",
        "description": "Top(SE pellet)에서 이온이 도달 가능한 AM 비율.\n\n100% = 모든 AM이 이온 경로에 연결.\nDead zone(빨강) = SE 미연결 → 반응 불가 AM.",
        "origin_tip": "Line+Fill plot → X: Configuration, Y: Ionic Active (%).\nGreen fill: Active zone (0 ~ line).\nRed fill: Dead zone (line ~ 100).\nY축: 95~100 확대 권장.",
    },
    "coverage": {
        "func": plot_coverage,
        "file": "coverage.png",
        "title": "AM Coverage",
        "description": "AM 표면의 SE 피복률.\n= (SE 접촉 면적) / (AM 자유 표면적) × 100%\n\nAM_P가 클수록 SE 접촉 면적↑ → Coverage↑.\nError bar = 입자 간 편차(std).",
        "origin_tip": "Grouped Bar + Error Bar.\nAM_P: Green #548235, AM_S: Light Green #A9D18E.\nCap size 4, Line width 1.2.\nX: Configuration, Y: Coverage (%).",
    },
}


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", required=True)
    parser.add_argument("-n", "--names", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-p", "--plots", nargs="+",
                        default=list(PLOT_REGISTRY.keys()) + ["four_panel"])
    args = parser.parse_args()

    if len(args.inputs) != len(args.names):
        print(f"ERROR: inputs ({len(args.inputs)}) != names ({len(args.names)})")
        sys.exit(1)

    all_data = []
    for path in args.inputs:
        with open(path, "r") as f:
            all_data.append(json.load(f))

    os.makedirs(args.output, exist_ok=True)
    plot_info = {}

    # Generate CSV for each plot type
    csv_map = {
        'porosity': [('Case', None), ('Porosity(%)', lambda d: _get(d, 'porosity'))],
        'am_se_interface': [('Case', None),
            ('AM_P-SE(μm²)', lambda d: _resolve_am_se(d)[0]),
            ('AM_S-SE(μm²)', lambda d: _resolve_am_se(d)[1]),
            ('Total(μm²)', lambda d: _get(d, 'area_AM전체_SE_total') or sum(_resolve_am_se(d)))],
        'se_se_tradeoff': [('Case', None),
            ('SE-SE Count', lambda d: _get(d, 'area_SE_SE_n')),
            ('SE-SE Mean Area(μm²)', lambda d: _get(d, 'area_SE_SE_mean'))],
        'se_se_total': [('Case', None),
            ('SE-SE Total(μm²)', lambda d: _get(d, 'area_SE_SE_total'))],
        'percolation_tortuosity': [('Case', None),
            ('Percolation(%)', lambda d: _get(d, 'percolation_pct')),
            ('Tortuosity', lambda d: _get(d, 'tortuosity_mean'))],
        'ionic_active': [('Case', None),
            ('Ionic Active AM(%)', lambda d: _get(d, 'ionic_active_pct'))],
        'coverage': [('Case', None),
            ('Coverage AM_P(%)', lambda d: _resolve_coverage(d)[0]),
            ('Coverage AM_P std', lambda d: _resolve_coverage(d)[1]),
            ('Coverage AM_S(%)', lambda d: _resolve_coverage(d)[2]),
            ('Coverage AM_S std', lambda d: _resolve_coverage(d)[3])],
    }

    import csv
    for pname, cols in csv_map.items():
        if pname in args.plots:
            csv_path = os.path.join(args.output, f"{pname}.csv")
            with open(csv_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow([c[0] for c in cols])
                for i, name in enumerate(args.names):
                    row = []
                    for col_name, fn in cols:
                        if fn is None:
                            row.append(name)
                        else:
                            row.append(round(fn(all_data[i]), 4))
                    writer.writerow(row)

    for plot_name in args.plots:
        if plot_name == "four_panel":
            outpath = plot_four_panel(all_data, args.names, args.output)
            plot_info["four_panel"] = {
                "file": "four_panel.png",
                "title": "Four-Panel Composite",
                "description": "Porosity, AM-SE Interface, SE-SE Trade-off,\nPercolation/Tortuosity 4개 핵심 지표를 2×2 패널로 종합 비교.\n\n논문 Figure로 적합.",
                "origin_tip": "Graph > Merge Multiple Graphs.\n2×2 layout, 각 패널 개별 설정.\n전체 크기: 14×10 inch.",
            }
            print(f"  [OK] four_panel -> {outpath}")
            continue

        if plot_name not in PLOT_REGISTRY:
            print(f"  [SKIP] Unknown: {plot_name}")
            continue

        entry = PLOT_REGISTRY[plot_name]
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
        entry["func"](all_data, args.names, ax=ax)
        outpath = _save(fig, args.output, entry["file"])

        csv_file = f"{plot_name}.csv" if plot_name in csv_map else None
        plot_info[plot_name] = {
            "file": entry["file"],
            "csv": csv_file,
            "title": entry["title"],
            "description": entry["description"],
            "origin_tip": entry["origin_tip"],
        }
        print(f"  [OK] {plot_name} -> {outpath}")

    info_path = os.path.join(args.output, "plot_info.json")
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(plot_info, f, indent=2, ensure_ascii=False)
    print(f"\nTotal: {len(plot_info)} plots")


if __name__ == "__main__":
    main()
