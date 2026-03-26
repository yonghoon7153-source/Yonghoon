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
import matplotlib.ticker as ticker
import numpy as np


# ─── Color palette ────────────────────────────────────────────────────────────
BLUE = "#4472C4"
RED = "#C00000"
GREEN = "#548235"
LIGHT_GREEN = "#A9D18E"
BLACK = "#000000"

DPI = 150
FIG_SINGLE = (6, 4)
FIG_FOUR = (12, 10)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _apply_style(ax, ylabel, names):
    """Apply common academic style to an axes."""
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=9, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=10)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, color="#CCCCCC")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _save(fig, outdir, fname):
    fig.tight_layout()
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def _get(data, key, default=0.0):
    """Safely retrieve a numeric value from a dict."""
    v = data.get(key, default)
    if v is None:
        return default
    return float(v)


# ─── Individual plot functions ────────────────────────────────────────────────

def plot_porosity(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "porosity") for d in all_data]
    ax.plot(xs, ys, marker="s", markersize=10, color=BLACK, linewidth=1.5,
            markerfacecolor=BLACK, markeredgecolor=BLACK, zorder=3)
    _apply_style(ax, "Porosity (%)", names)
    ax.set_title("Porosity", fontsize=12, fontweight="bold")
    if standalone:
        return _save(fig, "", "")
    return ax


def _resolve_am_se(d):
    """Standard 케이스에서 AM-SE를 P:S에 맞게 AM_P-SE 또는 AM_S-SE로 분배."""
    am_p = _get(d, "area_AM_P_SE_total")
    am_s = _get(d, "area_AM_S_SE_total")
    # Standard mode: area_AM_SE_total만 있고 AM_P/AM_S 분리 안 됨
    if am_p == 0 and am_s == 0:
        total = _get(d, "area_AM_SE_total") or _get(d, "area_AM전체_SE_total")
        ps = d.get("ps_ratio", "")
        if ps == "P only" or "AM_P" in str(d.get("plate_z_source", "")):
            return total, 0
        else:
            # S only 또는 기본 → AM_S
            return 0, total
    return am_p, am_s


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
    ax.legend(fontsize=9, frameon=False)

    _apply_style(ax, r"Interface Area ($\mu m^2$)", names)
    ax.set_title("AM-SE Interface Area", fontsize=12, fontweight="bold")
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

    ax.bar(xs, counts, width, color=BLUE, alpha=0.8, zorder=3, label="Contact Count")
    _apply_style(ax, "SE-SE Contact Count", names)
    ax.set_title("SE-SE Contact: Count vs Mean Area", fontsize=12, fontweight="bold")
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    ax2 = ax.twinx()
    ax2.plot(xs, means, marker="o", color=RED, linewidth=1.5, markersize=7,
             zorder=4, label="Mean Area")
    ax2.set_ylabel(r"Mean Contact Area ($\mu m^2$)", fontsize=10, color=RED)
    ax2.tick_params(axis="y", labelcolor=RED)
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="upper right", fontsize=8, frameon=False)

    if standalone:
        return _save(fig, "", "")
    return ax


def plot_se_se_total(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "area_SE_SE_total") for d in all_data]
    ax.plot(xs, ys, marker="o", color=BLUE, linewidth=1.5, markersize=8, zorder=3)

    # Annotate max
    if ys:
        idx_max = int(np.argmax(ys))
        ax.annotate(f"max: {ys[idx_max]:.1f}",
                    xy=(idx_max, ys[idx_max]),
                    xytext=(0, 12), textcoords="offset points",
                    fontsize=9, ha="center", color=RED, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))

    _apply_style(ax, r"SE-SE Total Contact Area ($\mu m^2$)", names)
    ax.set_title("SE-SE Total Contact Area", fontsize=12, fontweight="bold")
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

    ax.plot(xs, perc, marker="s", color=BLUE, linewidth=1.5, markersize=8,
            zorder=3, label="Percolation %")
    _apply_style(ax, "Percolation (%)", names)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.set_title("Percolation & Tortuosity", fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    ax2 = ax.twinx()
    ax2.plot(xs, tort, marker="^", color=RED, linewidth=1.5, markersize=8,
             zorder=4, label="Tortuosity")
    ax2.set_ylabel("Tortuosity (mean)", fontsize=10, color=RED)
    ax2.tick_params(axis="y", labelcolor=RED)
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="upper right", fontsize=8, frameon=False)

    if standalone:
        return _save(fig, "", "")
    return ax


def plot_ionic_active(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "ionic_active_pct") for d in all_data]

    ax.plot(xs, ys, marker="o", color=GREEN, linewidth=1.5, markersize=8, zorder=3)
    ax.fill_between(xs, ys, 100, color=RED, alpha=0.15, label="Dead zone")
    ax.fill_between(xs, 0, ys, color=GREEN, alpha=0.15, label="Active zone")
    ax.set_ylim(0, 105)
    _apply_style(ax, "Ionic Active AM (%)", names)
    ax.set_title("Ionic Active AM Fraction", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_coverage(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)

    xs = np.arange(len(names))
    width = 0.3

    mean_p = [_get(d, "coverage_AM_P_mean") for d in all_data]
    std_p = [_get(d, "coverage_AM_P_std") for d in all_data]
    mean_s = [_get(d, "coverage_AM_S_mean") for d in all_data]
    std_s = [_get(d, "coverage_AM_S_std") for d in all_data]

    ax.bar(xs - width / 2, mean_p, width, yerr=std_p, capsize=3,
           color=GREEN, label="AM_P", zorder=3, error_kw=dict(lw=1))
    ax.bar(xs + width / 2, mean_s, width, yerr=std_s, capsize=3,
           color=LIGHT_GREEN, label="AM_S", zorder=3, error_kw=dict(lw=1))

    _apply_style(ax, "Coverage (%)", names)
    ax.set_title("AM Coverage (P vs S)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    if standalone:
        return _save(fig, "", "")
    return ax


# ─── Four-panel composite ────────────────────────────────────────────────────

def plot_four_panel(all_data, names, outdir):
    fig, axes = plt.subplots(2, 2, figsize=FIG_FOUR)
    plot_porosity(all_data, names, ax=axes[0, 0])
    plot_am_se_interface(all_data, names, ax=axes[0, 1])
    plot_se_se_tradeoff(all_data, names, ax=axes[1, 0])
    plot_percolation_tortuosity(all_data, names, ax=axes[1, 1])
    fig.suptitle("DEM Analysis Comparison", fontsize=14, fontweight="bold", y=0.98)
    return _save(fig, outdir, "four_panel.png")


# ─── Plot dispatch table ─────────────────────────────────────────────────────

PLOT_REGISTRY = {
    "porosity": {
        "func": plot_porosity,
        "file": "porosity.png",
        "title": "Porosity",
        "description": "AM_P 비율 증가 → Porosity 감소 (V-shape). 7:3에서 최저.",
        "origin_tip": "Line+Symbol plot. X: P:S Configuration, Y: Porosity(%). Symbol: Square, Size 10. Line: B-Spline.",
    },
    "am_se_interface": {
        "func": plot_am_se_interface,
        "file": "am_se_interface.png",
        "title": "AM-SE Interface Area",
        "description": "AM_P-SE와 AM_S-SE 접촉 면적을 Stacked Bar로 비교. AM_P 비율 증가 시 AM_P-SE 면적 증가, AM_S-SE 감소.",
        "origin_tip": "Stacked Bar. X: Configuration, Y: Interface Area (um2). Colors: Dark Green (AM_P), Light Green (AM_S).",
    },
    "se_se_tradeoff": {
        "func": plot_se_se_tradeoff,
        "file": "se_se_tradeoff.png",
        "title": "SE-SE Contact Trade-off",
        "description": "SE-SE 접촉 개수(bar)와 평균 면적(line) 간 trade-off 관계. 개수 증가 시 평균 면적 감소 경향.",
        "origin_tip": "Double-Y plot. Left Y: Bar (Contact Count, Blue). Right Y: Line+Symbol (Mean Area, Red). X: Configuration.",
    },
    "se_se_total": {
        "func": plot_se_se_total,
        "file": "se_se_total.png",
        "title": "SE-SE Total Contact Area",
        "description": "SE-SE 전체 접촉 면적 비교. 최대값 지점에 annotation 표시.",
        "origin_tip": "Line+Symbol plot. X: Configuration, Y: Total Contact Area (um2). Annotate max point with arrow.",
    },
    "percolation_tortuosity": {
        "func": plot_percolation_tortuosity,
        "file": "percolation_tortuosity.png",
        "title": "Percolation & Tortuosity",
        "description": "이온 전도 경로 형성률(percolation %)과 경로 꼬임 정도(tortuosity)를 동시에 비교.",
        "origin_tip": "Double-Y plot. Left Y: Line (Percolation %, Blue, Square). Right Y: Line (Tortuosity, Red, Triangle). X: Configuration.",
    },
    "ionic_active": {
        "func": plot_ionic_active,
        "file": "ionic_active.png",
        "title": "Ionic Active AM Fraction",
        "description": "이온 전도 경로에 연결된 AM 비율. 100%에서 부족한 부분(red)은 dead zone.",
        "origin_tip": "Line+Fill plot. X: Configuration, Y: Ionic Active (%). Green fill = active, Red fill = dead zone above line.",
    },
    "coverage": {
        "func": plot_coverage,
        "file": "coverage.png",
        "title": "AM Coverage",
        "description": "AM_P와 AM_S의 SE 피복률(coverage) 비교. Error bar는 입자 간 편차(std).",
        "origin_tip": "Grouped Bar + Error Bar. X: Configuration, Y: Coverage (%). Green (AM_P), Light Green (AM_S). Cap size 3.",
    },
}


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate comparison plots from multiple DEM analysis cases."
    )
    parser.add_argument(
        "-i", "--inputs", nargs="+", required=True,
        help="Paths to full_metrics.json files",
    )
    parser.add_argument(
        "-n", "--names", nargs="+", required=True,
        help="Case names for x-axis labels",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory for PNG files",
    )
    parser.add_argument(
        "-p", "--plots", nargs="+",
        default=list(PLOT_REGISTRY.keys()) + ["four_panel"],
        help="Which plots to generate (space-separated)",
    )
    args = parser.parse_args()

    if len(args.inputs) != len(args.names):
        print(f"ERROR: Number of inputs ({len(args.inputs)}) != "
              f"number of names ({len(args.names)})")
        sys.exit(1)

    # Load data
    all_data = []
    for path in args.inputs:
        with open(path, "r") as f:
            all_data.append(json.load(f))

    os.makedirs(args.output, exist_ok=True)

    # Generate requested plots
    plot_info = {}

    for plot_name in args.plots:
        if plot_name == "four_panel":
            outpath = plot_four_panel(all_data, args.names, args.output)
            plot_info["four_panel"] = {
                "file": "four_panel.png",
                "title": "Four-Panel Composite",
                "description": "Porosity, AM-SE Interface, SE-SE Trade-off, Percolation/Tortuosity를 2x2로 한눈에 비교.",
                "origin_tip": "2x2 Multi-panel layout (12x10 inch). Each sub-panel follows its own style. Use Graph > Merge for Origin.",
            }
            print(f"  [OK] four_panel -> {outpath}")
            continue

        if plot_name not in PLOT_REGISTRY:
            print(f"  [SKIP] Unknown plot type: {plot_name}")
            continue

        entry = PLOT_REGISTRY[plot_name]
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
        entry["func"](all_data, args.names, ax=ax)
        outpath = _save(fig, args.output, entry["file"])

        plot_info[plot_name] = {
            "file": entry["file"],
            "title": entry["title"],
            "description": entry["description"],
            "origin_tip": entry["origin_tip"],
        }
        print(f"  [OK] {plot_name} -> {outpath}")

    # Write metadata
    info_path = os.path.join(args.output, "plot_info.json")
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(plot_info, f, indent=2, ensure_ascii=False)
    print(f"\nMetadata written to {info_path}")
    print(f"Total plots generated: {len(plot_info)}")


if __name__ == "__main__":
    main()
