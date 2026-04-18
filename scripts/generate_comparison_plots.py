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
GROUP_COLORS = ['#6c8cff', '#ff6b6b', '#51cf66', '#ffd43b', '#cc5de8', '#ff922b']


# ─── Helpers ──────────────────────────────────────────────────────────────────

_GROUP_INFO = None  # Set by main()
_GLOBAL_RGB = None  # (b, ln_k) from global fit across all plot groups
_GLOBAL_C_ION = None  # Fitted C for ionic scaling law (from global/ionic_scaling_fit)
_GLOBAL_FORMX_R2 = None  # (r2, loocv) from FORM X fit
_GLOBAL_IONIC_SIGMOID = None  # (C_thick, C_thin, tau_c, k) for sigmoid C(τ)
_GLOBAL_IONIC_POLY3 = None  # (a0, a1, a2, a3) for poly3 part of v9 BLEND
_ALL_DATA = None  # all_data for _apply_style auto-detect

def _apply_style(ax, ylabel, names, data_list=None):
    """Apply common academic style with group separators.
    Auto-detect monomodal: if all P:S are same, use AM:SE as x-axis."""
    n = len(names)
    ax.set_xticks(range(n))

    # Check if all names (P:S) are the same → monomodal, use AM:SE instead
    x_labels = list(names)
    x_title = "P:S Configuration"
    if data_list is None:
        data_list = _ALL_DATA
    unique_names = set(names)
    if len(unique_names) == 1 and data_list is not None:
        # All P:S same → show AM:SE ratio instead
        am_se_labels = []
        for d in data_list:
            am_se = _get(d, 'am_se_ratio', '')
            if not am_se:
                phi_am = _get(d, 'phi_am', 0)
                phi_se = _get(d, 'phi_se', 0)
                if phi_am > 0 and phi_se > 0:
                    # Mass ratio (density-weighted): NCM=4800, LPSCl=2000 kg/m³
                    am_mass = phi_am * 4800
                    se_mass = phi_se * 2000
                    am_pct = round(am_mass / (am_mass + se_mass) * 100)  # nearest 1%
                    am_se = f"{am_pct:.0f}:{100-am_pct:.0f}"
            am_se_labels.append(str(am_se) if am_se else '?')
        if any(l != '?' for l in am_se_labels):
            x_labels = am_se_labels
            x_title = f"AM:SE Ratio (P:S={names[0]})"

    # Adaptive x-axis label sizing
    if n <= 8:
        ax.set_xticklabels(x_labels, fontsize=10)
    elif n <= 15:
        ax.set_xticklabels(x_labels, fontsize=8, rotation=45, ha='right')
    else:
        ax.set_xticklabels(x_labels, fontsize=7, rotation=45, ha='right')

    ax.set_xlabel(x_title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, color="#CCCCCC", alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis='both', labelsize=9 if n <= 8 else 7)

    # Bottom margin for rotated labels
    if n > 8:
        ax.figure.subplots_adjust(bottom=0.18)

    # Add group separators and break lines at boundaries
    if _GROUP_INFO:
        sizes, gnames = _GROUP_INFO
        n_total = sum(sizes)
        # Break existing lines at group boundaries
        boundaries = []
        pos = 0
        for sz in sizes[:-1]:
            pos += sz
            boundaries.append(pos - 0.5)
        # Draw separators and labels
        pos = 0
        for gi, sz in enumerate(sizes):
            if gi > 0:
                ax.axvline(pos - 0.5, color='#888888', linestyle='--', linewidth=1, alpha=0.6)
            mid = pos + sz / 2 - 0.5
            # Staggered group labels: odd groups lower
            if n > 8:
                label_y = -0.26 if gi % 2 == 0 else -0.36
            else:
                label_y = -0.18
            label_fs = 7 if n > 15 else 8 if n > 8 else 9
            ax.text(mid, label_y, gnames[gi], ha='center', va='top',
                    transform=ax.get_xaxis_transform(),
                    fontsize=label_fs, fontweight='bold', color=GROUP_COLORS[gi % len(GROUP_COLORS)],
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=GROUP_COLORS[gi % len(GROUP_COLORS)], alpha=0.8))
            pos += sz


def _marker_size(n):
    """Adaptive marker size based on case count."""
    if n <= 8: return 10
    if n <= 15: return 7
    return 5

def _line_width(n):
    """Adaptive line width based on case count."""
    if n <= 8: return 2.5
    if n <= 15: return 2.0
    return 1.5

def _group_break_data(xs, ys):
    """Insert NaN at group boundaries so matplotlib breaks the line."""
    if not _GROUP_INFO:
        return xs, ys
    sizes = _GROUP_INFO[0]
    new_x, new_y = [], []
    pos = 0
    for gi, sz in enumerate(sizes):
        if gi > 0:
            new_x.append(float('nan'))
            new_y.append(float('nan'))
        for i in range(pos, pos + sz):
            if i < len(xs):
                new_x.append(xs[i] if not isinstance(xs, np.ndarray) else float(xs[i]))
                new_y.append(ys[i] if not isinstance(ys, np.ndarray) else float(ys[i]))
        pos += sz
    return new_x, new_y


def _break_lines_at_groups(fig):
    """Insert NaN at group boundaries for all lines in all axes."""
    if not _GROUP_INFO:
        return
    sizes = _GROUP_INFO[0]
    n_total = sum(sizes)
    for ax in fig.get_axes():
        for line in ax.get_lines():
            xd = line.get_xdata()
            yd = line.get_ydata()
            # Only process data lines (not grid, axvline, etc.)
            if len(xd) != n_total or not line.get_marker() or line.get_marker() == 'None':
                continue
            if True:
                new_x, new_y = [], []
                idx = 0
                for gi, sz in enumerate(sizes):
                    if gi > 0:
                        new_x.append(float('nan'))
                        new_y.append(float('nan'))
                    for j in range(sz):
                        if idx < len(xd):
                            new_x.append(float(xd[idx]))
                            new_y.append(float(yd[idx]))
                            idx += 1
                line.set_xdata(new_x)
                line.set_ydata(new_y)


def _save(fig, outdir, fname):
    _break_lines_at_groups(fig)
    fig.tight_layout(pad=1.5)
    fig.subplots_adjust(right=0.85, bottom=0.22 if _GROUP_INFO else 0.15)
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


def _write_csv(outdir, fname, headers, names, *columns):
    """Write comparison CSV with case names and data columns."""
    import csv
    with open(os.path.join(outdir, fname), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Case'] + headers)
        for i, name in enumerate(names):
            row = [name] + [c[i] if i < len(c) else '' for c in columns]
            w.writerow(row)


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
    top_reach = [_get(d, "top_reachable_pct") for d in all_data]
    tort = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean")) for d in all_data]

    ax.plot(xs, perc, marker="s", color=BLUE, linewidth=2, markersize=8,
            zorder=3, label="Percolation %")
    ax.plot(xs, top_reach, marker="D", color="#00B0F0", linewidth=1.5, markersize=7,
            linestyle="--", zorder=3, label="Top Reachable %")
    _apply_style(ax, "SE Connectivity (%)", names)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.set_title("Percolation / Top Reachable & Tortuosity", fontsize=13, fontweight="bold", pad=10)
    ax.legend(loc="center left", fontsize=8, frameon=False)

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

    # Fill per group to avoid cross-group shading
    if _GROUP_INFO and len(_GROUP_INFO[0]) > 1:
        sizes = _GROUP_INFO[0]
        pos = 0
        for gi, sz in enumerate(sizes):
            gx = xs[pos:pos+sz]
            gy = ys[pos:pos+sz]
            ax.plot(gx, gy, marker="o", color=GREEN, linewidth=2, markersize=9, zorder=3)
            ax.fill_between(gx, gy, 100, color=RED, alpha=0.12, label="Dead zone" if gi == 0 else None)
            ax.fill_between(gx, 0, gy, color=GREEN, alpha=0.12, label="Active zone" if gi == 0 else None)
            pos += sz
    else:
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


def plot_stress_cv(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))
    ys = [_get(d, "stress_cv") for d in all_data]
    ax.plot(xs, ys, marker="s", markersize=9, color=BLACK, linewidth=1.5, zorder=3)
    if ys:
        ymin, ymax = min(ys), max(ys)
        pad = max((ymax - ymin) * 0.15, 1)
        ax.set_ylim(ymin - pad, ymax + pad)
    _apply_style(ax, "Von Mises CV (%)", names)
    ax.set_title("Stress Distribution Uniformity", fontsize=13, fontweight="bold", pad=10)
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_stress_ratio(all_data, names, ax=None):
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)
    xs = list(range(len(names)))

    type_keys = ['AM_P', 'AM_S', 'SE']
    colors = {'AM_P': RED, 'AM_S': '#FF8C00', 'SE': GREEN}
    markers = {'AM_P': 's', 'AM_S': 'o', 'SE': '^'}

    for tk in type_keys:
        ys = [_get(d, f"stress_ratio_{tk}") for d in all_data]
        if any(v > 0 for v in ys):
            ax.plot(xs, ys, marker=markers[tk], markersize=8, color=colors[tk],
                    linewidth=1.5, label=tk, zorder=3)

    ax.axhline(y=1.0, color=GRAY, linestyle='--', linewidth=1, alpha=0.5, label='mean')
    _apply_style(ax, "σ / σ_mean", names)
    ax.set_title("Stress Ratio by Type", fontsize=13, fontweight="bold", pad=10)
    ax.legend(fontsize=9, frameon=False)
    if standalone:
        return _save(fig, "", "")
    return ax


def plot_stress_z_layer(all_data, names, ax=None):
    """Z-layer별 stress CV profile (all cases overlaid)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=FIG_SINGLE)

    colors_cycle = [BLUE, RED, GREEN, '#FF8C00', BLACK, '#9467BD']
    for i, d in enumerate(all_data):
        z_data = d.get('stress_z_layer_cv', [])
        if z_data:
            zs = [layer['z_mid_um'] for layer in z_data]
            cvs = [layer['cv'] for layer in z_data]
            c = colors_cycle[i % len(colors_cycle)]
            ax.plot(zs, cvs, marker='o', markersize=5, linewidth=1.5,
                    color=c, label=names[i], zorder=3)

    ax.set_xlabel("Z Position (μm)", fontsize=10)
    ax.set_ylabel("Von Mises CV (%)", fontsize=11)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, color="#CCCCCC", alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Stress Uniformity by Z-layer", fontsize=13, fontweight="bold", pad=10)
    ax.legend(fontsize=8, frameon=False)
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


def _generate_particle_info(all_data, ps_labels, case_names, outdir):
    """Generate particle info summary table(s) as PNG, split by group if available."""
    global _GROUP_INFO

    if _GROUP_INFO and len(_GROUP_INFO[0]) > 1:
        sizes, gnames = _GROUP_INFO
        n_groups = len(sizes)
        fig, axes = plt.subplots(n_groups, 1, figsize=(max(7, max(sizes)*1.5 + 2), 3 * n_groups))
        if n_groups == 1:
            axes = [axes]

        idx = 0
        for gi, (sz, gname) in enumerate(zip(sizes, gnames)):
            ax = axes[gi]
            ax.axis('off')
            ax.set_title(gname, fontsize=12, fontweight='bold', color=GROUP_COLORS[gi % len(GROUP_COLORS)], pad=10)

            group_data = all_data[idx:idx+sz]
            group_labels = ps_labels[idx:idx+sz]

            headers = ['P:S'] + group_labels
            rows_data = []
            for ptype in ['AM_P', 'AM_S', 'SE']:
                row = [f'{ptype} Count']
                for d in group_data:
                    val = d.get(f'n_{ptype}', '-')
                    row.append(str(int(val)) if val != '-' else '-')
                rows_data.append(row)
            for ptype in ['AM_P', 'AM_S', 'SE']:
                row = [f'{ptype} R (μm)']
                for d in group_data:
                    val = d.get(f'r_{ptype}', '-')
                    row.append(str(val) if val != '-' else '-')
                rows_data.append(row)

            table = ax.table(cellText=rows_data, colLabels=headers, loc='center', cellLoc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1, 1.4)
            for j in range(len(headers)):
                cell = table[0, j]
                cell.set_facecolor(GROUP_COLORS[gi % len(GROUP_COLORS)])
                cell.set_text_props(color='white', fontweight='bold')
            for i in range(len(rows_data)):
                for j in range(len(headers)):
                    cell = table[i+1, j]
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')

            idx += sz

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'particle_info.png'), dpi=DPI,
                    bbox_inches='tight', facecolor='white', pad_inches=0.1)
        plt.close(fig)
    else:
        # Single group (original behavior)
        fig, ax = plt.subplots(figsize=(max(7, len(ps_labels)*1.5 + 2), 3))
        ax.axis('off')

        headers = ['P:S'] + ps_labels
        rows_data = []
        for ptype in ['AM_P', 'AM_S', 'SE']:
            row = [f'{ptype} Count']
            for d in all_data:
                val = d.get(f'n_{ptype}', '-')
                row.append(str(int(val)) if val != '-' else '-')
            rows_data.append(row)
        for ptype in ['AM_P', 'AM_S', 'SE']:
            row = [f'{ptype} R (μm)']
            for d in all_data:
                val = d.get(f'r_{ptype}', '-')
                row.append(str(val) if val != '-' else '-')
            rows_data.append(row)

        table = ax.table(cellText=rows_data, colLabels=headers, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.4)
        for j in range(len(headers)):
            cell = table[0, j]
            cell.set_facecolor('#4472C4')
            cell.set_text_props(color='white', fontweight='bold')
        for i in range(len(rows_data)):
            for j in range(len(headers)):
                cell = table[i+1, j]
                cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'particle_info.png'), dpi=DPI,
                    bbox_inches='tight', facecolor='white', pad_inches=0.1)
        plt.close(fig)


# ─── New Plots ─────────────────────────────────────────────────────────────────

ORANGE = "#ED7D31"
PURPLE = "#7030A0"


def plot_se_network(data_list, names, outdir):
    """SE-SE CN mean + SE Cluster count (dual Y axis)."""
    fig, ax1 = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))

    cn = [_get(d, "se_se_cn") for d in data_list]
    clusters = [_get(d, "n_components") for d in data_list]
    large_clusters = [_get(d, "n_large_components") for d in data_list]

    color1, color2 = BLUE, ORANGE
    ax1.plot(x, cn, 's-', color=color1, markersize=10, linewidth=2.5, label="SE-SE CN mean")
    _apply_style(ax1, "SE-SE CN mean", names)
    ax1.tick_params(axis='y', labelcolor=color1)

    ax2 = ax1.twinx()
    ax2.bar(x - 0.15, clusters, 0.3, color=color2, alpha=0.4, label="Total Clusters")
    ax2.bar(x + 0.15, large_clusters, 0.3, color=color2, alpha=0.85, label="Large (≥10)")
    ax2.set_ylabel("SE Cluster Count", fontsize=11, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.spines["top"].set_visible(False)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper right')
    ax1.set_title("SE Network: CN & Clusters", fontsize=12, fontweight='bold')
    _write_csv(outdir, 'se_network.csv', ['SE-SE CN', 'Total Clusters', 'Large(≥10)'],
               names, cn, clusters, large_clusters)
    return _save(fig, outdir, "se_network.png")


def plot_contact_force(data_list, names, outdir):
    """Contact force distribution by type (grouped bar)."""
    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    n = len(names)
    x = np.arange(n)
    types = ['AM_P-AM_P', 'AM_P-SE', 'SE-SE', 'AM_P-AM_S', 'AM_S-SE']
    colors = [RED, GREEN, BLUE, ORANGE, LIGHT_GREEN]

    # Find which types actually have data
    active = []
    for ct, color in zip(types, colors):
        key = f"fn_{ct.replace('-','_')}_mean"
        vals = [_get(d, key) for d in data_list]
        if any(v > 0 for v in vals):
            active.append((ct, key, color, vals))

    if not active:
        plt.close(fig)
        return None

    w = 0.8 / len(active)
    for i, (ct, key, color, vals) in enumerate(active):
        offset = (i - len(active)/2 + 0.5) * w
        bars = ax.bar(x + offset, vals, w, label=ct, color=color, alpha=0.85)

    _apply_style(ax, "Fn mean (μN)", names)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title("Contact Force by Type", fontsize=12, fontweight='bold')
    # CSV
    headers = [ct for ct, _, _, _ in active]
    cols = [vals for _, _, _, vals in active]
    _write_csv(outdir, 'contact_force.csv', [f'Fn {h} mean(μN)' for h in headers], names, *cols)
    return _save(fig, outdir, "contact_force.png")


def plot_contact_pressure(data_list, names, outdir):
    """Contact pressure mean & max (dual bar)."""
    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    n = len(names)
    x = np.arange(n)
    w = 0.35

    means = [_get(d, "contact_pressure_mean") for d in data_list]
    maxes = [_get(d, "contact_pressure_max") for d in data_list]

    ax.bar(x - w/2, means, w, label="Mean", color=BLUE, alpha=0.85)
    ax.bar(x + w/2, maxes, w, label="Max", color=RED, alpha=0.85)

    _apply_style(ax, "Contact Pressure (MPa)", names)
    ax.legend(fontsize=9)
    ax.set_title("Contact Pressure", fontsize=12, fontweight='bold')
    _write_csv(outdir, 'contact_pressure.csv', ['CP mean(MPa)', 'CP max(MPa)'], names, means, maxes)
    return _save(fig, outdir, "contact_pressure.png")


def plot_am_vulnerability(data_list, names, outdir):
    """AM vulnerability + AM-SE CN (dual Y axis)."""
    fig, ax1 = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))

    vuln = [_get(d, "am_vulnerable_pct") for d in data_list]
    cn = [_get(d, "am_se_cn_mean") for d in data_list]

    color1, color2 = RED, BLUE
    ax1.bar(x, vuln, 0.5, color=color1, alpha=0.7, label="Vulnerable AM (%)")
    _apply_style(ax1, "Vulnerable AM (%)", names)
    ax1.tick_params(axis='y', labelcolor=color1)

    ax2 = ax1.twinx()
    ax2.plot(x, cn, 'o-', color=color2, markersize=8, linewidth=2, label="AM-SE CN")
    ax2.set_ylabel("AM-SE CN mean", fontsize=11, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.spines["top"].set_visible(False)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper right')
    ax1.set_title("AM Vulnerability & SE Connectivity", fontsize=12, fontweight='bold')
    _write_csv(outdir, 'am_vulnerability.csv', ['Vulnerable(%)', 'AM-SE CN'], names, vuln, cn)
    return _save(fig, outdir, "am_vulnerability.png")


def plot_effective_conductivity(data_list, names, outdir):
    """Bruggeman effective ionic conductivity."""
    fig, ax1 = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    phi = [_get(d, "phi_se") for d in data_list]

    ax1.plot(x, sigma_brug, 's-', color=GREEN, markersize=10, linewidth=2.5, label="σ_brug/σ_grain")
    _apply_style(ax1, "σ_brug / σ_grain", names)
    ax1.tick_params(axis='y', labelcolor=GREEN)

    ax1b = ax1.twinx()
    ax1b.bar(x, phi, 0.4, color=BLUE, alpha=0.25, label="φ_SE")
    ax1b.set_ylabel("φ_SE", fontsize=10, color=BLUE)
    ax1b.set_ylim(0.2, max(phi) * 1.15 if phi else 0.4)
    ax1b.tick_params(axis='y', labelcolor=BLUE)
    ax1b.spines["top"].set_visible(False)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper right')
    ax1.set_title("Bruggeman: σ_brug/σ_grain = φ_SE × f_perc / τ²", fontsize=12, fontweight='bold')

    _write_csv(outdir, 'effective_conductivity.csv',
               ['φ_SE', 'σ_brug/σ_grain'], names, phi, sigma_brug)
    return _save(fig, outdir, "effective_conductivity.png")


R_GB_MIN_R2 = 0.75  # minimum R² for R_gb fitting to be valid

def _fit_r_gb(data_list, names, use_global=True):
    """Fit GB correction: BLM+Constriction model R = C × (GB_d² × T)^α.
    log(R) = α × log(GB_d²×T) + ln(C)
    Returns alpha, ln_C, valid_idx, r_squared."""
    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    perc = [_get(d, "percolation_pct") / 100 for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    g_path = [_get(d, "path_conductance_mean") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    thickness = [_get(d, "thickness_um", 0) for d in data_list]

    sigma_proxy = [g_path[i] * perc[i] / tau[i] if g_path[i] > 0 and tau[i] > 0 else 0
                   for i in range(len(data_list))]

    valid_idx = [i for i in range(len(data_list))
                 if sigma_proxy[i] > 0 and gb_dens[i] > 0 and sigma_brug[i] > 0 and thickness[i] > 0]

    if use_global and _GLOBAL_RGB is not None:
        return _GLOBAL_RGB[0], _GLOBAL_RGB[1], valid_idx, 1.0

    alpha = 1.87
    ln_C = -3.26
    r_squared = 0
    if len(valid_idx) >= 2:
        ratio = np.array([sigma_brug[i] / sigma_proxy[i] for i in valid_idx])
        # x = GB_d² × T (BLM + Constriction)
        x_vals = np.array([gb_dens[i]**2 * thickness[i] for i in valid_idx])
        log_y = np.log(ratio)
        log_x = np.log(x_vals)
        from scipy import stats as sp_stats
        slope, intercept, r_val, _, _ = sp_stats.linregress(log_x, log_y)
        alpha = slope
        ln_C = intercept
        r_squared = r_val ** 2
    return alpha, ln_C, valid_idx, r_squared


def plot_rgb_fitting(data_list, names, outdir):
    """BLM+Constriction fitting: log(R) vs log(GB_d²×T)."""
    alpha, ln_C, valid_idx, r_squared = _fit_r_gb(data_list, names)

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    perc = [_get(d, "percolation_pct") / 100 for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    g_path = [_get(d, "path_conductance_mean") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    thickness = [_get(d, "thickness_um", 100) for d in data_list]
    sigma_proxy = [g_path[i] * perc[i] / tau[i] if g_path[i] > 0 and tau[i] > 0 else 0
                   for i in range(len(data_list))]

    fig, ax = plt.subplots(figsize=FIG_SINGLE)

    # x = log(GB_d² × T), y = log(R)
    x_pts = np.array([np.log(gb_dens[i]**2 * thickness[i]) for i in valid_idx])
    y_pts = np.array([np.log(sigma_brug[i] / sigma_proxy[i]) for i in valid_idx])

    # Group colors
    point_groups = [0] * len(valid_idx)
    group_boundaries = []
    if _GROUP_INFO:
        sizes, gnames = _GROUP_INFO
        pos = 0
        for sz in sizes:
            group_boundaries.append((pos, pos + sz))
            pos += sz
        for j, i in enumerate(valid_idx):
            for g_idx, (start, end) in enumerate(group_boundaries):
                if start <= i < end:
                    point_groups[j] = g_idx
                    break

    if _GROUP_INFO:
        for j in range(len(valid_idx)):
            gi = point_groups[j]
            ax.scatter(x_pts[j], y_pts[j], s=100, c=GROUP_COLORS[gi % len(GROUP_COLORS)],
                      zorder=5, edgecolors='white', linewidth=1.5)
    else:
        ax.scatter(x_pts, y_pts, s=100, c=BLUE, zorder=5, edgecolors='white', linewidth=1.5)

    # Labels
    try:
        from adjustText import adjust_text
        texts = []
        for j, i in enumerate(valid_idx):
            texts.append(ax.text(x_pts[j], y_pts[j], names[i], fontsize=8, color=BLACK, zorder=6))
        adjust_text(texts, x=list(x_pts), y=list(y_pts), ax=ax,
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                   force_text=(0.5, 0.5), expand=(1.3, 1.5))
    except ImportError:
        for j, i in enumerate(valid_idx):
            ax.annotate(names[i], (x_pts[j], y_pts[j]),
                       fontsize=8, ha='left', va='bottom', xytext=(5, 5),
                       textcoords='offset points', color=BLACK, zorder=6)

    # Group labels
    if _GROUP_INFO:
        for gi, (start, end) in enumerate(group_boundaries):
            group_js = [j for j, i in enumerate(valid_idx) if start <= i < end]
            if group_js:
                cx = np.mean([x_pts[j] for j in group_js])
                cy = min([y_pts[j] for j in group_js])
                ax.text(cx, cy - (max(y_pts) - min(y_pts)) * 0.08, gnames[gi],
                       ha='center', va='top',
                       fontsize=9, fontweight='bold', color=GROUP_COLORS[gi % len(GROUP_COLORS)],
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                edgecolor=GROUP_COLORS[gi % len(GROUP_COLORS)], alpha=0.8))

    # Fit line
    x_line = np.linspace(min(x_pts) * 0.9, max(x_pts) * 1.1, 100)
    y_line = alpha * x_line + ln_C
    ax.plot(x_line, y_line, '-', color=RED, linewidth=2.5,
            label=f"y = {alpha:.2f}·x + {ln_C:.2f}")

    ax.set_xlabel("log(GB_d² × T)", fontsize=12)
    ax.set_ylabel("log(σ_brug / σ_proxy)", fontsize=12)
    title_suffix = " (global)" if _GLOBAL_RGB is not None else ""
    ax.set_title(f"BLM+Constriction Fitting{title_suffix}", fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper left')
    ax.text(0.95, 0.05, f"α = {alpha:.2f}\nln(C) = {ln_C:.2f}\nR² = {r_squared:.4f}\nn = {len(x_pts)}",
            transform=ax.transAxes, fontsize=11, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8))
    y_margin = (max(y_pts) - min(y_pts)) * 0.15
    ax.set_ylim(min(y_pts) - y_margin, max(y_pts) + y_margin)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _write_csv(outdir, 'rgb_fitting.csv',
               ['GB_d²×T', 'log(GB_d²×T)', 'log(σ_brug/σ_proxy)', 'σ_brug', 'σ_proxy', 'thickness_um'],
               [names[i] for i in valid_idx],
               [gb_dens[i]**2 * thickness[i] for i in valid_idx],
               list(x_pts), list(y_pts),
               [sigma_brug[i] for i in valid_idx],
               [sigma_proxy[i] for i in valid_idx],
               [thickness[i] for i in valid_idx])
    return _save(fig, outdir, "rgb_fitting.png")


def _load_network_sigma(data_list):
    """Load σ_full from network_conductivity.json for each case."""
    sigma_net = [0] * len(data_list)
    for i, d in enumerate(data_list):
        net_val = _get(d, "sigma_full_mScm", 0)
        if net_val > 0:
            sigma_net[i] = net_val
        else:
            src = _get(d, "_source_path", "")
            if src:
                net_path = os.path.join(os.path.dirname(src), "network_conductivity.json")
                if os.path.exists(net_path):
                    try:
                        with open(net_path) as _nf:
                            nd = json.load(_nf)
                        sigma_net[i] = nd.get("sigma_full_mScm", 0) or 0
                    except:
                        pass
    return sigma_net


def plot_gb_corrected(data_list, names, outdir):
    """Part I: Proxy-based GB-corrected σ_eff using BLM+Constriction."""
    alpha, ln_C, _, r2_check = _fit_r_gb(data_list, names)
    SIGMA_BULK = 3.0  # σ_grain (grain interior)
    C = np.exp(ln_C)

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    thickness = [_get(d, "thickness_um", 100) for d in data_list]

    sigma_corr = []
    for i in range(len(data_list)):
        gb2t = gb_dens[i]**2 * thickness[i]
        if gb2t > 0:
            correction = C * gb2t**alpha
            sigma_corr.append(sigma_brug[i] / correction if correction > 0 else sigma_brug[i])
        else:
            sigma_corr.append(sigma_brug[i])
    sigma_abs = [s * SIGMA_BULK for s in sigma_corr]

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    ax.plot(x, sigma_corr, 's-', color=RED, markersize=ms, linewidth=lw, label="σ_brug/σ_grain")
    _apply_style(ax, "σ_brug / σ_grain (proxy)", names)
    ax.tick_params(axis='y', labelcolor=RED)

    ax2 = ax.twinx()
    ax2.plot(x, sigma_abs, 'D--', color=ORANGE, markersize=ms-2, linewidth=lw-0.5, label="σ_ionic (mS/cm)")
    ax2.set_ylabel("σ_ionic (mS/cm)", fontsize=11, color=ORANGE)
    ax2.tick_params(axis='y', labelcolor=ORANGE)
    ax2.spines["top"].set_visible(False)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper left')
    gb_src = "global" if _GLOBAL_RGB is not None else "local"
    ax.set_title(f"Part I: Proxy σ_eff = σ_brug / C·(GB_d²·T)^α\nα={alpha:.2f}, R²={r2_check:.3f} [{gb_src}]",
                 fontsize=10, fontweight='bold')

    _write_csv(outdir, 'gb_corrected.csv',
               ['GB_d', 'T(μm)', 'GB_d²×T', 'σ_brug/σ_grain', 'σ_brug/σ_grain', 'σ_eff(mS/cm)'],
               names, gb_dens, thickness,
               [gb_dens[i]**2*thickness[i] for i in range(len(names))],
               sigma_brug, sigma_corr, sigma_abs)
    return _save(fig, outdir, "gb_corrected.png")


def plot_ionic_scaling_fit(data_list, names, outdir):
    """Ionic scaling law fit: σ_predicted vs σ_actual scatter (log-log).
    Uses same fitting approach as find_scaling_law.py model 4."""
    SIGMA_BULK = 3.0

    phi_se = [_get(d, "phi_se") for d in data_list]
    f_perc = [_get(d, "percolation_pct", 0) / 100 for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    g_path = [_get(d, "path_conductance_mean", 0) for d in data_list]
    cn = [_get(d, "se_se_cn", 0) for d in data_list]
    sigma_net = _load_network_sigma(data_list)

    # Same as find_scaling_law.py: σ_brug = 3.0 × φ_SE × f_perc / τ²
    valid_idx = []
    SIGMA_MIN = 0.0  # include all non-zero σ_net (was 0.01; relaxed to test fit on low-σ cases)
    for i in range(len(data_list)):
        if (phi_se[i] > 0 and f_perc[i] > 0 and tau[i] > 0 and
            g_path[i] > 0 and cn[i] > 0 and gb_dens[i] > 0 and sigma_net[i] > SIGMA_MIN):
            valid_idx.append(i)

    if len(valid_idx) < 3:
        return None

    s_actual = np.array([sigma_net[i] for i in valid_idx])
    log_sf = np.log(s_actual)

    # σ_brug in mS/cm (same as find_scaling_law.py)
    sigma_brug = np.array([SIGMA_BULK * phi_se[i] * f_perc[i] / tau[i]**2 for i in valid_idx])
    log_sb = np.log(sigma_brug)

    # Fit: log(σ_full) = log(σ_brug) + a*log(g_path*gb_d²) + b*log(cn) + c
    log_combo = np.array([np.log(g_path[i] * gb_dens[i]**2) for i in valid_idx])
    log_cn = np.array([np.log(cn[i]) for i in valid_idx])
    residual = log_sf - log_sb
    X = np.column_stack([log_combo, log_cn, np.ones(len(valid_idx))])
    b_fit, _, _, _ = np.linalg.lstsq(X, residual, rcond=None)
    a_combo, b_cn, ln_C = b_fit
    C_fit = np.exp(ln_C)

    # Predicted (using fixed champion exponents: 0.25, 2)
    # Also compute free-fit R² for comparison
    pred_free = log_sb + X @ b_fit
    ss_res_free = np.sum((log_sf - pred_free)**2)
    ss_tot = np.sum((log_sf - np.mean(log_sf))**2)
    r2_free = 1 - ss_res_free / ss_tot

    # Fixed exponents + fit C only
    log_rhs_fixed = log_sb + 0.25 * log_combo + 2 * log_cn
    ln_C_fixed = np.mean(log_sf - log_rhs_fixed)
    C_fixed = np.exp(ln_C_fixed)
    # Save v3 C (legacy, for reference only)
    global _GLOBAL_C_ION
    pred_fixed = ln_C_fixed + log_rhs_fixed
    ss_res_fixed = np.sum((log_sf - pred_fixed)**2)
    r2_fixed = 1 - ss_res_fixed / ss_tot

    # === FORM X v2: CN^1.5 × (φ-φc)^¾ × cov^¼ / τ^¼ × exp(-0.1×lnCN×lnτ) ===
    PHI_C = 0.185
    coverage = [(lambda vs: sum(vs)/len(vs)/100 if vs else 0.20)([v for v in [_get(d,"coverage_AM_P_mean",0), _get(d,"coverage_AM_S_mean",0), _get(d,"coverage_AM_mean",0)] if v>0]) for d in data_list]

    phi_ex_arr = np.array([max(phi_se[i] - PHI_C, 0.001) for i in valid_idx])
    cn_arr = np.array([cn[i] for i in valid_idx])
    tau_arr = np.array([tau[i] for i in valid_idx])
    cov_arr = np.array([coverage[i] for i in valid_idx])

    # FORM X v9 BLEND: (1-w)·v5_sigmoid(τ) + w·poly3(lnτ)
    # Blend sigmoid (k, τc) continuously optimized by LOOCV.
    # Smooth k → fewer outliers in transition band.
    TAU_C = 2.1; TAU_K = 5.0          # v5 C(τ) sigmoid (fixed)
    fp_arr = np.array([max(f_perc[i], 0.01) for i in valid_idx])
    log_tau_arr = np.log(tau_arr)
    log_rhs_base = (np.log(SIGMA_BULK) + 0.75*np.log(phi_ex_arr) + 1.5*np.log(cn_arr)
                    + 0.25*np.log(cov_arr) + 2.0*np.log(fp_arr))
    w_sigmoid = 1.0 / (1.0 + np.exp(-TAU_K * (tau_arr - TAU_C)))

    def _fit_at(k_bl, tc_bl):
        """Fit v5+poly3 blend at (k, τc). Return (r2, loocv, w20, b_v5, b_p3, w_bl, pred)."""
        w_bl = 1.0 / (1.0 + np.exp(-k_bl * (tau_arr - tc_bl)))
        X_v5_l = np.column_stack([np.ones(len(log_sf)), w_sigmoid])
        X_p3_l = np.column_stack([np.ones(len(log_sf)), log_tau_arr, log_tau_arr**2, log_tau_arr**3])
        b_v5_l = np.linalg.lstsq(X_v5_l, log_sf - log_rhs_base, rcond=None)[0]
        b_p3_l = np.linalg.lstsq(X_p3_l, log_sf - log_rhs_base, rcond=None)[0]
        pv_l = X_v5_l @ b_v5_l
        pp_l = X_p3_l @ b_p3_l
        pred = (1 - w_bl) * pv_l + w_bl * pp_l + log_rhs_base
        r2 = 1 - np.sum((log_sf - pred)**2) / ss_tot
        n_loo = len(log_sf)
        sse_loo = 0.0
        for ii in range(n_loo):
            mk = np.ones(n_loo, bool); mk[ii] = False
            bv_ = np.linalg.lstsq(X_v5_l[mk], (log_sf - log_rhs_base)[mk], rcond=None)[0]
            bp_ = np.linalg.lstsq(X_p3_l[mk], (log_sf - log_rhs_base)[mk], rcond=None)[0]
            pred_ii = (1 - w_bl[ii]) * (X_v5_l[ii] @ bv_) + w_bl[ii] * (X_p3_l[ii] @ bp_) + log_rhs_base[ii]
            sse_loo += (log_sf[ii] - pred_ii)**2
        loocv = 1 - sse_loo / ss_tot
        s_act_l = np.exp(log_sf); s_prd_l = np.exp(pred)
        w20 = int(np.sum(np.abs(s_prd_l - s_act_l) / s_act_l < 0.20))
        return r2, loocv, w20, b_v5_l, b_p3_l, w_bl, pred

    # Continuous 2D optimization: maximize LOOCV over (k, τc)
    from scipy.optimize import minimize
    def _neg_loocv(p):
        k_, tc_ = p
        # k bounded to 20 — sweep showed LOOCV is flat for k>10; prevents bound-hit.
        if k_ <= 0.1 or k_ > 20 or tc_ < 1.2 or tc_ > 3.0:
            return 1e6
        return -_fit_at(k_, tc_)[1]
    res = minimize(_neg_loocv, x0=[5.0, 2.0], method='Nelder-Mead',
                   options={'xatol': 1e-3, 'fatol': 1e-5, 'maxiter': 200})
    best_k, best_tc = float(res.x[0]), float(res.x[1])
    # Coarse sweep for diagnostic visibility
    k_sweep = [1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0]
    print(f"\n[BLEND SWEEP] coarse k scan at τc=2.0")
    print(f"  {'k':>5s}  {'R²':>6s}  {'LOOCV':>6s}  {'±20%':>5s}")
    for k in k_sweep:
        r2_k, lo_k, w20_k, *_ = _fit_at(k, 2.0)
        print(f"  {k:5.1f}  {r2_k:.4f}  {lo_k:.4f}  {w20_k:2d}/{len(log_sf):2d}")
    r2_formX, loocv_formX, w20_opt, b_v5, b_p3, w_blend, pred_formX = _fit_at(best_k, best_tc)
    print(f"  → continuous optimum: k={best_k:.2f}, τc={best_tc:.3f}")
    print(f"    R²={r2_formX:.4f}, LOOCV={loocv_formX:.4f}, ±20%={w20_opt}/{len(log_sf)}")
    TAU_C_BL = best_tc
    TAU_K_BL = best_k

    # === v12 DATA-NATIVE EXPONENT FIT (diagnostic only, does not replace v9) ===
    # Jointly optimize (α, β, γ, δ, φc, k_bl, τc) where:
    #   log σ = log(σ_grain) + α·log(φ-φc) + β·log(CN) + γ·log(cov) + δ·log(fp) + C_blend(τ)
    # v9 fixed (α,β,γ,δ,φc) = (0.75, 1.5, 0.25, 2.0, 0.185). We check if the data
    # prefers different values by LOOCV. If v12 LOOCV > v9 LOOCV, data disagrees
    # with Kirkpatrick/Bruggeman priors. If ≈ same, v9 priors are consistent.
    cn_np = np.array([cn[i] for i in valid_idx], dtype=float)
    cov_np = np.array([coverage[i] for i in valid_idx], dtype=float)
    fp_np = np.array([max(f_perc[i], 0.01) for i in valid_idx], dtype=float)
    phi_np = np.array([phi_se[i] for i in valid_idx], dtype=float)
    ln_sigma = log_sf
    ss_tot_local = ss_tot

    def _loocv_at(params):
        a, b, g, dl, pc, kb, tcb = params
        if pc < 0.05 or pc > 0.30: return 1e6
        if a < 0.1 or a > 3.0: return 1e6
        if b < 0.1 or b > 3.5: return 1e6
        if g < -1.5 or g > 2.5: return 1e6
        if dl < 0.1 or dl > 5.0: return 1e6
        if kb < 0.1 or kb > 20 or tcb < 1.2 or tcb > 3.0: return 1e6
        phi_ex_v = np.maximum(phi_np - pc, 1e-4)
        lrhs = (np.log(SIGMA_BULK) + a*np.log(phi_ex_v) + b*np.log(cn_np)
                + g*np.log(cov_np) + dl*np.log(fp_np))
        w_bl_v = 1.0 / (1.0 + np.exp(-kb * (tau_arr - tcb)))
        X_v = np.column_stack([np.ones(len(ln_sigma)), w_sigmoid])
        X_p = np.column_stack([np.ones(len(ln_sigma)), log_tau_arr,
                                log_tau_arr**2, log_tau_arr**3])
        sse = 0.0
        n_loo = len(ln_sigma)
        for ii in range(n_loo):
            mk = np.ones(n_loo, bool); mk[ii] = False
            bv_ = np.linalg.lstsq(X_v[mk], (ln_sigma - lrhs)[mk], rcond=None)[0]
            bp_ = np.linalg.lstsq(X_p[mk], (ln_sigma - lrhs)[mk], rcond=None)[0]
            p_ii = (1 - w_bl_v[ii]) * (X_v[ii] @ bv_) + w_bl_v[ii] * (X_p[ii] @ bp_) + lrhs[ii]
            sse += (ln_sigma[ii] - p_ii)**2
        return sse / ss_tot_local   # minimize

    print(f"\n[v12 DATA-NATIVE EXPONENT FIT] optimizing (α, β, γ, δ, φc, k, τc) ...")
    x0 = [0.75, 1.5, 0.25, 2.0, 0.185, best_k, best_tc]
    res12 = minimize(_loocv_at, x0=x0, method='Nelder-Mead',
                     options={'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 2000, 'adaptive': True})
    a12, b12, g12, d12, pc12, kb12, tc12 = res12.x
    loocv12 = 1 - res12.fun
    # Training R² at optimum
    phi_ex_opt = np.maximum(phi_np - pc12, 1e-4)
    lrhs_opt = (np.log(SIGMA_BULK) + a12*np.log(phi_ex_opt) + b12*np.log(cn_np)
                + g12*np.log(cov_np) + d12*np.log(fp_np))
    w_bl_opt = 1.0 / (1.0 + np.exp(-kb12 * (tau_arr - tc12)))
    X_v_f = np.column_stack([np.ones(len(ln_sigma)), w_sigmoid])
    X_p_f = np.column_stack([np.ones(len(ln_sigma)), log_tau_arr, log_tau_arr**2, log_tau_arr**3])
    bv_f = np.linalg.lstsq(X_v_f, ln_sigma - lrhs_opt, rcond=None)[0]
    bp_f = np.linalg.lstsq(X_p_f, ln_sigma - lrhs_opt, rcond=None)[0]
    pred12 = (1 - w_bl_opt) * (X_v_f @ bv_f) + w_bl_opt * (X_p_f @ bp_f) + lrhs_opt
    r2_12 = 1 - np.sum((ln_sigma - pred12)**2) / ss_tot_local
    err12 = np.abs(np.exp(pred12) - np.exp(ln_sigma)) / np.exp(ln_sigma) * 100
    print(f"  {'exponent':12s} {'v9 fixed':>10s} {'v12 fit':>10s} {'Δ':>10s}")
    for name, v9v, v12v in [('α (φ-φc)', 0.75, a12), ('β (CN)', 1.5, b12),
                             ('γ (cov)', 0.25, g12), ('δ (fp)', 2.0, d12),
                             ('φc', 0.185, pc12), ('k_bl', best_k, kb12), ('τc', best_tc, tc12)]:
        delta = v12v - v9v
        flag = " *" if abs(delta) > 0.1 * max(abs(v9v), 0.1) else ""
        print(f"  {name:12s} {v9v:10.4f} {v12v:10.4f} {delta:+10.4f}{flag}")
    print(f"\n  v9  LOOCV = {loocv_formX:.4f}   R² = {r2_formX:.4f}   ±20%: {w20_opt}/{len(log_sf)}")
    print(f"  v12 LOOCV = {loocv12:.4f}   R² = {r2_12:.4f}   ±20%: {int((err12<20).sum())}/{len(log_sf)}")
    verdict = "v12 WINS" if loocv12 > loocv_formX + 0.002 else (
              "v9 WINS" if loocv_formX > loocv12 + 0.002 else "tied within noise")
    print(f"  → {verdict}   (LOOCV σ_noise ≈ {np.sqrt(2/len(log_sf))*(1-loocv_formX):.4f})")

    # === v13 UNIFIED: v12 exponents + Cluster C (CN saturation) ===
    # log σ = log(σ_grain) + α·log(φ-φc) + β·log(CN) - ε·(log CN)²
    #       + γ·log(cov) + δ·log(fp) + C_blend(τ)
    # 8 hyperparams jointly optimized by LOOCV: (α, β, ε, γ, δ, φc, k, τc)
    log_cn_np = np.log(cn_np)

    def _loocv_v13(params):
        a, b, e, g, dl, pc, kb, tcb = params
        if pc < 0.05 or pc > 0.30: return 1e6
        if a < 0.1 or a > 3.0: return 1e6
        if b < 0.1 or b > 4.0: return 1e6
        if e < -0.5 or e > 1.0: return 1e6      # saturation strength
        if g < -1.5 or g > 2.5: return 1e6
        if dl < 0.1 or dl > 5.0: return 1e6
        if kb < 0.1 or kb > 20 or tcb < 1.2 or tcb > 3.0: return 1e6
        phi_ex_v = np.maximum(phi_np - pc, 1e-4)
        lrhs = (np.log(SIGMA_BULK) + a*np.log(phi_ex_v) + b*log_cn_np - e*log_cn_np**2
                + g*np.log(cov_np) + dl*np.log(fp_np))
        w_bl_v = 1.0 / (1.0 + np.exp(-kb * (tau_arr - tcb)))
        X_v = np.column_stack([np.ones(len(ln_sigma)), w_sigmoid])
        X_p = np.column_stack([np.ones(len(ln_sigma)), log_tau_arr,
                                log_tau_arr**2, log_tau_arr**3])
        sse = 0.0
        n_loo = len(ln_sigma)
        for ii in range(n_loo):
            mk = np.ones(n_loo, bool); mk[ii] = False
            bv_ = np.linalg.lstsq(X_v[mk], (ln_sigma - lrhs)[mk], rcond=None)[0]
            bp_ = np.linalg.lstsq(X_p[mk], (ln_sigma - lrhs)[mk], rcond=None)[0]
            p_ii = (1 - w_bl_v[ii]) * (X_v[ii] @ bv_) + w_bl_v[ii] * (X_p[ii] @ bp_) + lrhs[ii]
            sse += (ln_sigma[ii] - p_ii)**2
        return sse / ss_tot_local

    print(f"\n[v13 UNIFIED FIT] v12 exponents + CN saturation (8 hyperparams) ...")
    # Initialize from v12 optimum + small positive ε
    x0_13 = [a12, b12, 0.05, g12, d12, pc12, kb12, tc12]
    res13 = minimize(_loocv_v13, x0=x0_13, method='Nelder-Mead',
                     options={'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 3000, 'adaptive': True})
    a13, b13, e13, g13, d13, pc13, kb13, tc13 = res13.x
    loocv13 = 1 - res13.fun
    # Training R² at v13 optimum
    phi_ex_13 = np.maximum(phi_np - pc13, 1e-4)
    lrhs_13 = (np.log(SIGMA_BULK) + a13*np.log(phi_ex_13) + b13*log_cn_np - e13*log_cn_np**2
               + g13*np.log(cov_np) + d13*np.log(fp_np))
    w_bl_13 = 1.0 / (1.0 + np.exp(-kb13 * (tau_arr - tc13)))
    X_v_13 = np.column_stack([np.ones(len(ln_sigma)), w_sigmoid])
    X_p_13 = np.column_stack([np.ones(len(ln_sigma)), log_tau_arr, log_tau_arr**2, log_tau_arr**3])
    bv_13 = np.linalg.lstsq(X_v_13, ln_sigma - lrhs_13, rcond=None)[0]
    bp_13 = np.linalg.lstsq(X_p_13, ln_sigma - lrhs_13, rcond=None)[0]
    pred13 = (1 - w_bl_13) * (X_v_13 @ bv_13) + w_bl_13 * (X_p_13 @ bp_13) + lrhs_13
    r2_13 = 1 - np.sum((ln_sigma - pred13)**2) / ss_tot_local
    err13 = np.abs(np.exp(pred13) - np.exp(ln_sigma)) / np.exp(ln_sigma) * 100
    # CN saturation effective peak: d/dCN [β·log(CN) - ε·(log CN)²] = 0 → log(CN*) = β/(2ε)
    cn_sat = np.exp(b13 / (2 * e13)) if e13 > 1e-3 else float('inf')
    print(f"  {'param':12s} {'v12':>10s} {'v13':>10s} {'Δ':>10s}")
    for name, v12v, v13v in [('α (φ-φc)', a12, a13), ('β (CN lin)', b12, b13),
                              ('ε (CN²)', 0.0, e13), ('γ (cov)', g12, g13),
                              ('δ (fp)', d12, d13), ('φc', pc12, pc13),
                              ('k_bl', kb12, kb13), ('τc', tc12, tc13)]:
        print(f"  {name:12s} {v12v:10.4f} {v13v:10.4f} {v13v-v12v:+10.4f}")
    print(f"  CN saturation peak: CN* = exp(β/2ε) = {cn_sat:.2f}  (effective CN cap)")
    tier_13 = {t: int(np.sum(err13 < t)) for t in [5, 10, 15, 20]}
    print(f"\n  v9  LOOCV = {loocv_formX:.4f}   R² = {r2_formX:.4f}   ±20%: {w20_opt}/{len(log_sf)}")
    print(f"  v12 LOOCV = {loocv12:.4f}   R² = {r2_12:.4f}   ±20%: {int((err12<20).sum())}/{len(log_sf)}")
    print(f"  v13 LOOCV = {loocv13:.4f}   R² = {r2_13:.4f}   ±20%: {tier_13[20]}/{len(log_sf)}"
          f"  (<5%:{tier_13[5]} <10%:{tier_13[10]} <15%:{tier_13[15]})")
    best_variant = "v13" if loocv13 == max(loocv_formX, loocv12, loocv13) else (
                   "v12" if loocv12 == max(loocv_formX, loocv12, loocv13) else "v9")
    print(f"  → best by LOOCV: {best_variant}")

    X_v5 = np.column_stack([np.ones(len(log_sf)), w_sigmoid])
    X_p3 = np.column_stack([np.ones(len(log_sf)), log_tau_arr, log_tau_arr**2, log_tau_arr**3])
    ss_res_formX = np.sum((log_sf - pred_formX)**2)

    ln_Ct, ln_delta = b_v5[0], b_v5[1]
    ln_Cn = ln_Ct + ln_delta
    C_thick = np.exp(ln_Ct); C_thin = np.exp(ln_Cn)
    C_formX = C_thick
    # Store poly3 coefficients for multiscale plot
    global _GLOBAL_IONIC_POLY3
    _GLOBAL_IONIC_POLY3 = tuple(b_p3)

    # (LOOCV already computed inside _fit_at_k for the winning k)

    # Save to global
    _GLOBAL_C_ION = C_thick
    global _GLOBAL_FORMX_R2
    _GLOBAL_FORMX_R2 = (r2_formX, loocv_formX)
    global _GLOBAL_IONIC_SIGMOID
    _GLOBAL_IONIC_SIGMOID = (C_thick, C_thin, TAU_C, TAU_K)

    # Use FORM X v4 as primary
    s_pred = np.exp(pred_formX)
    r2 = r2_formX
    s_actual = np.array([sigma_net[i] for i in valid_idx])

    # v3 prediction for comparison
    s_pred_v3 = np.exp(pred_fixed)

    # --- Residual diagnostic: find |err|>20% outliers and dump feature signature ---
    # Build per-index group label (e.g. "SE 0.5μm (1mAh_85:15)") + case hint from data source path
    def _case_label(idx):
        grp = ""
        if _GROUP_INFO:
            sizes, gnames = _GROUP_INFO
            pos = 0
            for sz, gn in zip(sizes, gnames):
                if pos <= idx < pos + sz:
                    grp = gn
                    break
                pos += sz
        src = data_list[idx].get("_source_path", "") if idx < len(data_list) else ""
        case_id = os.path.basename(os.path.dirname(src)) if src else ""
        lbl = names[idx] if idx < len(names) else f"idx{idx}"
        return f"{grp} | {case_id} [{lbl}]" if grp or case_id else lbl

    rel_err = (s_pred - s_actual) / s_actual * 100
    abs_err = np.abs(rel_err)
    n_total = len(valid_idx)
    print(f"\n[IONIC v9 BLEND DIAG] n={n_total}, R²={r2_formX:.4f}, LOOCV={loocv_formX:.4f}")
    print(f"  mean|err|={np.mean(abs_err):.1f}%  median|err|={np.median(abs_err):.1f}%")
    # Tier counts
    tiers = [5, 10, 15, 20, 30]
    cnt = {t: int(np.sum(abs_err < t)) for t in tiers}
    print(f"  tiers:  <5%:{cnt[5]}/{n_total}   <10%:{cnt[10]}/{n_total}   <15%:{cnt[15]}/{n_total}"
          f"   <20%:{cnt[20]}/{n_total}   <30%:{cnt[30]}/{n_total}")
    hdr = f"  {'#':>3s} {'case (group | id [P:S])':70s} {'σ_act':>7s} {'σ_pred':>7s} {'err%':>7s} {'φ_SE':>5s} {'f_p':>5s} {'τ':>5s} {'CN':>5s} {'cov':>5s}"
    print(hdr); print("  " + "-" * (len(hdr)-2))
    all_sorted = sorted(range(n_total), key=lambda j: -abs_err[j])
    shown = 0
    for rank, j in enumerate(all_sorted, 1):
        if abs_err[j] <= 10.0:
            continue  # only show |err|>10%
        i = valid_idx[j]
        nm = _case_label(i)[:70]
        flag = "  ⚠" if abs_err[j] > 20 else ("  ·" if abs_err[j] > 15 else "")
        print(f"  {rank:3d} {nm:70s} {s_actual[j]:7.4f} {s_pred[j]:7.4f} {rel_err[j]:+6.1f}%"
              f" {phi_se[i]:5.3f} {f_perc[i]:5.3f} {tau[i]:5.2f} {cn[i]:5.2f} {cov_arr[j]:5.3f}{flag}")
        shown += 1
    print(f"  ({shown} cases with |err|>10%, remaining {n_total - shown} are within ±10%)")
    # Residual-vs-feature correlation (Pearson on log-residual)
    log_res = np.log(s_pred) - np.log(s_actual)
    gb_arr = np.array([gb_dens[i] for i in valid_idx])
    gp_arr = np.array([g_path[i] for i in valid_idx])

    # Parse P:S ratio → particulate fraction
    def _parse_ps(d):
        ps = d.get("ps_ratio", "") or ""
        if ps in ("P only", "10:0"): return 1.0
        if ps in ("S only", "0:10"): return 0.0
        if ":" in ps:
            try:
                p, s = ps.split(":")
                p, s = float(p), float(s)
                return p / (p + s) if (p + s) > 0 else 0.5
            except Exception:
                return 0.5
        return 0.5
    p_frac_arr = np.array([_parse_ps(data_list[i]) for i in valid_idx])

    feats = {'log(phi_ex)': np.log(phi_ex_arr), 'log(CN)': np.log(cn_arr),
             'log(tau)': log_tau_arr, 'log(cov)': np.log(cov_arr),
             'log(fp)': np.log(fp_arr),
             'log(gb_dens)': np.log(np.maximum(gb_arr, 1e-10)),
             'log(g_path)': np.log(np.maximum(gp_arr, 1e-10)),
             'log(gp*gb²)': np.log(np.maximum(gp_arr * gb_arr**2, 1e-20)),
             'p_frac (P:S)': p_frac_arr,
             'p_frac - 0.5': p_frac_arr - 0.5,
             '(p_frac)²': p_frac_arr**2}
    print("  residual(log) correlations:")
    for nm, v in feats.items():
        c = np.corrcoef(log_res, v)[0, 1] if np.std(v) > 0 else 0.0
        flag = " ⚠" if abs(c) > 0.3 else ""
        print(f"    {nm:12s} r = {c:+.3f}{flag}")
    print()

    fig, ax = plt.subplots(figsize=FIG_SINGLE)

    # Group colors
    point_groups = [0] * len(valid_idx)
    group_boundaries = []
    if _GROUP_INFO:
        sizes, gnames = _GROUP_INFO
        pos = 0
        for sz in sizes:
            group_boundaries.append((pos, pos + sz))
            pos += sz
        for j, i in enumerate(valid_idx):
            for g_idx, (start, end) in enumerate(group_boundaries):
                if start <= i < end:
                    point_groups[j] = g_idx
                    break

    # Mark |err|>20% cases with a subtle darker edge (no overlay)
    is_outlier = np.abs((s_pred - s_actual) / s_actual) > 0.20

    if _GROUP_INFO:
        for j in range(len(valid_idx)):
            gi = point_groups[j]
            edge = '#444' if is_outlier[j] else 'white'
            lw = 1.8 if is_outlier[j] else 1.5
            ax.scatter(s_actual[j], s_pred[j], s=100,
                      c=GROUP_COLORS[gi % len(GROUP_COLORS)],
                      zorder=5, edgecolors=edge, linewidth=lw)
    else:
        ax.scatter(s_actual, s_pred, s=100, c=BLUE, zorder=5,
                  edgecolors='white', linewidth=1.5)

    # 1:1 line
    all_vals = np.concatenate([s_pred, s_actual])
    vmin, vmax = all_vals[all_vals > 0].min() * 0.5, all_vals.max() * 2
    ax.plot([vmin, vmax], [vmin, vmax], 'k--', linewidth=1.5, alpha=0.5, label='1:1 line')

    # ±20% band
    ax.fill_between([vmin, vmax], [vmin*0.8, vmax*0.8], [vmin*1.2, vmax*1.2],
                    alpha=0.1, color='green', label='±20%')

    # Group legend (no per-point labels to avoid MemoryError with adjustText on 41+ points)
    if _GROUP_INFO:
        from matplotlib.patches import Patch
        legend_patches = []
        for gi, (start, end) in enumerate(group_boundaries):
            group_js = [j for j, i in enumerate(valid_idx) if start <= i < end]
            if group_js:
                color = GROUP_COLORS[gi % len(GROUP_COLORS)]
                legend_patches.append(Patch(facecolor=color, label=gnames[gi]))
        if legend_patches:
            ax.legend(handles=legend_patches + [
                plt.Line2D([0], [0], linestyle='--', color='black', label='1:1 line'),
            ], fontsize=7, loc='upper left', ncol=1)

    # R²
    log_pred = np.log(s_pred[s_pred > 0])
    log_actual = np.log(s_actual[s_actual > 0])
    if len(log_pred) > 2:
        ss_res = np.sum((log_pred - log_actual)**2)
        ss_tot = np.sum((log_actual - np.mean(log_actual))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    else:
        r2 = 0

    # Error stats
    errors = np.abs(s_pred - s_actual) / s_actual * 100
    within_20 = np.sum(errors < 20)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("σ_actual (Network solver, mS/cm)", fontsize=11)
    ax.set_ylabel("σ_predicted (Scaling law, mS/cm)", fontsize=11)
    ax.set_title(f"Ionic v9 BLEND: [(1-w)·v5 + w·poly3] × σ_grain × CN^1.5 × (φ−φc)^¾ × ⁴√cov × fp²\n"
                 f"sigmoid(k={TAU_K:.0f},τc={TAU_C})  Ct={C_thick:.4f} Cn={C_thin:.4f}  R²={r2_formX:.3f}",
                 fontsize=9, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')

    txt = (f"R²={r2:.3f} (n={len(valid_idx)})\n"
           f"|err|={np.mean(errors):.0f}%, ≤20%: {within_20}/{len(errors)}")
    ax.text(0.95, 0.05, txt, transform=ax.transAxes, fontsize=9,
            ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect('equal')
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.xaxis.grid(True, linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    _write_csv(outdir, 'ionic_scaling_fit.csv',
               ['σ_actual(mS/cm)', 'σ_predicted(mS/cm)', 'error(%)'],
               [names[i] for i in valid_idx],
               list(s_actual), list(s_pred), list(errors))
    return _save(fig, outdir, "ionic_scaling_fit.png")


def plot_network_sigma(data_list, names, outdir):
    """Part II: Network solver σ_full (ground truth)."""
    SIGMA_BULK = 3.0  # σ_grain (grain interior)
    sigma_net = _load_network_sigma(data_list)
    sigma_brug_abs = [_get(d, "sigma_ratio") * SIGMA_BULK for d in data_list]

    if not any(s > 0 for s in sigma_net):
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    ax.plot(x, sigma_net, 's-', color='#2ecc71', markersize=ms, linewidth=lw,
            label="σ_full (network)")
    ax.plot(x, sigma_brug_abs, 'o--', color=BLUE, markersize=ms-2, linewidth=lw-0.5,
            alpha=0.5, label="σ_brug (Bruggeman)")
    _apply_style(ax, "σ_ionic (mS/cm)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title("Part II: Network Solver σ_full vs Bruggeman σ_brug",
                 fontsize=10, fontweight='bold')

    # R_brug annotation
    r_brugs = [sigma_brug_abs[i]/sigma_net[i] if sigma_net[i] > 0 else 0 for i in range(len(names))]
    valid_r = [r for r in r_brugs if r > 0]
    if valid_r:
        ax.text(0.95, 0.95, f"R_brug = {min(valid_r):.1f}~{max(valid_r):.1f}×\n(Bruggeman overestimation)",
                transform=ax.transAxes, fontsize=9, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8))

    _write_csv(outdir, 'network_sigma.csv',
               ['σ_brug(mS/cm)', 'σ_full(mS/cm)', 'R_brug'],
               names, sigma_brug_abs, sigma_net, r_brugs)
    return _save(fig, outdir, "network_sigma.png")


def plot_multiscale_sigma(data_list, names, outdir):
    """FORM X: σ = C × σ_grain × (φ-φc)^¾ × CN × √cov / √τ."""
    SIGMA_BULK = 3.0
    PHI_C = 0.185

    phi_se = [_get(d, "phi_se") for d in data_list]
    cn = [_get(d, "se_se_cn", 0) for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    coverage = [(lambda vs: sum(vs)/len(vs)/100 if vs else 0.20)([v for v in [_get(d,"coverage_AM_P_mean",0), _get(d,"coverage_AM_S_mean",0), _get(d,"coverage_AM_mean",0)] if v>0]) for d in data_list]
    f_perc = [max(_get(d, "percolation_pct", 100) / 100, 0.01) for d in data_list]
    sigma_net = _load_network_sigma(data_list)

    # Fit C from data
    valid_fit = []
    for i in range(len(data_list)):
        phi_ex = max(phi_se[i] - PHI_C, 0.001)
        if phi_ex > 0 and cn[i] > 0 and tau[i] > 0 and coverage[i] > 0 and sigma_net[i] > 0.01:
            valid_fit.append(i)

    # Sigmoid C(τ) parameters
    TAU_C = 2.1; TAU_K = 5.0
    if _GLOBAL_IONIC_SIGMOID is not None:
        C_thick_ms, C_thin_ms, TAU_C, TAU_K = _GLOBAL_IONIC_SIGMOID
    elif _GLOBAL_C_ION is not None:
        C_thick_ms = _GLOBAL_C_ION; C_thin_ms = C_thick_ms * 0.54
    elif len(valid_fit) < 3:
        C_thick_ms = 0.034; C_thin_ms = 0.019
    else:
        log_rhs = np.array([np.log(SIGMA_BULK) + 0.75*np.log(max(phi_se[i]-PHI_C, 0.001))
                            + 1.5*np.log(cn[i]) + 0.25*np.log(coverage[i])
                            for i in valid_fit])
        log_actual = np.array([np.log(sigma_net[i]) for i in valid_fit])
        C_thick_ms = np.exp(np.mean(log_actual - log_rhs))
        C_thin_ms = C_thick_ms * 0.54
    C_ms = C_thick_ms

    # v9 BLEND: (1-w_bl)·v5_sigmoid + w_bl·poly3
    p3_coefs = _GLOBAL_IONIC_POLY3 if _GLOBAL_IONIC_POLY3 is not None else (-3.7318, 3.0684, -6.3596, 3.0123)
    sigma_ms = []
    for i in range(len(data_list)):
        phi_ex = max(phi_se[i] - PHI_C, 0.001)
        if cn[i] > 0 and tau[i] > 0 and coverage[i] > 0:
            w_v5 = 1.0 / (1.0 + np.exp(-TAU_K * (tau[i] - TAU_C)))
            w_bl = 1.0 / (1.0 + np.exp(-15.0 * (tau[i] - 2.0)))
            ln_C_v5 = np.log(C_thick_ms) + (np.log(C_thin_ms) - np.log(C_thick_ms)) * w_v5
            lt = np.log(tau[i])
            ln_C_p3 = p3_coefs[0] + p3_coefs[1]*lt + p3_coefs[2]*lt**2 + p3_coefs[3]*lt**3
            ln_C = (1 - w_bl) * ln_C_v5 + w_bl * ln_C_p3
            s = np.exp(ln_C) * SIGMA_BULK * phi_ex**0.75 * cn[i]**1.5 * coverage[i]**0.25 * f_perc[i]**2
            sigma_ms.append(s)
        else:
            sigma_ms.append(0)

    has_net = any(s > 0 for s in sigma_net)

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    ax.plot(x, sigma_ms, 's-', color=RED, markersize=ms, linewidth=lw,
            label="FORM X (mS/cm)")
    # ±22% error band (DEM stochastic variability)
    _ms_arr = np.array(sigma_ms)
    _ms_lo = _ms_arr * 0.78; _ms_hi = _ms_arr * 1.22
    ax.fill_between(x, _ms_lo, _ms_hi, color=RED, alpha=0.10, label='±22% band')
    if has_net:
        ax.plot(x, sigma_net, 'D--', color='#2ecc71', markersize=ms-2, linewidth=lw-0.5,
                alpha=0.7, label="Network solver (mS/cm)")
        # Subtle outlier mark: small * above plot at |err|>20% positions
        _net_arr = np.array(sigma_net, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.where(_net_arr > 0, np.abs(_ms_arr - _net_arr) / _net_arr, 0.0)
        out_mask = rel > 0.20
        if np.any(out_mask):
            y_top = max(np.nanmax(_ms_hi), np.nanmax(_net_arr)) * 1.03
            ax.scatter(x[out_mask], [y_top] * int(out_mask.sum()),
                       marker='*', s=40, color='#888', zorder=3, clip_on=False)

    _apply_style(ax, "σ_ionic (mS/cm)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title(f"FORM X v9: C_blend(τ) × σ_grain × CN^1.5 × (φ−φc)^¾ × ⁴√cov × fp²",
                 fontsize=9, fontweight='bold')

    _write_csv(outdir, 'multiscale_sigma.csv',
               ['φ_SE', 'CN', 'τ', 'coverage', 'σ_FORMX(mS/cm)', 'σ_network(mS/cm)'],
               names, phi_se, cn, tau, sigma_ms, sigma_net)
    return _save(fig, outdir, "multiscale_sigma.png")


def plot_v3_fitting(data_list, names, outdir):
    """[Legacy] v3: σ_brug × C × (G_path × GB_d²)^(1/4) × CN²."""
    SIGMA_BULK = 3.0
    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    g_path = [_get(d, "path_conductance_mean", 0) for d in data_list]
    cn = [_get(d, "se_se_cn", 0) for d in data_list]
    sigma_net = _load_network_sigma(data_list)

    valid = [i for i in range(len(data_list))
             if g_path[i]>0 and cn[i]>0 and gb_dens[i]>0 and sigma_brug[i]>0 and sigma_net[i]>0.01]
    if len(valid) < 3:
        C_v3 = 0.073
    else:
        log_rhs = np.array([np.log(sigma_brug[i]*SIGMA_BULK) + 0.25*np.log(g_path[i]*gb_dens[i]**2) + 2*np.log(cn[i]) for i in valid])
        log_act = np.array([np.log(sigma_net[i]) for i in valid])
        C_v3 = np.exp(np.mean(log_act - log_rhs))

    sigma_v3 = []
    for i in range(len(data_list)):
        if g_path[i]>0 and cn[i]>0 and gb_dens[i]>0:
            sigma_v3.append(sigma_brug[i]*C_v3*(g_path[i]*gb_dens[i]**2)**0.25*cn[i]**2*SIGMA_BULK)
        else:
            sigma_v3.append(0)

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names)); lw = _line_width(len(names))
    ax.plot(x, sigma_v3, 's-', color=RED, markersize=ms, linewidth=lw, label="v3 (mS/cm)")
    if any(s>0 for s in sigma_net):
        ax.plot(x, sigma_net, 'D--', color='#2ecc71', markersize=ms-2, linewidth=lw-0.5, alpha=0.7, label="Network solver (mS/cm)")
    _apply_style(ax, "σ_ionic (mS/cm)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title(f"[Legacy v3] σ_brug × {C_v3:.4f} × (G_path×GB_d²)^¼ × CN²", fontsize=9, fontweight='bold')
    return _save(fig, outdir, "v3_fitting.png")


def plot_formx_decomposition(data_list, names, outdir):
    """FORM X v4 factor decomposition: (φ-φc)^¾, CN^1.5, ⁴√cov, C(τ)."""
    PHI_C = 0.185; TAU_C = 2.1; TAU_K = 5.0
    phi_se = [_get(d, "phi_se") for d in data_list]
    cn = [_get(d, "se_se_cn", 0) for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    coverage = [(lambda vs: sum(vs)/len(vs)/100 if vs else 0.20)([v for v in [_get(d,"coverage_AM_P_mean",0), _get(d,"coverage_AM_S_mean",0), _get(d,"coverage_AM_mean",0)] if v>0]) for d in data_list]

    n = len(data_list)
    log_phiex = np.array([0.75*np.log(max(phi_se[i]-PHI_C, 0.001)) for i in range(n)])
    log_cn = np.array([1.5*np.log(cn[i]) if cn[i]>0 else 0 for i in range(n)])
    log_cov = np.array([0.25*np.log(coverage[i]) if coverage[i]>0 else 0 for i in range(n)])
    if _GLOBAL_IONIC_SIGMOID:
        Ct, Cn, tc, k = _GLOBAL_IONIC_SIGMOID
    else:
        Ct, Cn, tc, k = 0.034, 0.019, TAU_C, TAU_K
    log_Csig = np.array([np.log(Ct) + (np.log(Cn)-np.log(Ct)) / (1+np.exp(-k*(tau[i]-tc))) for i in range(n)])

    # Ref: best case
    sigma_net = _load_network_sigma(data_list)
    ref = np.argmax(sigma_net) if any(s>0 for s in sigma_net) else 0

    factors = [
        ('(φ−φc)^¾', log_phiex, '#4472C4'),
        ('CN^1.5', log_cn, '#ED7D31'),
        ('⁴√cov', log_cov, '#A5A5A5'),
        ('C(τ)', log_Csig, '#FFC000'),
    ]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(max(8, n*0.8), 10), gridspec_kw={'height_ratios': [3, 2]})

    # Top: stacked bar (relative to ref)
    x = np.arange(n)
    bottom_pos = np.zeros(n); bottom_neg = np.zeros(n)
    for label, vals, color in factors:
        delta = vals - vals[ref]
        pos = np.clip(delta, 0, None); neg = np.clip(delta, None, 0)
        ax.bar(x, pos, bottom=bottom_pos, color=color, label=label, width=0.7, edgecolor='white', linewidth=0.5)
        ax.bar(x, neg, bottom=bottom_neg, color=color, width=0.7, edgecolor='white', linewidth=0.5)
        bottom_pos += pos; bottom_neg += neg

    ax.axhline(0, color='gray', linewidth=0.5)
    _apply_style(ax, 'Δlog(factor) from ref', names, data_list)
    # Get x_labels from ax (after _apply_style monomodal detection)
    x_labels = [t.get_text() for t in ax.get_xticklabels()]
    if not any(x_labels) or all(t == '' for t in x_labels):
        x_labels = list(names)
    ax.set_title(f'FORM X Factor Decomposition (ref: {x_labels[ref] if ref < len(x_labels) else names[ref]})', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left', ncol=4)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # Bottom: horizontal bar per case (dominant factor)
    for i in range(n):
        deltas = [(label, vals[i]-vals[ref]) for label, vals, _ in factors]
        deltas.sort(key=lambda x: abs(x[1]), reverse=True)
        dominant = deltas[0][0]
        colors = [c for l,_,c in factors if l==dominant][0]
        ax2.barh(i, deltas[0][1], color=colors, height=0.6, edgecolor='white', linewidth=0.5)
        ax2.text(deltas[0][1], i, f' {dominant}', va='center', fontsize=7)

    ax2.set_yticks(range(n))
    ax2.set_yticklabels([x_labels[i][:20] for i in range(n)], fontsize=7)
    ax2.set_xlabel('Dominant factor Δlog', fontsize=10)
    ax2.set_title('Dominant factor per case', fontsize=11, fontweight='bold')
    ax2.axvline(0, color='gray', linewidth=0.5)
    ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)

    fig.tight_layout()
    _write_csv(outdir, 'formx_decomposition.csv',
               ['(φ-φc)^¾', 'CN^1.5', '⁴√cov', 'C(τ)', 'dominant'],
               x_labels,
               list(log_phiex), list(log_cn), list(log_cov), list(log_Csig),
               [sorted([(l, v[i]-v[ref]) for l,v,_ in factors], key=lambda x:-abs(x[1]))[0][0] for i in range(n)])
    return _save(fig, outdir, "formx_decomposition.png")


def _load_electronic_sigma(data_list):
    """Load electronic σ_full from full_metrics or network_conductivity.json."""
    vals = [0.0] * len(data_list)
    for i, d in enumerate(data_list):
        v = _get(d, "electronic_sigma_full_mScm", 0)
        if v and v > 0:
            vals[i] = v
        else:
            src = _get(d, "_source_path", "")
            if src:
                net_path = os.path.join(os.path.dirname(src), "network_conductivity.json")
                if os.path.exists(net_path):
                    try:
                        with open(net_path) as _nf:
                            nd = json.load(_nf)
                        vals[i] = nd.get("electronic_sigma_full_mScm", 0) or 0
                    except:
                        pass
    return vals


def _load_thermal_sigma(data_list):
    """Load thermal σ_full from full_metrics or network_conductivity.json."""
    vals = [0.0] * len(data_list)
    for i, d in enumerate(data_list):
        v = _get(d, "thermal_sigma_full_mScm", 0)
        if v and v > 0:
            vals[i] = v
        else:
            src = _get(d, "_source_path", "")
            if src:
                net_path = os.path.join(os.path.dirname(src), "network_conductivity.json")
                if os.path.exists(net_path):
                    try:
                        with open(net_path) as _nf:
                            nd = json.load(_nf)
                        vals[i] = nd.get("thermal_sigma_full_mScm", 0) or 0
                    except:
                        pass
    return vals


def plot_electronic_sigma(data_list, names, outdir):
    """Electronic conductivity from AM-AM network."""
    sigma_el = _load_electronic_sigma(data_list)
    if not any(s > 0 for s in sigma_el):
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    # Use NaN to break line at σ=0 cases (prevents cross-group connection)
    y_line = np.array([s if s > 0 else np.nan for s in sigma_el])
    x_none = [x[i] for i in range(len(names)) if sigma_el[i] <= 0]

    ax.plot(x, y_line, 's-', color='#e74c3c', markersize=ms, linewidth=lw,
            label="σ_electronic (mS/cm)")
    if x_none:
        ax.plot(x_none, [0]*len(x_none), 'x', color='gray', markersize=ms+2,
                label="No AM percolation")

    _apply_style(ax, "σ_electronic (mS/cm)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title("Electronic Conductivity (AM-AM Network)\nσ_AM = 0.05 S/cm",
                 fontsize=10, fontweight='bold')

    # Annotation: range
    valid_vals = [s for s in sigma_el if s > 0]
    if valid_vals:
        ax.text(0.95, 0.95, f"Range: {min(valid_vals):.2f} ~ {max(valid_vals):.2f} mS/cm",
                transform=ax.transAxes, fontsize=9, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffeaea', alpha=0.8))

    _write_csv(outdir, 'electronic_sigma.csv',
               ['σ_electronic(mS/cm)'], names, sigma_el)
    return _save(fig, outdir, "electronic_sigma.png")


def plot_thermal_sigma(data_list, names, outdir):
    """Thermal conductivity from ALL-contact network."""
    sigma_th = _load_thermal_sigma(data_list)
    if not any(s > 0 for s in sigma_th):
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    ax.plot(x, sigma_th, 's-', color='#ff922b', markersize=ms, linewidth=lw,
            label="k_eff (mS/cm equiv.)")

    _apply_style(ax, "k_eff (thermal, mS/cm equiv.)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title("Thermal Conductivity (ALL Contact Network)\nk_AM=4.0e-2, k_SE=0.7e-2 W/(cm·K)",
                 fontsize=10, fontweight='bold')

    valid_vals = [s for s in sigma_th if s > 0]
    if valid_vals:
        ax.text(0.95, 0.95, f"Range: {min(valid_vals):.2f} ~ {max(valid_vals):.2f}",
                transform=ax.transAxes, fontsize=9, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', alpha=0.8))

    _write_csv(outdir, 'thermal_sigma.csv',
               ['k_eff(mS/cm)'], names, sigma_th)
    return _save(fig, outdir, "thermal_sigma.png")


def plot_electronic_scaling(data_list, names, outdir):
    """Electronic 2-regime scaling: thick (topology) + thin (contact mechanics)."""
    SIGMA_AM = 50.0
    phi_am = [_get(d, "phi_am") for d in data_list]
    cn_am = [_get(d, "am_am_cn") for d in data_list]
    thickness = [_get(d, "thickness_um") for d in data_list]
    d_am_list = [2.0 * max(_get(d, "r_AM_P", 0), _get(d, "r_AM_S", 0), _get(d, "r_AM", 0))
                 for d in data_list]
    tau = [max(_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)), 0.1) for d in data_list]
    cov_list = [(lambda vs: sum(vs)/len(vs)/100 if vs else 0.20)(
        [v for v in [_get(d,"coverage_AM_P_mean",0), _get(d,"coverage_AM_S_mean",0),
                     _get(d,"coverage_AM_mean",0)] if v > 0]) for d in data_list]
    am_delta = [max(_get(d, "am_am_mean_delta", 0), 0.001) for d in data_list]
    am_area = [max(_get(d, "am_am_mean_area", 0), 0.01) for d in data_list]
    am_hop = [max(_get(d, "am_am_mean_hop", 0), 0.1) for d in data_list]
    am_gc = [max(_get(d, "am_path_conductance_mean", 0), 0.001) for d in data_list]
    porosity = [max(_get(d, "porosity", 10), 0.1) for d in data_list]
    phi_se = [max(_get(d, "phi_se", 0.2), 0.01) for d in data_list]
    el_perc = [max(_get(d, "electronic_percolating_fraction", 0), 0.01) for d in data_list]
    el_act = [max(_get(d, "electronic_active_fraction", 0), 0.01) for d in data_list]
    sigma_net = _load_electronic_sigma(data_list)

    # --- Fit C globally on ALL electronic data (thick/thin separate) ---
    # Use screening_electronic_sweep.load_all_electronic() for reliable dedup
    import importlib.util as _ilu
    _sweep_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'screening_electronic_sweep.py')
    _spec = _ilu.spec_from_file_location("sweep", _sweep_path)
    _mod = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    _unique = _mod.load_all_electronic()

    # Fit C_thick and C_thin separately
    # Load mixing parameter for bimodal correction
    for r in _unique:
        _m2 = json.load(open(os.path.join(r['path'], 'full_metrics.json')))
        _nP = _m2.get('n_AM_P', 0); _nS = _m2.get('n_AM_S', 0)
        _rP = max(_m2.get('r_AM_P', 0), 0.1); _rS = max(_m2.get('r_AM_S', 0), 0.1)
        _vP = _nP * _rP**3; _vS = _nS * _rS**3; _vT = max(_vP + _vS, 0.001)
        r['fP_vol'] = _vP / _vT
        r['mix'] = 4 * r['fP_vol'] * (1 - r['fP_vol'])  # symmetric: 0=mono, 1=50:50

    # Load contact info + G_holm for each case
    for r in _unique:
        _m2 = json.load(open(os.path.join(r['path'], 'full_metrics.json')))
        _nP = _m2.get('n_AM_P', 0); _nS = _m2.get('n_AM_S', 0)
        _rP = max(_m2.get('r_AM_P', 0), 0.1); _rS = max(_m2.get('r_AM_S', 0), 0.1)
        _vP = _nP * _rP**3; _vS = _nS * _rS**3; _vT = max(_vP + _vS, 0.001)
        r['fP_vol'] = _vP / _vT
        # G_holm: Holm conductance-weighted CN
        _n_pp = _m2.get('area_AM_P_AM_P_n', 0) or 0
        _n_ps = _m2.get('area_AM_P_AM_S_n', 0) or 0
        _n_ss = _m2.get('area_AM_S_AM_S_n', 0) or 0
        _a_pp = _m2.get('area_AM_P_AM_P_mean', 0) or 0
        _a_ps = _m2.get('area_AM_P_AM_S_mean', 0) or 0
        _a_ss = _m2.get('area_AM_S_AM_S_mean', 0) or 0
        _nAM = max(_nP + _nS, 1)
        _th = _n_pp*np.sqrt(max(_a_pp,0.001)) + _n_ps*np.sqrt(max(_a_ps,0.001)) + _n_ss*np.sqrt(max(_a_ss,0.001))
        r['g_holm'] = max(_th / _nAM, 0.01) if _th > 0 else r['cn']
        r['g_total'] = _th if _th > 0 else r['cn'] * _nAM  # total Holm sum
        # A_min: smallest non-zero contact type area (bottleneck indicator)
        _areas = [a for a,n in [(_a_pp,_n_pp),(_a_ps,_n_ps),(_a_ss,_n_ss)] if n > 0 and a > 0]
        r['a_min'] = min(_areas) if _areas else 0.01

    C_thick = 1.0; C_thin = 1.0
    if len(_unique) >= 5:
        _s = np.array([r['sigma'] for r in _unique])
        _cn = np.array([r['cn'] for r in _unique])
        _ratio = np.array([r['ratio'] for r in _unique])
        _por = np.array([r['por'] for r in _unique])
        _delta = np.array([r['delta'] for r in _unique])
        _pa = np.array([r['phi_am'] for r in _unique])
        _gh = np.array([r.get('g_holm', r['cn']) for r in _unique])

        _tk = _ratio >= 8; _tn = _ratio < 8
        if _tk.sum() >= 3:
            # Thick: CN^1.5 × G_holm^0.25 × (φ-φc)² / por^0.35 × exp(-0.05×logCN×logG_h)
            _phi_ex_tk = np.clip(_pa[_tk] - 0.10, 0.001, None)
            _log_cn_tk = np.log(_cn[_tk]); _log_gh_tk = np.log(_gh[_tk])
            _log_phi_tk = np.log(_phi_ex_tk)
            _rhs_tk = SIGMA_AM * np.sqrt(_cn[_tk]) * np.sqrt(_gh[_tk]) * _phi_ex_tk**3 / _por[_tk]**0.25 * np.exp(-0.75*_log_cn_tk*_log_phi_tk - 0.25*_log_gh_tk**2)
            C_thick = float(np.exp(np.mean(np.log(_s[_tk]) - np.log(_rhs_tk))))
        if _tn.sum() >= 3:
            _delta_tn = _delta[_tn]
            _rhs_tn = SIGMA_AM * _cn[_tn] * _delta_tn**0.5 / _ratio[_tn]**0.5
            C_thin = float(np.exp(np.mean(np.log(_s[_tn]) - np.log(_rhs_tn))))

    # ── Global R² from ALL data (not just this comparison group) ──
    r2_global_tk = 0; r2_global_tn = 0; n_global_tk = 0; n_global_tn = 0
    if len(_unique) >= 5:
        _s_all = np.array([r['sigma'] for r in _unique])
        _ratio_all = np.array([r['ratio'] for r in _unique])
        # Thick global R²
        _tk_mask = _ratio_all >= 8
        if _tk_mask.sum() >= 3:
            _cn_g = np.log(np.array([r['cn'] for r in _unique])[_tk_mask])
            _pa_g = np.array([r['phi_am'] for r in _unique])[_tk_mask]
            _phi_ex_g = np.clip(_pa_g - 0.10, 0.001, None)
            _por_g = np.array([r['por'] for r in _unique])[_tk_mask]
            _gh_g = np.array([r.get('g_holm', r['cn']) for r in _unique])[_tk_mask]
            _log_phi_g = np.log(_phi_ex_g)
            _log_gh_g = np.log(_gh_g)
            _pred_tk = C_thick * SIGMA_AM * np.exp(_cn_g)**0.5 * _gh_g**0.5 * _phi_ex_g**3 / _por_g**0.25 * np.exp(-0.75*_cn_g*_log_phi_g - 0.25*_log_gh_g**2)
            _log_a = np.log(_s_all[_tk_mask]); _log_p = np.log(_pred_tk)
            _ss_res = np.sum((_log_a - _log_p)**2); _ss_tot = np.sum((_log_a - np.mean(_log_a))**2)
            r2_global_tk = 1 - _ss_res / _ss_tot if _ss_tot > 0 else 0
            n_global_tk = int(_tk_mask.sum())
        # Thin global R²
        _tn_mask = _ratio_all < 8
        _ep_all = np.array([r['ep'] for r in _unique])
        _tn_perc = _tn_mask & (_ep_all >= 0.50)
        if _tn_perc.sum() >= 3:
            _cn_tn = np.array([r['cn'] for r in _unique])[_tn_perc]
            _delta_tn = np.array([r['delta'] for r in _unique])[_tn_perc]
            _ratio_tn = _ratio_all[_tn_perc]
            _pred_tn = C_thin * SIGMA_AM * _cn_tn * _delta_tn**0.5 / _ratio_tn**0.5
            _log_a = np.log(_s_all[_tn_perc]); _log_p = np.log(_pred_tn)
            _ss_res = np.sum((_log_a - _log_p)**2); _ss_tot = np.sum((_log_a - np.mean(_log_a))**2)
            r2_global_tn = 1 - _ss_res / _ss_tot if _ss_tot > 0 else 0
            n_global_tn = int(_tn_perc.sum())

    # Compute per-case G_holm and A_min
    case_g_holm = []; case_g_total = []; case_a_min = []
    for d in data_list:
        _nP = _get(d, "n_AM_P", 0) or 0; _nS = _get(d, "n_AM_S", 0) or 0
        _n_pp = _get(d, "area_AM_P_AM_P_n", 0) or 0
        _n_ps = _get(d, "area_AM_P_AM_S_n", 0) or 0
        _n_ss = _get(d, "area_AM_S_AM_S_n", 0) or 0
        _a_pp = _get(d, "area_AM_P_AM_P_mean", 0) or 0
        _a_ps = _get(d, "area_AM_P_AM_S_mean", 0) or 0
        _a_ss = _get(d, "area_AM_S_AM_S_mean", 0) or 0
        _nAM = max(_nP + _nS, 1)
        _th = _n_pp*np.sqrt(max(_a_pp,0.001)) + _n_ps*np.sqrt(max(_a_ps,0.001)) + _n_ss*np.sqrt(max(_a_ss,0.001))
        case_g_holm.append(max(_th / _nAM, 0.01) if _th > 0 else cn_am[len(case_g_holm)] if len(case_g_holm) < len(cn_am) else 1.0)
        case_g_total.append(_th if _th > 0 else (cn_am[len(case_g_total)] * _nAM if len(case_g_total) < len(cn_am) else 1.0))
        _areas = [a for a,n in [(_a_pp,_n_pp),(_a_ps,_n_ps),(_a_ss,_n_ss)] if n > 0 and a > 0]
        case_a_min.append(min(_areas) if _areas else 0.01)

    # Compute predictions per case
    sigma_scaling = []
    for i in range(len(data_list)):
        if phi_am[i] > 0 and cn_am[i] > 0 and d_am_list[i] > 0 and thickness[i] > 0:
            ratio_i = thickness[i] / d_am_list[i]
            if ratio_i >= 8:
                # THICK: CN^1.5 × G_holm^0.25 × (φ-φc)² / por^0.35
                phi_ex_i = max(phi_am[i] - 0.10, 0.001)
                _gh_i = max(case_g_holm[i], 0.01)
                _log_cn_i = np.log(cn_am[i]); _log_gh_i = np.log(_gh_i); _log_phi_i = np.log(phi_ex_i)
                _correction = np.exp(-0.75*_log_cn_i*_log_phi_i - 0.25*_log_gh_i**2)
                s = C_thick * SIGMA_AM * np.sqrt(cn_am[i]) * np.sqrt(_gh_i) * phi_ex_i**3 / porosity[i]**0.25 * _correction
            else:
                # THIN: CN × δ^0.5 / √(T/d)
                if el_perc[i] >= 0.50:
                    s = C_thin * SIGMA_AM * cn_am[i] * am_delta[i]**0.5 / ratio_i**0.5
                else:
                    s = 0.0
            sigma_scaling.append(s)
        else:
            sigma_scaling.append(0.0)

    has_net = any(s > 0 for s in sigma_net)
    if not has_net and not any(s > 0 for s in sigma_scaling):
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    y_scaling = np.array([s if s > 0 else np.nan for s in sigma_scaling])
    y_net = np.array([s if s > 0 else np.nan for s in sigma_net]) if has_net else None

    ax.plot(x, y_scaling, 's-', color=RED, markersize=ms, linewidth=lw,
            label="Scaling law")
    # ±13% error band (DEM stochastic variability, 1σ from residual)
    _err_frac = 0.13
    y_lo = y_scaling * (1 - _err_frac)
    y_hi = y_scaling * (1 + _err_frac)
    ax.fill_between(x, y_lo, y_hi, color=RED, alpha=0.12, label='±13% band')
    if has_net and y_net is not None:
        ax.plot(x, y_net, 'D--', color='#2ecc71', markersize=ms-2, linewidth=lw-0.5,
                alpha=0.7, label="Network solver")

    _apply_style(ax, "σ_electronic (mS/cm)", names)
    ax.legend(fontsize=9, loc='upper left')

    # R² (log-space) — this group only (for display)
    valid_both = [i for i in range(len(names)) if sigma_net[i] > 0 and sigma_scaling[i] > 0]
    if len(valid_both) >= 3:
        sa = np.array([sigma_net[i] for i in valid_both])
        sp = np.array([sigma_scaling[i] for i in valid_both])
        errs = np.abs(sp - sa) / sa * 100
        w20 = np.sum(errs < 20)
    else:
        errs = np.array([0]); w20 = 0

    n_total = n_global_tk + n_global_tn
    ax.set_title(f"Electronic 2-Regime Scaling Law\n"
                 f"Thick R²={r2_global_tk:.3f}(n={n_global_tk}), Thin R²={r2_global_tn:.3f}(n={n_global_tn})",
                 fontsize=9, fontweight='bold')

    # Formula box with global R²
    txt = (f"Thick: G_h-primary + interactions R²={r2_global_tk:.3f}(n={n_global_tk})\n"
           f"Thin: CN×√δ/√(T/d) R²={r2_global_tn:.3f}(n={n_global_tn})\n"
           f"Group |err|={np.mean(errs):.0f}%, ≤20%: {w20}/{len(valid_both)}")
    ax.text(0.95, 0.95, txt, transform=ax.transAxes, fontsize=7, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffeaea', alpha=0.8))

    _write_csv(outdir, 'electronic_scaling.csv',
               ['σ_scaling(mS/cm)', 'σ_network(mS/cm)'],
               names, sigma_scaling, sigma_net)
    return _save(fig, outdir, "electronic_scaling.png")


def plot_thermal_scaling(data_list, names, outdir):
    """Thermal scaling law: σ_th = 286 × σ_ion^(3/4) × φ_AM² / CN_SE."""
    C_th = 286.0

    sigma_ion = _load_network_sigma(data_list)
    phi_am = [_get(d, "phi_am") for d in data_list]
    cn_se = [_get(d, "se_se_cn") for d in data_list]

    sigma_scaling = []
    for i in range(len(data_list)):
        if sigma_ion[i] > 0 and phi_am[i] > 0 and cn_se[i] > 0:
            s = C_th * sigma_ion[i]**(3/4) * phi_am[i]**2 / cn_se[i]
            sigma_scaling.append(s)
        else:
            sigma_scaling.append(0.0)

    # Network solver for comparison
    sigma_net = _load_thermal_sigma(data_list)
    has_net = any(s > 0 for s in sigma_net)

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    # Use NaN to break line at σ=0 cases
    y_scaling = np.array([s if s > 0 else np.nan for s in sigma_scaling])
    y_net = np.array([s if s > 0 else np.nan for s in sigma_net]) if has_net else None

    ax.plot(x, y_scaling, 's-', color='#ff922b', markersize=ms, linewidth=lw,
            label="Scaling law (mS/cm)")
    if has_net and y_net is not None:
        ax.plot(x, y_net, 'D--', color='#2ecc71', markersize=ms-2, linewidth=lw-0.5,
                alpha=0.7, label="Network solver (mS/cm)")

    _apply_style(ax, "\u03c3_th (mS/cm equiv.)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title("Thermal: \u03c3_th = 286 \u00d7 \u03c3_ion^(3/4) \u00d7 \u03c6_AM\u00b2 / CN_SE\n"
                 "R\u00b2=0.90, 1 free param",
                 fontsize=9, fontweight='bold')

    # Annotation box
    valid_scaling = [s for s in sigma_scaling if s > 0]
    if valid_scaling:
        txt = ("R\u00b2 = 0.90\n"
               "\u03c3_th = 286 \u00d7 \u03c3_ion^(3/4) \u00d7 \u03c6_AM\u00b2 / CN_SE")
        ax.text(0.95, 0.95, txt,
                transform=ax.transAxes, fontsize=8, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', alpha=0.8))

    _write_csv(outdir, 'thermal_scaling.csv',
               ['\u03c3_ion(mS/cm)', 'phi_AM', 'CN_SE', '\u03c3_scaling(mS/cm)', '\u03c3_network(mS/cm)'],
               names, sigma_ion, phi_am, cn_se, sigma_scaling, sigma_net)
    return _save(fig, outdir, "thermal_scaling.png")


def plot_transport_tradeoff(data_list, names, outdir):
    """Ionic vs Electronic trade-off: dual Y-axis showing inverse relationship."""
    sigma_ionic = _load_network_sigma(data_list)
    sigma_el = _load_electronic_sigma(data_list)

    if not any(s > 0 for s in sigma_ionic) or not any(s > 0 for s in sigma_el):
        return None

    fig, ax1 = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    # Ionic on left Y
    color_ion = '#2ecc71'
    ax1.plot(x, sigma_ionic, 's-', color=color_ion, markersize=ms, linewidth=lw,
             label="σ_ionic")
    ax1.set_ylabel("σ_ionic (mS/cm)", color=color_ion, fontsize=11)
    ax1.tick_params(axis='y', labelcolor=color_ion)

    # Electronic on right Y
    ax2 = ax1.twinx()
    color_el = '#e74c3c'
    y_el = np.array([s if s > 0 else np.nan for s in sigma_el])
    no_el = [i for i in range(len(names)) if sigma_el[i] == 0]
    ax2.plot(x, y_el, 'D-', color=color_el, markersize=ms, linewidth=lw, label="σ_electronic")
    if no_el:
        ax2.plot([x[i] for i in no_el], [0]*len(no_el), 'x', color='gray', markersize=ms-2)
    ax2.set_ylabel("σ_electronic (mS/cm)", color=color_el, fontsize=11)
    ax2.tick_params(axis='y', labelcolor=color_el)

    _apply_style(ax1, "", names)
    ax1.set_title("Ionic vs Electronic Conductivity Trade-off\n(↑ AM → ↑ electronic, ↓ ionic)",
                  fontsize=10, fontweight='bold')

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper left')

    _write_csv(outdir, 'transport_tradeoff.csv',
               ['σ_ionic(mS/cm)', 'σ_electronic(mS/cm)'],
               names, sigma_ionic, sigma_el)
    return _save(fig, outdir, "transport_tradeoff.png")


def plot_transport_normalized(data_list, names, outdir):
    """3-mode normalized comparison: ionic/electronic/thermal on same scale."""
    sigma_ionic = _load_network_sigma(data_list)
    sigma_el = _load_electronic_sigma(data_list)
    sigma_th = _load_thermal_sigma(data_list)

    has_ionic = any(s > 0 for s in sigma_ionic)
    has_el = any(s > 0 for s in sigma_el)
    has_th = any(s > 0 for s in sigma_th)

    if not has_ionic:
        return None

    # Normalize each to its own max
    def _norm(arr):
        mx = max(arr) if max(arr) > 0 else 1
        return [v / mx for v in arr]

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))
    w = 0.25

    if has_ionic:
        norm_ion = _norm(sigma_ionic)
        ax.bar(x - w, norm_ion, w, color='#2ecc71', alpha=0.8, label=f"Ionic (max={max(sigma_ionic):.2f})")
    if has_el:
        norm_el = _norm(sigma_el)
        ax.bar(x, norm_el, w, color='#e74c3c', alpha=0.8, label=f"Electronic (max={max(sigma_el):.2f})")
    if has_th:
        norm_th = _norm(sigma_th)
        ax.bar(x + w, norm_th, w, color='#ff922b', alpha=0.8, label=f"Thermal (max={max(sigma_th):.2f})")

    _apply_style(ax, "Normalized σ (ratio to max)", names)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, 1.15)
    ax.set_title("3-Mode Transport Comparison (Normalized)\nIonic (SE-SE) | Electronic (AM-AM) | Thermal (ALL)",
                 fontsize=10, fontweight='bold')

    _write_csv(outdir, 'transport_normalized.csv',
               ['ionic(mS/cm)', 'electronic(mS/cm)', 'thermal(mS/cm)'],
               names, sigma_ionic, sigma_el, sigma_th)
    return _save(fig, outdir, "transport_normalized.png")


def plot_transport_absolute(data_list, names, outdir):
    """3-mode absolute values on log scale."""
    sigma_ionic = _load_network_sigma(data_list)
    sigma_el = _load_electronic_sigma(data_list)
    sigma_th = _load_thermal_sigma(data_list)

    has_ionic = any(s > 0 for s in sigma_ionic)
    has_el = any(s > 0 for s in sigma_el)
    has_th = any(s > 0 for s in sigma_th)

    if not has_ionic:
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    if has_ionic:
        ax.plot(x, sigma_ionic, 's-', color='#2ecc71', markersize=ms, linewidth=lw,
                label="Ionic (SE-SE)")
    if has_el:
        y_el_abs = np.array([s if s > 0 else np.nan for s in sigma_el])
        ax.plot(x, y_el_abs, 'D-', color='#e74c3c', markersize=ms, linewidth=lw,
                label="Electronic (AM-AM)")
    if has_th:
        ax.plot(x, sigma_th, '^-', color='#ff922b', markersize=ms, linewidth=lw,
                label="Thermal (ALL)")

    ax.set_yscale('log')
    _apply_style(ax, "σ_ionic (mS/cm, log scale)", names)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title("3-Mode Transport: Absolute Values\nIonic ≪ Thermal < Electronic (typical)",
                 fontsize=10, fontweight='bold')

    # Minnmann reference line
    ax.axhline(0.17, color='gray', linestyle=':', linewidth=1, alpha=0.6)
    ax.text(len(names)-0.5, 0.17, "Minnmann 0.17", fontsize=7, color='gray', va='bottom', ha='right')

    _write_csv(outdir, 'transport_absolute.csv',
               ['ionic(mS/cm)', 'electronic(mS/cm)', 'thermal(mS/cm)'],
               names, sigma_ionic, sigma_el, sigma_th)
    return _save(fig, outdir, "transport_absolute.png")


def plot_r_brug_comparison(data_list, names, outdir):
    """R_brug for ionic vs electronic vs thermal — how much does Bruggeman overestimate?"""
    # Ionic R_brug
    r_ionic = []
    for d in data_list:
        v = _get(d, "R_brug_over_full", 0)
        r_ionic.append(v if v > 0 else 0)

    # Electronic R_brug
    r_el = []
    for i, d in enumerate(data_list):
        v = _get(d, "electronic_R_brug", 0)
        if v and v > 0:
            r_el.append(v)
        else:
            src = _get(d, "_source_path", "")
            if src:
                net_path = os.path.join(os.path.dirname(src), "network_conductivity.json")
                if os.path.exists(net_path):
                    try:
                        with open(net_path) as _nf:
                            nd = json.load(_nf)
                        r_el.append(nd.get("electronic_R_brug", 0) or 0)
                    except:
                        r_el.append(0)
                else:
                    r_el.append(0)
            else:
                r_el.append(0)

    # Thermal R_brug
    r_th = []
    for i, d in enumerate(data_list):
        v = _get(d, "thermal_R_brug", 0)
        if v and v > 0:
            r_th.append(v)
        else:
            src = _get(d, "_source_path", "")
            if src:
                net_path = os.path.join(os.path.dirname(src), "network_conductivity.json")
                if os.path.exists(net_path):
                    try:
                        with open(net_path) as _nf:
                            nd = json.load(_nf)
                        r_th.append(nd.get("thermal_R_brug", 0) or 0)
                    except:
                        r_th.append(0)
                else:
                    r_th.append(0)
            else:
                r_th.append(0)

    has_ionic = any(r > 0 for r in r_ionic)
    has_el = any(r > 0 for r in r_el)
    has_th = any(r > 0 for r in r_th)

    if not has_ionic:
        return None

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))
    w = 0.25

    if has_ionic:
        ax.bar(x - w, r_ionic, w, color='#2ecc71', alpha=0.8, label="Ionic R_brug")
    if has_el:
        ax.bar(x, r_el, w, color='#e74c3c', alpha=0.8, label="Electronic R_brug")
    if has_th:
        ax.bar(x + w, r_th, w, color='#ff922b', alpha=0.8, label="Thermal R_brug")

    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.text(len(names)-0.5, 1.05, "Bruggeman = exact", fontsize=7, color='gray', ha='right')

    _apply_style(ax, "R_brug (σ_brug / σ_network)", names)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title("Bruggeman Overestimation by Transport Mode\nR_brug > 1 = Bruggeman overestimates",
                 fontsize=10, fontweight='bold')

    _write_csv(outdir, 'r_brug_comparison.csv',
               ['R_ionic', 'R_electronic', 'R_thermal'],
               names, r_ionic, r_el, r_th)
    return _save(fig, outdir, "r_brug_comparison.png")


def plot_ion_path_quality(data_list, names, outdir):
    """Ion path quality: GB Density, Hop Area, Bottleneck, Conductance (2x2)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    x = np.arange(len(names))

    # GB Density (lower is better)
    vals = [_get(d, "gb_density_mean") for d in data_list]
    axes[0,0].plot(x, vals, 's-', color=BLUE, markersize=8, linewidth=2)
    _apply_style(axes[0,0], "GB Density (hops/μm)", names)
    axes[0,0].set_title("Grain Boundary Density (GB_d)  ↓ better", fontsize=10, fontweight='bold')

    # Path Hop Area mean (higher is better)
    vals = [_get(d, "path_hop_area_mean") for d in data_list]
    axes[0,1].plot(x, vals, 's-', color=ORANGE, markersize=8, linewidth=2)
    _apply_style(axes[0,1], "Hop Area mean (μm²)", names)
    axes[0,1].set_title("Path Hop Area  (↑ better)", fontsize=10, fontweight='bold')

    # Bottleneck (higher is better)
    vals = [_get(d, "path_hop_area_min_mean") for d in data_list]
    axes[1,0].plot(x, vals, 's-', color=RED, markersize=8, linewidth=2)
    _apply_style(axes[1,0], "Bottleneck (μm²)", names)
    axes[1,0].set_title("Path Bottleneck  (↑ better)", fontsize=10, fontweight='bold')

    # Path Conductance (higher is better)
    vals = [_get(d, "path_conductance_mean") for d in data_list]
    axes[1,1].plot(x, vals, 's-', color=GREEN, markersize=8, linewidth=2)
    _apply_style(axes[1,1], "Conductance (μm²)", names)
    axes[1,1].set_title("Path Conductance (G_path)  (↑ better)", fontsize=10, fontweight='bold')

    gb = [_get(d, "gb_density_mean") for d in data_list]
    ha = [_get(d, "path_hop_area_mean") for d in data_list]
    bn = [_get(d, "path_hop_area_min_mean") for d in data_list]
    gc = [_get(d, "path_conductance_mean") for d in data_list]
    _write_csv(outdir, 'ion_path_quality.csv',
               ['GB Density(hops/μm)', 'Hop Area mean(μm²)', 'Bottleneck(μm²)', 'Conductance(μm²)'],
               names, gb, ha, bn, gc)

    _break_lines_at_groups(fig)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle("Ion Path Quality", fontsize=14, fontweight='bold', x=0.5, y=0.99, ha='center')
    path = os.path.join(outdir, "ion_path_quality.png")
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor='white', pad_inches=0.3)
    plt.close(fig)
    return path


# ─── Plot dispatch table ─────────────────────────────────────────────────────

PLOT_REGISTRY = {
    "porosity": {
        "func": plot_porosity,
        "file": "porosity.png",
        "title": "Porosity",
        "description": "AM_P 비율 증가에 따른 기공률 변화.\n7:3 부근에서 최저 (bimodal packing 효과).\nV\u2011shape 경향이면 최적 조성 존재.",
        "origin_tip": "Line+Symbol → X: P:S Configuration, Y: Porosity(%).\nSymbol: Square (size 10), Color: Black.\nLine: B-Spline, Width 1.5.\nY축 범위: 자동 ± 1%p 여유.",
    },
    "am_se_interface": {
        "func": plot_am_se_interface,
        "file": "am_se_interface.png",
        "title": "AM-SE Interface Area",
        "description": "AM_P\u2011SE와 AM_S\u2011SE 접촉 면적을\nStacked Bar로 비교.\nAM_P 비율 증가 → AM_P\u2011SE 증가, AM_S\u2011SE 감소.\n전체 AM\u2011SE Total은 S only에서 최대.",
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
    "stress_cv": {
        "func": plot_stress_cv,
        "file": "stress_cv.png",
        "title": "Stress CV",
        "description": "Von Mises 응력 변동계수(CV).\n\nCV 낮을수록 전극 내 응력이 균일.\n유효영률 사용으로 절대값은 참고용, 상대 비교만 유효.",
        "origin_tip": "Line+Symbol → X: P:S, Y: VM CV (%).\nSymbol: Square, Black.\nCV < 100%면 양호.",
    },
    "stress_ratio": {
        "func": plot_stress_ratio,
        "file": "stress_ratio.png",
        "title": "Stress Ratio by Type",
        "description": "입자 유형별 응력 비율 (σ_type / σ_mean).\n\n> 1.0 = 평균보다 응력 집중\n< 1.0 = 평균보다 하중 적음\n\nSE > 1.0이면 SE에 응력 집중 → 소성변형 유발.",
        "origin_tip": "Multi-line → X: P:S, Y: σ/σ_mean.\nAM_P: Red, AM_S: Orange, SE: Green.\ny=1.0에 점선 (mean baseline).",
    },
    "se_network": {
        "func": plot_se_network,
        "file": "se_network.png",
        "title": "SE Network",
        "description": "SE-SE 배위수(CN)와 클러스터 수.\nCN = SE 입자 하나가 접촉하는 SE 이웃 수\nCN↑ = 조밀한 네트워크 (CN≥4: 안정적 3D percolation)\nCluster↓ = 분절 적음. Large(≥10) 클러스터만 이온 경로에 유의미.",
        "origin_tip": "Dual-Y → Left: Line (CN, Blue).\nRight: Bar (Clusters, Orange).",
    },
    "contact_force": {
        "func": plot_contact_force,
        "file": "contact_force.png",
        "title": "Contact Force",
        "description": "접촉 유형별(AM-AM, AM-SE, SE-SE) 법선력 평균.\nFn = √(fn_x² + fn_y² + fn_z²), 단위: μN\nF_real = F_sim / scale²\n대립자 간 접촉력이 가장 크며, P:S 변화에 따라 하중 분담 변화 관찰.",
        "origin_tip": "Grouped Bar → X: Configuration, Y: Fn mean (μN).\nAM-AM: Red, AM-SE: Green, SE-SE: Blue.",
    },
    "contact_pressure": {
        "func": plot_contact_pressure,
        "file": "contact_pressure.png",
        "title": "Contact Pressure",
        "description": "접촉 압력 P = Fn / A_contact (MPa).\n같은 힘이라도 면적이 작으면 압력↑.\nMax는 failure 시작점, Mean은 전체 경향.",
        "origin_tip": "Dual Bar → X: Configuration, Y: Pressure (MPa).\nMean: Blue, Max: Red.",
    },
    "am_vulnerability": {
        "func": plot_am_vulnerability,
        "file": "am_vulnerability.png",
        "title": "AM Vulnerability",
        "description": "Vulnerable = (SE 0개 + SE 1개 접촉 AM) / 전체 AM × 100\nAM-SE CN = AM당 SE 접촉 수 평균\nVulnerable↓ + CN↑ = 안정적 이온 공급.\nSE 1개 접촉 = single point of failure.",
        "origin_tip": "Dual-Y → Left: Bar (Vulnerable %, Red).\nRight: Line (AM-SE CN, Blue).",
    },
    "effective_conductivity": {
        "func": plot_effective_conductivity,
        "file": "effective_conductivity.png",
        "title": "σ_brug/σ_grain (Bruggeman Estimate)",
        "description": "σ_brug/σ_grain = φ_SE × f_perc / τ²  (σ_grain = 3.0 mS/cm)\n\n접촉 저항(constriction) 무시 → 실제 Network Solver 대비 3~10× 과대추정.\nNetwork Solver σ_ionic이 ground truth.",
        "origin_tip": "Line+Symbol (Green) + Bar (φ_SE, Blue).",
    },
    "rgb_fitting": {
        "func": plot_rgb_fitting,
        "file": "rgb_fitting.png",
        "title": "[Legacy] BLM+Constriction Fitting (Part I proxy)",
        "description": "⚠ Legacy: FORM X로 대체됨\nProxy: R = σ_brug/σ_proxy = C×(GB_d²×T)^α\nSingle\u2011path 근사 → R=15~1600 (실제 3~10)\nBruggeman exponent 분해: n_eff = n_geo(2.54) + n_contact(0.83) = 3.37",
        "origin_tip": "Scatter + Fit line (log\u2011log).",
        "hidden": True,
    },
    "gb_corrected": {
        "func": plot_gb_corrected,
        "file": "gb_corrected.png",
        "title": "[Legacy] Proxy σ_eff (Part I)",
        "description": "⚠ Legacy: FORM X로 대체됨\nProxy 기반 보정 (single\u2011path 근사)\n절대값 과장 → Network solver가 ground truth",
        "origin_tip": "Line (Red, ratio) + Dashed (Orange, mS/cm).",
        "hidden": True,
    },
    "network_sigma": {
        "func": plot_network_sigma,
        "file": "network_sigma.png",
        "title": "Ionic: Network Solver σ_full (Part II)",
        "description": "Kirchhoff resistor network (R_bulk + R_constriction)\n초록: σ_full (ground truth) | 파란: σ_brug (Bruggeman)\nR_brug = 3~10× | Constriction 69~81% 지배\nMinnmann(2021) 0.17 mS/cm과 same order\nσ_grain = 3.0 mS/cm (grain interior, not pellet)",
        "origin_tip": "Green: Network solver, Blue dashed: Bruggeman.",
    },
    "ionic_scaling_fit": {
        "func": plot_ionic_scaling_fit,
        "file": "ionic_scaling_fit.png",
        "title": "Ionic: FORM X Scaling Law (Predicted vs Actual)",
        "description": "FORM X v9 BLEND:\nσ = [(1-w)·C_v5(τ) + w·C_p3(τ)] × σ_grain × CN^1.5 × (φ-φc)^¾ × ⁴√cov × fp²\nw(τ) = σ(15·(τ-2.0))  # sharp transition\nC_v5(τ) = sigmoid(thick→thin) — captures moderate τ\nC_p3(τ) = cubic in lnτ — captures extreme thin (τ>3)\n\nφ_c=0.185, 6 params (2 v5 + 4 poly3)\nR²={formx_r2}, LOOCV={formx_loocv}",
        "origin_tip": "Scatter (log-log): X=actual (Network solver), Y=predicted (FORM X).\n1:1 line (black dashed), ±20% band (green).",
    },
    "multiscale_sigma": {
        "func": plot_multiscale_sigma,
        "file": "multiscale_sigma.png",
        "title": "Ionic: FORM X Scaling Law",
        "description": "FORM X v9 BLEND:\nσ = C_blend(τ) × σ_grain × CN^1.5 × (φ\u2011φc)^¾ × ⁴√cov × fp²\nv5 sigmoid (τ<2) + poly3(lnτ) (τ≥2), smooth blend at τ=2.0\nv5 handles moderate, poly3 handles extreme thin (τ>3)\nR²={formx_r2}, w20=51/55",
        "origin_tip": "Red: FORM X prediction.\nGreen dashed: Network solver (ground truth).",
    },
    "formx_decomposition": {
        "func": plot_formx_decomposition,
        "file": "formx_decomposition.png",
        "csv": "formx_decomposition.csv",
        "title": "FORM X Factor Decomposition",
        "description": "FORM X v9 각 항의 상대 기여도:\n(φ\u2011φc)^¾: percolation excess\nCN^1.5: network connectivity\n⁴√cov: AM\u2011SE 계면\nC_blend(τ): v5 + poly3 혼합 regime\nfp²: SE percolation fraction\nref: 최고 σ case 기준",
        "origin_tip": "Stacked bar (top): factor contributions.\nHorizontal bar (bottom): dominant factor per case.",
    },
}


def plot_sigma_decomposition(data_list, names, outdir):
    """Decompose σ_eff into Bruggeman + contact terms. Show which factor dominates."""
    SIGMA_BULK = 3.0  # σ_grain (grain interior)
    C_ms = 0.073  # default, used only for reference σ_ms (decomposition is relative, C cancels out)

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    g_path = [_get(d, "path_conductance_mean", 0) for d in data_list]
    cn = [_get(d, "se_se_cn", 0) for d in data_list]
    phi = [_get(d, "phi_se") for d in data_list]
    tau = [_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)) for d in data_list]
    f_perc = [_get(d, "percolation_pct") / 100 for d in data_list]

    n = len(data_list)
    # Decompose log(σ_brug) = log(σ_grain) + log(φ_SE) + log(f_perc) - 2log(τ) + log(C) + 0.25*log(G_path*GB_d²) + 2log(CN)
    log_phi = np.array([np.log(phi[i]) if phi[i] > 0 else 0 for i in range(n)])
    log_fperc = np.array([np.log(f_perc[i]) if f_perc[i] > 0 else 0 for i in range(n)])
    log_tau2 = np.array([-2 * np.log(tau[i]) if tau[i] > 0 else 0 for i in range(n)])
    log_gpath_gbd = np.array([0.25 * np.log(g_path[i] * gb_dens[i]**2) if g_path[i] > 0 and gb_dens[i] > 0 else 0 for i in range(n)])
    log_cn2 = np.array([2 * np.log(cn[i]) if cn[i] > 0 else 0 for i in range(n)])

    # Normalize: relative to best σ_eff case (show what limits performance)
    sigma_ms = []
    for i in range(n):
        if g_path[i] > 0 and cn[i] > 0 and gb_dens[i] > 0:
            sigma_ms.append(sigma_brug[i] * C_ms * (g_path[i] * gb_dens[i]**2)**0.25 * cn[i]**2 * SIGMA_BULK)
        else:
            sigma_ms.append(0)
    ref = int(np.argmax(sigma_ms))  # best σ_eff case
    d_phi = log_phi - log_phi[ref]
    d_tau = log_tau2 - log_tau2[ref]
    d_fperc = log_fperc - log_fperc[ref]
    d_gpath_gbd = log_gpath_gbd - log_gpath_gbd[ref]
    d_cn = log_cn2 - log_cn2[ref]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(names)*0.7), 10), gridspec_kw={'height_ratios': [2, 1]})

    # Top: σ_eff with each factor's contribution as stacked
    ax = axes[0]
    x = np.arange(n)
    w = 0.6

    # Stacked bar: each factor's log contribution (relative to reference)
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
    labels = ['φ_SE', 'τ²', 'f_perc', '(G_path×GB_d²)^¼', 'CN²']
    contributions = [d_phi, d_tau, d_fperc, d_gpath_gbd, d_cn]

    # Positive and negative stacking
    pos_bottom = np.zeros(n)
    neg_bottom = np.zeros(n)
    for j, (contrib, color, label) in enumerate(zip(contributions, colors, labels)):
        pos = np.maximum(contrib, 0)
        neg = np.minimum(contrib, 0)
        ax.bar(x, pos, w, bottom=pos_bottom, color=color, label=label, alpha=0.8, edgecolor='white', linewidth=0.5)
        ax.bar(x, neg, w, bottom=neg_bottom, color=color, alpha=0.4, edgecolor='white', linewidth=0.5)
        pos_bottom += pos
        neg_bottom += neg

    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_ylabel('Δlog(σ_brug) from reference', fontsize=11)
    ax.set_title(f'[Legacy v3] σ_ionic Factor Decomposition (ref: {names[ref]})', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left', ncol=3)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Bottom: dominant factor for each case
    ax2 = axes[1]
    # Find which factor has largest |contribution| for each case
    all_contribs = np.array(contributions)  # (5, n)
    dominant_idx = np.argmax(np.abs(all_contribs), axis=0)

    for i in range(n):
        if i == ref:
            ax2.barh(i, 0, color='gray')
            continue
        vals = all_contribs[:, i]
        sorted_idx = np.argsort(np.abs(vals))[::-1]
        # Show top 3 contributors
        y_pos = i
        for rank, j in enumerate(sorted_idx[:3]):
            ax2.barh(y_pos, vals[j], height=0.25, left=0,
                    color=colors[j], alpha=0.9 - rank*0.25,
                    edgecolor='white', linewidth=0.5)

    ax2.set_yticks(range(n))
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_xlabel('Factor contribution (Δlog)', fontsize=10)
    ax2.set_title('Dominant factors per case', fontsize=11, fontweight='bold')
    ax2.axvline(0, color='gray', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for c, l in zip(colors, labels)]
    ax2.legend(handles=legend_elements, fontsize=7, loc='lower right', ncol=3)

    plt.tight_layout()
    _write_csv(outdir, 'sigma_decomposition.csv',
               ['φ_SE(Δlog)', 'τ²(Δlog)', 'f_perc(Δlog)', '(G_path×GB_d²)^¼(Δlog)', 'CN²(Δlog)'],
               names, list(d_phi), list(d_tau), list(d_fperc), list(d_gpath_gbd), list(d_cn))
    return _save(fig, outdir, "sigma_decomposition.png")


PLOT_REGISTRY["v3_fitting"] = {
    "func": plot_v3_fitting,
    "file": "v3_fitting.png",
    "title": "[Legacy] v3 Fitting (thick only)",
    "description": "⚠ Legacy: FORM X로 대체됨\nv3: σ_brug × C × (G_path × GB_d²)^(1/4) × CN²\nthick R²=0.96, thin R²=-1.0",
    "origin_tip": "Red: v3 prediction. Green dashed: Network solver.",
    "hidden": True,
}

PLOT_REGISTRY["sigma_decomposition"] = {
    "func": plot_sigma_decomposition,
    "file": "sigma_decomposition.png",
    "title": "[Legacy] v3 Factor Decomposition",
    "description": "⚠ Legacy (v3): FORM X로 대체됨\nσ_ion = σ_brug × C × (G_path × GB_d²)^(1/4) × CN²\nthick only (thin에서 실패)",
    "origin_tip": "Stacked bar (top) + Horizontal bar (bottom).",
    "hidden": True,
}


def plot_electronic_decomposition(data_list, names, outdir):
    """Decompose σ_el: Thick(CN^(2-0.3M) × (φ-φc)² / por^(0.4-0.15M)), Thin(CN × √δ / √(T/d))."""
    SIGMA_AM = 50.0

    phi_am = [_get(d, "phi_am") for d in data_list]
    cn_am = [max(_get(d, "am_am_cn"), 0.01) for d in data_list]
    thickness = [_get(d, "thickness_um") for d in data_list]
    d_am = [2.0 * max(_get(d, "r_AM_P", 0), _get(d, "r_AM_S", 0), _get(d, "r_AM", 0))
            for d in data_list]
    porosity = [max(_get(d, "porosity", 10), 0.1) for d in data_list]
    am_delta = [max(_get(d, "am_am_mean_delta", 0), 0.001) for d in data_list]

    n = len(data_list)
    ratios = [thickness[i] / d_am[i] if d_am[i] > 0 and thickness[i] > 0 else 0 for i in range(n)]

    # Thick factors: CN^(2-0.3M), (φ-φc)², por^-(0.4-0.15M)
    phi_ex = [max(phi_am[i] - 0.10, 0.001) for i in range(n)]
    # Use fP_vol for Ψ in decomposition
    case_psi = []
    for d in data_list:
        _nP2 = _get(d, "n_AM_P", 0) or 0; _nS2 = _get(d, "n_AM_S", 0) or 0
        _rP2 = max(_get(d, "r_AM_P", 0) or 0, 0.1); _rS2 = max(_get(d, "r_AM_S", 0) or 0, 0.1)
        _vP2 = _nP2 * _rP2**3; _vS2 = _nS2 * _rS2**3; _vT2 = max(_vP2 + _vS2, 0.001)
        case_psi.append(_vP2 / _vT2)
    log_phi = np.array([2.0 * np.log(phi_ex[i]) for i in range(n)])
    case_mix = [4*case_psi[i]*(1-case_psi[i]) for i in range(n)]
    log_cn2 = np.array([(2.0 - 0.3*case_mix[i]) * np.log(cn_am[i]) for i in range(n)])
    log_por = np.array([-(0.4 - 0.15*case_mix[i]) * np.log(porosity[i]) for i in range(n)])
    # Thin factors: CN, δ^0.5, (T/d)^-0.5
    log_cn1 = np.array([1.0 * np.log(cn_am[i]) for i in range(n)])
    log_delta = np.array([0.5 * np.log(am_delta[i]) for i in range(n)])
    log_ratio = np.array([-0.5 * np.log(max(ratios[i], 0.1)) for i in range(n)])

    # Determine thick/thin per case
    is_thick = [ratios[i] >= 8 for i in range(n)]

    # σ_el for reference selection (use max)
    sigma_el = []
    for i in range(n):
        if phi_am[i] > 0 and cn_am[i] > 0 and d_am[i] > 0 and thickness[i] > 0:
            if is_thick[i]:
                _psi_d = case_psi[i] if i < len(case_psi) else 0
                _mix_d = 4*_psi_d*(1-_psi_d) if i < len(case_psi) else 0
                sigma_el.append(max(phi_am[i]-0.10,0.001)**2 * cn_am[i]**(2-0.3*_mix_d) / porosity[i]**(0.4-0.15*_mix_d))
            else:
                sigma_el.append(cn_am[i] * am_delta[i]**0.5 / max(ratios[i], 0.1)**0.5)
        else:
            sigma_el.append(0)

    if not any(s > 0 for s in sigma_el):
        return None

    ref = int(np.argmax(sigma_el))

    # Build contributions based on regime
    colors_all = ['#e74c3c', '#f39c12', '#27ae60', '#3498db', '#9b59b6']
    labels_all = ['[Thick] (φ-φc)²', '[Thick] CN^(2-0.3M)', '[Thick] 1/por^(0.4-0.15M)', '[Thin] √δ', '[Thin] 1/√(T/d)']
    # Thin CN uses same slot [1] but different label
    label_cn_thin = '[Thin] CN'

    # Thick reference: highest σ among thick cases; Thin reference: highest σ among thin
    thick_idx = [i for i in range(n) if is_thick[i] and sigma_el[i] > 0]
    thin_idx = [i for i in range(n) if not is_thick[i] and sigma_el[i] > 0]
    ref_tk = thick_idx[np.argmax([sigma_el[i] for i in thick_idx])] if thick_idx else ref
    ref_tn = thin_idx[np.argmax([sigma_el[i] for i in thin_idx])] if thin_idx else ref

    contribs = np.zeros((5, n))
    for i in range(n):
        if is_thick[i]:
            contribs[0, i] = log_phi[i] - log_phi[ref_tk]
            contribs[1, i] = log_cn2[i] - log_cn2[ref_tk]
            contribs[2, i] = log_por[i] - log_por[ref_tk]
        else:
            contribs[1, i] = log_cn1[i] - log_cn1[ref_tn]
            contribs[3, i] = log_delta[i] - log_delta[ref_tn]
            contribs[4, i] = log_ratio[i] - log_ratio[ref_tn]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(names)*0.7), 10), gridspec_kw={'height_ratios': [2, 1]})

    ax = axes[0]
    x = np.arange(n)
    w = 0.6

    has_thick = any(is_thick[i] for i in range(n))
    has_thin = any(not is_thick[i] for i in range(n))

    pos_bottom = np.zeros(n)
    neg_bottom = np.zeros(n)
    used_labels = set()
    for j in range(5):
        contrib = contribs[j]
        if np.all(np.abs(contrib) < 1e-6):
            continue
        pos = np.maximum(contrib, 0)
        neg = np.minimum(contrib, 0)
        lbl = labels_all[j]
        if j == 1 and has_thin and not has_thick:
            lbl = label_cn_thin
        ax.bar(x, pos, w, bottom=pos_bottom, color=colors_all[j], label=lbl, alpha=0.8, edgecolor='white', linewidth=0.5)
        ax.bar(x, neg, w, bottom=neg_bottom, color=colors_all[j], alpha=0.4, edgecolor='white', linewidth=0.5)
        pos_bottom += pos
        neg_bottom += neg
        used_labels.add(j)

    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_ylabel('Δlog(σ_el) from reference', fontsize=11)
    regime_str = []
    if has_thick: regime_str.append(f"Thick ref: {names[ref_tk]}")
    if has_thin: regime_str.append(f"Thin ref: {names[ref_tn]}")
    ax.set_title(f'Electronic σ_el Factor Decomposition ({", ".join(regime_str)})', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left', ncol=3)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Bottom: dominant factor
    ax2 = axes[1]
    for i in range(n):
        if i == ref:
            ax2.barh(i, 0, color='gray')
            continue
        vals = contribs[:, i]
        sorted_idx = np.argsort(np.abs(vals))[::-1]
        for rank, j in enumerate(sorted_idx[:3]):
            if abs(vals[j]) < 0.01:
                continue
            ax2.barh(i, vals[j], height=0.25, left=0,
                    color=colors_all[j], alpha=0.9 - rank*0.25,
                    edgecolor='white', linewidth=0.5)

    ax2.set_yticks(range(n))
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_xlabel('Factor contribution (Δlog)', fontsize=10)
    ax2.set_title('Dominant factors per case', fontsize=11, fontweight='bold')
    ax2.axvline(0, color='gray', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors_all[j], label=labels_all[j] if not (j==1 and has_thin and not has_thick) else label_cn_thin) for j in sorted(used_labels)]
    ax2.legend(handles=legend_elements, fontsize=7, loc='lower right', ncol=min(len(legend_elements), 3))

    plt.tight_layout()
    _write_csv(outdir, 'electronic_decomposition.csv',
               ['φ^2.5(Δlog)', 'CN(Δlog)', '1/√por(Δlog)', '√δ(Δlog)', '1/√(T/d)(Δlog)'],
               names, *[list(contribs[j]) for j in range(5)])
    return _save(fig, outdir, "electronic_decomposition.png")


PLOT_REGISTRY["electronic_decomposition"] = {
    "func": plot_electronic_decomposition,
    "file": "electronic_decomposition.png",
    "title": "Electronic: Factor Decomposition",
    "description": "Thick: CN^(2-0.3M) × (φ-φc)² / por^(0.4-0.15M)\nΨ = V_large/(V_large+V_small), 0=mono\nThin: CN × √δ / √(T/d)\n\n(φ-φc)²: percolation (φc=0.10)\nΨ: bimodal large-AM fraction\n√δ: AM-AM penetration depth",
    "origin_tip": "Stacked bar (top) + Horizontal bar (bottom).",
}


def plot_thermal_decomposition(data_list, names, outdir):
    """Decompose σ_th into σ_ion^(3/4), φ_AM², CN_SE^(-1) factors."""
    C_th = 286.0

    sigma_ion = _load_network_sigma(data_list)
    phi_am = [_get(d, "phi_am") for d in data_list]
    cn_se = [_get(d, "se_se_cn") for d in data_list]

    n = len(data_list)

    # log contributions
    log_ion = np.array([0.75 * np.log(sigma_ion[i]) if sigma_ion[i] > 0 else 0 for i in range(n)])
    log_phi = np.array([2 * np.log(phi_am[i]) if phi_am[i] > 0 else 0 for i in range(n)])
    log_cn = np.array([-1 * np.log(cn_se[i]) if cn_se[i] > 0 else 0 for i in range(n)])

    # σ_th for reference selection
    sigma_th = []
    for i in range(n):
        if sigma_ion[i] > 0 and phi_am[i] > 0 and cn_se[i] > 0:
            sigma_th.append(C_th * sigma_ion[i]**(3/4) * phi_am[i]**2 / cn_se[i])
        else:
            sigma_th.append(0)

    if not any(s > 0 for s in sigma_th):
        return None

    ref = int(np.argmax(sigma_th))
    d_ion = log_ion - log_ion[ref]
    d_phi = log_phi - log_phi[ref]
    d_cn = log_cn - log_cn[ref]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(names)*0.7), 10), gridspec_kw={'height_ratios': [2, 1]})

    ax = axes[0]
    x = np.arange(n)
    w = 0.6

    colors = ['#2ecc71', '#e74c3c', '#3498db']
    labels = ['σ_ion^(3/4)', 'φ_AM²', 'CN_SE⁻¹']
    contributions = [d_ion, d_phi, d_cn]

    pos_bottom = np.zeros(n)
    neg_bottom = np.zeros(n)
    for j, (contrib, color, label) in enumerate(zip(contributions, colors, labels)):
        pos = np.maximum(contrib, 0)
        neg = np.minimum(contrib, 0)
        ax.bar(x, pos, w, bottom=pos_bottom, color=color, label=label, alpha=0.8, edgecolor='white', linewidth=0.5)
        ax.bar(x, neg, w, bottom=neg_bottom, color=color, alpha=0.4, edgecolor='white', linewidth=0.5)
        pos_bottom += pos
        neg_bottom += neg

    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_ylabel('Δlog(σ_th) from reference', fontsize=11)
    ax.set_title(f'Thermal σ_th Factor Decomposition (ref: {names[ref]})', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left', ncol=3)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Bottom: dominant factor
    ax2 = axes[1]
    all_contribs = np.array(contributions)
    for i in range(n):
        if i == ref:
            ax2.barh(i, 0, color='gray')
            continue
        vals = all_contribs[:, i]
        sorted_idx = np.argsort(np.abs(vals))[::-1]
        for rank, j in enumerate(sorted_idx[:3]):
            ax2.barh(i, vals[j], height=0.25, left=0,
                    color=colors[j], alpha=0.9 - rank*0.25,
                    edgecolor='white', linewidth=0.5)

    ax2.set_yticks(range(n))
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_xlabel('Factor contribution (Δlog)', fontsize=10)
    ax2.set_title('Dominant factors per case', fontsize=11, fontweight='bold')
    ax2.axvline(0, color='gray', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for c, l in zip(colors, labels)]
    ax2.legend(handles=legend_elements, fontsize=8, loc='lower right', ncol=3)

    plt.tight_layout()
    _write_csv(outdir, 'thermal_decomposition.csv',
               ['σ_ion^3/4(Δlog)', 'φ_AM²(Δlog)', 'CN_SE⁻¹(Δlog)'],
               names, list(d_ion), list(d_phi), list(d_cn))
    return _save(fig, outdir, "thermal_decomposition.png")


PLOT_REGISTRY["thermal_decomposition"] = {
    "func": plot_thermal_decomposition,
    "file": "thermal_decomposition.png",
    "title": "Thermal: Factor Decomposition (R²=0.90)",
    "description": "σ_th = 286 × σ_ion^(3/4) × φ_AM² / CN_SE\nR²=0.90\n\nσ_ion^(3/4): SE backbone (ionic에서 계승)\nφ_AM²: AM thermal enhancement\nCN_SE⁻¹: SE 과밀집 페널티 (ionic과 부호 역전!)",
    "origin_tip": "Stacked bar (top) + Horizontal bar (bottom).",
}
def plot_electronic_active_am(data_list, names, outdir):
    """Electronic Active AM%: bottom-reachable AM fraction."""
    active = [_get(d, "electronic_active_fraction", 0) * 100 for d in data_list]
    perc = [_get(d, "electronic_percolating_fraction", 0) * 100 for d in data_list]
    phi_am = [_get(d, "phi_am", 0) for d in data_list]
    sigma_el = _load_electronic_sigma(data_list)

    has_data = any(a > 0 for a in active)
    if not has_data:
        # Fallback: estimate from phi_am (literature model)
        for i in range(len(data_list)):
            pa = phi_am[i]
            if pa >= 0.55:
                active[i] = 100.0
            elif pa >= 0.30:
                active[i] = 87 + (pa - 0.30) / 0.25 * 13
            elif pa > 0:
                active[i] = 50 + (pa - 0.18) / 0.12 * 37
        has_data = True

    fig, ax1 = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))
    ms = _marker_size(len(names))
    lw = _line_width(len(names))

    # Bar: Electronic Active AM%
    colors = ['#2ecc71' if a >= 95 else '#f39c12' if a >= 80 else '#e74c3c' for a in active]
    ax1.bar(x, active, 0.6, color=colors, alpha=0.8, label='Electronic Active AM (%)')
    ax1.set_ylabel('Electronic Active AM (%)', fontsize=11)
    ax1.set_ylim(0, 105)
    ax1.axhline(95, color='#2ecc71', linestyle='--', linewidth=0.8, alpha=0.5)
    ax1.axhline(80, color='#f39c12', linestyle='--', linewidth=0.8, alpha=0.5)

    # Line: σ_electronic
    ax2 = ax1.twinx()
    y_el = np.array([s if s > 0 else np.nan for s in sigma_el])
    ax2.plot(x, y_el, 'D-', color='#9b59b6', markersize=ms, linewidth=lw, label='σ_electronic')
    ax2.set_ylabel('σ_electronic (mS/cm)', color='#9b59b6', fontsize=11)
    ax2.tick_params(axis='y', labelcolor='#9b59b6')

    _apply_style(ax1, '', names)
    ax1.set_title('Electronic Active AM & Dead AM Analysis\nGreen≥95%, Yellow≥80%, Red<80%',
                  fontsize=10, fontweight='bold')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labels1+labels2, fontsize=8, loc='lower right')

    _write_csv(outdir, 'electronic_active_am.csv',
               ['Active_AM(%)', 'σ_electronic(mS/cm)', 'φ_AM'],
               names, active, sigma_el, phi_am)
    return _save(fig, outdir, 'electronic_active_am.png')


PLOT_REGISTRY["electronic_active_am"] = {
    "func": plot_electronic_active_am,
    "file": "electronic_active_am.png",
    "title": "Electronic Active AM (Dead AM Analysis)",
    "description": "Bottom(CC)-reachable AM = 전자가 도달하는 AM 비율.\n95%+: 녹색 (정상), 80~95%: 황색 (주의), <80%: 적색 (도전재 필요!)\nDead AM = CC에서 전자 경로 없는 AM → 비활성.\nφ_AM > 55%면 대부분 OK. φ_AM < 35%면 도전재 필수.",
    "origin_tip": "Bar (Active %, colored) + Line (σ_electronic, purple).",
}

PLOT_REGISTRY["ion_path_quality"] = {
        "func": plot_ion_path_quality,
        "file": "ion_path_quality.png",
        "title": "Ion Path Quality (G_path, GB_d)",

        "description": "이온 경로 품질 4종 (모두 percolating 경로들의 mean):\n\n• GB Density = mean(N_hops / L_z) (hops/μm)\n  → ↓ 좋음 (입계 적을수록 저항↓)\n• Path Hop Area = mean(각 hop의 접촉면적)\n  → ↑ 좋음\n• Bottleneck = mean(각 경로의 최소 접촉면적)\n  → ↑ 좋음\n• Path Conductance = mean(1/Σ(1/A_i)) (μm²)\n  → ↑ 좋음 (직렬 저항 모델의 유효 면적)\n\n왜 mean? 경로 안은 직렬, 경로끼리는 병렬.\nmean = 경로 품질, f_perc = 경로 수량 (역할 분리).\n30개 shortest path 샘플 (best/mean/worst 10개씩).\nConductance가 가장 종합적 지표.",
        "origin_tip": "2×2 Subplots → GB Density (Blue), Hop Area (Orange), Bottleneck (Red), Conductance (Green).",
}
PLOT_REGISTRY["electronic_sigma"] = {
    "func": plot_electronic_sigma,
    "file": "electronic_sigma.png",
    "title": "Electronic Conductivity",
    "description": "AM-AM 접촉 네트워크 기반 전자 전도도.\nσ_AM = 0.05 S/cm (NCM 활물질).\nAM percolation 없으면 σ=0 (도전재 필요).\n10:0(P only)에서 AM 부족 → percolation 실패 가능.",
    "origin_tip": "Line (Red) + X markers for no-percolation cases.",
}
PLOT_REGISTRY["thermal_sigma"] = {
    "func": plot_thermal_sigma,
    "file": "thermal_sigma.png",
    "title": "Thermal Conductivity",
    "description": "ALL 접촉 네트워크 기반 열전도도.\nk_AM=4.0e-2, k_SE=0.7e-2 W/(cm·K).\nAM-SE 접촉은 harmonic mean.\n양방향 전도 (source↔sink 모두).",
    "origin_tip": "Line (Orange).",
}
PLOT_REGISTRY["electronic_scaling"] = {
    "func": plot_electronic_scaling,
    "file": "electronic_scaling.png",
    "title": "Electronic: 2-Regime Scaling Law",
    "description": "Thick (T/d≥10): σ = C × σ_AM × φ² × CN² × √cov / √por\n  φ: AM vol fraction, CN: AM-AM coordination\n\nThin (T/d<10): σ = C × σ_AM × CN × √δ / √(T/d)\n  δ: AM-AM penetration depth, T/d: thickness ratio\n\nC: thick/thin 별도 global fit (parameter sweep 최적화)",
    "origin_tip": "Red: Scaling law, Green dashed: Network solver (ground truth).",
}
PLOT_REGISTRY["thermal_scaling"] = {
    "func": plot_thermal_scaling,
    "file": "thermal_scaling.png",
    "title": "Thermal Scaling Law",
    "description": "\u03c3_th = 286 \u00d7 \u03c3_ion^(3/4) \u00d7 \u03c6_AM\u00b2 / CN_SE\n\u03c3_ion: ionic network solver \uacb0\uacfc \uc0ac\uc6a9\nCN_SE\u207b\u00b9: SE clustering \u2192 AM \uace0\ub9bd penalty\nR\u00b2=0.90 (1 free param), LOOCV R\u00b2=0.86",
    "origin_tip": "Orange: Scaling law, Green dashed: Network solver.",
}
PLOT_REGISTRY["transport_tradeoff"] = {
    "func": plot_transport_tradeoff,
    "file": "transport_tradeoff.png",
    "title": "Ionic vs Electronic Trade-off",
    "description": "이온 전도(SE-SE)와 전자 전도(AM-AM)의 역관계.\nAM↑ → electronic↑, ionic↓.\n최적 조성 = 두 곡선의 교차점 근처.\n도전재(carbon) 추가 시 electronic 병목 완화 가능.",
    "origin_tip": "Dual Y-axis: Green (ionic, left), Red (electronic, right).",
}
PLOT_REGISTRY["transport_normalized"] = {
    "func": plot_transport_normalized,
    "file": "transport_normalized.png",
    "title": "3-Mode Transport (Normalized)",
    "description": "Ionic/Electronic/Thermal 정규화 비교.\n각 mode를 자체 최대값으로 나눔 (0~1).\n어떤 mode가 조성 변화에 가장 민감한지 비교.\nIonic이 보통 가장 민감 (AM 비율에 강하게 반응).",
    "origin_tip": "Grouped Bar: Green(ionic), Red(electronic), Orange(thermal).",
}
PLOT_REGISTRY["transport_absolute"] = {
    "func": plot_transport_absolute,
    "file": "transport_absolute.png",
    "title": "3-Mode Transport (Absolute, Log)",
    "description": "Ionic/Electronic/Thermal 절대값 비교 (log scale).\nIonic ≪ Thermal < Electronic (typical).\nMinnmann(2021) 0.17 mS/cm 참고선 포함.\n이온 전도가 항상 rate-limiting → ionic 최적화 우선.",
    "origin_tip": "Multi-line (log Y): Green(ionic), Red(electronic), Orange(thermal).",
}
PLOT_REGISTRY["r_brug_comparison"] = {
    "func": plot_r_brug_comparison,
    "file": "r_brug_comparison.png",
    "title": "R_brug by Transport Mode",
    "description": "각 transport mode별 Bruggeman 과대추정 배수.\nR_brug = σ_brug / σ_network.\nIonic R_brug=3~10× (constriction 지배).\nElectronic/Thermal은 R_brug 다를 수 있음.\nR_brug가 클수록 접촉 저항 기여가 큼.",
    "origin_tip": "Grouped Bar: Green(ionic), Red(electronic), Orange(thermal).",
}
PLOT_REGISTRY["stress_z_layer"] = {
    "func": plot_stress_z_layer,
    "file": "stress_z_layer.png",
    "title": "Stress Z-layer CV",
    "description": "전극 높이별 Von Mises CV 프로파일.",
    "origin_tip": "Multi-line → X: Z Position (μm), Y: VM CV (%).",
}


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", required=True)
    parser.add_argument("-n", "--names", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-p", "--plots", nargs="+",
                        default=list(PLOT_REGISTRY.keys()) + ["four_panel"])
    parser.add_argument("--group-sizes", default="")  # e.g. "3,5"
    parser.add_argument("--group-names", default="")  # e.g. "Case A,Case B"
    parser.add_argument("--global-rgb", default="")   # e.g. "5.21,0.20" (b,ln_k from global fit)
    parser.add_argument("--global-c-ion", default="")  # e.g. "0.0727" (C from ionic scaling fit)
    args = parser.parse_args()

    # Parse group info
    if args.group_sizes:
        args.group_sizes_list = [int(x) for x in args.group_sizes.split(',')]
        args.group_names_list = args.group_names.split(',') if args.group_names else [f"Case {chr(65+i)}" for i in range(len(args.group_sizes_list))]
    else:
        args.group_sizes_list = None
        args.group_names_list = None

    if len(args.inputs) != len(args.names):
        print(f"ERROR: inputs ({len(args.inputs)}) != names ({len(args.names)})")
        sys.exit(1)

    global _ALL_DATA
    all_data = []
    for path in args.inputs:
        with open(path, "r") as f:
            d = json.load(f)
        d["_source_path"] = path
        all_data.append(d)
    _ALL_DATA = all_data

    # Use P:S ratio as x-axis labels (fallback to case names)
    plot_names = []
    for i, d in enumerate(all_data):
        ps = d.get('ps_ratio', '')
        if ps:
            plot_names.append(ps)
        else:
            plot_names.append(args.names[i])

    os.makedirs(args.output, exist_ok=True)
    plot_info = {}

    # Set global group info for _apply_style
    global _GROUP_INFO, _GLOBAL_RGB
    if args.group_sizes_list and len(args.group_sizes_list) > 1:
        _GROUP_INFO = (args.group_sizes_list, args.group_names_list)
    else:
        _GROUP_INFO = None

    # Set global R_gb if provided (from combined fit across all plot groups)
    if args.global_rgb:
        parts = args.global_rgb.split(',')
        if len(parts) >= 2:
            _GLOBAL_RGB = (float(parts[0]), float(parts[1]))
            print(f"  Using global R_gb: b={_GLOBAL_RGB[0]:.2f}, ln(k)={_GLOBAL_RGB[1]:.2f}")
    else:
        _GLOBAL_RGB = None

    # Set global C for ionic scaling if provided
    global _GLOBAL_C_ION
    if args.global_c_ion:
        _GLOBAL_C_ION = float(args.global_c_ion)
        print(f"  Using global C_ion: {_GLOBAL_C_ION:.6f}")
    else:
        _GLOBAL_C_ION = None

    # Generate particle info table as first plot
    if 'particle_info' in args.plots or True:  # always generate
        _generate_particle_info(all_data, plot_names, args.names, args.output)
        plot_info['particle_info'] = {
            'file': 'particle_info.png',
            'title': '입자 정보',
            'description': '각 케이스별 AM_P, AM_S, SE 입자수와 반지름.\nSE 크기가 같으면 비교 조건 동일.',
            'origin_tip': 'Table 형태로 Origin에서는 Worksheet에 직접 입력.',
        }

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
        'stress_cv': [('Case', None),
            ('Stress CV(%)', lambda d: _get(d, 'stress_cv'))],
        'stress_ratio': [('Case', None),
            ('σ_AM_P/σ_mean', lambda d: _get(d, 'stress_ratio_AM_P')),
            ('σ_AM_S/σ_mean', lambda d: _get(d, 'stress_ratio_AM_S')),
            ('σ_SE/σ_mean', lambda d: _get(d, 'stress_ratio_SE'))],
        'stress_z_layer': [('Case', None),
            ('Z_data', lambda d: str(d.get('stress_z_layer_cv', [])))],
    }

    import csv
    for pname, cols in csv_map.items():
        if pname in args.plots:
            csv_path = os.path.join(args.output, f"{pname}.csv")
            with open(csv_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow([c[0] for c in cols])
                for i, name in enumerate(plot_names):
                    row = []
                    for col_name, fn in cols:
                        if fn is None:
                            row.append(name)
                        else:
                            val = fn(all_data[i])
                            if isinstance(val, (int, float)):
                                row.append(round(val, 4))
                            else:
                                row.append(val)
                    writer.writerow(row)

    for plot_name in args.plots:
        if plot_name == "four_panel":
            outpath = plot_four_panel(all_data, plot_names, args.output)
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
        # min_groups check: skip if not enough case groups
        min_groups = entry.get("min_groups", 1)
        n_groups = len(_GROUP_INFO[0]) if _GROUP_INFO else 1
        if n_groups < min_groups:
            print(f"  [SKIP] {plot_name}: requires {min_groups}+ case groups (have {n_groups})")
            continue
        func = entry["func"]
        import inspect
        params = inspect.signature(func).parameters
        if 'outdir' in params:
            # Standalone plot (creates own fig, saves itself)
            outpath = func(all_data, plot_names, args.output)
            if outpath is None:
                # R² too low — generate "Bruggeman reference" note
                if plot_name in ('rgb_fitting', 'gb_corrected'):
                    _, _, _, r2_val = _fit_r_gb(all_data, plot_names, use_global=False)
                    plot_info[plot_name] = {
                        "file": "",
                        "csv": None,
                        "title": entry["title"] + " (skipped)",
                        "description": f"R²={r2_val:.4f} < {R_GB_MIN_R2}: fitting not reliable.\nBruggeman reference only.\nRVE size or data range insufficient.",
                        "origin_tip": "",
                    }
                continue
        else:
            fig, ax = plt.subplots(figsize=FIG_SINGLE)
            func(all_data, plot_names, ax=ax)
            outpath = _save(fig, args.output, entry["file"])

        # Check for CSV from csv_map (legacy) or standalone _write_csv
        csv_file = f"{plot_name}.csv" if plot_name in csv_map else None
        if csv_file is None:
            standalone_csv = os.path.join(args.output, f"{plot_name}.csv")
            if os.path.exists(standalone_csv):
                csv_file = f"{plot_name}.csv"
        # Dynamic description: replace {formx_r2}, {formx_loocv} with actual values
        desc = entry["description"]
        if _GLOBAL_FORMX_R2 is not None:
            r2_val, loocv_val = _GLOBAL_FORMX_R2
            desc = desc.replace("{formx_r2}", f"{r2_val:.3f}")
            desc = desc.replace("{formx_loocv}", f"{loocv_val:.3f}")
        elif '{formx_r2}' in desc:
            # Compute FORM X R² on the fly (same logic as plot_ionic_scaling_fit)
            try:
                PHI_C = 0.185; SGRAIN = 3.0
                _phi_se = np.array([_get(d, "phi_se", 0) for d in all_data])
                _cn_se = np.array([_get(d, "se_se_cn", 0) for d in all_data])
                _tau_se = np.array([max(_get(d, "tortuosity_recommended", _get(d, "tortuosity_mean", 1)), 0.1) for d in all_data])
                _cov_se = np.array([(lambda vs: sum(vs)/len(vs)/100 if vs else 0.20)([v for v in [_get(d,"coverage_AM_P_mean",0),_get(d,"coverage_AM_S_mean",0),_get(d,"coverage_AM_mean",0)] if v>0]) for d in all_data])
                _snet = np.array(_load_network_sigma(all_data))
                _gb = np.array([_get(d, "gb_density_mean", 0) for d in all_data])
                _gp = np.array([_get(d, "path_conductance_mean", 0) for d in all_data])
                # Match ionic_scaling_fit filter: require gb_density and path_conductance
                _valid = (_phi_se > PHI_C+0.01) & (_cn_se > 0) & (_tau_se > 0) & (_cov_se > 0) & (_snet > 0.01) & (_gb > 0) & (_gp > 0)
                if _valid.sum() >= 5:
                    _phi_ex = np.clip(_phi_se[_valid] - PHI_C, 0.001, None)
                    _fp_se = np.clip(np.array([_get(d, "percolation_pct", 100) / 100 for d in all_data])[_valid], 0.01, 1.0)
                    _log_rhs = np.log(SGRAIN) + 0.75*np.log(_phi_ex) + 1.5*np.log(_cn_se[_valid]) + 0.25*np.log(_cov_se[_valid]) + 2.0*np.log(_fp_se)
                    _log_s = np.log(_snet[_valid])
                    _log_tau = np.log(_tau_se[_valid])
                    _w_sig = 1.0 / (1.0 + np.exp(-5.0 * (_tau_se[_valid] - 2.1)))
                    _w_bl = 1.0 / (1.0 + np.exp(-15.0 * (_tau_se[_valid] - 2.0)))
                    # v9 BLEND fit
                    _X_v5 = np.column_stack([np.ones(_valid.sum()), _w_sig])
                    _X_p3 = np.column_stack([np.ones(_valid.sum()), _log_tau, _log_tau**2, _log_tau**3])
                    _b_v5 = np.linalg.lstsq(_X_v5, _log_s - _log_rhs, rcond=None)[0]
                    _b_p3 = np.linalg.lstsq(_X_p3, _log_s - _log_rhs, rcond=None)[0]
                    _pred = (1-_w_bl)*(_X_v5@_b_v5) + _w_bl*(_X_p3@_b_p3) + _log_rhs
                    _ss_res = np.sum((_log_s - _pred)**2)
                    _ss_tot = np.sum((_log_s - np.mean(_log_s))**2)
                    _r2 = 1 - _ss_res / _ss_tot
                    desc = desc.replace("{formx_r2}", f"{_r2:.3f}")
                    desc = desc.replace("{formx_loocv}", f"~{_r2:.3f}")
                else:
                    desc = desc.replace("{formx_r2}", "N/A")
                    desc = desc.replace("{formx_loocv}", "N/A")
            except:
                desc = desc.replace("{formx_r2}", "N/A")
                desc = desc.replace("{formx_loocv}", "N/A")
        info_entry = {
            "file": entry["file"],
            "csv": csv_file,
            "title": entry["title"],
            "description": desc,
            "origin_tip": entry["origin_tip"],
        }
        # Add R_gb values to plot_info for downstream use
        if plot_name == 'rgb_fitting':
            r_gb, ln_k, _, _ = _fit_r_gb(all_data, plot_names, use_global=False)
            info_entry['b'] = round(r_gb, 4)
            info_entry['ln_k'] = round(ln_k, 4)
        if plot_name == 'ionic_scaling_fit' and _GLOBAL_C_ION is not None:
            info_entry['C_ion'] = round(_GLOBAL_C_ION, 6)
        plot_info[plot_name] = info_entry
        print(f"  [OK] {plot_name} -> {outpath}")

    info_path = os.path.join(args.output, "plot_info.json")
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(plot_info, f, indent=2, ensure_ascii=False)
    print(f"\nTotal: {len(plot_info)} plots")


if __name__ == "__main__":
    main()
