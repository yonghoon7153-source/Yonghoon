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

def _apply_style(ax, ylabel, names):
    """Apply common academic style with group separators."""
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=10)
    ax.set_xlabel("P:S Configuration", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, color="#CCCCCC", alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis='both', labelsize=9)
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
            ax.text(mid, -0.18, gnames[gi], ha='center', va='top',
                    transform=ax.get_xaxis_transform(),
                    fontsize=9, fontweight='bold', color=GROUP_COLORS[gi % len(GROUP_COLORS)],
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=GROUP_COLORS[gi % len(GROUP_COLORS)], alpha=0.8))
            pos += sz


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
    tort = [_get(d, "tortuosity_mean") for d in all_data]

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

    ax1.plot(x, sigma_brug, 's-', color=GREEN, markersize=10, linewidth=2.5, label="σ_eff/σ_bulk")
    _apply_style(ax1, "σ_eff / σ_bulk", names)
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
    ax1.set_title("Bruggeman: σ_eff/σ_bulk = φ_SE × f_perc / τ²", fontsize=12, fontweight='bold')

    _write_csv(outdir, 'effective_conductivity.csv',
               ['φ_SE', 'σ_eff/σ_bulk'], names, phi, sigma_brug)
    return _save(fig, outdir, "effective_conductivity.png")


def _fit_r_gb(data_list, names):
    """Fit R_gb from Bruggeman vs Path Conductance data. Returns r_gb, log_k, valid data."""
    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    perc = [_get(d, "percolation_pct") / 100 for d in data_list]
    tau = [_get(d, "tortuosity_mean", 1) for d in data_list]
    g_path = [_get(d, "path_conductance_mean") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]

    sigma_proxy = [g_path[i] * perc[i] / tau[i] if g_path[i] > 0 and tau[i] > 0 else 0
                   for i in range(len(data_list))]

    valid_idx = [i for i in range(len(data_list))
                 if sigma_proxy[i] > 0 and gb_dens[i] > 0 and sigma_brug[i] > 0]

    r_gb = 0.5
    ln_k = 0
    if len(valid_idx) >= 2:
        # Linear regression with intercept: log(σ_brug/σ_proxy) = b × GB_d + ln(k)
        # slope = b (grain boundary resistance), intercept = ln(k) (unit conversion)
        y_vals = np.array([sigma_brug[i] / sigma_proxy[i] for i in valid_idx])
        x_vals = np.array([gb_dens[i] for i in valid_idx])
        log_y = np.log(y_vals)
        from scipy import stats as sp_stats
        slope, intercept, _, _, _ = sp_stats.linregress(x_vals, log_y)
        r_gb = slope
        ln_k = intercept
    return r_gb, ln_k, valid_idx


def plot_rgb_fitting(data_list, names, outdir):
    """R_gb fitting scatter: log(σ_brug/σ_proxy) vs GB_density."""
    r_gb, log_k, valid_idx = _fit_r_gb(data_list, names)

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    perc = [_get(d, "percolation_pct") / 100 for d in data_list]
    tau = [_get(d, "tortuosity_mean", 1) for d in data_list]
    g_path = [_get(d, "path_conductance_mean") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    sigma_proxy = [g_path[i] * perc[i] / tau[i] if g_path[i] > 0 and tau[i] > 0 else 0
                   for i in range(len(data_list))]

    fig, ax = plt.subplots(figsize=FIG_SINGLE)

    x_pts = np.array([gb_dens[i] for i in valid_idx])
    y_pts = np.array([np.log(sigma_brug[i] / sigma_proxy[i]) for i in valid_idx])

    # Determine group index for each valid point
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

    # Scatter with group colors
    if _GROUP_INFO:
        for j in range(len(valid_idx)):
            gi = point_groups[j]
            ax.scatter(x_pts[j], y_pts[j], s=100, c=GROUP_COLORS[gi % len(GROUP_COLORS)],
                      zorder=5, edgecolors='white', linewidth=1.5)
    else:
        ax.scatter(x_pts, y_pts, s=100, c=BLUE, zorder=5, edgecolors='white', linewidth=1.5)

    # Individual point labels (P:S ratio) - adjustText for non-overlapping
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

    # Group labels at bottom center of each cluster
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

    # Linear fit with intercept: log(y) = b × GB_d + ln(k)
    x_line = np.linspace(min(x_pts) * 0.9, max(x_pts) * 1.15, 100)
    y_line = r_gb * x_line + log_k
    ax.plot(x_line, y_line, '-', color=RED, linewidth=2.5,
            label=f"y = {r_gb:.2f}·x + {log_k:.2f}")

    # R² for linear regression with intercept
    y_pred = r_gb * x_pts + log_k
    ss_res = np.sum((y_pts - y_pred) ** 2)
    ss_tot = np.sum((y_pts - np.mean(y_pts)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    ax.set_xlabel("GB Density (hops/μm)", fontsize=12)
    ax.set_ylabel("log(σ_brug / σ_proxy)", fontsize=12)
    ax.set_title("R_gb Fitting", fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper left')
    ax.text(0.95, 0.05, f"b = {r_gb:.2f}\nln(k) = {log_k:.2f}\nR² = {r_squared:.4f}\nn = {len(x_pts)}",
            transform=ax.transAxes, fontsize=11, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8))
    # Y축 하단 여백 확보 (그룹 라벨 공간)
    y_margin = (max(y_pts) - min(y_pts)) * 0.15
    ax.set_ylim(min(y_pts) - y_margin, max(y_pts) + y_margin)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _write_csv(outdir, 'rgb_fitting.csv',
               ['GB_density', 'log(σ_brug/σ_proxy)', 'σ_brug', 'σ_proxy'],
               [names[i] for i in valid_idx],
               list(x_pts), list(y_pts),
               [sigma_brug[i] for i in valid_idx],
               [sigma_proxy[i] for i in valid_idx])
    return _save(fig, outdir, "rgb_fitting.png")


def plot_gb_corrected(data_list, names, outdir):
    """GB-corrected effective conductivity using fitted R_gb."""
    r_gb, _, _ = _fit_r_gb(data_list, names)
    SIGMA_BULK = 1.3  # mS/cm, Li₆PS₅Cl cold-pressed

    sigma_brug = [_get(d, "sigma_ratio") for d in data_list]
    gb_dens = [_get(d, "gb_density_mean") for d in data_list]
    phi = [_get(d, "phi_se") for d in data_list]
    perc = [_get(d, "percolation_pct") / 100 for d in data_list]
    tau = [_get(d, "tortuosity_mean", 1) for d in data_list]

    # σ_eff_real = σ_brug × e^(-b × GB_d)
    # GB_d = 0: no GB → σ_eff_real = σ_brug (Bruggeman exact)
    sigma_corr = [sigma_brug[i] * np.exp(-r_gb * gb_dens[i]) if gb_dens[i] > 0 else sigma_brug[i]
                  for i in range(len(data_list))]

    fig, ax = plt.subplots(figsize=FIG_SINGLE)
    x = np.arange(len(names))

    ax.plot(x, sigma_corr, 's-', color=RED, markersize=10, linewidth=2.5, label="σ_eff_real/σ_bulk")
    _apply_style(ax, "σ_eff_real / σ_bulk", names)
    ax.legend(fontsize=9, loc='best')
    ax.set_title(f"GB-Corrected σ_eff/σ_bulk  (b={r_gb:.2f})",
                 fontsize=12, fontweight='bold')

    _write_csv(outdir, 'gb_corrected.csv',
               ['φ_SE', 'f_perc', 'τ', 'GB_density', 'R_gb',
                'σ_brug/σ_bulk', 'σ_corr/σ_bulk'],
               names, phi, perc, tau, gb_dens,
               [r_gb]*len(names), sigma_brug, sigma_corr)
    return _save(fig, outdir, "gb_corrected.png")


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
        "title": "Effective Conductivity (Bruggeman)",
        "description": "σ_eff/σ_bulk = φ_SE × f_perc / τ²\n\nSE 부피 분율(φ_SE), percolation(f_perc), tortuosity(τ) 기반 추정.\n입계(GB) 저항을 무시한 이론값 → 실제보다 과대평가 가능.",
        "origin_tip": "Line+Symbol (Green) + Bar (φ_SE, Blue).",
    },
    "rgb_fitting": {
        "func": plot_rgb_fitting,
        "file": "rgb_fitting.png",
        "title": "R_gb Fitting",
        "description": "log(σ_brug/σ_proxy) = b × GB_d + ln(k)\n\n선형 회귀로 b(기울기=입계 저항 계수)와 ln(k)(절편=단위 변환)를 분리.\nX축: GB Density, Y축: log(σ_brug/σ_proxy)\nR²가 높으면 입계 저항이 GB_density에 비례함을 확인.\n\nb는 한 번 구하면 같은 재료/조건에서 재사용 가능.",
        "origin_tip": "Scatter + Fit line (log Y scale).\nBlue dots: data, Red line: fitted.",
    },
    "gb_corrected": {
        "func": plot_gb_corrected,
        "file": "gb_corrected.png",
        "title": "GB-Corrected σ_eff",
        "description": "σ_eff_real/σ_bulk = φ_SE × f_perc / τ² × e^(-b × GB_d)\n\nb: R_gb fitting에서 결정된 입계 저항 계수\nGB_d=0 → Bruggeman 정확 | GB_d↑ → 지수적 감소\n\n⚠ 케이스 간 상대 비교용.\nDEM 구형 입자 + Hertz 접촉 한계로\n절대값(mS/cm)은 실측과 차이가 클 수 있음.",
        "origin_tip": "Line (Red, σ_corr/σ_bulk) + Dashed (Orange, mS/cm).",
    },
    "ion_path_quality": {
        "func": plot_ion_path_quality,
        "file": "ion_path_quality.png",
        "title": "Ion Path Quality (G_path, GB_d)",

        "description": "이온 경로 품질 4종 (모두 percolating 경로들의 mean):\n\n• GB Density = mean(N_hops / L_z) (hops/μm)\n  → ↓ 좋음 (입계 적을수록 저항↓)\n• Path Hop Area = mean(각 hop의 접촉면적)\n  → ↑ 좋음\n• Bottleneck = mean(각 경로의 최소 접촉면적)\n  → ↑ 좋음\n• Path Conductance = mean(1/Σ(1/A_i)) (μm²)\n  → ↑ 좋음 (직렬 저항 모델의 유효 면적)\n\n왜 mean? 경로 안은 직렬, 경로끼리는 병렬.\nmean = 경로 품질, f_perc = 경로 수량 (역할 분리).\n30개 shortest path 샘플 (best/mean/worst 10개씩).\nConductance가 가장 종합적 지표.",
        "origin_tip": "2×2 Subplots → GB Density (Blue), Hop Area (Orange), Bottleneck (Red), Conductance (Green).",
    },
    "stress_z_layer": {
        "func": plot_stress_z_layer,
        "file": "stress_z_layer.png",
        "title": "Stress Z-layer CV",
        "description": "전극 높이별 Von Mises CV 프로파일.\n\n균일 압축 → 전 층에서 CV 비슷.\n경계 효과 → 상/하단 CV 증가.\nAM_P↑ → force chain으로 Z 관통 → 균일.",
        "origin_tip": "Multi-line → X: Z Position (μm), Y: VM CV (%).\n각 케이스별 다른 색.\n경계(상/하단) 영역 음영 표시 권장.",
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
    parser.add_argument("--group-sizes", default="")  # e.g. "3,5"
    parser.add_argument("--group-names", default="")  # e.g. "Case A,Case B"
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

    all_data = []
    for path in args.inputs:
        with open(path, "r") as f:
            all_data.append(json.load(f))

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
    global _GROUP_INFO
    if args.group_sizes_list and len(args.group_sizes_list) > 1:
        _GROUP_INFO = (args.group_sizes_list, args.group_names_list)
    else:
        _GROUP_INFO = None

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
        func = entry["func"]
        import inspect
        params = inspect.signature(func).parameters
        if 'outdir' in params:
            # Standalone plot (creates own fig, saves itself)
            outpath = func(all_data, plot_names, args.output)
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
