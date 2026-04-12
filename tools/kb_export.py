#!/usr/bin/env python3
"""
Knowledge Base Export Tool — generate tables and formatted output for paper writing.

Usage:
    python tools/kb_export.py --table eos
    python tools/kb_export.py --table elastic
    python tools/kb_export.py --table electronic
    python tools/kb_export.py --table adhesion
    python tools/kb_export.py --table bonds
    python tools/kb_export.py --table all
    python tools/kb_export.py --format latex --table eos
    python tools/kb_export.py --format csv --table elastic
    python tools/kb_export.py --summary
    python tools/kb_export.py --comparison comp1 comp5
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "db"


def load_json(path: Path) -> dict:
    try:
        with open(path) as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        return {}


def load_all_comps() -> Dict[str, dict]:
    comps = {}
    for fp in sorted((DB / "compositions").glob("*.json")):
        data = load_json(fp)
        if data and "id" in data:
            comps[data["id"]] = data
    return comps


def load_props(name: str) -> dict:
    return load_json(DB / "properties" / f"{name}.json")


def fmt(v, digits=1):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


# ── Table generators ─────────────────────────────────────────────────

def table_eos(comps: dict, format: str = "markdown") -> str:
    prop = load_props("eos")
    rows = []
    if "results" in prop:
        for r in prop["results"]:
            rows.append({
                "id": r["id"],
                "formula": comps.get(r["id"], {}).get("formula", r["id"]),
                "B0": fmt(r.get("B0_GPa")),
                "B0p": fmt(r.get("B0_prime"), 2),
                "pts": str(r.get("pts", "—")),
                "R2": fmt(r.get("R2"), 6),
                "notes": r.get("notes", r.get("basin", "")),
            })

    if format == "latex":
        lines = [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\caption{Birch-Murnaghan EOS parameters for argyrodite compositions.}",
            r"\label{tab:eos}",
            r"\begin{tabular}{llccccl}",
            r"\hline",
            r"ID & Formula & $B_0$ (GPa) & $B_0'$ & Points & $R^2$ & Notes \\",
            r"\hline",
        ]
        for r in rows:
            lines.append(f"{r['id']} & {r['formula']} & {r['B0']} & {r['B0p']} & {r['pts']} & {r['R2']} & {r['notes']} \\\\")
        lines.extend([r"\hline", r"\end{tabular}", r"\end{table}"])
        return "\n".join(lines)

    elif format == "csv":
        lines = ["id,formula,B0_GPa,B0_prime,points,R2,notes"]
        for r in rows:
            lines.append(f"{r['id']},{r['formula']},{r['B0']},{r['B0p']},{r['pts']},{r['R2']},{r['notes']}")
        return "\n".join(lines)

    else:  # markdown
        lines = [
            "| ID | Formula | B0 (GPa) | B0' | Points | R^2 | Notes |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            lines.append(f"| {r['id']} | {r['formula']} | {r['B0']} | {r['B0p']} | {r['pts']} | {r['R2']} | {r['notes']} |")
        return "\n".join(lines)


def table_elastic(comps: dict, format: str = "markdown") -> str:
    prop = load_props("elastic")
    rows = []
    if "mlip_600K_snapshot" in prop and "results" in prop["mlip_600K_snapshot"]:
        for r in prop["mlip_600K_snapshot"]["results"]:
            rows.append({
                "id": r["id"],
                "formula": comps.get(r["id"], {}).get("formula", r["id"]),
                "E": fmt(r.get("E")),
                "E_std": fmt(r.get("E_std")),
                "K": fmt(r.get("K")),
                "G": fmt(r.get("G")),
                "expt": fmt(r.get("expt_E")) if r.get("expt_E") else "—",
            })

    if format == "latex":
        lines = [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\caption{MLIP 600K snapshot elastic constants for argyrodite compositions.}",
            r"\label{tab:elastic}",
            r"\begin{tabular}{llccccc}",
            r"\hline",
            r"ID & Formula & $E$ (GPa) & $\sigma_E$ & $K$ (GPa) & $G$ (GPa) & Expt $E$ \\",
            r"\hline",
        ]
        for r in rows:
            lines.append(f"{r['id']} & {r['formula']} & {r['E']} & {r['E_std']} & {r['K']} & {r['G']} & {r['expt']} \\\\")
        lines.extend([r"\hline", r"\end{tabular}", r"\end{table}"])
        return "\n".join(lines)

    elif format == "csv":
        lines = ["id,formula,E_GPa,E_std,K_GPa,G_GPa,expt_E"]
        for r in rows:
            lines.append(f"{r['id']},{r['formula']},{r['E']},{r['E_std']},{r['K']},{r['G']},{r['expt']}")
        return "\n".join(lines)

    else:
        lines = [
            "| ID | Formula | E (GPa) | E_std | K (GPa) | G (GPa) | Expt E |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            lines.append(f"| {r['id']} | {r['formula']} | {r['E']} | {r['E_std']} | {r['K']} | {r['G']} | {r['expt']} |")
        return "\n".join(lines)


def table_electronic(comps: dict, format: str = "markdown") -> str:
    prop = load_props("electronic")
    rows = []
    if "band_gaps" in prop:
        for r in prop["band_gaps"]:
            rows.append({
                "id": r["id"],
                "gap": fmt(r.get("gap_eV"), 2),
                "vbm": r.get("vbm", "—"),
                "cbm": r.get("cbm", "—"),
            })

    if format == "latex":
        lines = [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\caption{Electronic structure of argyrodite compositions.}",
            r"\label{tab:electronic}",
            r"\begin{tabular}{lcll}",
            r"\hline",
            r"ID & Gap (eV) & VBM & CBM \\",
            r"\hline",
        ]
        for r in rows:
            lines.append(f"{r['id']} & {r['gap']} & {r['vbm']} & {r['cbm']} \\\\")
        lines.extend([r"\hline", r"\end{tabular}", r"\end{table}"])
        return "\n".join(lines)
    else:
        lines = [
            "| ID | Gap (eV) | VBM | CBM |",
            "|---|---|---|---|",
        ]
        for r in rows:
            lines.append(f"| {r['id']} | {r['gap']} | {r['vbm']} | {r['cbm']} |")
        return "\n".join(lines)


def table_adhesion(comps: dict, format: str = "markdown") -> str:
    prop = load_props("adhesion")
    rows = []

    if "surface_energies" in prop:
        for r in prop["surface_energies"]:
            v1 = None
            if "adhesion_v1_random" in prop:
                for a in prop["adhesion_v1_random"]:
                    if a["id"] == r["id"]:
                        v1 = a
                        break
            rows.append({
                "id": r["id"],
                "formula": comps.get(r["id"], {}).get("formula", r["id"]),
                "gamma": fmt(r.get("gamma_SE"), 3),
                "wad": fmt(v1.get("Wad_Jm2"), 3) if v1 else "—",
                "expt_aJ": fmt(v1.get("expt_aJ")) if v1 and v1.get("expt_aJ") else "—",
            })

    if format == "latex":
        lines = [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\caption{Surface and adhesion energies of argyrodite compositions.}",
            r"\label{tab:adhesion}",
            r"\begin{tabular}{llccc}",
            r"\hline",
            r"ID & Formula & $\gamma_{SE}$ (J/m$^2$) & $W_{ad}$ (J/m$^2$) & Expt (aJ) \\",
            r"\hline",
        ]
        for r in rows:
            lines.append(f"{r['id']} & {r['formula']} & {r['gamma']} & {r['wad']} & {r['expt_aJ']} \\\\")
        lines.extend([r"\hline", r"\end{tabular}", r"\end{table}"])
        return "\n".join(lines)
    else:
        lines = [
            "| ID | Formula | gamma_SE (J/m2) | Wad (J/m2) | Expt (aJ) |",
            "|---|---|---|---|---|",
        ]
        for r in rows:
            lines.append(f"| {r['id']} | {r['formula']} | {r['gamma']} | {r['wad']} | {r['expt_aJ']} |")
        return "\n".join(lines)


def table_bonds(comps: dict, format: str = "markdown") -> str:
    prop = load_props("bonds")
    lines = [
        "| Composition | Li-Cl (A) | Li-Br (A) | Li-S (A) | P-S (A) |",
        "|---|---|---|---|---|",
    ]
    if "results" in prop:
        for cid, bonds in sorted(prop["results"].items()):
            lines.append(
                f"| {cid} | {fmt(bonds.get('Li-Cl'), 3)} | {fmt(bonds.get('Li-Br'), 3)} "
                f"| {fmt(bonds.get('Li-S'), 3)} | {fmt(bonds.get('P-S'), 3)} |"
            )
    return "\n".join(lines)


def generate_summary(comps: dict) -> str:
    """Generate a one-page summary of all results."""
    lines = [
        "# Argyrodite Mechanical Properties — Results Summary",
        "",
        "## Compositions",
        "| ID | Formula | Family | Atoms | Cell |",
        "|---|---|---|---|---|",
    ]
    for cid, c in sorted(comps.items()):
        lines.append(f"| {cid} | {c.get('formula','?')} | {c.get('family','?')} | {c.get('atoms','?')} | {c.get('cell_type','?')} |")

    lines.extend(["", "## Bulk Modulus (DFT EOS)", ""])
    lines.append(table_eos(comps))

    lines.extend(["", "## Young's Modulus (MLIP 600K Snapshot)", ""])
    lines.append(table_elastic(comps))

    lines.extend(["", "## Band Gap", ""])
    lines.append(table_electronic(comps))

    lines.extend(["", "## Surface & Adhesion Energy", ""])
    lines.append(table_adhesion(comps))

    lines.extend(["", "## Bond Lengths", ""])
    lines.append(table_bonds(comps))

    lines.extend([
        "",
        "## Key Trends",
        "- **B0**: Br↑ → B0↓ within each family. Vacancy → B0↓ across families.",
        "- **E (Young's)**: 29.1 > 28.6 > 27.3 > 26.4 > 25.8 GPa — perfectly monotonic with Br.",
        "- **Band gap**: 2.0-2.3 eV, Br slightly reduces. PS4 framework invariant.",
        "- **Adhesion**: Vacancy → chemical anchor → Wad↑ for Li5.4 family.",
        "- **Li ordering**: comp5 ΔC44=12.7 GPa (47%) — shear most sensitive.",
    ])

    return "\n".join(lines)


def comparison(comps: dict, ids: List[str]) -> str:
    """Side-by-side comparison of specific compositions."""
    lines = [f"# Comparison: {' vs '.join(ids)}", ""]

    for cid in ids:
        if cid not in comps:
            lines.append(f"## {cid}: NOT FOUND")
            continue
        c = comps[cid]
        lines.append(f"## {cid}: {c.get('formula', '?')}")
        lines.append(f"- Family: {c.get('family', '?')}, Atoms: {c.get('atoms', '?')}, Cell: {c.get('cell_type', '?')}")
        lines.append(f"- Vacancy: {'Yes' if c.get('has_vacancy') else 'No'}")

        if "eos" in c:
            eos = c["eos"]
            if isinstance(eos, dict) and "B0_GPa" in eos:
                lines.append(f"- B0: {eos['B0_GPa']} GPa")
            elif isinstance(eos, dict) and "basin_A" in eos:
                lines.append(f"- B0 (Basin A): {eos['basin_A'].get('B0_GPa', '?')} GPa")
                if "basin_B" in eos:
                    lines.append(f"- B0 (Basin B): {eos['basin_B'].get('B0_GPa', '?')} GPa")

        if "band_gap_eV" in c:
            bg = c["band_gap_eV"]
            if isinstance(bg, dict):
                for k, v in bg.items():
                    lines.append(f"- Band gap ({k}): {v} eV")
            else:
                lines.append(f"- Band gap: {bg} eV")

        if "elastic_mlip_600K_GPa" in c:
            el = c["elastic_mlip_600K_GPa"]
            if isinstance(el, dict):
                e = el.get("E_VRH", el.get("basin_B", {}).get("E_VRH", "?"))
                lines.append(f"- E (600K snap): {e} GPa")

        if "surface_energy_Jm2" in c:
            lines.append(f"- Surface energy: {c['surface_energy_Jm2']} J/m2")

        lines.append("")

    return "\n".join(lines)


TABLE_FUNCS = {
    "eos": table_eos,
    "elastic": table_elastic,
    "electronic": table_electronic,
    "adhesion": table_adhesion,
    "bonds": table_bonds,
}


def main():
    parser = argparse.ArgumentParser(description="Argyrodite Knowledge Base Export")
    parser.add_argument("--table", choices=list(TABLE_FUNCS.keys()) + ["all"],
                        help="Generate a specific table")
    parser.add_argument("--format", choices=["markdown", "latex", "csv"], default="markdown")
    parser.add_argument("--summary", action="store_true", help="Generate full summary")
    parser.add_argument("--comparison", nargs="+", metavar="COMP_ID",
                        help="Compare specific compositions")
    parser.add_argument("-o", "--output", help="Output file path")

    args = parser.parse_args()
    comps = load_all_comps()

    output = ""

    if args.summary:
        output = generate_summary(comps)
    elif args.comparison:
        output = comparison(comps, args.comparison)
    elif args.table:
        if args.table == "all":
            parts = []
            for name, func in TABLE_FUNCS.items():
                parts.append(f"## {name.upper()}\n")
                parts.append(func(comps, args.format))
                parts.append("")
            output = "\n".join(parts)
        else:
            output = TABLE_FUNCS[args.table](comps, args.format)
    else:
        parser.print_help()
        return

    if args.output:
        Path(args.output).write_text(output)
        print(f"Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
