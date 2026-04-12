#!/usr/bin/env python3
"""
Knowledge Base Indexer — build and maintain the search index.

Features:
  1. Validates all JSON files in db/
  2. Builds a flat index of all data points for fast search
  3. Checks cross-references between compositions, properties, and literature
  4. Reports missing data and inconsistencies

Usage:
    python tools/kb_index.py --rebuild       # Full rebuild
    python tools/kb_index.py --validate      # Validate only
    python tools/kb_index.py --check-refs    # Check cross-references
    python tools/kb_index.py --summary       # Print database summary
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Any, Dict, List, Tuple
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "db"
KB = ROOT / "kb"
INDEX_FILE = DB / "_index.json"


def load_json_safe(path: Path) -> Tuple[dict, List[str]]:
    errors = []
    try:
        with open(path) as f:
            data = json.load(f)
        return data, errors
    except json.JSONDecodeError as e:
        errors.append(f"JSON parse error in {path}: {e}")
        return {}, errors
    except FileNotFoundError:
        errors.append(f"File not found: {path}")
        return {}, errors


def flatten(obj: Any, prefix: str = "") -> List[Tuple[str, Any]]:
    items = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            items.extend(flatten(v, f"{prefix}.{k}" if prefix else k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            items.extend(flatten(v, f"{prefix}[{i}]"))
    else:
        items.append((prefix, obj))
    return items


def validate_compositions() -> Tuple[int, List[str]]:
    """Validate all composition JSON files."""
    errors = []
    count = 0
    required_fields = ["id", "formula", "family"]

    for fp in sorted((DB / "compositions").glob("*.json")):
        data, parse_errors = load_json_safe(fp)
        errors.extend(parse_errors)
        if not data:
            continue
        count += 1

        for field in required_fields:
            if field not in data:
                errors.append(f"{fp.name}: missing required field '{field}'")

        # Sanity checks
        if "eos" in data:
            eos = data["eos"]
            if isinstance(eos, dict):
                for basin_key in ["B0_GPa", "basin_A", "basin_B"]:
                    if basin_key in eos:
                        b0_data = eos[basin_key] if basin_key != "B0_GPa" else eos
                        if isinstance(b0_data, dict) and "B0_GPa" in b0_data:
                            b0 = b0_data["B0_GPa"]
                            if not (5 < b0 < 200):
                                errors.append(f"{fp.name}: B0={b0} GPa out of range (5-200)")

    return count, errors


def validate_properties() -> Tuple[int, List[str]]:
    errors = []
    count = 0
    for fp in sorted((DB / "properties").glob("*.json")):
        data, parse_errors = load_json_safe(fp)
        errors.extend(parse_errors)
        if data:
            count += 1
            if "property" not in data:
                errors.append(f"{fp.name}: missing 'property' field")
    return count, errors


def validate_literature() -> Tuple[int, List[str]]:
    errors = []
    count = 0
    data, parse_errors = load_json_safe(DB / "literature" / "refs.json")
    errors.extend(parse_errors)
    if data and "references" in data:
        for ref in data["references"]:
            count += 1
            for field in ["id", "year"]:
                if field not in ref:
                    errors.append(f"refs.json: reference missing '{field}' — {ref.get('id', '?')}")
    return count, errors


def check_cross_references() -> List[str]:
    """Check that all referenced ref IDs exist in literature."""
    issues = []

    # Load all ref IDs
    data, _ = load_json_safe(DB / "literature" / "refs.json")
    ref_ids = set()
    if data and "references" in data:
        for r in data["references"]:
            ref_ids.add(r.get("id", ""))

    # Check composition files
    for fp in sorted((DB / "compositions").glob("*.json")):
        comp_data, _ = load_json_safe(fp)
        if not comp_data:
            continue
        if "references" in comp_data:
            for rid in comp_data["references"]:
                if rid not in ref_ids:
                    issues.append(f"{fp.name}: references unknown ref '{rid}'")

    return issues


def build_index() -> dict:
    """Build flat search index from all data."""
    index = {
        "built": datetime.now().isoformat(),
        "compositions": {},
        "properties": {},
        "literature_count": 0,
        "kb_articles": [],
        "data_points": [],
    }

    # Index compositions
    for fp in sorted((DB / "compositions").glob("*.json")):
        data, _ = load_json_safe(fp)
        if data and "id" in data:
            cid = data["id"]
            index["compositions"][cid] = {
                "formula": data.get("formula", "?"),
                "family": data.get("family", "?"),
                "atoms": data.get("atoms", "?"),
                "cell_type": data.get("cell_type", "?"),
                "has_eos": "eos" in data,
                "has_elastic": any("elastic" in k for k in data.keys()),
                "has_band_gap": "band_gap_eV" in data,
                "has_adhesion": any("adhesion" in k for k in data.keys()),
            }
            # Flatten all data points
            for path, val in flatten(data):
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    index["data_points"].append({
                        "comp": cid,
                        "path": path,
                        "value": val,
                    })

    # Index properties
    for fp in sorted((DB / "properties").glob("*.json")):
        data, _ = load_json_safe(fp)
        if data and "property" in data:
            index["properties"][data["property"]] = {
                "n_results": len(data.get("results", [])),
                "has_trends": "trends" in data,
            }

    # Index literature
    data, _ = load_json_safe(DB / "literature" / "refs.json")
    if data:
        index["literature_count"] = len(data.get("references", []))

    # Index KB articles
    if KB.exists():
        for md in sorted(KB.rglob("*.md")):
            text = md.read_text(errors="replace")
            index["kb_articles"].append({
                "path": str(md.relative_to(ROOT)),
                "title": text.splitlines()[0].lstrip("# ") if text.strip() else md.stem,
                "words": len(text.split()),
            })

    return index


def print_summary(index: dict):
    print("\n" + "=" * 72)
    print("  ARGYRODITE KNOWLEDGE BASE — SUMMARY")
    print("=" * 72)

    print(f"\n  Built: {index.get('built', '?')}")

    # Compositions
    print(f"\n  Compositions ({len(index['compositions'])}):")
    for cid, info in sorted(index["compositions"].items()):
        flags = []
        if info.get("has_eos"): flags.append("EOS")
        if info.get("has_elastic"): flags.append("Elastic")
        if info.get("has_band_gap"): flags.append("Gap")
        if info.get("has_adhesion"): flags.append("Wad")
        print(f"    {cid:<10} {info['formula']:<30} [{', '.join(flags)}]")

    # Properties
    print(f"\n  Property Tables ({len(index['properties'])}):")
    for pid, info in sorted(index["properties"].items()):
        print(f"    {pid:<25} {info['n_results']} entries")

    # Literature
    print(f"\n  Literature References: {index['literature_count']}")

    # KB Articles
    print(f"\n  Knowledge Base Articles ({len(index['kb_articles'])}):")
    for art in index["kb_articles"]:
        print(f"    {art['path']:<50} ({art['words']} words)")

    # Data points
    n_dp = len(index.get("data_points", []))
    print(f"\n  Total numeric data points: {n_dp}")

    print()


def main():
    parser = argparse.ArgumentParser(description="Argyrodite Knowledge Base Indexer")
    parser.add_argument("--rebuild", action="store_true", help="Rebuild search index")
    parser.add_argument("--validate", action="store_true", help="Validate all data files")
    parser.add_argument("--check-refs", action="store_true", help="Check cross-references")
    parser.add_argument("--summary", action="store_true", help="Print database summary")

    args = parser.parse_args()

    if not any([args.rebuild, args.validate, args.check_refs, args.summary]):
        args.summary = True  # Default action

    all_errors = []

    if args.validate or args.rebuild:
        print("Validating data files...")
        n_comp, comp_errs = validate_compositions()
        n_prop, prop_errs = validate_properties()
        n_ref, ref_errs = validate_literature()
        all_errors.extend(comp_errs + prop_errs + ref_errs)
        print(f"  Compositions: {n_comp} files, {len(comp_errs)} errors")
        print(f"  Properties:   {n_prop} files, {len(prop_errs)} errors")
        print(f"  Literature:   {n_ref} refs, {len(ref_errs)} errors")

    if args.check_refs or args.rebuild:
        print("Checking cross-references...")
        xref_issues = check_cross_references()
        all_errors.extend(xref_issues)
        print(f"  Cross-ref issues: {len(xref_issues)}")

    if args.rebuild:
        print("Building index...")
        index = build_index()
        with open(INDEX_FILE, "w") as f:
            json.dump(index, f, indent=2, ensure_ascii=False)
        print(f"  Index written to {INDEX_FILE}")
        print_summary(index)

    if args.summary and not args.rebuild:
        if INDEX_FILE.exists():
            with open(INDEX_FILE) as f:
                index = json.load(f)
            print_summary(index)
        else:
            print("No index found. Run with --rebuild first.")
            index = build_index()
            print_summary(index)

    if all_errors:
        print(f"\n{'!' * 72}")
        print(f"  ERRORS ({len(all_errors)}):")
        for e in all_errors:
            print(f"    - {e}")
        print()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
