#!/usr/bin/env python3
"""
Knowledge Base Search Tool — query the argyrodite research database.

Searches across:
  1. db/compositions/*.json  — structured computational results
  2. db/properties/*.json    — cross-composition property tables
  3. db/literature/refs.json — literature references
  4. kb/**/*.md              — knowledge base articles
  5. db/pipelines/*.json     — pipeline definitions

Usage:
    python tools/kb_search.py "comp5 bulk modulus"
    python tools/kb_search.py "Br effect elastic"
    python tools/kb_search.py "basin analysis"
    python tools/kb_search.py "adhesion vacancy"
    python tools/kb_search.py --comp comp5 --prop elastic
    python tools/kb_search.py --ref pustorino2025
    python tools/kb_search.py --list-comps
    python tools/kb_search.py --list-props
"""

import json
import os
import re
import sys
import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "db"
KB = ROOT / "kb"

# ── Composition aliases ──────────────────────────────────────────────
COMP_ALIASES = {
    "comp1": "comp1", "li6ps5cl": "comp1", "lpscl": "comp1",
    "comp2": "comp2", "li6ps5cl0.5br0.5": "comp2",
    "comp3": "comp3", "li5.4ps4.4cl1.0br0.6": "comp3",
    "comp4": "comp4", "li5.4ps4.4cl0.8br0.8": "comp4",
    "comp5": "comp5", "li5.4ps4.4cl0.6br1.0": "comp5",
    "modelc": "modelc", "model_c": "modelc", "model c": "modelc",
    "li5.4ps4.4cl1.6": "modelc",
}

# ── Property keywords ────────────────────────────────────────────────
PROP_KEYWORDS = {
    "eos": ["eos", "bulk modulus", "b0", "birch", "murnaghan", "equation of state", "volume"],
    "elastic": ["elastic", "young", "shear", "c11", "c12", "c44", "cij", "modulus", "stiffness"],
    "electronic": ["band gap", "dos", "pdos", "electronic", "bader", "charge"],
    "bonds": ["bond", "bond length", "li-cl", "li-br", "li-s", "p-s", "distance"],
    "adhesion": ["adhesion", "wad", "surface energy", "interface", "ncm", "gamma"],
}

TOPIC_KEYWORDS = {
    "basin": ["basin", "li ordering", "ordering sensitivity", "multi-basin"],
    "vacancy": ["vacancy", "compositional vacancy", "li5.4", "halogen-rich"],
    "br_effect": ["br effect", "br substitution", "halogen substitution", "cl to br"],
    "pipeline": ["pipeline", "annealing", "screening", "enumerate"],
    "intrinsic_extrinsic": ["intrinsic", "extrinsic", "grain boundary", "gb pinning"],
}


def load_json(path: Path) -> Optional[dict]:
    try:
        with open(path) as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        return None


def load_all_compositions() -> Dict[str, dict]:
    comps = {}
    for fp in sorted((DB / "compositions").glob("*.json")):
        data = load_json(fp)
        if data and "id" in data:
            comps[data["id"]] = data
    return comps


def load_all_properties() -> Dict[str, dict]:
    props = {}
    for fp in sorted((DB / "properties").glob("*.json")):
        data = load_json(fp)
        if data and "property" in data:
            props[data["property"]] = data
    return props


def load_literature() -> Dict[str, dict]:
    refs = {}
    data = load_json(DB / "literature" / "refs.json")
    if data and "references" in data:
        for r in data["references"]:
            refs[r["id"]] = r
    return refs


def load_pipelines() -> Dict[str, dict]:
    pipes = {}
    for fp in sorted((DB / "pipelines").glob("*.json")):
        data = load_json(fp)
        if data and "id" in data:
            pipes[data["id"]] = data
    return pipes


def search_kb_markdown(terms: List[str], limit: int = 5) -> List[dict]:
    """Full-text search across all kb/ markdown files."""
    results = []
    if not KB.exists():
        return results
    pattern = re.compile("|".join(re.escape(t) for t in terms), re.IGNORECASE)
    for md in KB.rglob("*.md"):
        text = md.read_text(errors="replace")
        matches = list(pattern.finditer(text))
        if matches:
            lines = text.splitlines()
            snippets = []
            seen_lines = set()
            for m in matches[:5]:
                ln = text[:m.start()].count("\n")
                if ln in seen_lines:
                    continue
                seen_lines.add(ln)
                start = max(0, ln - 1)
                end = min(len(lines), ln + 3)
                snippets.append((ln + 1, "\n".join(lines[start:end])))
            results.append({
                "path": str(md.relative_to(ROOT)),
                "n_matches": len(matches),
                "snippets": snippets[:3],
            })
    return sorted(results, key=lambda r: r["n_matches"], reverse=True)[:limit]


def flatten_json(obj: Any, prefix: str = "") -> List[Tuple[str, Any]]:
    """Flatten nested dict/list into (dotpath, value) pairs."""
    items = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            new_key = f"{prefix}.{k}" if prefix else k
            items.extend(flatten_json(v, new_key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            items.extend(flatten_json(v, f"{prefix}[{i}]"))
    else:
        items.append((prefix, obj))
    return items


def search_in_data(data: dict, terms: List[str]) -> List[Tuple[str, Any]]:
    """Search flattened data for term matches in keys or string values."""
    flat = flatten_json(data)
    results = []
    for path, val in flat:
        score = 0
        text = f"{path} {val}" if isinstance(val, str) else path
        text_lower = text.lower()
        for t in terms:
            if t.lower() in text_lower:
                score += 1
        if score > 0:
            results.append((path, val, score))
    return sorted(results, key=lambda x: -x[2])


def detect_comp(query: str) -> Optional[str]:
    q = query.lower()
    for alias, cid in COMP_ALIASES.items():
        if alias in q:
            return cid
    return None


def detect_prop(query: str) -> Optional[str]:
    q = query.lower()
    for prop, kws in PROP_KEYWORDS.items():
        for kw in kws:
            if kw in q:
                return prop
    return None


def detect_topic(query: str) -> Optional[str]:
    q = query.lower()
    for topic, kws in TOPIC_KEYWORDS.items():
        for kw in kws:
            if kw in q:
                return topic
    return None


def print_section(title: str, char: str = "─"):
    w = 72
    print(f"\n{char * 3} {title} {char * max(1, w - len(title) - 5)}")


def format_value(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, dict):
        return json.dumps(v, indent=2, ensure_ascii=False)
    if isinstance(v, list) and len(v) <= 6:
        return str(v)
    return str(v)


def cmd_search(args):
    query = " ".join(args.terms)
    terms = [t for t in query.split() if len(t) > 1]

    comp_id = args.comp or detect_comp(query)
    prop_id = args.prop or detect_prop(query)
    topic = detect_topic(query)

    comps = load_all_compositions()
    props = load_all_properties()
    refs = load_literature()
    pipes = load_pipelines()

    print(f"{'=' * 72}")
    print(f"  Query: {query}")
    if comp_id:
        print(f"  Composition: {comp_id}")
    if prop_id:
        print(f"  Property: {prop_id}")
    if topic:
        print(f"  Topic: {topic}")
    print(f"{'=' * 72}")

    found = False

    # ── 1. Composition data ──────────────────────────────────────────
    target_comps = {comp_id: comps[comp_id]} if comp_id and comp_id in comps else comps

    if prop_id and prop_id in props:
        print_section(f"Property Table: {prop_id}")
        prop_data = props[prop_id]
        # Print results table if it exists
        if "results" in prop_data:
            for entry in prop_data["results"]:
                eid = entry.get("id", "?")
                if comp_id and eid != comp_id:
                    continue
                parts = [f"  {eid}:"]
                for k, v in entry.items():
                    if k != "id":
                        parts.append(f"{k}={format_value(v)}")
                print("  ".join(parts))
                found = True
        # Print trends
        if "trends" in prop_data:
            print("\n  Trends:")
            for k, v in prop_data["trends"].items():
                print(f"    {k}: {v}")
                found = True

    if comp_id and comp_id in comps:
        print_section(f"Composition: {comps[comp_id].get('formula', comp_id)}")
        c = comps[comp_id]
        # Show targeted property if specified
        if prop_id and prop_id in c:
            hits = search_in_data({prop_id: c[prop_id]}, terms)
            for path, val, _ in hits[:15]:
                print(f"  {path} = {format_value(val)}")
                found = True
        elif prop_id:
            # Try alternate key names
            for key in c:
                if prop_id in key.lower():
                    hits = search_in_data({key: c[key]}, terms)
                    for path, val, _ in hits[:10]:
                        print(f"  {path} = {format_value(val)}")
                        found = True
        else:
            # General search
            hits = search_in_data(c, terms)
            for path, val, score in hits[:20]:
                print(f"  [{score}] {path} = {format_value(val)}")
                found = True

    # ── 2. Literature ────────────────────────────────────────────────
    if args.ref_id:
        if args.ref_id in refs:
            print_section(f"Reference: {args.ref_id}")
            r = refs[args.ref_id]
            print(f"  {r.get('authors', '?')} ({r.get('year', '?')})")
            print(f"  {r.get('title', '?')}")
            print(f"  {r.get('journal', '?')} {r.get('volume', '')} {r.get('pages', '')}")
            print(f"  Key: {r.get('key_content', '')}")
            found = True
    else:
        ref_hits = []
        for rid, r in refs.items():
            score = 0
            searchable = f"{rid} {r.get('title','')} {r.get('key_content','')} {' '.join(r.get('tags',[]))}"
            for t in terms:
                if t.lower() in searchable.lower():
                    score += 1
            if score > 0:
                ref_hits.append((rid, r, score))
        if ref_hits:
            ref_hits.sort(key=lambda x: -x[2])
            print_section("Literature Matches")
            for rid, r, score in ref_hits[:5]:
                print(f"  [{score}] {rid}: {r.get('authors','')} ({r.get('year','')}) — {r.get('title','')[:80]}")
                print(f"       Key: {r.get('key_content','')[:100]}")
                found = True

    # ── 3. Knowledge Base Articles ───────────────────────────────────
    kb_hits = search_kb_markdown(terms)
    if kb_hits:
        print_section("Knowledge Base Articles")
        for hit in kb_hits:
            print(f"\n  [{hit['n_matches']} matches] {hit['path']}")
            for ln, snippet in hit["snippets"]:
                for line in snippet.splitlines()[:3]:
                    print(f"    L{ln}: {line.strip()[:90]}")
            found = True

    # ── 4. Pipeline data ─────────────────────────────────────────────
    if topic == "pipeline" or any(t in query.lower() for t in ["pipeline", "annealing", "screening"]):
        print_section("Pipeline Information")
        for pid, p in pipes.items():
            hits = search_in_data(p, terms)
            if hits:
                print(f"\n  {pid}: {p.get('name', '')}")
                for path, val, _ in hits[:8]:
                    print(f"    {path} = {format_value(val)}")
                found = True

    if not found:
        print("\n  No results found.")
        print(f"  Available compositions: {', '.join(sorted(comps.keys()))}")
        print(f"  Available properties: {', '.join(sorted(props.keys()))}")
        print(f"  Try: python tools/kb_search.py --list-comps")

    print()


def cmd_list_comps(args):
    comps = load_all_compositions()
    print(f"\n{'ID':<10} {'Formula':<30} {'Family':<8} {'Atoms':<6} {'Cell':<10}")
    print("─" * 70)
    for cid, c in sorted(comps.items()):
        print(f"{cid:<10} {c.get('formula','?'):<30} {c.get('family','?'):<8} {c.get('atoms','?'):<6} {c.get('cell_type','?'):<10}")
    print()


def cmd_list_props(args):
    props = load_all_properties()
    print("\nAvailable property tables:")
    for pid, p in sorted(props.items()):
        desc = p.get("description", "")[:60]
        n = len(p.get("results", []))
        print(f"  {pid:<25} ({n} entries) — {desc}")
    print()


def cmd_ref(args):
    refs = load_literature()
    if args.ref_id:
        if args.ref_id in refs:
            r = refs[args.ref_id]
            print(f"\n  ID: {args.ref_id}")
            print(f"  Authors: {r.get('authors', '?')}")
            print(f"  Title: {r.get('title', '?')}")
            print(f"  Journal: {r.get('journal', '?')} {r.get('volume', '')} {r.get('pages', '')} ({r.get('year', '')})")
            print(f"  Key Content: {r.get('key_content', '')}")
            print(f"  Tags: {', '.join(r.get('tags', []))}")
            if "data" in r:
                print(f"  Data: {json.dumps(r['data'], indent=4)}")
        else:
            print(f"  Reference '{args.ref_id}' not found.")
            print(f"  Available: {', '.join(sorted(refs.keys()))}")
    else:
        print(f"\n{'ID':<20} {'Year':<6} {'Authors':<30} {'Tags'}")
        print("─" * 90)
        for rid, r in sorted(refs.items(), key=lambda x: x[1].get("year", 0)):
            tags = ", ".join(r.get("tags", [])[:3])
            print(f"{rid:<20} {r.get('year','?'):<6} {str(r.get('authors','?'))[:30]:<30} {tags}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Argyrodite Knowledge Base Search",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python tools/kb_search.py "comp5 bulk modulus"
  python tools/kb_search.py "Br effect on elastic"
  python tools/kb_search.py --comp comp1 --prop eos
  python tools/kb_search.py --ref pustorino2025
  python tools/kb_search.py --list-comps
  python tools/kb_search.py --list-props
  python tools/kb_search.py --list-refs
"""
    )
    parser.add_argument("terms", nargs="*", help="Search terms")
    parser.add_argument("--comp", help="Filter by composition ID")
    parser.add_argument("--prop", help="Filter by property (eos, elastic, electronic, bonds, adhesion)")
    parser.add_argument("--ref", dest="ref_id", help="Look up a specific reference")
    parser.add_argument("--list-comps", action="store_true", help="List all compositions")
    parser.add_argument("--list-props", action="store_true", help="List all property tables")
    parser.add_argument("--list-refs", action="store_true", help="List all references")

    args = parser.parse_args()

    if args.list_comps:
        cmd_list_comps(args)
    elif args.list_props:
        cmd_list_props(args)
    elif args.list_refs:
        cmd_ref(type("Args", (), {"ref_id": None})())
    elif args.ref_id:
        cmd_ref(args)
    elif args.terms or args.comp or args.prop:
        if not args.terms:
            args.terms = []
        cmd_search(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
