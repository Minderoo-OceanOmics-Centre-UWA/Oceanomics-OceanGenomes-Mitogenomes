#!/usr/bin/env python3
"""
Generate assets/taxonomy/mito_origin_anchors.json from the curated reference DBs.

WHAT THIS MEASURES
------------------
For every record in assets/refdb/<group>/, which gene sits at position 1 -- i.e.
where the submitter chose to linearise the circle. Tallied per taxonomic order,
per class and per group, that is the deposition convention for each lineage, and
it is what the pipeline should re-origin a published mitogenome to.

WHY IT EXISTS
-------------
bin/mitos_to_emma.py used to re-origin EVERY invertebrate to tRNA-Met, with a
docstring claiming that matched how NCBI coral mitogenomes are deposited.
Measured against these same databases, trnM is the deposited origin for Porifera
0%, Annelida 0%, Ctenophora 0%, Echinodermata 0.7%, Mollusca 0.8% and Arthropoda
1.9%. The claim was true only for where it came from: Scleractinia, at 63.8%.

That last number is why the table is keyed on ORDER FIRST. The Anthozoa class
aggregate (rrnL 33.9%, cox1 28.5%, trnM 17.6%) has no majority and hides four
orders that each have a clear and different convention:

    Scleractinia     n=58  TM   63.8%     <- stony corals keep tRNA-Met
    Malacalcyonacea  n=75  RNR2 69.3%
    Zoantharia       n=29  CO1  58.6%
    Scleralcyonacea  n=31  CO1  51.6%

A class-level table would have rotated every already-submitted stony coral off
its existing origin, changing its ENA sequence checksum for no reason at all.

RESOLUTION LADDER
-----------------
A taxon's own plurality when it clears BOTH --min-fraction and --min-records,
else the next level up: order -> class -> group (phylum) -> --default-anchor.
Unlike the genetic-code table, an unmapped taxon never aborts; it takes the
default. A wrong anchor rotates a circle (cosmetic, reversible, and the adapter
falls back when the gene is not annotated), whereas no anchor is a hard failure.

Usage:
    build_origin_anchor_table.py --taxdump-dir /path/to/taxdump
    build_origin_anchor_table.py --taxdump-dir ... --check   # CI: is the asset stale?
"""

import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mito_gene_order import is_valid_anchor, refdb_gene_to_emma  # noqa: E402
from taxdump_lineage import TaxdumpLineage  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_REFDB = REPO_ROOT / "assets" / "refdb"
DEFAULT_OUT = REPO_ROOT / "assets" / "taxonomy" / "mito_origin_anchors.json"
GROOVY_SETS = REPO_ROOT / "lib" / "InvertTaxonGroups.groovy"

# Groovy set name -> assets/refdb/<group> directory, mirroring
# InvertTaxonGroups.seedDbGroup(). Parsed rather than duplicated so a class added
# on the Groovy side cannot silently go missing from the anchor table.
SET_TO_GROUP = {
    "CNIDARIA_CLASSES": "anthozoa",
    "PORIFERA_CLASSES": "porifera",
    "MOLLUSCA_CLASSES": "mollusca",
    "ARTHROPODA_CLASSES": "arthropoda",
    "ECHINODERMATA_CLASSES": "echinodermata",
    "CTENOPHORA_CLASSES": "ctenophora",
    "TUNICATA_CLASSES": "tunicata",
    "ANNELIDA_CLASSES": "annelida",
}


def parse_groovy_class_groups(path):
    """class name (lower) -> refdb group, read straight out of InvertTaxonGroups.groovy.

    Parsed rather than re-typed here: this map and seedDbGroup() must agree, and a
    second hand-maintained copy is exactly the drift the genetic-code asset was
    created to remove.
    """
    text = Path(path).read_text()
    out = {}
    for set_name, group in SET_TO_GROUP.items():
        m = re.search(
            rf"static\s+final\s+Set<String>\s+{set_name}\s*=\s*\[(.*?)\]\s*as\s+Set",
            text, re.S)
        if not m:
            sys.exit(f"[anchors] could not parse {set_name} from {path}")
        for name in re.findall(r"'([^']+)'", m.group(1)):
            out[name.strip().lower()] = group
    if not out:
        sys.exit(f"[anchors] parsed no classes from {path}")
    return out


def first_feature_per_accession(features_tsv, unmapped_counter):
    """accession -> EMMA key of the feature with the lowest start coordinate.

    `parts` is Biopython 0-based "start-end:strand", comma-joined for a multi-exon
    feature; the first element is the 5'-most piece in the record's own order.
    """
    lowest = {}
    with open(features_tsv) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            m = re.match(r"(\d+)-", (row.get("parts") or "").split(",")[0])
            if not m:
                continue
            start = int(m.group(1))
            acc = row["accession"]
            if acc not in lowest or start < lowest[acc][0]:
                key = refdb_gene_to_emma(row.get("type"), row.get("gene"), row.get("product"))
                lowest[acc] = (start, key, (row.get("type"), row.get("gene"), row.get("product")))
    out = {}
    for acc, (_start, key, triple) in lowest.items():
        if key is None:
            unmapped_counter[triple] += 1
            continue
        out[acc] = key
    return out


def decide(counter, min_fraction, min_records):
    """(anchor, n, fraction, top) for a tally, or (None, n, fraction, top) if it fails a bar."""
    n = sum(counter.values())
    if n == 0:
        return None, 0, 0.0, None
    top, hits = counter.most_common(1)[0]
    fraction = hits / n
    if n < min_records or fraction < min_fraction or not is_valid_anchor(top):
        return None, n, fraction, top
    return top, n, fraction, top


def row(anchor, level, n, fraction, top, counter, reason=None, group=None):
    entry = {"anchor": anchor, "level": level, "n": n, "fraction": round(fraction, 3)}
    if group:
        entry["group"] = group
    if top and top != anchor:
        entry["top"] = top
    if counter is not None:
        entry["counts"] = dict(counter.most_common(3))
    if reason:
        entry["reason"] = reason
    return entry


def why(n, fraction, top, min_records, min_fraction):
    bits = []
    if n < min_records:
        bits.append(f"n={n} below min_records={min_records}")
    elif fraction < min_fraction:
        bits.append(f"plurality {top} {fraction:.1%} below min_fraction={min_fraction:.0%}")
    if top and not is_valid_anchor(top):
        bits.append(f"'{top}' is not an addressable MITOS feature key "
                    f"(the reference DBs do not disambiguate tRNA-Leu/tRNA-Ser copies)")
    return "; ".join(bits) or "no majority"


def build(refdb_root, taxdump_dir, min_fraction, min_records, default_anchor, max_unmapped):
    class_to_group = parse_groovy_class_groups(GROOVY_SETS)
    known_classes = set(class_to_group)

    taxdump = None
    if taxdump_dir:
        taxdump = TaxdumpLineage(taxdump_dir)
        if not taxdump.available:
            sys.exit(f"[anchors] --taxdump-dir {taxdump_dir} has no names.dmp/nodes.dmp")

    by_order = collections.defaultdict(collections.Counter)
    by_class = collections.defaultdict(collections.Counter)
    by_group = collections.defaultdict(collections.Counter)
    order_group = {}
    unmapped = collections.Counter()
    total = 0
    order_cache = {}

    for group_dir in sorted(Path(refdb_root).iterdir()):
        group = group_dir.name
        features = group_dir / f"{group}_mito_refdb.features.tsv"
        manifest = group_dir / f"{group}_mito_refdb.manifest.tsv"
        if not (features.is_file() and manifest.is_file()):
            continue
        with open(manifest) as handle:
            meta = {r["accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
        firsts = first_feature_per_accession(features, unmapped)
        total += len(firsts) + sum(unmapped.values())
        for acc, key in firsts.items():
            by_group[group][key] += 1
            organism = meta.get(acc, {}).get("organism", "")
            resolved = {}
            if taxdump:
                if organism not in order_cache:
                    try:
                        got = taxdump.lineage_for_name(organism) or {}
                    except Exception:
                        got = {}
                    order_cache[organism] = ((got.get("order") or "").strip(),
                                             (got.get("class") or "").strip())
                resolved = {"order": order_cache[organism][0],
                            "class": order_cache[organism][1]}

            # Rank-resolved class from the taxdump wins. The lineage fallback takes
            # the LAST matching token, not the first: the Groovy class sets carry the
            # phylum name too (so a sample identified only to phylum still routes), and
            # an NCBI lineage runs root -> tip, so the first match is 'Mollusca' where
            # the class is 'Gastropoda'. Taking hit[0] silently collapsed all 850
            # molluscs into one 'mollusca' class row.
            taxon_class = resolved.get("class", "")
            if taxon_class.lower() not in known_classes:
                lineage = (meta.get(acc, {}).get("lineage") or "")
                hit = [t.strip() for t in lineage.split(";")
                       if t.strip().lower() in known_classes]
                taxon_class = hit[-1] if hit else ""
            if taxon_class:
                by_class[taxon_class.lower()][key] += 1

            taxon_order = resolved.get("order", "")
            if taxon_order:
                by_order[taxon_order.lower()][key] += 1
                order_group[taxon_order.lower()] = group

    n_unmapped = sum(unmapped.values())
    denom = total or 1
    if n_unmapped / denom > max_unmapped:
        print(f"[anchors] ERROR: {n_unmapped}/{denom} ({n_unmapped/denom:.1%}) first features "
              f"could not be mapped to an EMMA key, above --max-unmapped {max_unmapped:.1%}.",
              file=sys.stderr)
        for triple, count in unmapped.most_common(20):
            print(f"    {count:>5}  type={triple[0]!r} gene={triple[1]!r} product={triple[2]!r}",
                  file=sys.stderr)
        sys.exit(1)
    if n_unmapped:
        print(f"[anchors] {n_unmapped}/{denom} ({n_unmapped/denom:.2%}) first features unmapped "
              f"(under the {max_unmapped:.0%} bar):", file=sys.stderr)
        for triple, count in unmapped.most_common(10):
            print(f"    {count:>5}  type={triple[0]!r} gene={triple[1]!r} product={triple[2]!r}",
                  file=sys.stderr)

    groups_out = {}
    for group, counter in sorted(by_group.items()):
        anchor, n, fraction, top = decide(counter, min_fraction, min_records)
        if anchor:
            groups_out[group] = row(anchor, "group", n, fraction, top, counter)
        else:
            groups_out[group] = row(default_anchor, "default", n, fraction, top, counter,
                                    reason=why(n, fraction, top, min_records, min_fraction))

    classes_out = {}
    for taxon_class, counter in sorted(by_class.items()):
        group = class_to_group.get(taxon_class)
        anchor, n, fraction, top = decide(counter, min_fraction, min_records)
        if anchor:
            classes_out[taxon_class] = row(anchor, "class", n, fraction, top, None, group=group)
            continue
        reason = why(n, fraction, top, min_records, min_fraction)
        group_row = groups_out.get(group)
        if group_row and group_row["level"] == "group":
            classes_out[taxon_class] = row(
                group_row["anchor"], "group", n, fraction, top, None,
                reason=f"{reason}; inherits the {group} group anchor", group=group)
        else:
            classes_out[taxon_class] = row(default_anchor, "default", n, fraction, top, None,
                                           reason=reason, group=group)

    orders_out = {}
    for taxon_order, counter in sorted(by_order.items()):
        anchor, n, fraction, top = decide(counter, min_fraction, min_records)
        if not anchor:
            # Below a bar: say nothing at order level and let the class/group
            # levels answer, rather than pinning the order to an inherited value
            # that would go stale the moment its class tally changes.
            continue
        orders_out[taxon_order] = row(anchor, "order", n, fraction, top, counter,
                                      group=order_group.get(taxon_order))

    levels = (["order"] if taxdump else []) + ["class", "group", "default"]
    if not taxdump:
        print("[anchors] WARNING: no --taxdump-dir, so the ORDER level was skipped. The table "
              "will resolve Scleractinia from the Anthozoa class aggregate (no majority) and "
              "rotate stony corals off tRNA-Met. Do not ship this.", file=sys.stderr)

    return {
        "_comment": [
            "Taxon -> the gene the published invertebrate mitogenome is re-origined to,",
            "measured from this repo's own curated RefSeq databases under assets/refdb/.",
            "",
            "GENERATED. Do not hand-edit: bin/build_origin_anchor_table.py rewrites this",
            "file wholesale from each group's .features.tsv and .manifest.tsv. Regenerate",
            "after any database rebuild. This is the one difference from the sibling asset",
            "mito_genetic_codes.json, which IS curated by hand.",
            "",
            "Read by lib/InvertTaxonGroups.groovy (originAnchor(), which the annotation",
            "subworkflow passes to MITOS2 and CORAL_ANNOTATION_FIX as --origin-gene) and by",
            "bin/mito_gene_order.py (load_origin_anchors(), used by the generator's own",
            "validation and the unit tests) -- two readers for the same reason the genetic",
            "code table has two: a Groovy copy and a Python copy would drift.",
            "",
            "This REPLACES an unconditional re-origin to tRNA-Met whose docstring claimed",
            "trnM matched the NCBI coral convention. Measured here, trnM is the deposited",
            "origin for Porifera 0%, Annelida 0%, Ctenophora 0%, Echinodermata 0.7%,",
            "Mollusca 0.8%, Arthropoda 1.9%. It IS the convention for Scleractinia (63.8%),",
            "which is where the rule came from -- hence the order level below.",
            "",
            "Anchors are MITOS/EMMA feature keys as produced by mitos_to_emma.map_gene_name:",
            "CO1, CO3, RNR1, RNR2, TM, TF, ... They are looked up directly in that script's",
            "features dict, so the two vocabularies cannot drift.",
            "",
            "LADDER (see 'policy.levels'): a taxon's own plurality when it clears both",
            "min_fraction and min_records; else the next level up; else default_anchor.",
            "A taxon absent from every level resolves to the default -- unlike the genetic",
            "code table, an unmapped taxon never aborts, because a wrong rotation is a",
            "presentation choice while a missing one is a hard failure.",
        ],
        "policy": {
            "min_fraction": min_fraction,
            "min_records": min_records,
            "default_anchor": default_anchor,
            "levels": levels,
            "basis": (
                f"min_fraction {min_fraction} = a plurality that is also a majority, so the "
                f"anchor is what more than half of published records for the lineage actually "
                f"use. min_records {min_records} stops a tiny tally deciding: Remipedia is n=1 "
                f"at TM=100%, Monoplacophora n=2, Nuda n=2, Ostracoda n=4. Below either bar the "
                f"taxon inherits the next level up, which is what keeps Homoscleromorpha (n=13) "
                f"and Hexactinellida (n=4) on the poriferan rrnL instead of dropping to cox1. "
                f"tRNA-Leu and tRNA-Ser can win a tally but are rejected as anchors: the "
                f"reference DBs record them as bare 'tRNA-Leu'/'tRNA-Ser' while MITOS emits "
                f"TL1/TL2/TS1/TS2, so the winning key would not be addressable."
            ),
        },
        "orders": orders_out,
        "classes": classes_out,
        "groups": groups_out,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdb-root", type=Path, default=DEFAULT_REFDB)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--taxdump-dir", default=None,
                    help="NCBI taxdump directory (names.dmp + nodes.dmp). Required for the "
                         "ORDER level; without it Scleractinia cannot be resolved.")
    ap.add_argument("--min-fraction", type=float, default=0.50)
    ap.add_argument("--min-records", type=int, default=20)
    ap.add_argument("--default-anchor", default="CO1")
    ap.add_argument("--max-unmapped", type=float, default=0.02)
    ap.add_argument("--check", action="store_true",
                    help="Exit 1 if the file on disk differs from what would be generated.")
    args = ap.parse_args()

    table = build(args.refdb_root, args.taxdump_dir, args.min_fraction,
                  args.min_records, args.default_anchor, args.max_unmapped)
    rendered = json.dumps(table, indent=2, sort_keys=False) + "\n"

    if args.check:
        if not args.out.exists():
            print(f"[anchors] {args.out} does not exist", file=sys.stderr)
            return 1
        if args.out.read_text() != rendered:
            print(f"[anchors] {args.out} is STALE -- rerun build_origin_anchor_table.py",
                  file=sys.stderr)
            return 1
        print(f"[anchors] {args.out} is up to date")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(rendered)
    print(f"[anchors] wrote {args.out}: {len(table['orders'])} orders, "
          f"{len(table['classes'])} classes, {len(table['groups'])} groups")
    return 0


if __name__ == "__main__":
    sys.exit(main())
