#!/usr/bin/env python3
"""Decide whether an EMMA annotation is a candidate for the tRNA rescue.

EMMA's tRNA model periodically misses one or a few tRNAs on a mitogenome that is
otherwise complete and correctly ordered -- the sequence is in the assembly, but
the feature never reaches the GFF/TBL. Those assemblies then fail
``annotation_stats.py`` (``passed=no`` when the shortfall exceeds
``annotation_trna_tolerance``) and are held out of QC/ENA.

This gate reads the EMMA GFF and emits a one-line decision the annotation
subworkflow branches on:

    FIX\\t<targets>   only tRNAs are missing, all 13 PCGs + both rRNAs are present,
                     and the present genes are already in canonical order
    PASS\\t-          anything else -- the rescue is not applicable

``<targets>`` is a comma-list of REF tRNA names (e.g. ``TP`` or ``TW,TA,TN``) in
canonical order. Only FIX assemblies go to TRNA_RESCUE; PASS assemblies flow
through untouched.

Any parse failure degrades to ``PASS\\t-`` (exit 0) so a malformed GFF can never
drop a sample from the run -- it just misses the rescue.

Usage:
    trna_rescue_gate.py --gff <emma.gff> --out <prefix>.trna_rescue_qc.txt
"""

import argparse
import sys
import traceback
from pathlib import Path

# The standard vertebrate gene order and its PCG / rRNA / tRNA partitions, shared
# with annotation_stats.py and the ND4L/ATP8 gate via bin/mito_gene_order.py so
# the gates and the QC step cannot drift on what "present and in order" means.
from mito_gene_order import REF_GENES, TRNA_GENES, RRNA_GENES, PCG_GENES


def parse_gff_attributes(attr_str):
    return dict(
        item.split("=", 1)
        for item in attr_str.strip().split(";")
        if "=" in item
    )


def genes_by_coord(gff_path):
    """Return REF gene names present in the GFF, ordered by genomic start.

    A gene written as more than one `gene` line (an origin-spanning feature is
    split into two) is kept once, at its LOWEST start. Taking the first line seen
    instead put such a gene at whichever half the file happened to list first,
    which could make the order check fail on an annotation that is in fact
    correctly ordered -- and a failed order check silently suppresses the rescue.
    """
    starts = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            attrs = parse_gff_attributes(parts[8])
            name = attrs.get("Name")
            if not name:
                continue
            gene = name.replace("MT-", "")
            start = int(parts[3])
            if gene not in starts or start < starts[gene]:
                starts[gene] = start
    return [g for g, _ in sorted(starts.items(), key=lambda kv: kv[1])]


def decide(gff_path):
    present = genes_by_coord(gff_path)
    present_set = set(present)
    missing = [g for g in REF_GENES if g not in present_set]

    if not missing:
        return "PASS", "-"
    if any(g not in TRNA_GENES for g in missing):
        return "PASS", "-"  # a PCG / rRNA is missing -> not our case, leave it held

    # The whole conserved core must be present; a tRNA rescue on a collapsed
    # assembly would just paper over a real defect.
    if not RRNA_GENES.issubset(present_set) or not PCG_GENES.issubset(present_set):
        return "PASS", "-"

    # The genes that ARE present must already be in REF order (equivalently:
    # order_correct would be "yes"). Inserting each missing tRNA at its true
    # coordinate then leaves the whole set ordered, and every missing tRNA's
    # REF neighbours are present so its insertion gap is well defined.
    ref_subset = [g for g in REF_GENES if g in present_set]
    if present != ref_subset:
        return "PASS", "-"

    targets = ",".join(g for g in REF_GENES if g in missing)
    return "FIX", targets


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    try:
        state, targets = decide(args.gff)
    except Exception:  # noqa: BLE001 - never drop a sample on a parse error
        traceback.print_exc(file=sys.stderr)
        state, targets = "PASS", "-"

    args.out.write_text(f"{state}\t{targets}\n")
    print(f"[trna_rescue_gate] {args.gff.name}: {state} {targets}", file=sys.stderr)


if __name__ == "__main__":
    main()
