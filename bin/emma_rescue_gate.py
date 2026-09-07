#!/usr/bin/env python3
"""Decide whether an EMMA annotation is a candidate for the ND4L / ATP8 rescue.

EMMA's ``rationalise_matches!`` step discards a short protein-coding gene when its
computed circular overlap with a longer neighbour exceeds half the shorter
feature's length. In practice this drops ND4L (vs ND4) and ATP8 (vs ATP6) from
mitogenomes that are otherwise complete and correctly ordered -- the gene is
present in the assembly, EMMA even reports the match, but it never reaches the
GFF/TBL. Those assemblies then fail ``annotation_stats.py`` (``passed=no``) and
are held out of QC/ENA.

This gate reads the EMMA GFF and emits a one-line decision the annotation
subworkflow branches on:

    FIX\\t<targets>   only ND4L and/or ATP8 are missing, both flanks of each are
                     present, and every other REF gene is present and in order
    PASS\\t-          anything else -- the rescue is not applicable

``<targets>`` is a comma-list drawn from {ND4L, ATP8}. Only FIX assemblies go to
EMMA_GENE_RESCUE; PASS assemblies flow through untouched.

Any parse failure degrades to ``PASS\\t-`` (exit 0) so a malformed GFF can never
drop a sample from the run -- it just misses the rescue.

Usage:
    emma_rescue_gate.py --gff <emma.gff> --out <prefix>.emma_rescue_qc.txt
"""

import argparse
import sys
import traceback
from pathlib import Path

# The standard vertebrate gene order, shared with annotation_stats.py and the
# tRNA gate via bin/mito_gene_order.py so the gates and the QC step cannot drift
# on what "present and in order" means -- which covers HOW the GFF is read as
# well as what the reference order is.
from mito_gene_order import REF_GENES, genes_by_coord

# Genes this rescue can recover, and the REF-order neighbours each one needs
# present for the flanking-coordinate window to be well defined.
RESCUABLE = {
    "ND4L": ("TR", "ND4"),
    "ATP8": ("TK", "ATP6"),
}


def decide(gff_path):
    present = genes_by_coord(gff_path)
    present_set = set(present)
    missing = [g for g in REF_GENES if g not in present_set]

    if not missing:
        return "PASS", "-"
    if any(g not in RESCUABLE for g in missing):
        return "PASS", "-"  # a non-rescuable gene is missing -> not our case

    # Every missing gene must have both REF-order neighbours present, else the
    # flanking window is undefined.
    for g in missing:
        left, right = RESCUABLE[g]
        if left not in present_set or right not in present_set:
            return "PASS", "-"

    # The genes that ARE present must already be in REF order (equivalently:
    # order_correct would be "yes"). Inserting the missing gene at its true
    # coordinate then leaves the whole set ordered.
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
    print(f"[emma_rescue_gate] {args.gff.name}: {state} {targets}", file=sys.stderr)


if __name__ == "__main__":
    main()
