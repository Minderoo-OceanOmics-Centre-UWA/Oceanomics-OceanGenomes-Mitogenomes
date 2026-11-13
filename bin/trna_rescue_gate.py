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
# the gates and the QC step cannot drift on what "present and in order" means --
# which covers HOW the GFF is read as well as what the reference order is.
from mito_gene_order import (
    REF_GENES,
    TRNA_GENES,
    RRNA_GENES,
    PCG_GENES,
    add_taxon_arguments,
    genes_by_coord,
    matching_order_for,
    taxon_from_args,
)


def decide(gff_path, taxon=None):
    present = genes_by_coord(gff_path)
    present_set = set(present)

    # See emma_rescue_gate.decide: the accepted order for this taxon, canonical
    # unless a curated variant rule applies.
    order_ref, _rule = matching_order_for(present, taxon)

    missing = [g for g in REF_GENES if g not in present_set]

    if not missing:
        return "PASS", "-"
    if any(g not in TRNA_GENES for g in missing):
        return "PASS", "-"  # a PCG / rRNA is missing -> not our case, leave it held

    # The whole conserved core must be present; a tRNA rescue on a collapsed
    # assembly would just paper over a real defect.
    if not RRNA_GENES.issubset(present_set) or not PCG_GENES.issubset(present_set):
        return "PASS", "-"

    # The genes that ARE present must already be in an accepted order (equivalently:
    # order_correct would be "yes" or "variant"). Inserting each missing tRNA at its
    # true coordinate then leaves the whole set ordered, and every missing tRNA's
    # neighbours in that order are present so its insertion gap is well defined.
    if order_ref is None:
        return "PASS", "-"

    targets = ",".join(g for g in order_ref if g in missing)
    return "FIX", targets


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    # FUTURE: see emma_rescue_gate.main -- no LCA exists yet at this point in the
    # pipeline, so this gate matches on meta taxonomy only.
    add_taxon_arguments(ap)
    args = ap.parse_args()

    try:
        state, targets = decide(args.gff, taxon_from_args(args))
    except Exception:  # noqa: BLE001 - never drop a sample on a parse error
        traceback.print_exc(file=sys.stderr)
        state, targets = "PASS", "-"

    args.out.write_text(f"{state}\t{targets}\n")
    print(f"[trna_rescue_gate] {args.gff.name}: {state} {targets}", file=sys.stderr)


if __name__ == "__main__":
    main()
