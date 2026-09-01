#!/usr/bin/env python3
"""
Decide whether a MITOS2 invertebrate annotation is deficient enough to route
through the fixer (coral_fix_bed.py, on the FIX branch).

The fixer is cheap but not free, and most invertebrates in a batch are annotated
correctly by MITOS2; we only want to repair the ones that are actually broken.
This gate inspects the EMMA-contract annotation MITOS2 produced and emits a
one-line decision:

    FIX\t<reason>     # deficient -> send to the fixer
    PASS\t<reason>    # complete  -> keep MITOS2 output unchanged

An annotation is deficient when any of:

  * (code-generic) a protein-coding gene translates without a valid start codon,
    without a terminal stop, or with an internal stop -- read from the spliced
    per-gene CDS nucleotides MITOS wrote to annotation/cds/<gene>.<prefix>.fa and
    checked against the sample's own NCBI translation table. This is what NCBI
    table2asn enforces (SEQ_FEAT.StartCodon / SEQ_FEAT.NoStop / internal stop),
    and the check the old gate lacked -- it let mis-placed boundaries (e.g. a
    MITOS ND1 off by 3 codons) go straight to table2asn and fail terminally.

  * (cnidarian only, code 4) a conserved core gene is missing from the GFF -- the
    cnidarian core is the 13 PCGs + RNR1 (12S) + RNR2 (16S); MITOS routinely
    drops RNR2 on divergent corals.

  * (cnidarian only, code 4) nad5 is present but truncated -- the group-I-intron
    split means MITOS often annotates only one exon, leaving the translated
    protein far short of the ~600 aa full length.

  * (cnidarian only, code 4) cox1 is present but truncated, for the same reason:
    some scleractinians carry a second group I intron in cox1, and MITOS then
    reports only one of its two exons. That leaves a ~290 aa CO1. Such a CDS
    usually trips the no-stop check too, but only by luck of where the exon ends,
    so the length is checked explicitly instead. A short CO1 deserves its own signal
    regardless: CO1 is the barcode the LCA calls species from, and half a CO1
    still BLASTs to plausible-looking neighbours.

The core-gene-presence, ND5 and CO1 heuristics are Anthozoa-specific and only run
for cnidarian (code 4) samples; the per-PCG ORF check runs for every invertebrate
lineage. Always exits 0; any parse error degrades to PASS.

Usage:
    annotation_qc_gate.py --gff annotation/PREFIX.gff --proteins annotation/proteins \
        --cds annotation/cds --genetic-code 4 --out PREFIX.coral_qc.txt \\
        [--min-nd5-aa 540] [--min-co1-aa 450]
"""

import argparse
import re
import sys
import traceback
from pathlib import Path

from orf_utils import classify_cds

# EMMA-contract names of the 13 protein-coding genes.
PCGS = [
    "ATP6", "ATP8", "CO1", "CO2", "CO3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
]

CNIDARIAN_CORE = PCGS + ["RNR1", "RNR2"]

# The cnidarian-specific heuristics (core-gene presence, ND5 intron truncation)
# only make sense for the Coelenterate translation table.
CNIDARIAN_CODE = 4


def genes_in_gff(gff_path):
    genes = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            m = re.search(r"Name=MT-([^;]+)", cols[8])
            if m:
                genes.add(m.group(1))
    return genes


def protein_seq(prot_dir, gene):
    """Return the translated protein string for a gene, or None if absent."""
    hits = sorted(Path(prot_dir).glob(f"MT-{gene}.*.fa*"))
    if not hits:
        return None
    seq = []
    with open(hits[0]) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).rstrip("*")


def cds_seq(cds_dir, gene):
    """Return the spliced CDS nucleotide string for a gene, or None if absent."""
    if not cds_dir:
        return None
    hits = sorted(Path(cds_dir).glob(f"{gene}.*.fa*"))
    if not hits:
        return None
    seq = []
    with open(hits[0]) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def evaluate(gff, proteins, cds, genetic_code, min_nd5_aa, min_co1_aa,
             is_cnidarian):
    """Return (decision, reason). decision is 'FIX' or 'PASS'."""
    reasons = []
    genes = genes_in_gff(gff)

    if is_cnidarian:
        missing = [g for g in CNIDARIAN_CORE if g not in genes]
        if missing:
            reasons.append("missing=" + ",".join(missing))

        nd5 = protein_seq(proteins, "ND5")
        if nd5 is not None:
            nd5_aa = len(nd5)
            if nd5_aa < min_nd5_aa:
                reasons.append(f"ND5_trunc={nd5_aa}aa<{min_nd5_aa}")
            elif not nd5.startswith("M"):
                reasons.append(
                    f"ND5_no_start_M(first={nd5[:1] or '-'},{nd5_aa}aa)")

        co1 = protein_seq(proteins, "CO1")
        if co1 is not None:
            co1_aa = len(co1.rstrip("*"))
            if co1_aa < min_co1_aa:
                reasons.append(f"CO1_trunc={co1_aa}aa<{min_co1_aa}")

    # Code-generic per-PCG ORF check: every invertebrate lineage. A gene absent
    # from the GFF but with no CDS file is left to the (cnidarian) missing-core
    # check above; a non-cnidarian invert with a genuinely absent PCG still gets
    # routed to FIX via the missing CDS only if a CDS file exists but is broken.
    for gene in PCGS:
        nt = cds_seq(cds, gene)
        if not nt:
            continue
        info = classify_cds(nt, genetic_code)
        if not info["start_ok"]:
            reasons.append(f"{gene}_no_start")
        if not info["stop_ok"]:
            reasons.append(f"{gene}_no_stop")
        if info["internal_stops"]:
            reasons.append(f"{gene}_internal_stop={info['internal_stops']}")

    decision = "FIX" if reasons else "PASS"
    reason = ";".join(reasons) if reasons else "core complete; PCGs clean ORFs"
    return decision, reason


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--cds", required=True,
                    help="Dir of spliced per-gene CDS FASTAs (annotation/cds)")
    ap.add_argument("--genetic-code", type=int, required=True,
                    help="NCBI mitochondrial translation table for this sample "
                         "(from meta.genetic_code); no default -- the caller "
                         "always resolves it.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-nd5-aa", type=int, default=540,
                    help="ND5 translation shorter than this is treated as "
                         "truncated (cnidarian only).")
    ap.add_argument("--min-co1-aa", type=int, default=450,
                    help="CO1 translation shorter than this is treated as "
                         "truncated (cnidarian only). Anthozoan cox1 is ~515 aa; "
                         "one exon of an intron-split cox1 is ~290.")
    args = ap.parse_args()

    is_cnidarian = int(args.genetic_code) == CNIDARIAN_CODE
    try:
        decision, reason = evaluate(
            args.gff, args.proteins, args.cds, int(args.genetic_code),
            args.min_nd5_aa, args.min_co1_aa, is_cnidarian)
    except Exception:  # noqa: BLE001 - never drop a sample on a parse error
        traceback.print_exc(file=sys.stderr)
        decision, reason = "PASS", "gate parse error -> PASS"

    Path(args.out).write_text(f"{decision}\t{reason}\n")
    print(f"[annotation_qc_gate] {decision}: {reason}")


if __name__ == "__main__":
    main()
