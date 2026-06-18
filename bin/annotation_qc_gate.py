#!/usr/bin/env python3
"""
Decide whether a MITOS2 anthozoan annotation is deficient enough to route through
the coral fixer (coral_fix_bed.py).

The fixer is cheap but not free, and most anthozoans in a batch are annotated
correctly by MITOS2; we only want to repair the ones that are actually broken.
This gate inspects the EMMA-contract annotation MITOS2 produced and emits a
one-line decision:

    FIX\t<reason>     # deficient -> send to coral fixer
    PASS\t<reason>    # complete  -> keep MITOS2 output unchanged

A coral annotation is deficient when either:
  * a conserved core gene is missing from the GFF -- the cnidarian core is the
    13 PCGs + RNR1 (12S) + RNR2 (16S); MITOS routinely drops RNR2 on divergent
    corals; or
  * nad5 is present but truncated -- the group-I-intron split means MITOS often
    annotates only one exon, leaving the translated protein far short of the
    ~600 aa full length (the intact spliced ND5).

Reads the GFF for gene presence and the proteins/ dir for the ND5 length, so it
works off exactly the files mitos_to_emma.py wrote. Always exits 0.

Usage:
    annotation_qc_gate.py --gff annotation/PREFIX.gff --proteins annotation/proteins \
        --out PREFIX.coral_qc.txt [--min-nd5-aa 500]
"""

import argparse
import re
from pathlib import Path

CNIDARIAN_CORE = [
    "ATP6", "ATP8", "CO1", "CO2", "CO3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "RNR1", "RNR2",
]


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
    for line in open(hits[0]):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).rstrip("*")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-nd5-aa", type=int, default=540,
                    help="ND5 translation shorter than this is treated as truncated "
                         "(intact anthozoan ND5 is ~580-605 aa). MITOS often calls a "
                         "partial first exon that still clears ~540 aa, so the M-start "
                         "check below is the more reliable truncation signal.")
    args = ap.parse_args()

    reasons = []
    genes = genes_in_gff(args.gff)
    missing = [g for g in CNIDARIAN_CORE if g not in genes]
    if missing:
        reasons.append("missing=" + ",".join(missing))

    nd5 = protein_seq(args.proteins, "ND5")
    if nd5 is not None:
        nd5_aa = len(nd5)
        if nd5_aa < args.min_nd5_aa:
            reasons.append(f"ND5_trunc={nd5_aa}aa<{args.min_nd5_aa}")
        # A complete anthozoan ND5 starts on Met; when MITOS truncates the
        # intron-split first exon the translation loses that start codon while
        # often staying just above the length floor (e.g. 540 aa, no M). Treat a
        # non-Met start as a deficient first exon so the fixer rebuilds it.
        elif not nd5.startswith("M"):
            reasons.append(f"ND5_no_start_M(first={nd5[:1] or '-'},{nd5_aa}aa)")

    decision = "FIX" if reasons else "PASS"
    reason = ";".join(reasons) if reasons else "core complete; ND5 full-length"
    Path(args.out).write_text(f"{decision}\t{reason}\n")
    print(f"[annotation_qc_gate] {decision}: {reason}")


if __name__ == "__main__":
    main()
