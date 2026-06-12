#!/usr/bin/env python3
"""
Adapt MITOS2 annotation output to the EMMA output contract used by the rest of
the OceanGenomes mitogenome pipeline.

MITOS2 emits result.bed / result.gff / result.faa / result.fas with its own gene
nomenclature (cox1, cob, rrnS, rrnL, trnF(gaa), ...). Downstream modules
(annotation_stats.py, extract_genes_gff.py, process_files.py, the BLAST->LCA mix)
were written against EMMA, which uses MT-<GENE> naming (CO1, CYTB, RNR1, RNR2,
TF, ...), an EMMA-style GFF (##sequence-region + gene/CDS/tRNA/rRNA lines with
Name= attributes) and per-gene FASTAs under cds/ and proteins/.

This script reads MITOS2's result.bed (stable, well-defined columns) plus the
assembly FASTA and writes an EMMA-shaped directory:

    <outdir>/<prefix>.fa                      (the assembly genome)
    <outdir>/<prefix>.gff                     (EMMA-style GFF)
    <outdir>/cds/<GENE>.<prefix>.fa           (per-gene nucleotide FASTA)
    <outdir>/proteins/MT-<GENE>.<prefix>.fa   (per-PCG translated FASTA)

Coordinates and sequences are derived from the BED so the result is independent
of MITOS2's GFF/FASTA header formatting, which varies between releases.

Usage:
    mitos_to_emma.py --bed result.bed --genome asm.fa \
        --prefix OG1.ilmn.250131.getorg1770.mitos2110 --outdir emma --code 5
"""

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# MITOS2 protein-coding gene name -> EMMA gene name
PCG_MAP = {
    "cox1": "CO1", "cox2": "CO2", "cox3": "CO3",
    "cob": "CYTB",
    "atp6": "ATP6", "atp8": "ATP8",
    "nad1": "ND1", "nad2": "ND2", "nad3": "ND3",
    "nad4": "ND4", "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6",
}

# MITOS2 rRNA name -> EMMA gene name
RRNA_MAP = {"rrnS": "RNR1", "rrnL": "RNR2"}

# Products for nicer GenBank output (downstream falls back to gene name otherwise)
PRODUCT = {
    "CO1": "cytochrome c oxidase subunit I",
    "CO2": "cytochrome c oxidase subunit II",
    "CO3": "cytochrome c oxidase subunit III",
    "CYTB": "cytochrome b",
    "ATP6": "ATP synthase F0 subunit 6",
    "ATP8": "ATP synthase F0 subunit 8",
    "ND1": "NADH dehydrogenase subunit 1",
    "ND2": "NADH dehydrogenase subunit 2",
    "ND3": "NADH dehydrogenase subunit 3",
    "ND4": "NADH dehydrogenase subunit 4",
    "ND4L": "NADH dehydrogenase subunit 4L",
    "ND5": "NADH dehydrogenase subunit 5",
    "ND6": "NADH dehydrogenase subunit 6",
    "RNR1": "12S ribosomal RNA",
    "RNR2": "16S ribosomal RNA",
}


def map_gene_name(raw):
    """Return (emma_name, feature_type) for a MITOS2 BED feature name.

    feature_type is one of 'CDS', 'rRNA', 'tRNA'. Returns (None, None) for names
    we don't recognise (kept out of the EMMA contract rather than guessed).
    """
    # Strip MITOS fragment suffixes (e.g. cox1_0, nad5_1) and anticodon parens.
    name = re.sub(r"_\d+$", "", raw.strip())
    bare = re.sub(r"\(.*?\)$", "", name)

    if bare in PCG_MAP:
        return PCG_MAP[bare], "CDS"
    if bare in RRNA_MAP:
        return RRNA_MAP[bare], "rRNA"
    if bare.startswith("trn"):
        # trnF -> TF, trnL2 -> TL2, trnS1 -> TS1
        return "T" + bare[3:].upper(), "tRNA"
    return None, None


def parse_bed(bed_path):
    """Parse MITOS2 result.bed, merging fragmented features by gene name.

    Returns dict keyed by emma_name -> (chrom, start1, end, strand, ftype),
    with 1-based inclusive start.
    """
    features = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                cols = line.split()
            if len(cols) < 6:
                continue
            chrom, start, end, raw_name, _score, strand = cols[:6]
            emma_name, ftype = map_gene_name(raw_name)
            if emma_name is None:
                continue
            start1, end = int(start) + 1, int(end)
            if emma_name in features:
                c, s0, e0, st, ft = features[emma_name]
                features[emma_name] = (c, min(s0, start1), max(e0, end), st, ft)
            else:
                features[emma_name] = (chrom, start1, end, strand, ftype)
    return features


def write_gff(gff_path, features, chrom, seq_len, species):
    ordered = sorted(features.items(), key=lambda kv: kv[1][1])
    with open(gff_path, "w") as out:
        out.write("##gff-version 3\n")
        out.write(f"##sequence-region {chrom} 1 {seq_len}\n")
        if species:
            out.write(f"##organism {species}\n")
        for emma_name, (_c, start, end, strand, ftype) in ordered:
            product = PRODUCT.get(emma_name, emma_name)
            attrs_gene = f"ID=gene-{emma_name};Name=MT-{emma_name};Product={product}"
            attrs_feat = f"ID={emma_name};Parent=gene-{emma_name};Name=MT-{emma_name};Product={product}"
            out.write(f"{chrom}\tmitos\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs_gene}\n")
            out.write(f"{chrom}\tmitos\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t{attrs_feat}\n")


def gene_nt_seq(genome, chrom, start, end, strand):
    seq = genome[chrom].seq[start - 1:end]
    return seq.reverse_complement() if strand == "-" else seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, type=Path, help="MITOS2 result.bed")
    ap.add_argument("--genome", required=True, type=Path, help="Assembly FASTA")
    ap.add_argument("--prefix", required=True, help="Output prefix (id.tech.date.assembler.mitos<ver>)")
    ap.add_argument("--outdir", required=True, type=Path, help="EMMA-style output dir")
    ap.add_argument("--code", type=int, default=5, help="Mito genetic code (default 5, invertebrate)")
    ap.add_argument("--species", default="", help="Species name for ##organism")
    args = ap.parse_args()

    if not args.bed.exists():
        sys.exit(f"BED not found: {args.bed}")
    if not args.genome.exists():
        sys.exit(f"Genome not found: {args.genome}")

    genome = SeqIO.to_dict(SeqIO.parse(str(args.genome), "fasta"))
    if not genome:
        sys.exit(f"No sequences read from {args.genome}")

    features = parse_bed(args.bed)
    if not features:
        print(f"WARNING: no recognised features parsed from {args.bed}", file=sys.stderr)

    # Anchor the GFF header on the contig MITOS annotated (fall back to first).
    chrom = next(iter(features.values()))[0] if features else next(iter(genome))
    if chrom not in genome:
        # MITOS may rename the contig; fall back to the single assembly record.
        chrom = next(iter(genome))
        features = {n: (chrom, s, e, st, ft) for n, (_c, s, e, st, ft) in features.items()}
    seq_len = len(genome[chrom].seq)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cds_dir = args.outdir / "cds"
    prot_dir = args.outdir / "proteins"
    cds_dir.mkdir(exist_ok=True)
    prot_dir.mkdir(exist_ok=True)

    # Genome FASTA under the EMMA prefix
    SeqIO.write(list(genome.values()), str(args.outdir / f"{args.prefix}.fa"), "fasta")

    write_gff(args.outdir / f"{args.prefix}.gff", features, chrom, seq_len, args.species)

    for emma_name, (c, start, end, strand, ftype) in features.items():
        nt = gene_nt_seq(genome, c, start, end, strand)
        header = f">{emma_name}.{args.prefix}"
        with open(cds_dir / f"{emma_name}.{args.prefix}.fa", "w") as out:
            out.write(f"{header}\n{nt}\n")
        if ftype == "CDS":
            aa = Seq(str(nt)).translate(table=args.code)
            with open(prot_dir / f"MT-{emma_name}.{args.prefix}.fa", "w") as out:
                out.write(f">MT-{emma_name}.{args.prefix}\n{aa}\n")

    print(f"[mitos_to_emma] wrote {len(features)} features to {args.outdir}")


if __name__ == "__main__":
    main()
