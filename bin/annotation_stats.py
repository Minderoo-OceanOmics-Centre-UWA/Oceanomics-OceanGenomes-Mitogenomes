#!/usr/bin/env python3
"""
Using the .gff, counts the length of the mitogenome and all of the genes.
Using the .gff, reports on if there are missing or extra genes and if theyre in the correct order.
Counts the length of each mitochondrial protein-coding gene
from translated FASTA (.faa/.fasta) files in a `proteins/` directory.

Usage:
    singularity run $SING/psycopg2:0.1.sif python emma_stats.py /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/emma/*.gff /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/emma/

"""

import os
import csv
import sys
from pathlib import Path
from Bio import SeqIO


# Reference gene order (tRNA, rRNA, CDS)
REF_GENES = [
    "TF", "RNR1", "TV", "RNR2", "TL2", "ND1", "TI", "TQ",
    "TM", "ND2", "TW", "TA", "TN", "TC", "TY", "CO1", "TS2",
    "TD", "CO2", "TK", "ATP8", "ATP6", "CO3", "TG", "ND3",
    "TR", "ND4L", "ND4", "TH", "TS1", "TL1", "ND5", "ND6",
    "TE", "CYTB", "TT", "TP"
]

# Protein-coding genes to pull from .faa/.fa files
PROT_GENES = [
    "ATP6", "ATP8", "CO1", "CO2", "CO3",
    "CYTB", "ND1", "ND2", "ND3", "ND4",
    "ND4L", "ND5", "ND6"
]

def parse_gff_attributes(attr_str):
    return dict(
        item.split("=", 1)
        for item in attr_str.strip().split(";")
        if "=" in item
    )

def extract_total_length(gff_path):
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) == 4:
                    return abs(int(parts[3]) - int(parts[2])) + 1
    return None

def get_annotation_name(gff_path):
    """Extracts annotation name from the GFF file basename (no extension)."""
    return Path(gff_path).stem

def process_gff(gff_path, annotation_name):
    parts = annotation_name.split(".")
    if len(parts) != 5:
        print(f"⚠️ Warning: Unexpected annotation_name format: {annotation_name}")
        og_id = tech = seq_date = code = annotation = ""
    else:
        og_id, tech, seq_date, code, annotation = parts

    gene_entries = []  # list of (gene_name, start, end, strand)
    seen = set()
    gene_lengths = {gene: "" for gene in REF_GENES}
    total_length = extract_total_length(gff_path)

    num_cds = 0
    num_trna = 0
    num_rrna = 0

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            feature_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            attributes = parse_gff_attributes(parts[8])

            if feature_type == "cds":
                num_cds += 1
            elif feature_type == "trna":
                num_trna += 1
            elif feature_type == "rrna":
                num_rrna += 1

            if parts[2] == "gene" and "Name" in attributes:
                gene_name = attributes["Name"].replace("MT-", "")
                if gene_name not in seen:
                    seen.add(gene_name)
                    gene_entries.append((gene_name, start, end, parts[6]))  # add strand too
                    if gene_name in gene_lengths:
                        gene_lengths[gene_name] = str(abs(end - start) + 1)

    gene_entries.sort(key=lambda x: x[1])  # sort by start
    found_by_coord = [g[0] for g in gene_entries]

    missing = [g for g in REF_GENES if g not in found_by_coord]
    extra = [g for g in found_by_coord if g not in REF_GENES]
    
    ref_subset = [g for g in REF_GENES if g in found_by_coord]
    order_ok = (found_by_coord == ref_subset)
    
    passed = len(missing) == 0 and order_ok

    gff_summary = {
        "og_id": og_id,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "annotation": annotation,
        "missing_genes": ";".join(missing) if missing else "no",
        "extra_genes": ";".join(extra) if extra else "no",
        "order_correct": "yes" if order_ok else "no",
        "passed": "yes" if passed else "no",
        "total_length": total_length if total_length is not None else "NA",
        "num_cds": num_cds,
        "num_trna": num_trna,
        "num_rrna": num_rrna
    }
    gff_summary.update(gene_lengths)
    return gff_summary

def process_protein_lengths(prot_dir, annotation_name):
    missing = []
    prot_lengths = {f"{gene}_trans": "" for gene in PROT_GENES}
    for gene in PROT_GENES:
        pattern = f"MT-{gene}.{annotation_name}.fa*"
        candidates = sorted(prot_dir.glob(pattern))
        if not candidates:
            prot_lengths[f"{gene}_trans"] = None
            missing.append(gene)
        else:
            faa_path = candidates[0]
            rec = next(SeqIO.parse(faa_path, "fasta"))
            prot_lengths[f"{gene}_trans"] = len(rec.seq)
    if missing:
        print(f"⚠️  Missing translated genes in {annotation_name}: {', '.join(missing)}")
    else:
        print(f"✅ All translated genes present in {annotation_name}")
    return prot_lengths

def main(gff_path):
    if not os.path.isfile(gff_path):
        sys.exit(f"❌ GFF file not found: {gff_path}")

    annotation_name = get_annotation_name(gff_path)
    gff_summary = process_gff(gff_path, annotation_name)
    prot_lengths = process_protein_lengths(prot_dir, annotation_name)

    combined = {**gff_summary, **prot_lengths}

    og_id = get_annotation_name(gff_path).split(".")[0]
    output_path = f"{og_id}.annotation_stats.csv"
    with open(output_path, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(combined.keys()))
        writer.writeheader()
        writer.writerow(combined)

    print(f"✅ Combined summary written to:\n  {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python annotation_stats.py path/to/file.gff path/to/proteinsdir")
        sys.exit(1)
    # config_file = sys.argv[1]
    gff_path = sys.argv[1]
    prot_dir = Path(sys.argv[2])

    main(gff_path)
