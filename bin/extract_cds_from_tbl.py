#!/usr/bin/env python3
# extract_cds_from_tbl.py
import argparse
import re
import sys
from pathlib import Path
from Bio import SeqIO

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract CDS FASTAs from a mitochondrial genome using a GenBank .tbl feature table (single sample)."
    )
    p.add_argument("--fasta",    required=True, type=Path, help="Genome FASTA")
    p.add_argument("--tbl",      required=True, type=Path, help="Matching .tbl feature table")
    p.add_argument("--outdir",   required=True, type=Path, help="Output directory (CDSs under outdir/cds)")
    p.add_argument("--assembly", default=None,             help="Assembly/sample ID; default=FASTA stem")
    return p.parse_args()

def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    args = parse_args()

    if not args.fasta.exists():
        sys.exit(f"FASTA not found: {args.fasta}")
    if not args.tbl.exists():
        sys.exit(f"TBL not found: {args.tbl}")

    assembly = args.assembly or args.fasta.stem
    out_cds  = args.outdir / "cds"
    out_cds.mkdir(parents=True, exist_ok=True)

    # Load genome & organism
    genome_record = SeqIO.read(str(args.fasta), "fasta")

    m = re.search(r"\[organism=([^\]]+)\]", genome_record.description or "")
    if not m:
        raise ValueError("Species not found in FASTA header. Expected [organism=...] tag.")
    species = m.group(1)

    # Read tbl
    lines = args.tbl.read_text().splitlines()

    # State across lines
    gene_name = "unknown_gene"

    # Iterate lines; emulate your original simple TBL layout
    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")

        # Capture gene name to apply to the next CDS block
        if len(parts) == 2 and parts[0] == "gene":
            gene_name = parts[1]
            continue

        # CDS rows look like: start \t end \t CDS
        if len(parts) >= 3 and parts[2] == "CDS":
            try:
                start, end = int(parts[0]), int(parts[1])
            except ValueError:
                # If coordinates are not plain ints (e.g. joins), fall back with a helpful error
                raise ValueError(
                    f"Unsupported TBL coordinate line at row {i+1}: '{line}'. "
                    "This light parser expects simple 'start\\tend\\tCDS' rows."
                )

            # Look ahead a few lines for qualifiers
            product_name = "unknown_product"
            gene_id      = "unknown_protein_id"
            for j in range(i + 1, min(i + 6, len(lines))):
                nxt = lines[j].strip()
                if "\tproduct\t" in nxt:
                    product_name = nxt.split("\t")[-1].strip()
                elif "\tprotein_id\t" in nxt:
                    raw_gene_id = nxt.split("\t")[-1].strip()   # e.g. gnl|Emma|7ccb…
                    gene_id = raw_gene_id.split("|")[-1]        # keep the last token
                # break early if we have both
                if product_name != "unknown_product" and gene_id != "unknown_protein_id":
                    break

            # Extract sequence with strand logic
            if start <= end:  # '+' strand
                cds_seq   = genome_record.seq[start - 1:end]
                coord_str = f"{start}-{end}"
                direction = "+"
            else:             # '-' strand
                cds_seq   = genome_record.seq[end - 1:start].reverse_complement()
                coord_str = f"complement({end}-{start})"
                direction = "-"

            # Build header consistent with your format
            header = (
                f">{assembly}|{coord_str}|{direction}|MT-{gene_name} "
                f"[organism={species}] [mgcode=2] "
                f"[geneid={gene_id}] [gene-coordinates={coord_str}({direction})] "
                f"{species} mitochondrially encoded {product_name}"
            )

            # Write one FASTA per CDS
            out_fa = out_cds / f"{gene_name}.{assembly}.fa"
            with out_fa.open("w") as fh:
                fh.write(f"{header}\n{wrap(str(cds_seq))}\n")

            print(f"[extract_cds_from_tbl] Extracted {gene_name} -> {out_fa}")

    print("[extract_cds_from_tbl] CDS extraction complete.")

if __name__ == "__main__":
    main()
