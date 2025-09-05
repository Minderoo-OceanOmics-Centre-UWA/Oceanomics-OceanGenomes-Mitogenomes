#!/usr/bin/env python3
# translate_genes.py
from pathlib import Path
import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser(
        description="Translate mitochondrial CDS FASTAs (table 2)."
    )
    p.add_argument("--input", required=True, type=Path,
                   help="Either a genes directory (containing *.fa) OR a single multi-FASTA (*.genes.fa)")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Output directory (proteins will be placed under outdir/proteins)")
    p.add_argument("--table", type=int, default=2,
                   help="NCBI translation table (default=2 vertebrate mitochondrial)")
    return p.parse_args()

def translate_file(src_fa: Path, dst_fa: Path, table: int):
    with src_fa.open() as fin, dst_fa.open("w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            prot = Seq(record.seq).translate(table=table)
            wrapped = "\n".join(str(prot)[i:i+70] for i in range(0, len(prot), 70))
            fout.write(f">{record.description}\n{wrapped}\n")

def main():
    args = parse_args()
    if not args.input.exists():
        sys.exit(f"Input not found: {args.input}")
    prot_dir = args.outdir / "proteins"
    prot_dir.mkdir(parents=True, exist_ok=True)

    if args.input.is_dir():
        # translate every .fa/.fasta in the directory
        for fa in sorted(args.input.glob("*.fa")) + sorted(args.input.glob("*.fasta")):
            # keep name, append .protein.fa
            out = prot_dir / (fa.name.replace(".fasta", ".protein.fa").replace(".fa", ".protein.fa"))
            translate_file(fa, out, args.table)
            print(f"[translate_genes] {fa} -> {out}")
    else:
        # single file; produce <name>.protein.fa (e.g., sample.genes.fa -> sample.genes.protein.fa)
        base = args.input.name
        out_name = base.replace(".fasta", ".protein.fa").replace(".fa", ".protein.fa")
        out = prot_dir / out_name
        translate_file(args.input, out, args.table)
        print(f"[translate_genes] {args.input} -> {out}")

if __name__ == "__main__":
    main()
