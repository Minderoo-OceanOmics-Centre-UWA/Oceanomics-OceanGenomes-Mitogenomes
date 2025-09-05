#!/usr/bin/env python3
# extract_genes.py
from pathlib import Path
import argparse, re, sys
from Bio import SeqIO

# Name→default Product for older GFFs (unchanged from your script)
NAME2DEFAULT_PRODUCT = {
    "RNR1": "12S rRNA",
    "RNR2": "16S rRNA",
    "TH"  : "tRNA-His(GUG)",
    "TW"  : "tRNA-Trp(UCA)",
    "TA"  : "tRNA-Ala(UGC)",
    "TN"  : "tRNA-Asn(GUU)",
    "TC"  : "tRNA-Cys(GCA)",
    "TY"  : "tRNA-Tyr(GUA)",
    "TS1" : "tRNA-Ser(GCU)",
    "TL1" : "tRNA-Leu(UAG)",
    "TS2" : "tRNA-Ser(UGA)",
    "TD"  : "tRNA-Asp(GUC)",
    "TK"  : "tRNA-Lys(UUU)",
    "TG"  : "tRNA-Gly(UCC)",
    "TR"  : "tRNA-Arg(UCG)",
    "TE"  : "tRNA-Glu(UUC)",
    "TT"  : "tRNA-Thr(UGU)",
    "TP"  : "tRNA-Pro(UGG)",
    "TF"  : "tRNA-Phe(GAA)",
    "TV"  : "tRNA-Val(UAC)",
    "TL2" : "tRNA-Leu(UAA)",
    "TI"  : "tRNA-Ile(GAU)",
    "TQ"  : "tRNA-Gln(UUG)",
    "TM"  : "tRNA-Met(CAU)",
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract mitochondrial genes from FASTA+GFF for a single sample."
    )
    p.add_argument("--fasta", required=True, type=Path, help="Genome FASTA")
    p.add_argument("--gff",   required=True, type=Path, help="Matching GFF")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("--assembly", default=None,
                   help="Assembly/sample ID; default=FASTA stem")
    return p.parse_args()

def main():
    args = parse_args()
    if not args.fasta.exists():
        sys.exit(f"FASTA not found: {args.fasta}")
    if not args.gff.exists():
        sys.exit(f"GFF not found: {args.gff}")

    assembly = args.assembly or args.fasta.stem
    args.outdir.mkdir(parents=True, exist_ok=True)
    gene_dir = args.outdir / "genes"
    gene_dir.mkdir(exist_ok=True)

    genome_record = SeqIO.read(str(args.fasta), "fasta")

    genes = {}               # name -> {start,end,strand,id,product}
    genes_with_cds = set()
    name2product = {}
    RNA_TYPES = {"tRNA", "rRNA"}
    species = "Unknown organism"

    with args.gff.open() as gff:
        for ln in gff:
            if ln.startswith("##organism"):
                species = ln.split("##organism", 1)[1].strip() or species
                continue
            if ln.startswith("#"):
                continue
            cols = ln.rstrip().split("\t")
            if len(cols) < 9:
                continue
            ftype, start, end, strand, attrs = cols[2], int(cols[3]), int(cols[4]), cols[6], cols[8]
            attr = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
            name = attr.get("Name")
            if not name:
                continue
            if "Product" in attr:
                name2product[name] = attr["Product"]

            if ftype == "gene":
                genes[name] = dict(
                    start=start, end=end, strand=strand,
                    id=attr.get("ID", "unk_id").split("|")[-1], product=""
                )
            elif ftype == "CDS":
                genes_with_cds.add(name)
                # ensure gene exists even if 'gene' line is missing
                genes.setdefault(name, dict(
                    start=start, end=end, strand=strand,
                    id=attr.get("ID", "unk_id").split("|")[-1], product=""
                ))
            elif ftype in RNA_TYPES:
                genes.setdefault(name, dict(
                    start=start, end=end, strand=strand,
                    id=attr.get("ID", "unk_id").split("|")[-1], product=""
                ))

    # fill products
    for name, info in genes.items():
        raw = name2product.get(name) or NAME2DEFAULT_PRODUCT.get(name) or name
        clean = re.sub(r"\([^)]*\)$", "", raw).strip()
        info["product"] = clean

    ordered = sorted(genes.items(), key=lambda kv: kv[1]["start"])
    concat_path = gene_dir / f"{assembly}.genes.fa"
    with concat_path.open("w") as concat_fh:
        for name, info in ordered:
            s, e, strand, gene_id, product = info["start"], info["end"], info["strand"], info["id"], info["product"]
            seq = (genome_record.seq[s-1:e] if strand == "+"
                   else genome_record.seq[s-1:e].reverse_complement())
            coord = f"{s}-{e}"
            header = (
                f">{assembly}|{coord}|{strand}|{name} "
                f"[organism={species}] [mgcode=2] [topology=linear] "
                f"[geneid={gene_id}] [gene-coordinates={coord}({strand})] "
                f"{species} mitochondrially encoded {product}"
            )
            wrapped = "\n".join(str(seq)[i:i+70] for i in range(0, len(seq), 70))
            concat_fh.write(f"{header}\n{wrapped}\n")

            if name in genes_with_cds:
                single_path = gene_dir / f"{name}.{assembly}.fa"
                with single_path.open("w") as out_fa:
                    out_fa.write(f"{header}\n{wrapped}\n")

    print(f"[extract_genes] wrote concatenated: {concat_path}")
    print(f"[extract_genes] wrote per-CDS FASTAs into: {gene_dir}")

if __name__ == "__main__":
    main()
