#!/usr/bin/env python3
"""Build a curated coral (Anthozoa) mitogenome reference database for the pipeline.

The annotation reference is otherwise resolved from each sample's species *label*
(MITOHIFI_FINDMITOREFERENCE -> NCBI by name), which silently yields a wrong-family
reference when the label is wrong or coarse. This builds a label-free database the
pipeline can pick a reference from by *sequence* similarity to the assembly.

Source: NCBI RefSeq Anthozoa complete mitogenomes ONLY, and only records that are
*fully annotated* -- they must carry the features the coral fixer transfers and a
complete gene set, so every reference in the DB is usable:

    * 16S rRNA (rrnL) feature,
    * 12S rRNA (rrnS) feature,
    * a nad5 CDS, and
    * at least --min-cds protein-coding genes (default 13, the cnidarian PCG set).

This RefSeq-only build is freely shareable. In-house assemblies are NOT included
by default (licensing); add them later for a private, more complete build with
--extra-gb (one or more GenBank files / directories), which are filtered by the
same completeness bar and merged in.

Outputs (into --out-dir):
    coral_mito_refdb.gb          combined GenBank (canonical artifact; has features)
    coral_mito_refdb.fasta       all genome sequences (selection + GetOrganelle seed -s)
    coral_mito_refdb.label.fasta GetOrganelle label/gene database (--genes), headers
                                 '>gene type - accession--Organism'
    coral_mito_refdb.manifest.tsv  accession, organism, family, length, n_cds, n_rrna, source
    coral_mito_refdb.nucl.*       BLAST nucleotide DB (if makeblastdb is available)

The .fasta (whole genomes) doubles as the GetOrganelle seed and the .label.fasta as
its custom label database -- both used ONLY for invertebrate (coral) samples.

Runs in the MITOS2 BioContainer (biopython + blast). NCBI access needs --email.

Usage:
    build_coral_reference_db.py --email you@uwa.edu.au --out-dir coral_refdb
    build_coral_reference_db.py --email you@uwa.edu.au --out-dir coral_refdb \
        --extra-gb my_curated_gb/        # private build with in-house records
    build_coral_reference_db.py --out-dir coral_refdb --no-download \
        --extra-gb my_curated_gb/        # offline; in-house records only
"""
import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO

# Anthozoa NCBI taxid; mitochondrial complete RefSeq genomes only.
DEFAULT_TAXON = "txid6101"  # Anthozoa
SEARCH_TEMPLATE = ('{taxon}[Organism:exp] AND refseq[filter] AND '
                   'mitochondrion[filter] AND "complete genome"[Title]')


def label(feat):
    q = feat.qualifiers
    return " ".join(q.get("gene", []) + q.get("product", [])).upper()


def record_is_complete(rec, min_cds):
    """A record qualifies only if fully annotated: 16S + 12S rRNA, a nad5 CDS, and
    >= min_cds protein-coding genes. Returns (ok, n_cds, n_rrna, reason)."""
    n_cds = sum(1 for f in rec.features if f.type == "CDS")
    rrna = [label(f) for f in rec.features if f.type == "rRNA"]
    has_16s = any("16S" in l or "RRNL" in l or "LARGE" in l for l in rrna)
    has_12s = any("12S" in l or "RRNS" in l or "SMALL" in l for l in rrna)
    has_nad5 = any(f.type == "CDS" and ("ND5" in label(f) or "NAD5" in label(f)
                   or "SUBUNIT 5" in label(f)) for f in rec.features)
    if not has_16s:
        return False, n_cds, len(rrna), "no 16S rRNA"
    if not has_12s:
        return False, n_cds, len(rrna), "no 12S rRNA"
    if not has_nad5:
        return False, n_cds, len(rrna), "no nad5 CDS"
    if n_cds < min_cds:
        return False, n_cds, len(rrna), f"only {n_cds} CDS (<{min_cds})"
    return True, n_cds, len(rrna), "ok"


def family_of(rec):
    lineage = rec.annotations.get("taxonomy", []) or []
    return next((t for t in lineage if t.endswith("idae")), "")


def write_label_db(records, path):
    """Write a GetOrganelle label database from the kept records' gene features.

    GetOrganelle's --genes database labels assembly-graph contigs by organelle gene
    during disentangling; for divergent animal/coral mitogenomes a custom one helps
    where the built-in animal_mt labels do not. Header format mirrors GetOrganelle's
    get_annotated_regions_from_gb.py: '>gene type - accession--Organism_no_spaces'.
    """
    from Bio.SeqIO import write as _write
    from Bio.SeqRecord import SeqRecord
    out = []
    for rec in records:
        acc = rec.id
        org = (rec.annotations.get("organism", "") or "").replace(" ", "_")
        for feat in rec.features:
            if feat.type not in ("CDS", "rRNA", "tRNA"):
                continue
            names = feat.qualifiers.get("gene") or feat.qualifiers.get("product")
            if not names:
                continue
            gene = names[0].replace(" ", "_")
            try:
                seq = feat.extract(rec.seq)
            except Exception:
                continue
            if len(seq) < 30:
                continue
            sr = SeqRecord(seq, id=f"{gene}", description="")
            sr.description = f"{feat.type} - {acc}--{org}"
            out.append(sr)
    with open(path, "w") as fh:
        for sr in out:
            fh.write(f">{sr.id} {sr.description}\n{str(sr.seq)}\n")
    return len(out)


def fetch_refseq(taxon, email, api_key, batch, retmax):
    """esearch + batched efetch of GenBank records from NCBI nucleotide."""
    from Bio import Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    query = SEARCH_TEMPLATE.format(taxon=taxon)
    print(f"[build_db] esearch: {query}", file=sys.stderr)
    with Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y") as h:
        res = Entrez.read(h)
    ids = res["IdList"]
    print(f"[build_db] {len(ids)} RefSeq accessions found", file=sys.stderr)
    webenv, qkey = res["WebEnv"], res["QueryKey"]
    out = []
    for start in range(0, len(ids), batch):
        for attempt in range(4):
            try:
                with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                   retstart=start, retmax=batch,
                                   webenv=webenv, query_key=qkey) as h:
                    out.extend(SeqIO.parse(h, "genbank"))
                break
            except Exception as exc:  # transient NCBI error -> back off and retry
                print(f"[build_db] efetch {start} attempt {attempt+1} failed: {exc}",
                      file=sys.stderr)
                time.sleep(3 * (attempt + 1))
        print(f"[build_db] fetched {min(start+batch, len(ids))}/{len(ids)}", file=sys.stderr)
        time.sleep(0.34)  # stay under NCBI rate limits
    return out


def load_extra(paths):
    recs = []
    for p in paths:
        p = Path(p)
        files = sorted(p.glob("*.gb")) + sorted(p.glob("*.gbk")) + sorted(p.glob("*.gbf")) \
            if p.is_dir() else [p]
        for f in files:
            for rec in SeqIO.parse(str(f), "genbank"):
                recs.append((rec, f"local:{f.name}"))
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--email", help="NCBI Entrez email (required unless --no-download)")
    ap.add_argument("--api-key", help="NCBI API key (optional; raises rate limit)")
    ap.add_argument("--taxon", default=DEFAULT_TAXON,
                    help="NCBI taxon for the search (default txid6101 = Anthozoa)")
    ap.add_argument("--min-cds", type=int, default=13,
                    help="Minimum protein-coding genes for a record to qualify (default 13).")
    ap.add_argument("--extra-gb", nargs="*", default=[],
                    help="In-house GenBank files/dirs to merge (private build). "
                         "Filtered by the same completeness bar. Off by default.")
    ap.add_argument("--no-download", action="store_true",
                    help="Skip NCBI; build from --extra-gb only.")
    ap.add_argument("--retmax", type=int, default=2000)
    ap.add_argument("--batch", type=int, default=200)
    args = ap.parse_args()

    if not args.no_download and not args.email:
        ap.error("--email is required for NCBI download (or pass --no-download)")

    candidates = []  # (rec, source)
    if not args.no_download:
        candidates += [(r, "refseq") for r in
                       fetch_refseq(args.taxon, args.email, args.api_key, args.batch, args.retmax)]
    if args.extra_gb:
        candidates += load_extra(args.extra_gb)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    kept, manifest, seen = [], [], set()
    n_drop = 0
    for rec, source in candidates:
        acc = rec.id.split(".")[0]
        if acc in seen:
            continue
        ok, n_cds, n_rrna, reason = record_is_complete(rec, args.min_cds)
        if not ok:
            n_drop += 1
            continue
        seen.add(acc)
        kept.append(rec)
        manifest.append((rec.id, rec.annotations.get("organism", ""), family_of(rec),
                         len(rec.seq), n_cds, n_rrna, source))

    if not kept:
        sys.exit("[build_db] no records passed the completeness filter; nothing written")

    gb = args.out_dir / "coral_mito_refdb.gb"
    fa = args.out_dir / "coral_mito_refdb.fasta"
    lab = args.out_dir / "coral_mito_refdb.label.fasta"
    mf = args.out_dir / "coral_mito_refdb.manifest.tsv"
    SeqIO.write(kept, str(gb), "genbank")
    SeqIO.write(kept, str(fa), "fasta")
    n_label = write_label_db(kept, lab)
    with open(mf, "w") as fh:
        fh.write("accession\torganism\tfamily\tlength_bp\tn_cds\tn_rrna\tsource\n")
        for row in sorted(manifest, key=lambda r: (r[2], r[1])):
            fh.write("\t".join(map(str, row)) + "\n")

    # BLAST nucleotide DB for the (future) sequence-based selection step.
    if shutil.which("makeblastdb"):
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                        "-out", str(args.out_dir / "coral_mito_refdb.nucl"),
                        "-title", "coral_mito_refdb"], check=True)
    else:
        print("[build_db] makeblastdb not found; skipped BLAST DB (run it later on the fasta)",
              file=sys.stderr)

    fams = sorted({m[2] for m in manifest if m[2]})
    print(f"[build_db] kept {len(kept)} records, dropped {n_drop} (incomplete).")
    print(f"[build_db] families ({len(fams)}): {', '.join(fams)}")
    print(f"[build_db] label db: {n_label} gene sequences")
    print(f"[build_db] wrote {gb.name}, {fa.name}, {lab.name}, {mf.name} to {args.out_dir}/")


if __name__ == "__main__":
    main()
