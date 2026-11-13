#!/usr/bin/env python3
"""
Rebuild a precomputed annotation bundle's cds/ and proteins/ from its GFF + FASTA, and
write the BLAST/LCA regions (CO1, 12S = RNR1, 16S = RNR2).

Used when --skip_mitogenome_annotation is given with --blast_precomputed_annotations:
the annotators did not run, so nothing produced this run's cds/, proteins/ or BLAST
inputs. A file glob cannot carry the bundle's cds/ and proteins/ directories either, and
without them ANNOTATION_STATS cannot run the PCG ORF / integrity check -- the QC gate then
records `not_evaluated` and passes a CO1 with no stop codon. So both are rebuilt here,
per gene, with the GFF as the source of truth: a hand-curated bundle can have an edited
GFF and a cds/ that was never regenerated, and BLASTing or QC-ing that would silently
check the pre-curation sequence.

Per protein-coding / rRNA gene, the GFF-derived sequence is compared with the bundle's
own cds/ file, found beside the GFF (following the staging symlink):

    match    identical.
    trimmed  cds/ is the GFF sequence minus up to MAX_TERMINAL_TRIM trailing bases. EMMA
             writes cds/ this way (it drops the stop codon and the last base or two of
             each rRNA), so it is consistent with the GFF, and the bundle's cds/ and
             proteins/ files are reused as they are -- results then match the run that
             produced the bundle.
    stale    anything else: the GFF was curated and cds/ not regenerated. The GFF wins and
             the protein is re-translated under --genetic-code.
    no_cds   no cds/ file to compare with; GFF-derived, protein translated.

Output (same names and headers as the annotators, so getAnnotationName(), the BLAST/LCA
filenames and annotation_stats.py all read it unchanged):

    <outdir>/<stem>.fa, <stem>.gff, ...    the bundle's top-level files, copied
    <outdir>/cds/<GENE>.<stem>.fa          header: ><stem>
    <outdir>/proteins/MT-<GENE>.<stem>.fa  protein-coding genes only
    <regions>/{CO1,RNR1,RNR2}.<stem>.fa    BLAST inputs; a region the GFF lacks is not written
    <regions>/<stem>.bundle_rebuild.tsv    one row per gene: status, length, cds_check, notes

Input files are given positionally (the staged bundle); the GFF and the FASTA with the
same stem are picked out of them.
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# Sibling import: Nextflow bind-mounts the whole bin/ dir onto PATH, so this
# resolves via sys.path[0] (same pattern as extract_genes_gff.py / orf_utils).
from mitos_to_emma import EXTRA_PCG_MAP, PCG_MAP, RRNA_MAP, gene_nt_seq, map_gene_name

REGIONS = ("CO1", "RNR1", "RNR2")

# EMMA's cds/ drops at most a stop codon's worth of trailing bases (see module docstring).
MAX_TERMINAL_TRIM = 3

PCG_NAMES = set(PCG_MAP.values()) | set(EXTRA_PCG_MAP.values())
RRNA_NAMES = set(RRNA_MAP.values())

# Leaf feature types that carry a gene's sequence. gene/mRNA rows span the whole locus
# (including any intron), so a gene row is only a fallback for a gene with no leaf row.
LEAF_TYPES = {"CDS", "rRNA"}
FALLBACK_TYPES = {"gene"}

# Names seen in GFFs that neither "MT-" stripping nor map_gene_name resolves.
ALIASES = {
    "COX1": "CO1", "COI": "CO1", "COXI": "CO1",
    "12S": "RNR1", "RRNS": "RNR1", "S-RRNA": "RNR1", "12S_RRNA": "RNR1",
    "16S": "RNR2", "RRNL": "RNR2", "L-RRNA": "RNR2", "16S_RRNA": "RNR2",
}


def gene_for(name):
    """Map a GFF Name= value to (EMMA gene name, 'CDS'|'rRNA'), or (None, None) for
    tRNAs, origins and anything unrecognised."""
    raw = (name or "").strip()
    bare = raw[3:] if raw.upper().startswith("MT-") else raw
    up = ALIASES.get(bare.upper(), bare.upper())
    if up in PCG_NAMES:
        return up, "CDS"
    if up in RRNA_NAMES:
        return up, "rRNA"
    emma_name, ftype, _frag, _ac = map_gene_name(bare)
    if ftype in ("CDS", "rRNA"):
        return emma_name, ftype
    return None, None


def parse_attrs(col):
    return dict(a.split("=", 1) for a in col.strip().split(";") if "=" in a)


def pick_bundle(files):
    """Return (stem, gff, fasta) from the staged bundle files."""
    paths = [Path(f) for f in files]
    gffs = [p for p in paths if p.suffix in (".gff", ".gff3") and p.is_file()]
    if not gffs:
        sys.exit("ERROR: no .gff in the annotation bundle: " + ", ".join(p.name for p in paths))
    if len(gffs) > 1:
        sys.exit("ERROR: more than one GFF in one assembly's annotation bundle ("
                 + ", ".join(sorted(p.name for p in gffs)) + "). Narrow "
                 "--precomputed_mitogenome_annotation_results to one annotation per assembly.")
    gff = gffs[0]
    stem = gff.stem
    fastas = [p for p in paths if p.is_file() and p.suffix in (".fa", ".fasta", ".fna")
              and p.stem == stem]
    if not fastas:
        sys.exit(f"ERROR: no {stem}.fa beside {gff.name} in the annotation bundle")
    return stem, gff, fastas[0]


def read_features(gff):
    """gene -> (kind, rows [(seqid, start, end, strand, attrs)] in file order), preferring
    leaf rows (CDS/rRNA) and falling back to gene rows only for a gene with no leaf.
    Genes are returned in the order they first appear in the GFF."""
    leaf, fallback, kinds, order = {}, {}, {}, []
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in LEAF_TYPES and ftype not in FALLBACK_TYPES:
                continue
            attrs = parse_attrs(cols[8])
            gene, kind = gene_for(attrs.get("Name"))
            if not gene:
                continue
            if gene not in kinds:
                kinds[gene] = kind
                order.append(gene)
            row = (cols[0], int(cols[3]), int(cols[4]), cols[6], attrs)
            (leaf if ftype in LEAF_TYPES else fallback).setdefault(gene, []).append(row)
    return {g: (kinds[g], leaf.get(g) or fallback[g]) for g in order}


def first_locus(rows):
    """Keep the rows of the first locus only. Exons of one gene share a Parent; a
    second Parent for the same gene is a duplicate copy, which the annotators also
    never emit twice into cds/."""
    parent = rows[0][4].get("Parent")
    if parent is None:
        return rows[:1], len(rows) > 1
    keep = [r for r in rows if r[4].get("Parent") == parent]
    return keep, len(keep) != len(rows)


def row_seq(genome, row):
    seqid, start, end, strand, attrs = row
    # Origin-spanning features are written up to the sequence end with the wrapped end
    # in the attributes (mitos_to_emma._gff_span); restore start > end so gene_nt_seq
    # joins tail + head across the circular join.
    if "origin-spanning" in attrs.get("Note", "") and "wrapped_end" in attrs:
        end = int(attrs["wrapped_end"])
    return gene_nt_seq(genome, seqid, start, end, strand)


def bundle_file(gff, subdir, names):
    """The bundle's own <subdir>/<name> file, or None. Looked for beside the GFF's real
    location, since the bundle is staged as symlinks and cds/ / proteins/ are not globbed."""
    for base in (Path(os.path.realpath(gff)).parent, Path(gff).parent):
        for name in names:
            if (base / subdir / name).is_file():
                return base / subdir / name
    return None


def fasta_seq(path):
    return "".join(str(r.seq) for r in SeqIO.parse(str(path), "fasta"))


def cds_check(extracted, published):
    ext, pub = extracted.upper(), published.upper()
    if pub == ext:
        return "match"
    if pub and ext.startswith(pub) and len(ext) - len(pub) <= MAX_TERMINAL_TRIM:
        return "trimmed"
    return "stale"


def translate(nt, code, stem, gene):
    """Translate the way mitos_to_emma.py does: the whole GFF span, terminator included,
    so a re-translated protein is identical to the one the MITOS path would publish."""
    if code is None:
        sys.exit(f"ERROR: {stem}: {gene} has to be re-translated but no --genetic-code was given")
    nt = nt[: len(nt) - len(nt) % 3]
    return str(Seq(nt).translate(table=code))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--outdir", required=True, type=Path, help="Rebuilt bundle directory")
    p.add_argument("--regions-dir", required=True, type=Path, help="BLAST region FASTAs + status TSV")
    p.add_argument("--genetic-code", type=int, default=None,
                   help="Mitochondrial translation table (meta.genetic_code); needed only "
                        "when a protein has to be re-translated")
    p.add_argument("files", nargs="+", help="Staged annotation bundle files")
    args = p.parse_args()

    stem, gff, fasta = pick_bundle(args.files)
    genome = SeqIO.to_dict(SeqIO.parse(str(fasta), "fasta"))
    if not genome:
        sys.exit(f"ERROR: {fasta.name} has no FASTA records")
    # A single-record FASTA is matched regardless of the GFF's seqid (curation tools
    # often rename the record); a multi-record one must match by seqid.
    only = next(iter(genome.values())) if len(genome) == 1 else None

    cds_dir, prot_dir = args.outdir / "cds", args.outdir / "proteins"
    for d in (cds_dir, prot_dir, args.regions_dir):
        d.mkdir(parents=True, exist_ok=True)
    for f in args.files:
        f = Path(f)
        if f.is_file():
            shutil.copyfile(f, args.outdir / f.name)

    features = read_features(gff)
    status_rows = []
    for gene, (kind, rows) in features.items():
        rows, duplicate = first_locus(rows)
        seq = ""
        for row in rows:
            g = genome if row[0] in genome else ({row[0]: only} if only is not None else None)
            if g is None:
                sys.exit(f"ERROR: {gff.name} references seqid '{row[0]}', not in {fasta.name}")
            seq += str(row_seq(g, row))
        if not seq:
            status_rows.append((gene, kind, "missing", 0, "", "empty_sequence"))
            continue

        cds = bundle_file(gff, "cds", (f"{gene}.{stem}.fa", f"MT-{gene}.{stem}.fa"))
        check = "no_cds" if cds is None else cds_check(seq, fasta_seq(cds))
        reuse = check in ("match", "trimmed")
        if reuse:
            seq = fasta_seq(cds)
        with open(cds_dir / f"{gene}.{stem}.fa", "w") as out:
            out.write(f">{stem}\n{seq}\n")

        notes = []
        if kind == "CDS":
            prot = bundle_file(gff, "proteins", (f"MT-{gene}.{stem}.fa",)) if reuse else None
            aa = fasta_seq(prot) if prot else translate(seq, args.genetic_code, stem, gene)
            if not prot:
                notes.append("protein_translated")
            with open(prot_dir / f"MT-{gene}.{stem}.fa", "w") as out:
                out.write(f">{stem}\n{aa}\n")
        if gene in REGIONS:
            with open(args.regions_dir / f"{gene}.{stem}.fa", "w") as out:
                out.write(f">{stem}\n{seq}\n")
        if len(rows) > 1:
            notes.append(f"exons={len(rows)}")
        if duplicate:
            notes.append("duplicate_locus_ignored")
        status_rows.append((gene, kind, "found", len(seq), check, ";".join(notes)))

    found = {r[0] for r in status_rows if r[2] == "found"}
    for region in REGIONS:
        if region not in features:
            status_rows.append((region, "CDS" if region == "CO1" else "rRNA", "missing", 0, "", ""))

    with open(args.regions_dir / f"{stem}.bundle_rebuild.tsv", "w") as out:
        out.write("annotation\tgene\tkind\tstatus\tlength\tcds_check\tnotes\n")
        for row in status_rows:
            out.write(stem + "\t" + "\t".join(str(v) for v in row) + "\n")

    stale = [r[0] for r in status_rows if r[4] == "stale"]
    if stale:
        print(f"WARNING: {stem}: bundle cds/ differs from the GFF for {', '.join(stale)}; "
              "using the GFF-derived sequence (cds/ was likely not regenerated after curation)",
              file=sys.stderr)
    n_regions = sum(1 for r in REGIONS if r in found)
    print(f"[rebuild_annotation_bundle] {stem}: {len(found)} genes rebuilt, "
          f"{n_regions}/3 BLAST regions")


if __name__ == "__main__":
    main()
