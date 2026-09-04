#!/usr/bin/env python3
"""Build curated per-phylum invertebrate mitogenome reference databases.

Each group's database is used twice by the pipeline, and both uses are label-free:

  * GetOrganelle reseed -- <group>_mito_refdb.fasta is the seed (-s) and
    <group>_mito_refdb.label.fasta the custom gene/label database (--genes) for
    an invertebrate whose first-pass assembly failed. Until now every
    invertebrate, mollusc and sea star included, was reseeded from the Anthozoa
    database, which is too divergent to help.
  * Annotation reference (Cnidaria only so far) -- SELECT_CORAL_REFERENCE BLASTs
    the assembly against <group>_mito_refdb.gb and picks a reference by sequence
    similarity, instead of resolving one from the sample's species *label*
    (MITOHIFI_FINDMITOREFERENCE -> NCBI by name), which silently yields a
    wrong-family reference when the label is wrong or coarse.

Source: NCBI RefSeq complete mitogenomes for the group's taxon, keeping only
records annotated well enough to be usable as a reference. The completeness bar
is per group (GROUPS below), because it is not the same across phyla: coral
references must carry the features the coral fixer transfers (both rRNAs and a
nad5 CDS) plus the full 13-PCG cnidarian set, while ctenophore mitogenomes are
genuinely reduced -- no atp6, no tRNAs, ~10 PCGs, rRNAs often unannotated -- so
the coral bar would reject every valid ctenophore record.

RefSeq-only builds are freely shareable. In-house assemblies are NOT included by
default (licensing); add them for a private, more complete build with --extra-gb
(one or more GenBank files / directories), filtered by the same bar and merged in.

Outputs (into --out-dir, default assets/refdb/<group>/):
    <group>_mito_refdb.gb          combined GenBank (canonical artifact; has features)
    <group>_mito_refdb.fasta       all genome sequences (selection + GetOrganelle seed -s)
    <group>_mito_refdb.label.fasta GetOrganelle label/gene database (--genes), headers
                                   '>gene type - accession--Organism'
    <group>_mito_refdb.manifest.tsv  group, accession, organism, family, length, n_cds, n_rrna, source
    <group>_mito_refdb.nucl.*      BLAST nucleotide DB (if makeblastdb is available)

Runs in the MITOS2 BioContainer (biopython + blast). NCBI access needs --email.

Usage:
    build_invert_reference_db.py --group mollusca --email you@uwa.edu.au
    build_invert_reference_db.py --all --email you@uwa.edu.au
    build_invert_reference_db.py --group anthozoa --email you@uwa.edu.au \
        --extra-gb my_curated_gb/        # private build with in-house records
    build_invert_reference_db.py --group anthozoa --no-download \
        --extra-gb my_curated_gb/        # offline; in-house records only
"""
import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO

SEARCH_TEMPLATE = ('({organism}) AND refseq[filter] AND '
                   'mitochondrion[filter] AND "complete genome"[Title]')

# Per-group build spec. `organism` is the Entrez organism expression for the
# group -- usually one taxid, but Arthropoda has to subtract the terrestrial
# radiations (Insecta/Hexapoda, Arachnida, Myriapoda) that otherwise swamp the
# search: they are >95% of arthropod RefSeq mitogenomes, none of them is
# anything OceanOmics sequences, and including them buries the marine
# crustaceans and pycnogonids the seed is actually for. min_cds is the group's
# expected protein-coding gene count; require_rrna demands BOTH rRNAs be
# annotated; require_nad5 is the coral fixer's extra requirement (it transfers
# nad5 across the giant group I intron). Keys match
# InvertTaxonGroups.seedDbGroup() in lib/ -- a group renamed here must be
# renamed there, or the reseed will look for a database that does not exist.
GROUPS = {
    # group:          organism expression                    min_cds  rrna   nad5
    "anthozoa":       ("txid6101[Organism:exp]",             13,      True,  True),
    "porifera":       ("txid6040[Organism:exp]",             13,      True,  True),
    "mollusca":       ("txid6447[Organism:exp]",             12,      True,  False),
    "arthropoda":     ("txid6657[Organism:exp] "
                       "NOT txid6960[Organism:exp] "     # Hexapoda (insects, springtails)
                       "NOT txid6854[Organism:exp] "     # Arachnida
                       "NOT txid61985[Organism:exp]",    # Myriapoda
                                                             13,      True,  False),
    "echinodermata":  ("txid7586[Organism:exp]",             13,      True,  False),
    "ctenophora":     ("txid10197[Organism:exp]",            10,      False, False),
    "tunicata":       ("txid7712[Organism:exp]",             12,      True,  False),
    "annelida":       ("txid6340[Organism:exp]",             12,      True,  False),
}

# rRNA product/gene synonyms across invertebrate annotation conventions. The
# coral-era matcher only knew 16S/RRNL/LARGE and 12S/RRNS/SMALL, which drops
# correctly annotated mollusc, arthropod and annelid records that use the
# l-rRNA / rnl / MT-RNR2 spellings instead.
LARGE_RRNA = ("16S", "RRNL", "LARGE", "L-RRNA", "LRRNA", "RNL", "MT-RNR2", "RRN16")
SMALL_RRNA = ("12S", "RRNS", "SMALL", "S-RRNA", "SRRNA", "RNS", "MT-RNR1", "RRN12")


def label(feat):
    q = feat.qualifiers
    return " ".join(q.get("gene", []) + q.get("product", [])).upper()


def record_is_complete(rec, min_cds, require_rrna, require_nad5):
    """A record qualifies only if annotated to the group's bar. Returns
    (ok, n_cds, n_rrna, reason)."""
    n_cds = sum(1 for f in rec.features if f.type == "CDS")
    rrna = [label(f) for f in rec.features if f.type == "rRNA"]
    has_16s = any(any(k in l for k in LARGE_RRNA) for l in rrna)
    has_12s = any(any(k in l for k in SMALL_RRNA) for l in rrna)
    has_nad5 = any(f.type == "CDS" and ("ND5" in label(f) or "NAD5" in label(f)
                   or "SUBUNIT 5" in label(f)) for f in rec.features)
    if require_rrna and not has_16s:
        return False, n_cds, len(rrna), "no 16S rRNA"
    if require_rrna and not has_12s:
        return False, n_cds, len(rrna), "no 12S rRNA"
    if require_nad5 and not has_nad5:
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
    during disentangling; for divergent animal mitogenomes a custom one helps where
    the built-in animal_mt labels do not. Header format mirrors GetOrganelle's
    get_annotated_regions_from_gb.py: '>gene type - accession--Organism_no_spaces'.
    """
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


def fetch_refseq(organism, email, api_key, batch, retmax):
    """esearch + batched efetch of GenBank records from NCBI nucleotide."""
    from Bio import Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    query = SEARCH_TEMPLATE.format(organism=organism)
    print(f"[build_db] esearch: {query}", file=sys.stderr)
    with Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y") as h:
        res = Entrez.read(h)
    ids = res["IdList"]
    print(f"[build_db] {len(ids)} RefSeq accessions found "
          f"(of {res['Count']} matching)", file=sys.stderr)
    if int(res["Count"]) > len(ids):
        # Silently building from the first retmax hits would ship a database that
        # is a arbitrary slice of the group rather than the group.
        sys.exit(f"[build_db] search returned {res['Count']} records but retmax is "
                 f"{retmax}; raise --retmax or narrow the group's organism expression")
    if not ids:
        return []
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


def default_out_dir(group):
    """assets/refdb/<group>/ relative to the repo this script lives in."""
    return Path(__file__).resolve().parent.parent / "assets" / "refdb" / group


def build_group(group, args):
    organism, min_cds, require_rrna, require_nad5 = GROUPS[group]
    organism = args.taxon or organism
    min_cds = args.min_cds if args.min_cds is not None else min_cds
    out_dir = args.out_dir or default_out_dir(group)

    print(f"[build_db] === {group}: {organism}, min_cds={min_cds}, "
          f"require_rrna={require_rrna}, require_nad5={require_nad5} ===", file=sys.stderr)

    candidates = []  # (rec, source)
    if not args.no_download:
        candidates += [(r, "refseq") for r in
                       fetch_refseq(organism, args.email, args.api_key, args.batch, args.retmax)]
    if args.extra_gb:
        candidates += load_extra(args.extra_gb)

    out_dir.mkdir(parents=True, exist_ok=True)
    kept, manifest, seen = [], [], set()
    n_drop = 0
    for rec, source in candidates:
        acc = rec.id.split(".")[0]
        if acc in seen:
            continue
        ok, n_cds, n_rrna, reason = record_is_complete(rec, min_cds, require_rrna, require_nad5)
        if not ok:
            n_drop += 1
            continue
        seen.add(acc)
        kept.append(rec)
        manifest.append((group, rec.id, rec.annotations.get("organism", ""), family_of(rec),
                         len(rec.seq), n_cds, n_rrna, source))

    if not kept:
        sys.exit(f"[build_db] {group}: no records passed the completeness filter; nothing written")

    base = out_dir / f"{group}_mito_refdb"
    gb, fa = base.with_suffix(".gb"), base.with_suffix(".fasta")
    lab = Path(f"{base}.label.fasta")
    mf = Path(f"{base}.manifest.tsv")
    SeqIO.write(kept, str(gb), "genbank")
    SeqIO.write(kept, str(fa), "fasta")
    n_label = write_label_db(kept, lab)
    with open(mf, "w") as fh:
        fh.write("group\taccession\torganism\tfamily\tlength_bp\tn_cds\tn_rrna\tsource\n")
        for row in sorted(manifest, key=lambda r: (r[3], r[2])):
            fh.write("\t".join(map(str, row)) + "\n")

    # BLAST nucleotide DB for the sequence-based reference selection step.
    if shutil.which("makeblastdb"):
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                        "-out", f"{base}.nucl", "-title", f"{group}_mito_refdb"], check=True)
    else:
        print("[build_db] makeblastdb not found; skipped BLAST DB (run it later on the fasta)",
              file=sys.stderr)

    fams = sorted({m[3] for m in manifest if m[3]})
    print(f"[build_db] {group}: kept {len(kept)} records, dropped {n_drop} (incomplete).")
    print(f"[build_db] {group}: families ({len(fams)}): {', '.join(fams)}")
    print(f"[build_db] {group}: label db: {n_label} gene sequences")
    print(f"[build_db] {group}: wrote {gb.name}, {fa.name}, {lab.name}, {mf.name} to {out_dir}/")
    return len(kept)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--group", choices=sorted(GROUPS),
                     help="Taxon group to build (see GROUPS in this script).")
    grp.add_argument("--all", action="store_true", help="Build every group in turn.")
    ap.add_argument("--out-dir", type=Path,
                    help="Output directory (default assets/refdb/<group>/). "
                         "Not allowed with --all, which writes each group to its own default.")
    ap.add_argument("--email", help="NCBI Entrez email (required unless --no-download)")
    ap.add_argument("--api-key", help="NCBI API key (optional; raises rate limit)")
    ap.add_argument("--taxon", metavar="EXPR",
                    help="Override the group's Entrez organism expression, e.g. "
                         "'txid6681[Organism:exp]'")
    ap.add_argument("--min-cds", type=int,
                    help="Override the group's minimum protein-coding gene count.")
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
    if args.all and args.out_dir:
        ap.error("--out-dir cannot be combined with --all")
    if args.all and (args.taxon or args.min_cds is not None):
        ap.error("--taxon/--min-cds are per-group overrides; use them with --group")

    groups = sorted(GROUPS) if args.all else [args.group]
    totals = {g: build_group(g, args) for g in groups}
    if len(totals) > 1:
        print("[build_db] summary: " + ", ".join(f"{g}={n}" for g, n in totals.items()))


if __name__ == "__main__":
    main()
