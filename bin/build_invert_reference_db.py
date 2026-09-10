#!/usr/bin/env python3
"""Build curated per-phylum invertebrate mitogenome reference databases.

Each group's database is used twice by the pipeline, and both uses are label-free.
Both go through SELECT_REFERENCE_DB, which BLASTs the sample's own assembly
against <group>_mito_refdb.fasta -- the second of the two narrowing stages that
resolve a seed (stage 1 being InvertTaxonGroups.seedDbGroup: class -> group):

  * GetOrganelle reseed -- the top-n matching records become the seed (-s) and
    their genes the label database (--genes) for an invertebrate whose first-pass
    assembly failed. Handing GetOrganelle the WHOLE group instead recruits reads
    from across the phylum: it took a coral reseed from 2 scaffolds to 12, and
    mollusca (850 genomes) and arthropoda (647) are several times worse than the
    anthozoa (221) that broke.
  * Annotation reference -- the single best record becomes the reference for
    CORAL_ANNOTATION_FIX and the QC checks, instead of one resolved from the
    sample's species *label* (MITOHIFI_FINDMITOREFERENCE -> NCBI by name), which
    silently yields a wrong-family reference when the label is wrong or coarse.

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

Outputs (into --out-dir, default assets/refdb/<group>/). The four TRACKED files
are what the pipeline reads; the .gb is a build intermediate, gitignored for
every group because the eight of them are 84 MB raw / 26 MB compressed and each
rebuild rewrites all of it into history. bin/refdb_record.py rebuilds an
equivalent record from the tracked four, so nothing needs the .gb at runtime.
    <group>_mito_refdb.fasta       all genome sequences (BLAST subject + GetOrganelle seed -s)
    <group>_mito_refdb.label.fasta GetOrganelle label/gene database (--genes), headers
                                   '>gene type - accession--Organism'
    <group>_mito_refdb.manifest.tsv  group, accession, organism, family, length, n_cds,
                                   n_rrna, source, lineage
    <group>_mito_refdb.features.tsv  accession, type, gene, product, exon coordinates
    <group>_mito_refdb.gb          combined GenBank (build intermediate; gitignored)
    <group>_mito_refdb.nucl.*      BLAST nucleotide DB (if makeblastdb is available)

Runs in the MITOS2 BioContainer (biopython + blast). NCBI access needs --email.

Usage:
    build_invert_reference_db.py --group mollusca --email you@uwa.edu.au
    build_invert_reference_db.py --all --email you@uwa.edu.au
    build_invert_reference_db.py --group anthozoa --email you@uwa.edu.au \
        --extra-gb my_curated_gb/        # private build with in-house records
    build_invert_reference_db.py --group anthozoa --no-download \
        --extra-gb my_curated_gb/        # offline; in-house records only
    build_invert_reference_db.py --all --refresh-derived
        # regenerate the tracked files from each group's existing .gb, no NCBI --
        # a schema refresh that cannot change which records the database holds
"""
import argparse
import collections
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO

SEARCH_TEMPLATE = ('({organism}) {refseq}AND '
                   'mitochondrion[filter] AND "complete genome"[Title]')
REFSEQ_TERM = 'AND refseq[filter] '

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
class GroupSpec(NamedTuple):
    organism: str
    min_cds: int
    require_rrna: bool
    require_nad5: bool
    # RefSeq is a curated subset, not a completeness bar -- the records it omits are
    # ordinary INSDC submissions that pass the same record_is_complete() check, and
    # restricting to it hides most of a phylum wherever it is the binding constraint.
    # It is what held ctenophora to 4 records out of 35, excluding every Platyctenida
    # genome -- including two Tjalfiella mitogenomes, the genus of a panel sample that
    # was written off as having no reference at all. It stays per group rather than
    # becoming a global constant because a group can legitimately want the curated
    # subset (one record per genome, no dedup needed), and because flipping one is
    # then a reviewable one-line change with its own rebuild.
    refseq_only: bool = True
    # Records to keep per organism (see cap_per_organism). Only bites once refseq_only
    # is off, which is when a heavily-resequenced species can otherwise supply most of
    # a seed panel. Kept per group rather than left to the CLI so a rebuild is
    # reproducible from --group alone, with no flag to remember: getting this wrong
    # silently changes what the database holds. 1 for the widened groups, because a
    # second isolate of the same species adds no taxonomic coverage, only bulk --
    # and at INSDC scale the tracked assets are plain git blobs rewritten wholesale on
    # every rebuild. 2 for ctenophora, which is what its shipped build used.
    max_per_organism: int = 2


# Every group searches all of INSDC (refseq_only=False) and keeps one record per
# organism. RefSeq-only was the original setting and was lifted group by group once
# it was established that no submitted mitogenome was ever built from these
# databases: assets/refdb/, select_reference_db.py and the origin anchor table all
# postdate v2.0.0, which is what produced the corals now in ENA. See
# assets/refdb/README.md for the per-group before/after counts.
GROUPS = {
    # group:          organism expression                    min_cds  rrna   nad5
    "anthozoa":       GroupSpec("txid6101[Organism:exp]",    13,      True,  True,
                                refseq_only=False, max_per_organism=1),
    "porifera":       GroupSpec("txid6040[Organism:exp]",    13,      True,  True,
                                refseq_only=False, max_per_organism=1),
    "mollusca":       GroupSpec("txid6447[Organism:exp]",    12,      True,  False,
                                refseq_only=False, max_per_organism=1),
    "arthropoda":     GroupSpec("txid6657[Organism:exp] "
                                "NOT txid6960[Organism:exp] "   # Hexapoda (insects, springtails)
                                "NOT txid6854[Organism:exp] "   # Arachnida
                                "NOT txid61985[Organism:exp]",  # Myriapoda
                                                          13,      True,  False,
                                refseq_only=False, max_per_organism=1),
    "echinodermata":  GroupSpec("txid7586[Organism:exp]",    13,      True,  False,
                                refseq_only=False, max_per_organism=1),
    # The group the restriction was lifted for first, and the only one that keeps a
    # cap of 2: at 16 records it is small enough that a second isolate of a species
    # is worth more than the bulk it costs, and 2 is what its shipped build used.
    "ctenophora":     GroupSpec("txid10197[Organism:exp]",   10,      False, False,
                                refseq_only=False),
    "tunicata":       GroupSpec("txid7712[Organism:exp]",    12,      True,  False,
                                refseq_only=False, max_per_organism=1),
    "annelida":       GroupSpec("txid6340[Organism:exp]",    12,      True,  False,
                                refseq_only=False, max_per_organism=1),
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


# A RefSeq accession carries an underscore (NC_045305.1); an INSDC one never does
# (MG655622.1). Used to prefer the curated copy of a duplicated genome.
def is_refseq(rec):
    return "_" in str(rec.id).split(".")[0]


def source_of(rec):
    return "refseq" if is_refseq(rec) else "genbank"


# RefSeq records name the INSDC record they were derived from, in the COMMENT:
#   "PROVISIONAL REFSEQ: ... The reference sequence is identical to MG655622.1."
_REFSEQ_TWIN = re.compile(r"reference sequence is identical to\s+([A-Z0-9_.]+)", re.I)


def refseq_twin(rec):
    """The INSDC accession a RefSeq record duplicates, without version, or ''."""
    if not is_refseq(rec):
        return ""
    comment = " ".join(str(rec.annotations.get("comment", "")).split())
    m = _REFSEQ_TWIN.search(comment)
    return m.group(1).split(".")[0] if m else ""


def drop_insdc_twins(entries):
    """Drop INSDC records that a RefSeq record in the same set already duplicates.

    Only reachable once refseq_only is off: with the filter on, a group holds NC_
    records exclusively and there is nothing to pair. Without this, ctenophora keeps
    both NC_038065 and its source MG655622 (Beroe forskalii), and NC_045864 next to
    MN544300/MN544301 (Hormiphora californensis) -- near-identical genomes that
    would fill a top-n seed panel with copies of one organism.
    """
    twins = {t for rec, *_ in entries if (t := refseq_twin(rec))}
    if not twins:
        return entries, 0
    kept = [e for e in entries if is_refseq(e[0]) or e[0].id.split(".")[0] not in twins]
    return kept, len(entries) - len(kept)


def cap_per_organism(entries, cap):
    """Keep at most `cap` records per organism, best first.

    The widened ctenophore search returns nine Vallicula multiformis isolates
    (PX922690-PX922698). All nine sit in one family, so neither a top-n seed panel
    nor select_fallback_seed's family balancing can dilute them: the cap has to
    happen here. Ranked RefSeq first, then longest, then most CDS -- all three are
    deterministic, so a rebuild is reproducible.
    """
    if cap <= 0:
        return entries, 0
    by_org = collections.defaultdict(list)
    for e in entries:
        by_org[(e[0].annotations.get("organism", "") or "").strip().lower()].append(e)
    kept = []
    for org in sorted(by_org):
        ranked = sorted(by_org[org],
                        key=lambda e: (is_refseq(e[0]), len(e[0].seq), e[2]),
                        reverse=True)
        kept.extend(ranked[:cap])
    kept.sort(key=lambda e: str(e[0].id))
    return kept, len(entries) - len(kept)


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


def lineage_of(rec):
    """Full NCBI taxonomy lineage as a ';'-joined string (manifest column)."""
    return ";".join(t.strip() for t in (rec.annotations.get("taxonomy", []) or []) if t.strip())


def write_features_tsv(records, path):
    """Write the feature coordinates the pipeline reads back out of a reference.

    This is what replaces the combined <group>_mito_refdb.gb as a tracked file:
    84 MB raw / 26 MB compressed across the eight groups becomes ~5 MB / ~1 MB,
    in a text table that diffs between rebuilds instead of an opaque blob. See
    bin/refdb_record.py for the reader that turns it back into a SeqRecord.

    The label database cannot serve this purpose: it stores each feature's
    SPLICED sequence (feat.extract), whereas coral_fix_bed.ref_features needs
    one entry per exon of the group-I-intron-split nad5. Coordinates keep the
    exon structure; sequence comes from the .fasta.
    """
    from refdb_record import FEATURE_TYPES, FEATURES_COLUMNS, format_parts
    n = 0
    with open(path, "w") as fh:
        fh.write("\t".join(FEATURES_COLUMNS) + "\n")
        for rec in records:
            for feat in rec.features:
                if feat.type not in FEATURE_TYPES:
                    continue
                gene = (feat.qualifiers.get("gene") or [""])[0]
                product = (feat.qualifiers.get("product") or [""])[0]
                if not gene and not product:
                    continue
                try:
                    parts = format_parts(feat.location)
                except Exception:
                    continue
                fh.write("\t".join([rec.id, feat.type,
                                    gene.replace("\t", " "),
                                    product.replace("\t", " "),
                                    parts]) + "\n")
                n += 1
    return n


MANIFEST_COLUMNS = ["group", "accession", "organism", "family", "length_bp",
                    "n_cds", "n_rrna", "source", "lineage"]


def write_manifest(manifest, path):
    """Provenance for the kept records, sorted by family then organism.

    'lineage' is the full NCBI taxonomy; reference_divergence_check.py grades a
    reference on genus/family/order and could only get family from the old
    columns. Readers must key on the header, not on position: anthozoa's
    manifest predates the per-group builder and has no 'group' column.
    """
    with open(path, "w") as fh:
        fh.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for row in sorted(manifest, key=lambda r: (r[3], r[2])):
            fh.write("\t".join(map(str, row)) + "\n")


def refresh_derived(group, out_dir):
    """Regenerate the derived artifacts from an existing local .gb, no NCBI.

    The .gb is a build intermediate, so the fasta/label/manifest/features
    quartet can always be rebuilt from it without re-querying NCBI -- which
    matters because a fresh download would also change WHICH records are in the
    database (anthozoa currently yields 278 rather than the shipped 221), and a
    schema refresh must not smuggle in a content change.
    """
    base = out_dir / f"{group}_mito_refdb"
    gb = base.with_suffix(".gb")
    if not gb.exists():
        sys.exit(f"[build_db] {group}: {gb} not found; --refresh-derived needs the "
                 f"build intermediate. Rebuild the group from NCBI instead.")
    kept = list(SeqIO.parse(str(gb), "genbank"))
    if not kept:
        sys.exit(f"[build_db] {group}: no records in {gb}")

    spec = GROUPS[group]
    manifest = []
    for rec in kept:
        _ok, n_cds, n_rrna, _reason = record_is_complete(
            rec, spec.min_cds, spec.require_rrna, spec.require_nad5)
        # source per record, not a hardcoded "refseq": a group built with
        # refseq_only=False holds both kinds and the manifest has to say which.
        manifest.append((group, rec.id, rec.annotations.get("organism", ""), family_of(rec),
                         len(rec.seq), n_cds, n_rrna, source_of(rec), lineage_of(rec)))

    SeqIO.write(kept, str(base.with_suffix(".fasta")), "fasta")
    n_label = write_label_db(kept, Path(f"{base}.label.fasta"))
    n_feat = write_features_tsv(kept, Path(f"{base}.features.tsv"))
    write_manifest(manifest, Path(f"{base}.manifest.tsv"))
    print(f"[build_db] {group}: refreshed derived artifacts from {gb.name} "
          f"({len(kept)} records, {n_label} label seqs, {n_feat} feature rows)")
    return len(kept)


def fetch_records(organism, email, api_key, batch, retmax, refseq_only=True):
    """esearch + batched efetch of GenBank records from NCBI nucleotide."""
    from Bio import Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    query = SEARCH_TEMPLATE.format(organism=organism,
                                   refseq=REFSEQ_TERM if refseq_only else "")
    print(f"[build_db] esearch: {query}", file=sys.stderr)
    with Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y") as h:
        res = Entrez.read(h)
    ids = res["IdList"]
    print(f"[build_db] {len(ids)} accessions found "
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
    spec = GROUPS[group]
    organism = args.taxon or spec.organism
    min_cds = args.min_cds if args.min_cds is not None else spec.min_cds
    refseq_only = spec.refseq_only if args.refseq_only is None else args.refseq_only
    max_per_organism = (spec.max_per_organism if args.max_per_organism is None
                        else args.max_per_organism)
    out_dir = args.out_dir or default_out_dir(group)

    print(f"[build_db] === {group}: {organism}, min_cds={min_cds}, "
          f"require_rrna={spec.require_rrna}, require_nad5={spec.require_nad5}, "
          f"refseq_only={refseq_only}, max_per_organism={max_per_organism} ===",
          file=sys.stderr)

    candidates = []  # (rec, source)
    if not args.no_download:
        candidates += [(r, source_of(r)) for r in
                       fetch_records(organism, args.email, args.api_key, args.batch,
                                     args.retmax, refseq_only=refseq_only)]
    if args.extra_gb:
        candidates += load_extra(args.extra_gb)

    out_dir.mkdir(parents=True, exist_ok=True)
    entries, seen = [], set()   # (rec, source, n_cds, n_rrna)
    n_drop = 0
    for rec, source in candidates:
        acc = rec.id.split(".")[0]
        if acc in seen:
            continue
        ok, n_cds, n_rrna, reason = record_is_complete(
            rec, min_cds, spec.require_rrna, spec.require_nad5)
        if not ok:
            n_drop += 1
            continue
        seen.add(acc)
        entries.append((rec, source, n_cds, n_rrna))

    # Both stages are no-ops for a refseq_only group (one curated record per genome,
    # rarely more than one isolate per organism) and only bite once the filter is
    # lifted, which is when a single organism can otherwise supply most of the panel.
    entries, n_twin = drop_insdc_twins(entries)
    entries, n_dup = cap_per_organism(entries, max_per_organism)

    kept = [e[0] for e in entries]
    manifest = [(group, rec.id, rec.annotations.get("organism", ""), family_of(rec),
                 len(rec.seq), n_cds, n_rrna, source, lineage_of(rec))
                for rec, source, n_cds, n_rrna in entries]

    if not kept:
        sys.exit(f"[build_db] {group}: no records passed the completeness filter; nothing written")

    base = out_dir / f"{group}_mito_refdb"
    gb, fa = base.with_suffix(".gb"), base.with_suffix(".fasta")
    lab = Path(f"{base}.label.fasta")
    mf = Path(f"{base}.manifest.tsv")
    ft = Path(f"{base}.features.tsv")
    # The .gb is a build-time intermediate, not a tracked artifact: it is
    # written so a build can be inspected and so --refresh-derived can re-read
    # it, but .gitignore excludes it for every group. What the pipeline reads is
    # the fasta + label + manifest + features quartet below.
    SeqIO.write(kept, str(gb), "genbank")
    SeqIO.write(kept, str(fa), "fasta")
    n_label = write_label_db(kept, lab)
    n_feat = write_features_tsv(kept, ft)
    write_manifest(manifest, mf)

    # BLAST nucleotide DB for the sequence-based reference selection step.
    if shutil.which("makeblastdb"):
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                        "-out", f"{base}.nucl", "-title", f"{group}_mito_refdb"], check=True)
    else:
        print("[build_db] makeblastdb not found; skipped BLAST DB (run it later on the fasta)",
              file=sys.stderr)

    fams = sorted({m[3] for m in manifest if m[3]})
    print(f"[build_db] {group}: kept {len(kept)} records, dropped {n_drop} (incomplete), "
          f"{n_twin} (INSDC twin of a RefSeq record), "
          f"{n_dup} (over the {max_per_organism}-per-organism cap).")
    print(f"[build_db] {group}: families ({len(fams)}): {', '.join(fams)}")
    print(f"[build_db] {group}: label db: {n_label} gene sequences")
    print(f"[build_db] {group}: features: {n_feat} rows")
    print(f"[build_db] {group}: wrote {fa.name}, {lab.name}, {mf.name}, {ft.name} to {out_dir}/ "
          f"({gb.name} is a build intermediate and is gitignored)")
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
    ap.add_argument("--refresh-derived", action="store_true",
                    help="Regenerate fasta/label/manifest/features from the group's "
                         "existing local .gb without contacting NCBI. Use after a schema "
                         "change, so the artifacts are refreshed without also changing "
                         "which records the database contains.")
    ap.add_argument("--refseq-only", dest="refseq_only", action="store_true", default=None,
                    help="Restrict the NCBI search to RefSeq, overriding the group's "
                         "refseq_only setting.")
    ap.add_argument("--no-refseq-filter", dest="refseq_only", action="store_false",
                    help="Search all of INSDC, not just RefSeq. Widens every group "
                         "2.5-8.8x and changes which reference is picked for existing "
                         "samples, so prefer setting refseq_only per group.")
    ap.add_argument("--max-per-organism", type=int, default=None, metavar="N",
                    help="Keep at most N records per organism (0 = no cap), overriding "
                         "the group's max_per_organism. Stops one heavily-resequenced "
                         "species filling the seed panel.")
    # An unfiltered group can exceed the old 2000: mollusca matches ~3169 and
    # arthropoda ~2391. fetch_records hard-exits rather than shipping the first
    # retmax hits, so too low a value aborts the build instead of corrupting it.
    ap.add_argument("--retmax", type=int, default=6000)
    ap.add_argument("--batch", type=int, default=200)
    args = ap.parse_args()

    if not args.no_download and not args.email and not args.refresh_derived:
        ap.error("--email is required for NCBI download (or pass --no-download)")
    if args.all and args.out_dir:
        ap.error("--out-dir cannot be combined with --all")
    if args.all and (args.taxon or args.min_cds is not None or args.refseq_only is not None
                     or args.max_per_organism is not None):
        ap.error("--taxon/--min-cds/--refseq-only/--max-per-organism are per-group "
                 "overrides; use them with --group")

    groups = sorted(GROUPS) if args.all else [args.group]
    if args.refresh_derived:
        totals = {g: refresh_derived(g, args.out_dir or default_out_dir(g)) for g in groups}
    else:
        totals = {g: build_group(g, args) for g in groups}
    if len(totals) > 1:
        print("[build_db] summary: " + ", ".join(f"{g}={n}" for g, n in totals.items()))

    # The published-origin anchor table is derived from these databases, so any change
    # to their contents makes it stale. Regenerating it needs the taxdump (the ORDER
    # level cannot be resolved without it), which this script does not otherwise take,
    # so print the exact command rather than guessing a path. Getting this wrong is not
    # cosmetic: a table without the order level resolves Scleractinia from the Anthozoa
    # class aggregate and rotates every submitted stony coral off tRNA-Met.
    if not args.refresh_derived:
        print("[build_db] NOTE: the databases changed, so "
              "assets/taxonomy/mito_origin_anchors.json is now stale. Regenerate it with:")
        print("[build_db]   bin/build_origin_anchor_table.py --taxdump-dir <taxdump>")
        print("[build_db] and check it in the same commit "
              "(tests/unit/test_origin_anchor_table.py enforces this).")


if __name__ == "__main__":
    main()
