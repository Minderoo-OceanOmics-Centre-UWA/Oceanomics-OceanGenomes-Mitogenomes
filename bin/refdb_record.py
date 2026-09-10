#!/usr/bin/env python3
"""Materialize a reference GenBank record from the tracked reference-DB artifacts.

The per-group databases under assets/refdb/<group>/ used to ship a combined
<group>_mito_refdb.gb -- 84 MB raw / 26 MB compressed across the eight groups,
rewritten wholesale into git history on every rebuild. Only anthozoa's was
tracked, which is what confined the sequence-based reference selection
(SELECT_REFERENCE_DB) to corals.

Everything the pipeline actually reads out of a reference GenBank is available
from three much smaller tracked files, so the .gb is now a build-time
intermediate and a record is rebuilt on demand instead:

    <group>_mito_refdb.fasta         ORIGIN sequence + length
    <group>_mito_refdb.manifest.tsv  organism, family, full taxonomy lineage
    <group>_mito_refdb.features.tsv  feature type, gene/product, exon coordinates

The rebuilt record is written as GenBank, so every downstream consumer keeps its
current interface (--ref-gb / --reference-gb) and none of them needed changing:
coral_fix_bed.py (rrnL / nad5 exons / cox1), reference_divergence_check.py
(organism + lineage), reference_relevance_check.py, check_getorganelle.py and
join_scaffolds_by_reference.py (sequence + length). That also keeps them working
unchanged for vertebrates, where the reference really is a GenBank file
downloaded by findMitoReference.

Exon structure is the reason a coordinates table is required rather than reusing
<group>_mito_refdb.label.fasta: the label DB stores each feature's SPLICED
sequence (feat.extract), while coral_fix_bed.ref_features needs one entry per
exon of the group-I-intron-split nad5.

Used as a library by select_reference_db.py and select_fallback_seed.py, which
also share the coverage and label-subsetting helpers below; the CLI is for
spot-checking and for the round-trip test that gates dropping the .gb files.

Usage:
    refdb_record.py --refdb-dir assets/refdb/anthozoa --group anthozoa \
        --accession NC_027611.1 --out-gb ref.gb
"""
import argparse
import csv
import sys
from pathlib import Path

FEATURES_COLUMNS = ["accession", "type", "gene", "product", "parts"]

# Feature types worth carrying. Everything the consumers look at is a gene
# feature; source/misc_feature rows would triple the table for no reader.
FEATURE_TYPES = ("CDS", "rRNA", "tRNA")


def refdb_paths(refdb_dir, group):
    """The three tracked artifacts for a group, as (fasta, manifest, features)."""
    base = Path(refdb_dir) / f"{group}_mito_refdb"
    return (Path(f"{base}.fasta"), Path(f"{base}.manifest.tsv"),
            Path(f"{base}.features.tsv"))


def format_parts(location):
    """Serialise a feature location as 'start-end:strand[,start-end:strand...]'.

    Coordinates are Biopython's own 0-based half-open ints and parts stay in
    transcript order, so parse_parts round-trips them exactly. Strand is stored
    per part rather than per feature: it is almost always uniform, but a
    trans-spliced feature is not, and per-part costs nothing.
    """
    out = []
    for part in location.parts:
        strand = 0 if part.strand is None else int(part.strand)
        out.append(f"{int(part.start)}-{int(part.end)}:{strand}")
    return ",".join(out)


def parse_parts(spec):
    """Inverse of format_parts. Returns a Simple/CompoundLocation."""
    from Bio.SeqFeature import CompoundLocation
    try:
        from Bio.SeqFeature import SimpleLocation as Location
    except ImportError:  # biopython < 1.81
        from Bio.SeqFeature import FeatureLocation as Location

    locs = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        span, _, strand = chunk.partition(":")
        start, _, end = span.partition("-")
        strand = int(strand) if strand else 0
        locs.append(Location(int(start), int(end), strand=(strand or None)))
    if not locs:
        raise ValueError(f"no parsable location in {spec!r}")
    return locs[0] if len(locs) == 1 else CompoundLocation(locs)


# Shared by the two reference selectors -- select_reference_db.py, which ranks a
# group database against a sample's own first-pass assembly, and
# select_fallback_seed.py, which ranks a taxonomy-bounded shortlist against its
# reads when there is no assembly to rank against. They live here rather than in
# select_reference_db.py so that the fallback selector's pure-logic half stays
# importable without Biopython, the same reason Bio is imported lazily below.


def union_coverage(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered, clo, chi = 0, *intervals[0]
    for lo, hi in intervals[1:]:
        if lo > chi + 1:
            covered += chi - clo + 1
            clo, chi = lo, hi
        else:
            chi = max(chi, hi)
    return covered + (chi - clo + 1)


def label_accession(header):
    """Accession out of a label-DB header: '>GENE type - ACCESSION--Organism'.

    Matched on the accession alone, never the gene name: the gene vocabulary
    differs by group (COX1 / cox2 / cytochrome_c_oxidase_subunit_I / URF1), the
    accession does not.
    """
    _, _, tail = header.partition(" - ")
    acc, _, _ = tail.partition("--")
    return acc.strip()


def subset_label_db(label_path, accessions, out_path):
    """Copy the label records belonging to `accessions`, preserving order."""
    wanted = set(accessions)
    n = 0
    with open(label_path) as fh, open(out_path, "w") as out:
        keep = False
        for line in fh:
            if line.startswith(">"):
                keep = label_accession(line[1:].rstrip()) in wanted
                n += 1 if keep else 0
            if keep:
                out.write(line)
    return n


def read_manifest(path):
    """accession -> {organism, family, lineage[list], length_bp, ...}.

    Keyed by the header row, never by position: anthozoa's manifest predates the
    per-group builder and has no 'group' column, so the other seven are one
    column wider.
    """
    rows = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if not acc:
                continue
            lineage = [t.strip() for t in (row.get("lineage") or "").split(";") if t.strip()]
            rows[acc] = {
                "organism": (row.get("organism") or "").strip(),
                "family": (row.get("family") or "").strip(),
                "lineage": lineage,
                "length_bp": int(row["length_bp"]) if (row.get("length_bp") or "").strip() else None,
            }
    return rows


def read_features(path, accessions=None):
    """accession -> [(type, gene, product, parts_spec), ...] in file order.

    `accessions`, when given, restricts the scan so selecting one reference out
    of mollusca does not materialise all 31,718 of its feature rows.
    """
    wanted = set(accessions) if accessions is not None else None
    feats = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if not acc or (wanted is not None and acc not in wanted):
                continue
            feats.setdefault(acc, []).append((
                (row.get("type") or "").strip(),
                (row.get("gene") or "").strip(),
                (row.get("product") or "").strip(),
                (row.get("parts") or "").strip(),
            ))
    return feats


def read_sequence(fasta_path, accession):
    """The record's sequence, matched on the accession token of the FASTA header."""
    from Bio import SeqIO
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        if rec.id == accession or rec.id.split(".")[0] == accession.split(".")[0]:
            return rec.seq
    return None


def build_record(accession, seq, manifest_row, feature_rows):
    """Reassemble a SeqRecord equivalent to the source GenBank record.

    Equivalent for every field the pipeline reads -- sequence, id, organism,
    taxonomy, and gene features with their exon structure -- not a byte-for-byte
    reproduction of the original file.
    """
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(seq, id=accession, name=accession.split(".")[0],
                    description=f"{manifest_row.get('organism', '')} mitochondrion, complete genome")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    rec.annotations["organism"] = manifest_row.get("organism", "")
    rec.annotations["source"] = manifest_row.get("organism", "")
    rec.annotations["taxonomy"] = list(manifest_row.get("lineage", []))

    for ftype, gene, product, parts in feature_rows:
        if not parts:
            continue
        quals = {}
        if gene:
            quals["gene"] = [gene]
        if product:
            quals["product"] = [product]
        rec.features.append(SeqFeature(parse_parts(parts), type=ftype, qualifiers=quals))
    return rec


def materialize(refdb_dir, group, accession):
    """Rebuild one record from a group's tracked artifacts, or None if absent."""
    fasta, manifest, features = refdb_paths(refdb_dir, group)
    seq = read_sequence(fasta, accession)
    if seq is None:
        return None
    rows = read_manifest(manifest)
    feats = read_features(features, accessions=[accession])
    return build_record(accession, seq, rows.get(accession, {}), feats.get(accession, []))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdb-dir", required=True, type=Path)
    ap.add_argument("--group", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--out-gb", required=True, type=Path)
    args = ap.parse_args()

    from Bio import SeqIO
    rec = materialize(args.refdb_dir, args.group, args.accession)
    if rec is None:
        sys.exit(f"[refdb_record] {args.accession} not found in {args.refdb_dir}")
    SeqIO.write([rec], str(args.out_gb), "genbank")
    print(f"[refdb_record] wrote {args.out_gb} "
          f"({len(rec.seq)} bp, {len(rec.features)} features)", file=sys.stderr)


if __name__ == "__main__":
    main()
