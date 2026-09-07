#!/usr/bin/env python3
"""
Using the .gff, counts the length of the mitogenome and all of the genes.
Using the .gff, reports on if there are missing or extra genes and if theyre in the correct order.
Counts the length of each mitochondrial protein-coding gene
from translated FASTA (.faa/.fasta) files in a `proteins/` directory.

Usage:
    singularity run $SING/psycopg2:0.1.sif python emma_stats.py /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/annotation/*.gff /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/annotation/

"""

import os
import csv
import sys
import argparse
from pathlib import Path


# Reference gene order (tRNA, rRNA, CDS) — the standard vertebrate set. Shared
# with the rescue gates via bin/mito_gene_order.py so the gates and this QC step
# cannot drift on what "present and in order" means.
from mito_gene_order import (
    REF_GENES,
    TRNA_GENES as TRNA_NAMES,
    accepted_orders_for,
    gene_entries_by_coord,
    parse_gff_attributes,
)

# Protein-coding genes to pull from .faa/.fa files
PROT_GENES = [
    "ATP6", "ATP8", "CO1", "CO2", "CO3",
    "CYTB", "ND1", "ND2", "ND3", "ND4",
    "ND4L", "ND5", "ND6"
]

# Taxonomic classes with legitimately reduced/atypical mt tRNA complements:
# Cnidaria follows the well-known pattern of 13 PCGs + 2 rRNAs but only ~2 mt
# tRNAs (trnM/trnW), the rest nuclear-encoded and imported. Porifera groups
# with it here too -- sponge mt tRNA counts are highly variable across
# lineages (2-27 genes), including documented tRNA-Phe loss in some clades --
# so the same relaxed expectation applies. The vertebrate 22-tRNA expectation
# would wrongly fail both groups, so for them completeness is judged on the
# conserved protein-coding + rRNA core only and gene order is not evaluated.
REDUCED_TRNA_CLASSES = {
    "anthozoa", "hydrozoa", "scyphozoa", "cubozoa",
    "staurozoa", "myxozoa", "polypodiozoa", "cnidaria",
    "demospongiae", "calcarea", "hexactinellida", "homoscleromorpha", "porifera",
}
# The conserved protein-coding + rRNA core. Named for the lineage that motivated it,
# but it is the completeness set for EVERY non-vertebrate translation table now (see
# completeness_profile): 13 PCGs + 2 rRNAs is the shared expectation, and gene order
# is not evaluated because no non-vertebrate reference order is defined here.
CNIDARIAN_CORE = PROT_GENES + ["RNR1", "RNR2"]

# The 22 vertebrate mt tRNAs (the T* entries of REF_GENES) and the conserved
# protein-coding + rRNA core. A vertebrate mitogenome that carries the whole core
# in the right order but is short a small number of tRNAs is an annotation
# limitation (EMMA's tRNA model misses divergent copies), not an assembly defect,
# so it is allowed to pass with the shortfall recorded in trna_advisory. See
# --trna-tolerance and the matching TRNA_TOLERANCE in mitogenome_assembly_summary.py.
VERT_CORE = set(PROT_GENES) | {"RNR1", "RNR2"}            # 13 PCG + 2 rRNA
DEFAULT_TRNA_TOLERANCE = 2

# Interior intergenic gap, in bp, at or above which annotation_gaps reports a span.
#
# 50 and not 100. A single missed mitochondrial tRNA leaves a 55-75 bp hole, so a
# 100 bp threshold sits ABOVE the exact signal this diagnostic exists to catch --
# it is what distinguishes a real transposition (a tRNA moved) from a missed call
# (a tRNA-shaped hole left at both the origin and the destination of the apparent
# move). The cost of reaching that low is the OriL, which runs 30-51 bp and so
# straddles the threshold: a passing assembly reporting one ~51 bp TN->TC span is
# expected, not a regression, and is exactly what "purely advisory" is for. Do not
# raise the default to silence it.
DEFAULT_GAP_THRESHOLD = 50

# Adjacent pairs whose intergenic span is legitimately large and must never be
# reported. TP->TF is the vertebrate control region, routinely ~1 kb.
GAP_EXEMPT_PAIRS = {("TP", "TF")}


def has_reduced_trna_expectation(class_name):
    return (class_name or "").strip().lower() in REDUCED_TRNA_CLASSES


# The vertebrate profile is the only one with a defined reference gene order and a
# 22-tRNA expectation, and it applies to exactly one translation table: code 2.
VERTEBRATE_CODE = 2


def completeness_profile(genetic_code=None, class_name=""):
    """Return 'vertebrate' or 'core' -- which completeness profile to judge by.

    Keyed on the resolved mitochondrial genetic code, because that is what the
    rest of the pipeline already resolves per sample (meta.genetic_code) and it is
    the property that actually determines whether the vertebrate 37-gene profile
    applies. The class string is only a fallback for callers that have no code:
    keying on the class list alone meant a code-9 echinoderm or code-5 mollusc was
    judged against vertebrate gene order and the 22-tRNA count, and failed for
    being what it is.

    Every non-vertebrate table gets the conserved PCG+rRNA core with gene order
    reported NA, since no non-vertebrate reference order is defined here.
    """
    if genetic_code is not None:
        return "vertebrate" if int(genetic_code) == VERTEBRATE_CODE else "core"
    return "core" if has_reduced_trna_expectation(class_name) else "vertebrate"


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

def annotation_gaps(gene_entries, threshold=DEFAULT_GAP_THRESHOLD):
    """Interior intergenic spans at or above `threshold`, as a formatted list.

    Each is rendered 'ND4:11768-TS1:11982(213)'. The coordinate wrap (last gene
    back to first) is skipped, because it is not an interior gap; so is any
    consecutive pair in GAP_EXEMPT_PAIRS.

    PURELY ADVISORY. This must never influence `passed`. Its job is to make a
    held assembly explain itself -- a tRNA-sized hole sitting immediately beside
    an apparently transposed tRNA points at the tRNA call, not at the gene order,
    and that distinction is otherwise only visible by reading the GFF by hand.

    Runs for EVERY completeness profile, not only the vertebrate one. A
    core-profile assembly is judged on gene presence alone, so without this a
    coral with all 15 core genes and a kilobase of unannotated sequence between
    two of them passes with nothing recorded at all.
    """
    gaps = []
    for (left_name, _ls, left_end, _lstr), (right_name, right_start, _re, _rstr) in zip(
            gene_entries, gene_entries[1:]):
        if (left_name, right_name) in GAP_EXEMPT_PAIRS:
            continue
        size = right_start - left_end - 1
        if size >= threshold:
            gaps.append(f"{left_name}:{left_end}-{right_name}:{right_start}({size})")
    return gaps


def order_deviation(found_by_coord, ref_subset):
    """How far out of order an annotation is: genes outside the longest common
    subsequence of the observed order against the reference order.

    order_correct=no collapses "one tRNA out of place" and "half the genome is
    inverted" into the same value, which makes a held pile impossible to triage.
    This splits it: a low non-zero deviation says "look at this one, it is nearly
    right", a large one says "this assembly is broken". A GROUP of samples sharing
    a low non-zero deviation is also what a missing ORDER_VARIANTS row looks like,
    so this is the cheapest signal for spotting the next clade variant.

    PURELY ADVISORY, like annotation_gaps. difflib is stdlib, so this adds no
    dependency to the container.
    """
    import difflib
    matcher = difflib.SequenceMatcher(a=ref_subset, b=found_by_coord, autojunk=False)
    matched = sum(block.size for block in matcher.get_matching_blocks())
    return len(found_by_coord) - matched


def process_gff(gff_path, annotation_name, class_name="",
                trna_tolerance=DEFAULT_TRNA_TOLERANCE, genetic_code=None,
                taxon=None, gap_threshold=DEFAULT_GAP_THRESHOLD):
    parts = annotation_name.split(".")
    if len(parts) != 5:
        print(f"⚠️ Warning: Unexpected annotation_name format: {annotation_name}")
        og_id = tech = seq_date = code = annotation = ""
    else:
        og_id, tech, seq_date, code, annotation = parts

    # Gene entries come from the SHARED reader, which dedups a gene written as more
    # than one `gene` line at its lowest start. This used to be an inline first-seen
    # dedup, which on an origin-spanning gene put it at whichever half the file
    # listed first -- so the same GFF could be judged in-order by the rescue gates
    # (which already deduped by lowest start) and out-of-order here.
    gene_entries = gene_entries_by_coord(gff_path)
    gene_lengths = {gene: "" for gene in REF_GENES}
    for gene_name, start, end, _strand in gene_entries:
        if gene_name in gene_lengths:
            gene_lengths[gene_name] = str(abs(end - start) + 1)
    total_length = extract_total_length(gff_path)

    # Count features by the gene they belong to, not by feature line, so an
    # intron-split gene written as multiple exon lines (e.g. a coral nad5 with a
    # group I intron -> two CDS lines sharing one Parent) counts once.
    cds_genes = set()
    trna_genes = set()
    rrna_genes = set()

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

            # Group by Parent (the gene) when present, else fall back to Name.
            gene_key = attributes.get("Parent") or attributes.get("Name") or f"{start}-{end}"
            if feature_type == "cds":
                cds_genes.add(gene_key)
            elif feature_type == "trna":
                trna_genes.add(gene_key)
            elif feature_type == "rrna":
                rrna_genes.add(gene_key)

    found_by_coord = [entry[0] for entry in gene_entries]

    trna_advisory = []
    # "no" is the established nothing-to-report value in this file (missing_genes,
    # trna_advisory, extra_genes all use it), so the new columns follow suit rather
    # than overloading NULL or an empty string.
    order_variant = "no"
    deviation = "no"
    gaps = annotation_gaps(gene_entries, gap_threshold)
    profile = completeness_profile(genetic_code, class_name)
    if profile == "core":
        # Judge completeness on the conserved protein-coding + rRNA core only;
        # cnidarians and sponges legitimately lack most tRNAs, and no non-vertebrate
        # lineage follows the vertebrate gene order, so order is reported as NA rather
        # than failed.
        missing = [g for g in CNIDARIAN_CORE if g not in found_by_coord]
        extra = [g for g in found_by_coord if g not in REF_GENES]
        order_correct = "NA"
        passed = len(missing) == 0
    else:
        # Every order this taxon may legitimately show, canonical FIRST. A curated
        # variant rule ADDS an accepted order rather than replacing the canonical
        # one, so a canonical member of a rule-carrying taxon keeps passing and
        # does not acquire an order_variant value merely because its family has a
        # rule. Checking canonical first is what makes that hold.
        candidates = accepted_orders_for(taxon)
        order_ref, matched_rule = candidates[0]
        order_ok = False
        for candidate_order, rule_id in candidates:
            candidate_subset = [g for g in candidate_order if g in found_by_coord]
            if found_by_coord == candidate_subset:
                order_ref, matched_rule, order_ok = candidate_order, rule_id, True
                break

        # missing is computed against the order that matched, so the reported
        # sequence of missing genes stays consistent with the order being judged.
        # missing and extra are set-membership tests, so a variant can only change
        # the ORDER in which missing genes are listed, never which are missing.
        missing = [g for g in order_ref if g not in found_by_coord]
        extra = [g for g in found_by_coord if g not in REF_GENES]
        ref_subset = [g for g in order_ref if g in found_by_coord]
        if order_ok and matched_rule:
            order_correct = "variant"
            order_variant = matched_rule
        else:
            order_correct = "yes" if order_ok else "no"
        passed = len(missing) == 0 and order_ok
        if not order_ok:
            # Deviation against the CLOSEST accepted order, so a taxon with a rule
            # is scored against whichever ordering it is nearer to rather than
            # being penalised for the rule existing.
            deviation = str(min(
                order_deviation(found_by_coord,
                                [g for g in candidate_order if g in found_by_coord])
                for candidate_order, _rule in candidates))

        # Tolerate a small tRNA-only shortfall on an otherwise complete, correctly
        # ordered mitogenome. All 13 PCGs + both rRNAs must be present, gene order
        # correct, and every missing gene a tRNA, with at most trna_tolerance of
        # them. missing_genes still lists them; trna_advisory records which were
        # waived so the pass stays auditable and downstream can flag it.
        if (not passed and order_ok and missing
                and all(g in TRNA_NAMES for g in missing)
                and VERT_CORE.issubset(found_by_coord)
                and len(missing) <= trna_tolerance):
            trna_advisory = list(missing)  # already in REF order
            passed = True

    gff_summary = {
        "og_id": og_id,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "annotation": annotation,
        "missing_genes": ";".join(missing) if missing else "no",
        "trna_advisory": ";".join(trna_advisory) if trna_advisory else "no",
        "extra_genes": ";".join(extra) if extra else "no",
        "order_correct": order_correct,
        # Which curated rule accepted a non-canonical order, so a relaxed gate stays
        # auditable after the fact. Same precedent as trna_advisory.
        "order_variant": order_variant,
        # Advisory only: neither of these may affect `passed`.
        "order_deviation": deviation,
        "annotation_gaps": ";".join(gaps) if gaps else "no",
        "passed": "yes" if passed else "no",
        # Which profile this verdict was reached under. mitogenome_assembly_summary.py
        # reads it so the run-level report cannot judge a core-profile assembly
        # against the vertebrate 37-gene expectation this step deliberately did not
        # apply -- the two would otherwise disagree on the same annotation.
        "completeness_profile": profile,
        "total_length": total_length if total_length is not None else "NA",
        "num_cds": len(cds_genes),
        "num_trna": len(trna_genes),
        "num_rrna": len(rrna_genes)
    }
    gff_summary.update(gene_lengths)
    return gff_summary

def process_protein_lengths(prot_dir, annotation_name):
    from Bio import SeqIO  # local: keeps process_gff importable without Biopython
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

def main(gff_path, prot_dir, class_name="", trna_tolerance=DEFAULT_TRNA_TOLERANCE,
         genetic_code=None, taxon=None, gap_threshold=DEFAULT_GAP_THRESHOLD):
    if not os.path.isfile(gff_path):
        sys.exit(f"❌ GFF file not found: {gff_path}")

    annotation_name = get_annotation_name(gff_path)
    gff_summary = process_gff(gff_path, annotation_name, class_name, trna_tolerance,
                              genetic_code, taxon=taxon, gap_threshold=gap_threshold)
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
    ap = argparse.ArgumentParser(
        description="Summarise mitogenome annotation stats from a GFF + proteins dir."
    )
    ap.add_argument("gff", help="path to the annotation .gff")
    ap.add_argument("proteins", help="path to the proteins/ directory")
    ap.add_argument("--class", dest="class_name", default="",
                    help="taxonomic class (e.g. Anthozoa). Fallback profile selector "
                         "for callers with no --genetic-code: cnidarian classes are "
                         "judged on the PCG+rRNA core only, everything else gets the "
                         "vertebrate 22-tRNA profile. --genetic-code wins when both "
                         "are given.")
    ap.add_argument("--genetic-code", dest="genetic_code", type=int, default=None,
                    help="NCBI mitochondrial translation table (from meta.genetic_code). "
                         "Selects the completeness profile: code 2 is judged against the "
                         "vertebrate 37-gene set and gene order, every other table "
                         "against the conserved PCG+rRNA core with order reported NA.")
    ap.add_argument("--trna-tolerance", dest="trna_tolerance", type=int,
                    default=DEFAULT_TRNA_TOLERANCE,
                    help="max tRNAs a non-cnidarian mitogenome may be missing and "
                         "still pass, provided all 13 PCGs + 2 rRNAs are present and "
                         "gene order is correct. Tolerated tRNAs are recorded in the "
                         "trna_advisory column; missing_genes stays truthful. "
                         "0 requires a complete 37-gene annotation. Default %(default)s.")
    ap.add_argument("--family", default="",
                    help="taxonomic family (meta.family). Used ONLY to look up a "
                         "curated gene-order variant; it never affects completeness.")
    ap.add_argument("--order", dest="taxon_order", default="",
                    help="taxonomic order (meta.order). See --family.")
    ap.add_argument("--genus", default="",
                    help="genus, derived as the first whitespace token of "
                         "meta.nominal_species_id -- there is no separate genus field, "
                         "and this is how reference_divergence_check.py derives it too. "
                         "See --family.")
    ap.add_argument("--gap-threshold", dest="gap_threshold", type=int,
                    default=DEFAULT_GAP_THRESHOLD,
                    help="interior intergenic span, in bp, at or above which a gap is "
                         "reported in annotation_gaps. Purely advisory: it never "
                         "affects passed. Default %(default)s, deliberately below the "
                         "55-75 bp hole a single missed mt-tRNA leaves.")
    args = ap.parse_args()

    # Rank -> value, in the shape ref_order_for expects. Unresolved values are
    # normalised away inside the resolver, so passing '' here is safe.
    taxon = {
        "genus": args.genus,
        "family": args.family,
        "order": args.taxon_order,
        "class": args.class_name,
    }

    main(args.gff, Path(args.proteins), args.class_name, args.trna_tolerance,
         args.genetic_code, taxon=taxon, gap_threshold=args.gap_threshold)
