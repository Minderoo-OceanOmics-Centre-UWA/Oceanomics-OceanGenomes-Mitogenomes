#!/usr/bin/env python3
"""
A NON-VALIDATING view of what the BLAST/LCA evidence actually says.

Why this exists
---------------
In NOVA_260909_LA, 65 of 68 assemblies were held and 47 of those had a clean
annotation and no stated reason beyond `species=unknown`. All 68 validation
records read, identically:

    "action": "skip", "reason": "not validated (no Found_in_blast_YN=Yes)"

That single sentence was covering three completely different situations:

  * the specimen was never identified to species, so the species-level validator
    was inapplicable (2 of 80 samples carry a species-level claim at all);
  * the assembly matched something at high identity that CONTRADICTS the label
    (OG2617, nominal Zoanthus, hits three Parazoanthus species at 98.9% -- a
    different family, i.e. a specimen labelling problem, not an assembly problem);
  * nothing close enough exists in the reference databases (60 of 68 assemblies
    produced no LCA row at all).

Those need three different actions from the person holding the jar. This file
tells them apart.

What this is NOT
----------------
It is not a second validator, and it has no write path to anything
authoritative. It reads the existing blast_combined / lca_combined outputs and
writes one TSV. It never sets Found_in_blast_YN, never touches lca_validation,
never feeds the proceed_qc validation term and never reaches submission metadata.
The >=98% identity BLAST/LCA route remains the sole authority for validating a
nominal identification, unchanged.

It also surfaces evidence the run currently discards outright. OG3060
(Scleralcyonacea sp.) has three hits at 98.6 / 98.3 / 98.0% to Chrysogorgia
tricaulis, Ideogorgia capensis and Heliopora coerulea -- all octocorals,
consistent with its nominal order -- which fail the query-span threshold
(354 nt aligned of a 1788 nt query) and vanish without trace. Sub-threshold hits
are reported here WITH their span, explicitly flagged non-qualifying.
"""

import argparse
import csv
import sys
from pathlib import Path

from nominal_rank import parse_rank

# Columns of blast_combined.*.tsv that this reads. The file is headerless.
BLAST_SUBJECT_SPECIES = 4
BLAST_PIDENT = 6
BLAST_ALN_LEN = 7
BLAST_QUERY_LEN = 8

# The pipeline's own qualifying thresholds, mirrored here for REPORTING only.
# Changing them changes what this file calls "qualifying"; it cannot change what
# the validator accepts, which is deliberate.
QUALIFYING_PIDENT = 98.0
QUALIFYING_QUERY_COVER = 0.60

OUTCOMES = (
    "nominal_id_validated",
    "nominal_rank_consistent",
    "high_identity_evidence_discordant",
    "no_species_level_claim",
    "no_qualifying_high_identity_hit",
    "no_annotated_validation_marker",
    "taxonomy_unresolved",
)

FIELDS = [
    "assembly", "og_id", "nominal_label", "nominal_rank", "nominal_taxon",
    "validation_outcome", "nominal_id_validated",
    "deepest_supported_rank", "deepest_supported_taxon",
    "markers_with_qualifying_hits", "best_hit_species", "best_hit_pident",
    "best_hit_query_cover", "best_hit_qualifying", "agrees_with_nominal_at_rank",
    "sub_threshold_hits",
]

RANK_COLUMNS = ("phylum", "class", "order", "family", "genus")


def read_blast(path):
    """Best hit per subject species: (species, pident, query_cover)."""
    hits = []
    if not path or not Path(path).is_file():
        return hits
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= BLAST_QUERY_LEN:
                continue
            try:
                pident = float(f[BLAST_PIDENT])
                aln = float(f[BLAST_ALN_LEN])
                qlen = float(f[BLAST_QUERY_LEN])
            except (ValueError, IndexError):
                continue
            cover = aln / qlen if qlen else 0.0
            hits.append((f[BLAST_SUBJECT_SPECIES].strip(), pident, cover))
    hits.sort(key=lambda h: (-h[1], -h[2]))
    return hits


def read_lca(path):
    """Rows of lca_combined.*.tsv as dicts (may be empty: 60/68 in this batch)."""
    if not path or not Path(path).is_file():
        return []
    with open(path) as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r]


def deepest_supported(lca_rows):
    """Deepest rank the LCA rows agree on, as (rank, taxon)."""
    for rank in reversed(RANK_COLUMNS):
        values = {(r.get(rank) or "").strip() for r in lca_rows}
        values.discard("")
        values.discard("dropped")
        if len(values) == 1:
            return rank, values.pop()
    return "", ""


def agrees_at(nominal_rank, nominal_taxon, lca_rows):
    """The rank at which the supported taxon agrees with the nominal label, or ''."""
    if not nominal_taxon:
        return ""
    want = nominal_taxon.split()[0].lower()
    for rank in RANK_COLUMNS:
        for row in lca_rows:
            value = (row.get(rank) or "").strip().lower()
            if value and value == want:
                return rank
    # The LCA rank columns are often 'dropped' for a divergent invertebrate, so
    # fall back to the species list. Report the rank the NOMINAL LABEL claims,
    # not "species": finding 'Zoanthus' inside 'Zoanthus sansibaricus' confirms
    # the genus, and calling that a species-level agreement would overstate it.
    want_full = nominal_taxon.strip().lower()
    for row in lca_rows:
        listed = (row.get("species_in_LCA") or "").lower()
        if want_full and want_full in listed:
            return "species" if nominal_rank in ("species", "tentative_species") else nominal_rank
        if want in [s.strip().split()[0] for s in listed.split(",") if s.strip()]:
            return nominal_rank if nominal_rank != "species" else "genus"
    return ""


def classify(nominal_label, lineage, hits, lca_rows, validated):
    rank, taxon = parse_rank(nominal_label, lineage)
    qualifying = [h for h in hits
                  if h[1] >= QUALIFYING_PIDENT and h[2] >= QUALIFYING_QUERY_COVER]
    agree_rank = agrees_at(rank, taxon, lca_rows)

    # Rank is tested BEFORE hit availability, because it is the structural
    # blocker: a label above species can never satisfy the species validator no
    # matter what the databases contain, and that is a different instruction
    # ("identify the specimen further") from "no close reference exists".
    if validated:
        outcome = "nominal_id_validated"
    elif rank == "unresolved":
        outcome = "taxonomy_unresolved"
    elif rank not in ("species", "tentative_species"):
        # The label claims a rank above species. Where the evidence is consistent
        # with that claim, say so -- reported, never validating: the >=98% species
        # route remains the only thing that validates a nominal ID.
        outcome = "nominal_rank_consistent" if agree_rank else "no_species_level_claim"
    elif not hits:
        outcome = "no_annotated_validation_marker"
    elif qualifying and not agree_rank:
        outcome = "high_identity_evidence_discordant"
    else:
        outcome = "no_qualifying_high_identity_hit"
    return rank, taxon, outcome, qualifying, agree_rank


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, help="mt_assembly_prefix")
    ap.add_argument("--og-id", required=True)
    ap.add_argument("--nominal", default="", help="nominal_species_id as recorded")
    ap.add_argument("--blast", default=None, help="blast_combined.<prefix>.tsv")
    ap.add_argument("--lca", default=None, help="lca_combined.<prefix>.tsv")
    ap.add_argument("--class", dest="class_name", default="")
    ap.add_argument("--family", default="")
    ap.add_argument("--order", default="")
    ap.add_argument("--validated", action="store_true",
                    help="the authoritative validator set Found_in_blast_YN=Yes. "
                         "Read only so this report can agree with it; this script "
                         "never sets it.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lineage = {"class": args.class_name, "family": args.family, "order": args.order}
    hits = read_blast(args.blast)
    lca_rows = read_lca(args.lca)
    rank, taxon, outcome, qualifying, agree_rank = classify(
        args.nominal, lineage, hits, lca_rows, args.validated)

    d_rank, d_taxon = deepest_supported(lca_rows)
    best = hits[0] if hits else ("", "", "")
    sub = [f"{h[0]}:{h[1]:.1f}%/{h[2]:.2f}cov" for h in hits[:5]
           if not (h[1] >= QUALIFYING_PIDENT and h[2] >= QUALIFYING_QUERY_COVER)]

    row = {
        "assembly": args.assembly,
        "og_id": args.og_id,
        "nominal_label": args.nominal,
        "nominal_rank": rank,
        "nominal_taxon": taxon,
        "validation_outcome": outcome,
        # Explicit, so no consumer has to infer it from the outcome string.
        "nominal_id_validated": "true" if args.validated else "false",
        "deepest_supported_rank": d_rank or "none",
        "deepest_supported_taxon": d_taxon or "none",
        "markers_with_qualifying_hits": str(len(qualifying)),
        "best_hit_species": best[0] or "none",
        "best_hit_pident": f"{best[1]:.3f}" if best[0] else "none",
        "best_hit_query_cover": f"{best[2]:.3f}" if best[0] else "none",
        "best_hit_qualifying": "true" if qualifying else "false",
        "agrees_with_nominal_at_rank": agree_rank or "none",
        # Evidence the validator correctly rejects but a curator still wants:
        # OG3060's 98.6% octocoral hits at 20% query cover, for instance.
        "sub_threshold_hits": ";".join(sub) if sub else "none",
    }
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerow(row)
    print(f"[taxonomic_evidence] {args.assembly}: {outcome} "
          f"(nominal_rank={rank}, nominal_id_validated={row['nominal_id_validated']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
