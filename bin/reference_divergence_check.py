#!/usr/bin/env python3
"""Flag, *before* assembly, when the mitogenome reference is not a close relative
of the sample.

MITOHIFI_FINDMITOREFERENCE resolves a reference by walking the sample's NCBI
lineage (species -> genus -> family -> ...) and downloading the first complete
mitogenome it finds. For a species with no congeneric record on NCBI -- common
for deep-sea / poorly-sampled taxa -- it silently falls back to a different genus
(or worse). MitoHiFi then recruits reads by mapping to that divergent reference
and discards reads longer than it (the NUMT size filter); the most divergent gene
blocks recruit poorly and collapse out of the assembly. The result looks "clean"
(even coverage, plausible length) but is missing protein-coding genes, with no
obvious cause in the logs.

Empirically (OG2102, Rouleina attrita): even a *same-family* reference
(Alepocephalus agassizii, Alepocephalidae) was too divergent and produced a
7/13-PCG collapse. So the review trigger here is "not congeneric": a congeneric
reference is trusted; anything coarser is flagged for review. Family / order /
class only grade the severity of the message, they do not relax the trigger.

This check is pre-assembly and needs only the reference GenBank plus the sample's
own taxon name (and, if available, its class/family/order from the samplesheet).
It compares taxonomy, never sequence, so it is independent of -- and complementary
to -- the post-assembly REFERENCE_RELEVANCE BLAST check.

Emits one line:  TIER \t <details>  where TIER is one of
CONGENERIC | CONFAMILIAL | DIFFERENT_FAMILY | NON_CONGENERIC | UNKNOWN
and always exits 0, so a missing/odd reference can never break the run.

Grading uses genus and (when supplied) family only. Class is deliberately not
used for matching: the OceanOmics species table reports the modern NCBI class
(e.g. Actinopteri) while GenBank lineages still carry the older name
(Actinopterygii), so a class-string match gives false CROSS_CLASS calls. Genus
and family names are consistent between the two sources, so the comparison stays
robust; class, if supplied, only enriches the human-readable detail.

Usage:
    reference_divergence_check.py --reference-gb ref.gb \
        --sample-species "Rouleina attrita" --out PREFIX.reference_divergence.txt \
        [--sample-class Actinopteri] [--sample-family Alepocephalidae]
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO


def first_token(name):
    """Genus = first whitespace-delimited token of a taxon name, else ''."""
    return (name or "").strip().split(" ")[0] if (name or "").strip() else ""


def parse_reference(gb_path):
    """Return (organism, genus, family, lineage[list]) from a reference GenBank."""
    recs = list(SeqIO.parse(str(gb_path), "genbank"))
    if not recs:
        raise ValueError("no records in reference GenBank")
    rec = recs[0]
    organism = rec.annotations.get("organism", rec.id) or ""
    lineage = [t.strip() for t in (rec.annotations.get("taxonomy", []) or []) if t.strip()]
    genus = first_token(organism)
    family = next((t for t in lineage if t.endswith("idae")), "")
    return organism, genus, family, lineage


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-gb", required=True, type=Path)
    ap.add_argument("--sample-species", required=True,
                    help="Sample nominal species name (genus taken as first token).")
    ap.add_argument("--out", required=True, type=Path)
    # Optional sample taxonomy. --sample-family upgrades the tier from
    # NON_CONGENERIC to CONFAMILIAL / DIFFERENT_FAMILY when the samplesheet
    # supplies it; --sample-class is used only in the detail text. Both degrade
    # gracefully when absent.
    ap.add_argument("--sample-class", default="")
    ap.add_argument("--sample-family", default="")
    args = ap.parse_args()

    def finish(tier, msg):
        args.out.write_text(f"{tier}\t{msg}\n")
        print(f"[reference_divergence] {tier}: {msg}", file=sys.stderr)
        sys.exit(0)

    if not args.reference_gb.exists() or args.reference_gb.stat().st_size == 0:
        finish("UNKNOWN", f"reference missing/empty: {args.reference_gb.name}")

    try:
        ref_org, ref_genus, ref_family, ref_lineage = parse_reference(args.reference_gb)
    except Exception as exc:  # noqa: BLE001 - never break the run on a parse error
        finish("UNKNOWN", f"could not parse reference {args.reference_gb.name}: {exc}")

    sample_genus = first_token(args.sample_species)
    if not sample_genus or not ref_genus:
        finish("UNKNOWN",
               f"insufficient taxonomy: sample='{args.sample_species}' "
               f"ref_organism='{ref_org}'")

    ref_desc = ref_org + (f" [{ref_family}]" if ref_family else "")
    sample_cls = f" (sample class {args.sample_class})" if args.sample_class else ""
    detail = f"sample={args.sample_species}{sample_cls} ref={ref_desc}"

    # Congeneric: the reference is in the sample's own genus -> trusted.
    if sample_genus.lower() == ref_genus.lower():
        finish("CONGENERIC", f"genus={ref_genus}; {detail}")

    # Not congeneric. Every non-congeneric tier is review-worthy (a same-family
    # reference still collapsed OG2102); family only grades how far the NCBI
    # fallback walked, to tell the curator how surprising the result is.
    if args.sample_family and ref_family:
        if args.sample_family.lower() == ref_family.lower():
            finish("CONFAMILIAL",
                   f"family={ref_family}, sample genus {sample_genus} != ref genus "
                   f"{ref_genus}; {detail}")
        finish("DIFFERENT_FAMILY",
               f"sample family {args.sample_family} != ref family {ref_family} "
               f"(reference fell back past family level); {detail}")

    # No family supplied: we can robustly say only that it is non-congeneric.
    finish("NON_CONGENERIC",
           f"sample genus {sample_genus} != ref genus {ref_genus}, "
           f"family not supplied for finer grading; {detail}")


if __name__ == "__main__":
    main()
