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
CROSS_ORDER | CONGENERIC | CONFAMILIAL | DIFFERENT_FAMILY | NON_CONGENERIC | UNKNOWN
and always exits 0, so a missing/odd reference can never break the run.

CROSS_ORDER is the most severe tier: the NCBI fallback walked past order level and
grabbed a reference from a different order entirely (e.g. a beachsalmon sample vs a
flatfish reference). At that distance MitoHiFi's read recruitment maps almost
nothing and the assembly fails outright, so a CROSS_ORDER call is the pipeline's
signal to skip the reference-based path and use the reference-free HiFi assembler
instead. It is detected only when a sample order is supplied (from the samplesheet)
and the reference GenBank lineage carries an order; without those it degrades to
the existing genus/family grading with no change in behaviour.

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


def first_token(name):
    """Genus = first whitespace-delimited token of a taxon name, else ''."""
    return (name or "").strip().split(" ")[0] if (name or "").strip() else ""


def order_from_lineage(lineage):
    """Order = the lineage token ending in 'iformes' (standard for fish and most
    actinopterygian/vertebrate orders), else ''."""
    return next((t for t in lineage if t.lower().endswith("iformes")), "")


def parse_reference(gb_path):
    """Return (organism, genus, family, order, lineage[list]) from a reference GenBank."""
    # Import Bio lazily so the taxonomy-grading logic (classify_divergence) can be
    # imported and unit-tested without biopython/numpy installed.
    from Bio import SeqIO
    recs = list(SeqIO.parse(str(gb_path), "genbank"))
    if not recs:
        raise ValueError("no records in reference GenBank")
    rec = recs[0]
    organism = rec.annotations.get("organism", rec.id) or ""
    lineage = [t.strip() for t in (rec.annotations.get("taxonomy", []) or []) if t.strip()]
    genus = first_token(organism)
    family = next((t for t in lineage if t.endswith("idae")), "")
    order = order_from_lineage(lineage)
    return organism, genus, family, order, lineage


def classify_divergence(sample_species, ref_org, ref_genus, ref_family, ref_order,
                        sample_family="", sample_order="", sample_class=""):
    """Grade the reference against the sample from taxonomy alone. Pure (no I/O),
    so it is directly unit-testable without biopython. Returns (tier, message)."""
    sample_genus = first_token(sample_species)
    if not sample_genus or not ref_genus:
        return "UNKNOWN", (f"insufficient taxonomy: sample='{sample_species}' "
                           f"ref_organism='{ref_org}'")

    ref_desc = ref_org + (f" [{ref_family}{('; ' + ref_order) if ref_order else ''}]"
                          if ref_family or ref_order else "")
    sample_cls = f" (sample class {sample_class})" if sample_class else ""
    detail = f"sample={sample_species}{sample_cls} ref={ref_desc}"

    # Congeneric: the reference is in the sample's own genus -> trusted.
    if sample_genus.lower() == ref_genus.lower():
        return "CONGENERIC", f"genus={ref_genus}; {detail}"

    # Cross-order: the reference is from a different order entirely. This is the
    # worst case (recruitment will fail; the assembly collapses to nothing) and is
    # the signal to switch to the reference-free assembler. Detected only when both
    # orders are known; otherwise fall through to genus/family grading unchanged.
    if sample_order and ref_order and sample_order.lower() != ref_order.lower():
        return "CROSS_ORDER", (f"sample order {sample_order} != ref order {ref_order} "
                               f"(reference fell back past order level); {detail}")

    # Not congeneric. Every non-congeneric tier is review-worthy (a same-family
    # reference still collapsed OG2102); family only grades how far the NCBI
    # fallback walked, to tell the curator how surprising the result is.
    if sample_family and ref_family:
        if sample_family.lower() == ref_family.lower():
            return "CONFAMILIAL", (f"family={ref_family}, sample genus {sample_genus} "
                                   f"!= ref genus {ref_genus}; {detail}")
        return "DIFFERENT_FAMILY", (f"sample family {sample_family} != ref family "
                                    f"{ref_family} (reference fell back past family level); {detail}")

    # No family supplied: we can robustly say only that it is non-congeneric.
    return "NON_CONGENERIC", (f"sample genus {sample_genus} != ref genus {ref_genus}, "
                              f"family not supplied for finer grading; {detail}")


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
    # --sample-order enables the CROSS_ORDER tier (reference from a different order
    # -> reference-based recruitment will fail; route to reference-free assembly).
    ap.add_argument("--sample-order", default="")
    args = ap.parse_args()

    def finish(tier, msg):
        args.out.write_text(f"{tier}\t{msg}\n")
        print(f"[reference_divergence] {tier}: {msg}", file=sys.stderr)
        sys.exit(0)

    if not args.reference_gb.exists() or args.reference_gb.stat().st_size == 0:
        finish("UNKNOWN", f"reference missing/empty: {args.reference_gb.name}")

    try:
        ref_org, ref_genus, ref_family, ref_order, _lineage = parse_reference(args.reference_gb)
    except Exception as exc:  # noqa: BLE001 - never break the run on a parse error
        finish("UNKNOWN", f"could not parse reference {args.reference_gb.name}: {exc}")

    tier, msg = classify_divergence(
        args.sample_species, ref_org, ref_genus, ref_family, ref_order,
        sample_family=args.sample_family, sample_order=args.sample_order,
        sample_class=args.sample_class,
    )
    finish(tier, msg)


if __name__ == "__main__":
    main()
