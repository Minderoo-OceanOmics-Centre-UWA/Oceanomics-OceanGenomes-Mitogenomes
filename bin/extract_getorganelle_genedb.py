#!/usr/bin/env python3
"""Build a GetOrganelle custom label (gene) database from a reference GenBank.

GetOrganelle's own ``get_annotated_regions_from_gb.py`` would do this, but the
GetOrganelle BioConda container ships without Biopython, so that utility prints
"Python package biopython not found!" and silently produces nothing. This
standalone extractor runs in the MitoHiFi container (which has Biopython) and
emits FASTA records in the exact header convention GetOrganelle's bundled
``animal_mt`` label database uses and its ``--genes`` parser expects:

    >cox1 gene - NC_072168

i.e. ``>GENE_NAME gene - ACCESSION``. The gene name is taken from the feature's
``gene`` qualifier (falling back to ``product``/``locus_tag``); the nucleotide
sequence is extracted with Biopython so compound (spliced/wrapping) CDS
locations and strand are handled correctly.
"""

import argparse
import re
import sys

from Bio import SeqIO


def gene_name(feature):
    """Pick a usable gene label from a feature's qualifiers."""
    for key in ("gene", "product", "locus_tag"):
        vals = feature.qualifiers.get(key)
        if vals:
            # Normalise to a single whitespace-free token: GetOrganelle splits
            # the header on ' - ', so the label must not contain that delimiter.
            name = re.sub(r"\s+", "_", vals[0].strip())
            name = name.replace(" - ", "_")
            if name:
                return name
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("genbank", help="Reference GenBank file (one or more records).")
    ap.add_argument("out_fasta", help="Output label-database FASTA.")
    ap.add_argument(
        "--types",
        default="CDS,rRNA",
        help="Comma-separated feature types to extract (default: CDS,rRNA).",
    )
    args = ap.parse_args()

    wanted = {t.strip() for t in args.types.split(",") if t.strip()}
    written = 0
    seen = set()

    with open(args.out_fasta, "w") as out:
        for record in SeqIO.parse(args.genbank, "genbank"):
            # Drop the version suffix to match GetOrganelle's accession style
            # (>... - AB478570), though keeping it would also be harmless.
            accession = record.id.split(".")[0] if record.id else record.name
            for feature in record.features:
                if feature.type not in wanted:
                    continue
                name = gene_name(feature)
                if not name:
                    sys.stderr.write(
                        "WARN: %s feature at %s has no gene/product/locus_tag; skipped\n"
                        % (feature.type, feature.location)
                    )
                    continue
                try:
                    seq = str(feature.extract(record.seq))
                except Exception as exc:  # noqa: BLE001 - keep going on odd features
                    sys.stderr.write(
                        "WARN: could not extract %s %s: %s\n" % (name, feature.location, exc)
                    )
                    continue
                if not seq:
                    continue
                # Keep headers unique within a single GenBank so makeblastdb does
                # not choke on duplicate sequence IDs, without changing the gene
                # label (the first token, which GetOrganelle uses to classify):
                # only the trailing accession field is disambiguated.
                acc = accession
                dup = 2
                while ("%s gene - %s" % (name, acc)) in seen:
                    acc = "%s_%d" % (accession, dup)
                    dup += 1
                header = "%s gene - %s" % (name, acc)
                seen.add(header)
                out.write(">%s\n%s\n" % (header, seq))
                written += 1

    sys.stderr.write("Extracted %d gene records into %s\n" % (written, args.out_fasta))


if __name__ == "__main__":
    main()
