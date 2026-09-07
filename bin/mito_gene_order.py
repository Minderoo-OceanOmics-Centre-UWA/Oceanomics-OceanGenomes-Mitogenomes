#!/usr/bin/env python3
"""
Shared vertebrate mitochondrial gene-order and tRNA-naming tables.

Pure stdlib: no Biopython, no BLAST. Nextflow bind-mounts the whole bin/ dir onto
PATH, so a sibling import resolves via sys.path[0] (same pattern as
bin/orf_utils.py and bin/geo_loc_name_utils.py). This module must import cleanly
in the stdlib-only tylerpeirce/psycopg2:0.1 container.

REF_GENES was previously copy-pasted into annotation_stats.py,
emma_rescue_gate.py, trna_rescue_gate.py and rescue_trna.py, each carrying a
"keep this in sync" comment that had already drifted (three of them said "change
both"/"change all three" while there were four copies). The gates and the QC step
must agree byte-for-byte on what "present and in order" means, so the list lives
here and nowhere else.

Agreeing on the LIST turned out not to be enough: they also have to agree on how
a GFF is read into an ordered gene list, and they did not. parse_gff_attributes
had three copies and genes_by_coord had two byte-identical ones, while
annotation_stats.py kept a third, different implementation inline that deduped a
repeated gene by FIRST LINE SEEN rather than by lowest start. On a gene split
across the origin those two readings disagree, so the same GFF could be judged
in-order by a rescue gate and out-of-order by the QC step. Both readers now live
here, for the same reason REF_GENES does.
"""

# Reference gene order (tRNA, rRNA, CDS) -- the standard vertebrate set.
REF_GENES = [
    "TF", "RNR1", "TV", "RNR2", "TL2", "ND1", "TI", "TQ",
    "TM", "ND2", "TW", "TA", "TN", "TC", "TY", "CO1", "TS2",
    "TD", "CO2", "TK", "ATP8", "ATP6", "CO3", "TG", "ND3",
    "TR", "ND4L", "ND4", "TH", "TS1", "TL1", "ND5", "ND6",
    "TE", "CYTB", "TT", "TP",
]

TRNA_GENES = {g for g in REF_GENES if g.startswith("T")}   # the 22 mt tRNAs
RRNA_GENES = {"RNR1", "RNR2"}
PCG_GENES = {g for g in REF_GENES if g not in TRNA_GENES | RRNA_GENES}  # 13 PCGs

# EMMA-style tRNA gene suffix (the part after the leading 'T', e.g. "W", "S1",
# "L2") -> 3-letter amino acid, for tRNA /product strings ("tRNA-Trp(UCA)").
TRNA_AA = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L1": "Leu", "L2": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S1": "Ser", "S2": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}

# Canonical anticodon (DNA form, as tRNAscan-SE reports it) for each REF tRNA.
# The RNA form of these is what appears in the /product string, so this table and
# TRNA_PRODUCT below are two views of the same fact.
TRNA_ANTICODON = {
    "TF": "GAA", "TV": "TAC", "TL2": "TAA", "TL1": "TAG", "TI": "GAT",
    "TQ": "TTG", "TM": "CAT", "TW": "TCA", "TA": "TGC", "TN": "GTT",
    "TC": "GCA", "TY": "GTA", "TS2": "TGA", "TS1": "GCT", "TD": "GTC",
    "TK": "TTT", "TG": "TCC", "TR": "TCG", "TH": "GTG", "TE": "TTC",
    "TT": "TGT", "TP": "TGG",
}


def trna_product(emma_name, anticodon):
    """Return an EMMA-style tRNA /product string, e.g. 'tRNA-Trp(UCA)'.

    ``anticodon`` may be given in DNA (TCA) or RNA (UCA) form; it is transcribed
    to RNA for the product string. Falls back to a bare gene-name product if the
    suffix or anticodon isn't recognised -- a cosmetic field should never fail
    the run.
    """
    suffix = emma_name[1:] if emma_name.startswith("T") else emma_name
    aa = TRNA_AA.get(suffix)
    if aa and anticodon:
        return f"tRNA-{aa}({anticodon.upper().replace('T', 'U')})"
    return f"tRNA-{emma_name}"


# REF tRNA name -> canonical /product string, e.g. 'TW' -> 'tRNA-Trp(UCA)'.
TRNA_PRODUCT = {
    name: trna_product(name, ac) for name, ac in TRNA_ANTICODON.items()
}


# ----------------------------------------------------------------- GFF reading
#
# Shared by annotation_stats.py, emma_rescue_gate.py and trna_rescue_gate.py. See
# the module docstring: the gates and the QC step have to agree on how a GFF
# becomes an ordered gene list, not just on what the reference order is.


def parse_gff_attributes(attr_str):
    """Parse a GFF9 attribute column into a dict, ignoring malformed entries."""
    return dict(
        item.split("=", 1)
        for item in attr_str.strip().split(";")
        if "=" in item
    )


def gene_entries_by_coord(gff_path):
    """Return [(name, start, end, strand)] for the GFF's `gene` lines, by start.

    A gene written as more than one `gene` line -- which is how an origin-spanning
    feature is represented -- is kept ONCE, at its LOWEST start. Keeping the first
    line seen instead puts such a gene at whichever half the file happened to list
    first, which can make the order check fail on an annotation that is in fact
    correctly ordered; and a failed order check both fails the QC gate and silently
    suppresses the rescue gates.

    Lowest-start is deterministic where first-seen was not, so for the usual
    coordinate-sorted GFF nothing changes, and for a split gene the result no
    longer depends on file order. The returned start/end are the LOWER fragment's,
    which is what makes the ordering meaningful; a caller wanting the true span of
    an origin-spanning gene has to reconstruct it from the GFF itself.

    The `MT-` prefix EMMA writes is stripped, so names match REF_GENES.
    """
    best = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            attrs = parse_gff_attributes(parts[8])
            name = attrs.get("Name")
            if not name:
                continue
            gene = name.replace("MT-", "")
            start = int(parts[3])
            if gene not in best or start < best[gene][1]:
                best[gene] = (gene, start, int(parts[4]), parts[6])
    return sorted(best.values(), key=lambda entry: entry[1])


def genes_by_coord(gff_path):
    """Gene names present in the GFF, ordered by genomic start."""
    return [entry[0] for entry in gene_entries_by_coord(gff_path)]
