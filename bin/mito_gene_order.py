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
