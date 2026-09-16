#!/usr/bin/env python3
"""
Gene-name translation between the pipeline's internal EMMA contract and the
names INSDC/ENA actually accept.

Pure stdlib, no imports outside this file, so it loads in every container the
annotation chain uses (the psycopg2 image, the MITOS image, the EMMA image).

Why this exists
---------------
Internally every gene is an uppercase bare name (CO1, ND4L, RNR2, ATP9, MTMUTS)
so it can key the `features` dict, the per-gene FASTA filenames and the
annotation-stats columns uniformly. `mitos_to_emma.py` writes it into the .tbl as
`MT-<NAME>`, and `process_files.py` strips the `MT-` on the way to table2asn.

That round-trip is lossless for every vertebrate gene, because their INSDC names
are uppercase anyway. It is NOT lossless for the non-vertebrate genes: the
octocoral mismatch-repair gene is `mtMutS` in INSDC nomenclature, and `MTMUTS` is
not a recognised gene name. The case has to be restored at the submission
boundary, and it has to be restored from a map rather than a special case, so the
next non-vertebrate gene with a mixed-case INSDC name has somewhere to go.
"""

# Internal uppercase key -> the /gene value INSDC expects.
# Only genes whose INSDC spelling is NOT simply the uppercase form belong here;
# anything absent passes through unchanged, which is the right answer for CO1,
# ND4L, ATP9 and every vertebrate gene.
SUBMISSION_GENE_NAME = {
    # Octocoral (and some demosponge) mitochondrial mismatch-repair gene. INSDC
    # spelling is mixed-case mtMutS; MITOS calls it msh1. See NC_ octocoral
    # reference records, which uniformly use /gene="mtMutS".
    "MTMUTS": "mtMutS",
}

# Products for the non-vertebrate genes, kept beside the names they belong to.
SUBMISSION_PRODUCT = {
    "ATP9": "ATP synthase F0 subunit 9",
    "MTMUTS": "mismatch repair protein mtMutS",
}


def submission_gene_name(internal_name):
    """The /gene value to write for an internal gene key.

    Accepts the name with or without the MT- prefix and returns it prefix-free,
    because the only caller that matters (process_files.py) is stripping MT- at
    the same point it needs the mapping.
    """
    bare = bare_gene_name(internal_name)
    return SUBMISSION_GENE_NAME.get(bare.upper(), bare)


def bare_gene_name(value):
    """Strip the MT- prefix from a gene name, a qualifier line, or a list of them.

    Tolerant on purpose: process_files.py holds qualifier lines as raw tab-joined
    .tbl text, and the caller should not have to parse them to ask this question.
    Returns '' when no gene name is present.
    """
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        for key in ("gene", "product"):
            for q in value:
                cols = str(q).split("\t")
                if len(cols) >= 2 and cols[-2] == key:
                    return bare_gene_name(cols[-1])
        return ""
    text = str(value).strip()
    if "\t" in text:
        cols = text.split("\t")
        text = cols[-1].strip()
    if text.upper().startswith("MT-"):
        text = text[3:]
    return text
