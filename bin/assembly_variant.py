#!/usr/bin/env python3
"""Which assembly variants are real molecules, and which are scaffolding artefacts.

A multi-scaffold assembly can be published in a form that is NOT a molecule anyone
sequenced:

  <prefix>_concat  -- every contig glued end to end by SANITISE_FASTA
  <prefix>_rgj     -- contigs welded in reference order with 100bp N spacers by
                      GETORGANELLE_JOIN (off by default, see skip_getorganelle_join)

Both exist for one reason: so a partial mitogenome can still be annotated and get a
species ID. Neither is submittable, because the junctions are inferred rather than
observed and the sequence across them does not exist in the organism.

`<prefix>_collapsed` is deliberately NOT in this set. A collapse rewrites a confirmed
head-to-tail concatemer back to its monomer, which IS the real molecule.

Nothing in the pipeline used to distinguish these: the QC gate is content-based and
submission_ready is gate arithmetic, so OG2936's _rgj reached submission_ready=true in
NOVA_260909_LA and would have been submitted to ENA.

The Nextflow half of the same rule is isRescueMolecule() in
subworkflows/local/ena_submission_prep/main.nf -- a file-scope def rather than a lib/ class,
because lib/*.groovy resolves against the launched script's directory and a
nextflow_pipeline test whose script lives under tests/ cannot see it.
tests/unit/test_assembly_variant.py asserts the two suffix lists have not drifted apart.
"""

# Suffixes on the assembly CODE (the 4th dot-token of a prefix, e.g. "getorg1770_rgj")
# that mark an artificial molecule.
RESCUE_SUFFIXES = ("_rgj", "_concat")


def is_rescue_molecule(assembly_code_or_prefix: str) -> bool:
    """True when this assembly is a scaffolding artefact, not a real molecule.

    Accepts either a bare code ("getorg1770_rgj") or a full prefix / SeqID
    ("OG1.ilmn.260909.getorg1770_rgj", "OG1.ilmn.260909.getorg1770_rgj.mitos2110"),
    because callers hold different amounts of the name.
    """
    if not assembly_code_or_prefix:
        return False
    # A SeqID's annotation version is a 5th token, so test every token rather than
    # guessing which one is the code.
    return any(
        token.endswith(suffix)
        for token in str(assembly_code_or_prefix).split(".")
        for suffix in RESCUE_SUFFIXES
    )
