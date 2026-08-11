#!/usr/bin/env python3
"""
Shared species-name helpers.

Kept dependency-free so any script in bin/ can import it: Nextflow bind-mounts
the whole bin/ directory into the task container and puts it on PATH, so a
sibling import resolves via sys.path[0].
"""
import re


# 'Genus sp', 'Genus spp', 'Genus spp.' -- the qualifier must be the entire
# remainder of the name. Anything with trailing cruft ('Centrodraco sp 2',
# 'Nesogobius sp. `groove cheek`', 'Synodus macrops cf') is deliberately left
# alone rather than guessed at.
OPEN_NOMENCLATURE_RE = re.compile(r"^(?P<genus>[A-Za-z][A-Za-z-]*)\s+spp?\.?$")


def normalise_open_nomenclature(name):
    """
    Normalise an open-nomenclature species name to the ENA-submittable form.

    ENA/NCBI only recognise 'Genus sp.' for an undescribed species. 'Genus sp'
    (no period) and 'Genus spp.' are not taxa, and webin-cli rejects the
    submission with 'Organism is not Submittable' -- the /organism= value in the
    EMBL flatfile is taken verbatim from the nominal species name.

        'Chaunax sp'    -> 'Chaunax sp.'
        'Blachea spp.'  -> 'Blachea sp.'
        'Zenion spp'    -> 'Zenion sp.'

    Idempotent, and a no-op for real binomials and for anything that doesn't
    match the bare 'Genus <qualifier>' shape.
    """
    if not name:
        return name
    s = str(name).strip()
    m = OPEN_NOMENCLATURE_RE.match(s)
    if not m:
        return s
    return f"{m.group('genus')} sp."


def genus_of(name):
    """Genus = first whitespace-delimited token of a taxon name, else ''.

    Shared by the reference-divergence (taxonomy) and reference-relevance (BLAST)
    checks so both grade "is this reference congeneric with the sample?" the same
    way. Genus names are consistent between the OceanOmics species table and
    GenBank lineages, which is why the comparison is done at genus rather than at
    class -- see reference_divergence_check.classify_divergence.
    """
    return (name or "").strip().split(" ")[0] if (name or "").strip() else ""


def same_genus(a, b):
    """True when two taxon names share a genus (case-insensitive)."""
    ga, gb = genus_of(a).lower(), genus_of(b).lower()
    return bool(ga) and ga == gb
