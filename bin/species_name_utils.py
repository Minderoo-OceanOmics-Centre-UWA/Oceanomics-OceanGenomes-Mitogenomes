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

# A single alphabetic token, i.e. a name asserting no species at all.
BARE_GENUS_RE = re.compile(r"^[A-Za-z][A-Za-z-]*$")

# Qualifiers that mark a tentative determination. The name still asserts a
# species, just an uncertain one, which is why they are classified apart from
# 'sp.' below rather than lumped with it.
UNCERTAINTY_TOKENS = {"cf", "cf.", "aff", "aff.", "nr", "nr.", "?"}

# Values that mean "the taxonomy did not resolve". All of these occur in real
# samplesheets -- see isUnresolvedTaxon in subworkflows/local/prepare_samplesheet
# -- and none of them may ever be normalised into something submittable-looking
# or matched at any rank.
UNRESOLVED_TAXON_TOKENS = {"", "unknown", "na", "n/a", "none", "null", "tbc"}

# Family-name suffixes under the ICZN. This is a HEURISTIC, and it is safe here
# only because of the direction it can fail in: a false positive downgrades a
# match to a coarser rank or refuses to normalise, and can never invent a match
# that the evidence does not support.
#
# Note it under-fires on names asserting an ORDER or a CLASS, which invertebrate
# field labels carry far more often than fish ones ('Scleractinia', 'Actiniaria',
# '-ida', '-acea'). Those classify as unparseable and the sample stays held,
# which is the safe direction; it is recorded here so the behaviour is a known
# limitation rather than a surprise.
FAMILY_SUFFIXES = ("idae", "inae")


def _looks_like_family(token):
    return token.lower().endswith(FAMILY_SUFFIXES)


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
        'Serrivomer'    -> 'Serrivomer sp.'

    The bare-genus case is a correctness fix, not a convenience: ENA recognises
    'Serrivomer sp.' and does not recognise 'Serrivomer', and a sample submitted
    under the bare genus is rejected at Webin with 'Organism is not Submittable'
    after passing every gate upstream.

    Two guards on that case. A single token that looks like a FAMILY
    ('-idae'/'-inae') must not become 'Ophidiidae sp.' here -- that decision
    belongs to the family-rank handling in the validator, not to the normaliser.
    And the unresolved-taxon sentinels pass through untouched: 'unknown sp.' is a
    worse organism string than 'unknown', because it looks submittable.

    Idempotent, and a no-op for real binomials and for anything that doesn't
    match the bare 'Genus <qualifier>' shape.
    """
    if not name:
        return name
    s = str(name).strip()
    if s.lower() in UNRESOLVED_TAXON_TOKENS:
        return s
    m = OPEN_NOMENCLATURE_RE.match(s)
    if m:
        return f"{m.group('genus')} sp."
    if BARE_GENUS_RE.match(s) and not _looks_like_family(s):
        return f"{s} sp."
    return s


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


def parse_nominal(name):
    """Classify a nominal species ID by the RANK it actually asserts.

    Returns (rank, value), where rank is one of 'species', 'species_uncertain',
    'genus', 'family', or None when nothing can be asserted.

        'Notacanthus abbotti'   -> ('species', 'Notacanthus abbotti')
        'Diaphus sp.'           -> ('genus', 'Diaphus')
        'Diaphus sp 1'          -> ('genus', 'Diaphus')
        'Serrivomer'            -> ('genus', 'Serrivomer')
        'Ophidiidae'            -> ('family', 'Ophidiidae')
        'Squalus notocaudatus?' -> ('species_uncertain', 'Squalus')
        'unknown' / '' / None   -> (None, None)

    This is where the GUESS lives, deliberately separate from
    normalise_open_nomenclature. That function's regex refuses trailing cruft
    ('Genus sp 2') on purpose, and its comment about leaving such names alone is
    still right FOR NORMALISATION -- other callers depend on it. Loosening the
    regex to serve this would change their behaviour too.

    'species_uncertain' is reported as a genus-level value on purpose. A 'cf.' or
    '?' label asserts a tentative species, so resolving it at genus DROPS an
    uncertain claim rather than asserting one, which is safe in the direction
    that matters. The caller records that the downgrade happened.
    """
    if name is None:
        return (None, None)
    s = str(name).strip()
    if not s or s.lower() in UNRESOLVED_TAXON_TOKENS:
        return (None, None)

    # Drop parenthetical field notes ('Diaphus sp 2 (short jaw group)',
    # 'Spectrunculus grandis (Brown color)'), which are annotations on the name
    # rather than part of it, and are a documented cause of downstream breakage.
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    if not s:
        return (None, None)

    # A trailing '?' marks the whole determination as uncertain.
    uncertain = s.endswith("?")
    tokens = [t for t in s.rstrip("?").strip().split() if t]
    if not tokens:
        return (None, None)

    lowered = [t.lower() for t in tokens]
    if any(t in UNCERTAINTY_TOKENS for t in lowered):
        uncertain = True
        tokens = [t for t, low in zip(tokens, lowered) if low not in UNCERTAINTY_TOKENS]
    if not tokens:
        return (None, None)

    head = tokens[0]
    if not BARE_GENUS_RE.match(head):
        return (None, None)

    # A family-suffixed first token settles the rank before anything else is
    # considered. A family name cannot be the genus half of a binomial, so
    # whatever follows it is a field note rather than an epithet -- 'Macrouridae
    # black' asserts a family and a colour, not a species, and reading it as one
    # is what sent a reference lookup searching for a species that does not exist.
    if _looks_like_family(head):
        return ("family", head)

    genus = head

    # 'Genus sp', 'Genus spp.', 'Genus sp 1', 'Genus sp. 2' -- an explicit refusal
    # to name a species, so the claim is genus-level however much cruft trails it.
    if len(tokens) > 1 and re.match(r"^spp?\.?$", tokens[1], re.IGNORECASE):
        return ("genus", genus)

    if len(tokens) == 1:
        return ("genus", genus)

    epithet = tokens[1]
    if not re.match(r"^[A-Za-z][A-Za-z-]*$", epithet):
        # Not a usable epithet ('Genus 12'); fall back to the rank the first token
        # supports rather than asserting a species.
        return ("genus", genus)

    if uncertain:
        return ("species_uncertain", genus)
    return ("species", f"{genus} {epithet}")
