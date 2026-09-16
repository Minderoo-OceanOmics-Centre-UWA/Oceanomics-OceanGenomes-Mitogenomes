#!/usr/bin/env python3
"""
What taxonomic rank does a specimen's nominal_species_id actually claim?

The pipeline's only route to proceed_qc=true is validating a SPECIES-level
nominal claim through the >=98% BLAST/LCA path. In NOVA_260909_LA that test was
inapplicable to almost the whole batch: of 80 samples, 3 carry a binomial, 47 are
'Genus sp.', 12 name a family, 10 a class and 8 an order. Those 77 were reported
as "not validated (no Found_in_blast_YN=Yes)" -- indistinguishable from a genuine
contradiction, and 47 of them were held with a clean annotation and no stated
reason beyond species=unknown.

Knowing the rank does not let anything new through. It lets the pipeline say
"this specimen was never identified to species" instead of "this specimen failed
species validation", which are very different instructions to the person holding
the jar.

Pure stdlib.
"""

import re

# Ranks the parser can return, ordered most to least specific.
RANKS = ("species", "tentative_species", "genus", "family", "order", "class",
         "phylum", "higher", "unresolved")

# 'higher' means: above species, but the label alone does not say which rank.
# Latinised endings cannot separate class from order -- Asteroidea (class) and
# Scleralcyonacea (order) are both -acea/-oidea, Hexactinellida (class) and
# Tetractinellida (order) are both -ida. Guessing produces confident wrong labels
# in run reports, so the parser says 'higher' and lets the caller pin the rank
# from the resolved lineage, which the pipeline already has in
# samplesheet/taxonomy_resolution.tsv.

# Explicit rank in parentheses, as the sample database writes it:
#   'Tetractinellida (Order)', 'Cladorhizidae (Family)', 'Hexactinellida (Class)'
_PARENTHETICAL = re.compile(r"\(\s*(family|order|class|phylum|genus)\s*\)\s*$", re.I)

# Open nomenclature: 'Genus sp.', 'Genus spp.', 'Genus sp', 'Genus cf. species',
# 'Genus aff. species', 'Munidopsis cf, aries' (yes, with a comma).
_OPEN_NOMEN = re.compile(r"^\s*([A-Z][a-z]+)\s+(?:sp{1,2}\.?|cf\.?,?|aff\.?,?)(\s|$)", re.I)

# Latinised family/superfamily endings. -aceae is botanical but harmless here.
_FAMILY_SUFFIX = ("idae", "aceae")
# Common latinised higher-rank endings, used only when nothing more specific fits.
_HIGHER_SUFFIX = ("iformes", "ida", "ales", "acea", "oidea", "ozoa", "phora", "poda")

# Labels that are not taxonomy at all. These must not be fuzzy-matched to a taxon:
# 'Larval fish' is a life stage, not a name, and guessing one would be worse than
# admitting the specimen is unidentified.
_NON_TAXONOMIC = re.compile(
    r"^\s*(?:(?:larval|juvenile|unknown|unidentified|n/?a|none|tbd)\b|[?\-\s]+$)", re.I)


def _clean(label):
    text = str(label or "").strip()
    # Drop a trailing parenthetical that is NOT a rank marker, e.g.
    # 'Halieutaea sp (coccinea?)' -> 'Halieutaea sp'
    if not _PARENTHETICAL.search(text):
        text = re.sub(r"\s*\([^)]*\)\s*$", "", text).strip()
    return text


def parse_rank(label, lineage=None):
    """Return (rank, taxon) for a nominal_species_id.

    `taxon` is the name the claim is about at that rank: the binomial for a
    species, the genus for 'Genus sp.', the family/order/class name otherwise.

    `lineage` is an optional mapping of resolved taxonomy
    ({'genus':.., 'family':.., 'order':.., 'class':.., 'phylum':..}). When the
    name extracted from the label matches one of those, that rank is used
    verbatim -- it is a resolved fact rather than a guess at a Latin suffix.
    """
    text = _clean(label)
    if not text or _NON_TAXONOMIC.match(text):
        return "unresolved", ""

    m = _PARENTHETICAL.search(text)
    if m:
        rank = m.group(1).lower()
        name = _PARENTHETICAL.sub("", text).strip()
        return rank, name

    m = _OPEN_NOMEN.match(text)
    if m:
        name = m.group(1)
        # 'Genus cf. species' still names a species, tentatively; 'Genus sp.' does not.
        tail = text[m.end(1):].strip().lower()
        if tail.startswith(("cf", "aff")) and len(text.split()) >= 3:
            return "tentative_species", text
        # 'X sp.' does NOT mean X is a genus. This batch writes 'Asteroidea sp.'
        # (a class), 'Porifera sp.' (a phylum) and 'Scleralcyonacea sp.' (an
        # order) in exactly this form, so the rank has to be resolved, not assumed.
        return _rank_of_name(name, lineage), name

    words = text.split()
    if len(words) == 1:
        return _rank_of_name(words[0], lineage), words[0]

    if len(words) >= 2 and words[1][:1].islower() and words[1].isalpha():
        return "species", f"{words[0]} {words[1]}"

    return "unresolved", text


def _rank_of_name(name, lineage=None):
    """Rank of a bare taxon name: resolved lineage first, suffix as a last resort."""
    low = str(name).strip().lower()
    for rank in ("genus", "family", "order", "class", "phylum"):
        value = str((lineage or {}).get(rank) or "").strip().lower()
        if value and value == low:
            return rank
    if low.endswith(_FAMILY_SUFFIX):
        return "family"
    if low.endswith(_HIGHER_SUFFIX):
        return "higher"
    # A bare capitalised word matching nothing resolved and carrying no latinised
    # higher-rank ending reads as a genus. It is the conservative answer: the one
    # rank below species that still cannot satisfy the species validator.
    return "genus"


def is_species_level(label):
    """True only when the label makes a claim the >=98% species validator can test."""
    return parse_rank(label)[0] in ("species", "tentative_species")
