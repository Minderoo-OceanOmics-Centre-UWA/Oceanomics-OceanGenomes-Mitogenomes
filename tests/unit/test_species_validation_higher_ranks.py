"""Nominal labels above family: resolved to their real rank, reported, not released.

NOVA_260909_LA regression. species_validation.py classified the nominal label with
species_name_utils.parse_nominal, whose vocabulary stops at family, so a class or
order label became a genus claim: 'Asteroidea sp.' read as ('genus', 'Asteroidea'),
a genus that does not exist. It could never match, and the sample was reported as a
failed species validation rather than as a specimen identified only to class. That
covered 39 of the 80 labels in that samplesheet.

Two invariants are asserted here, because they pull in opposite directions:
  - the rank is now resolved correctly, against the sample's own lineage rather than
    a Latin suffix (no suffix separates Asteroidea from Scleralcyonacea);
  - nothing new passes the QC gate. bin/nominal_rank.py states the principle:
    knowing the rank "does not let anything new through". Found_in_blast_YN stays
    "No" above family; only validated_rank changes, from 'unmatched' to
    '<rank>_unvalidated'.
"""
import sys
import unittest
from collections import namedtuple
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from species_validation import match_at_rank, parse_lineage_arg, resolve_nominal  # noqa: E402

BlastHits = namedtuple("BlastHits", "species genera")
NO_HITS = BlastHits(species=set(), genera=set())

ASTEROID_LINEAGE = {"class": "Asteroidea", "order": "Velatida", "family": "Pterasteridae"}
SPONGE_LINEAGE = {"class": "Demospongiae", "order": "Tetractinellida"}


def lca_row(**kw):
    row = {"genus": "", "family": "", "order": "", "class": "", "phylum": ""}
    row.update(kw)
    return row


class ResolveNominalTests(unittest.TestCase):
    def test_class_label_with_explicit_parenthetical(self):
        self.assertEqual(resolve_nominal("Asteroidea (Class)", ASTEROID_LINEAGE),
                         ("class", "Asteroidea"))

    def test_order_label_with_explicit_parenthetical(self):
        self.assertEqual(resolve_nominal("Tetractinellida (Order)", SPONGE_LINEAGE),
                         ("order", "Tetractinellida"))

    def test_class_label_wearing_open_nomenclature(self):
        """'Asteroidea sp.' is a class, not a genus. This is the actual defect."""
        self.assertEqual(resolve_nominal("Asteroidea sp.", ASTEROID_LINEAGE),
                         ("class", "Asteroidea"))

    def test_binomial_and_family_are_untouched(self):
        self.assertEqual(resolve_nominal("Ariomma indicum", {})[0], "species")
        self.assertEqual(resolve_nominal("Scalpellidae", {}), ("family", "Scalpellidae"))

    def test_genus_sp_still_resolves_to_its_genus(self):
        self.assertEqual(resolve_nominal("Chaceon sp.", {"class": "Malacostraca"}),
                         ("genus", "Chaceon"))

    def test_tentative_species_yields_the_genus_not_the_whole_label(self):
        """match_at_rank compares this value against the LCA genus.

        parse_rank returns the full label for a tentative species; handing that
        through unchanged stopped OG2982 ('Munidopsis cf, aries') matching at all.
        """
        rank, value = resolve_nominal("Munidopsis cf, aries", {"family": "Munidopsidae"})
        self.assertEqual(rank, "species_uncertain")
        self.assertEqual(value, "Munidopsis")


class MatchAtRankHigherTests(unittest.TestCase):
    def test_class_agreement_is_recorded_but_not_validated(self):
        found, rank = match_at_rank(("class", "Asteroidea"), NO_HITS,
                                    lca_row(**{"class": "Asteroidea"}))
        self.assertEqual(found, "No")
        self.assertEqual(rank, "class_unvalidated")

    def test_class_contradiction_stays_unmatched(self):
        """A genuine contradiction must remain distinguishable from "never claimed"."""
        found, rank = match_at_rank(("class", "Asteroidea"), NO_HITS,
                                    lca_row(**{"class": "Malacostraca"}))
        self.assertEqual((found, rank), ("No", "unmatched"))

    def test_order_agreement_is_recorded(self):
        found, rank = match_at_rank(("order", "Tetractinellida"), NO_HITS,
                                    lca_row(order="Tetractinellida"))
        self.assertEqual((found, rank), ("No", "order_unvalidated"))

    def test_higher_tries_order_then_class_then_phylum(self):
        _, rank = match_at_rank(("higher", "Porifera"), NO_HITS, lca_row(phylum="Porifera"))
        self.assertEqual(rank, "phylum_unvalidated")

    def test_tentative_species_still_validates_at_genus(self):
        hits = BlastHits(species={"munidopsis aries"}, genera={"munidopsis"})
        found, rank = match_at_rank(("species_uncertain", "Munidopsis"), hits,
                                    lca_row(genus="Munidopsis"))
        self.assertEqual((found, rank), ("Yes", "genus_downgraded"))


class LineageArgTests(unittest.TestCase):
    def test_parses_pairs_and_ignores_blanks(self):
        self.assertEqual(parse_lineage_arg("class=Asteroidea,order=,family=Pterasteridae"),
                         {"class": "Asteroidea", "family": "Pterasteridae"})

    def test_empty_is_an_empty_mapping(self):
        self.assertEqual(parse_lineage_arg(None), {})
        self.assertEqual(parse_lineage_arg(""), {})

    def test_without_lineage_a_higher_label_is_never_a_confident_class(self):
        """No lineage: the label must not be upgraded past 'higher' by suffix alone."""
        self.assertIn(resolve_nominal("Tetractinellida", {})[0], ("higher", "genus"))


if __name__ == "__main__":
    unittest.main()
