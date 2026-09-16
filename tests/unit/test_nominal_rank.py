"""What rank does a specimen label actually claim?

The pipeline's only route to proceed_qc=true validates a SPECIES-level claim. In
NOVA_260909_LA that test was inapplicable to nearly the whole batch: 2 of 80
samples carry a species-level claim, the rest name a genus, family, order or
class. Those were all reported as "not validated", which reads as "the assembly
contradicts the label" when it actually means "nobody identified this specimen to
species". Getting the rank right is what separates those two messages.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from nominal_rank import is_species_level, parse_rank


class ParseRankTests(unittest.TestCase):
    def test_binomial_is_species(self):
        self.assertEqual(parse_rank("Ariomma indicum"), ("species", "Ariomma indicum"))

    def test_explicit_parenthetical_rank_wins(self):
        for label, rank, taxon in [
            ("Tetractinellida (Order)", "order", "Tetractinellida"),
            ("Hexactinellida (Class)", "class", "Hexactinellida"),
            ("Cladorhizidae (Family)", "family", "Cladorhizidae"),
        ]:
            self.assertEqual(parse_rank(label), (rank, taxon), label)

    def test_latinised_family_suffix(self):
        self.assertEqual(parse_rank("Brisingidae"), ("family", "Brisingidae"))
        self.assertEqual(parse_rank("Synaphobranchidae"), ("family", "Synaphobranchidae"))

    def test_open_nomenclature_is_not_species(self):
        self.assertEqual(parse_rank("Zoanthus sp."), ("genus", "Zoanthus"))
        self.assertFalse(is_species_level("Zoanthus sp."))

    def test_cf_is_tentative_species(self):
        rank, _ = parse_rank("Munidopsis cf, aries")
        self.assertEqual(rank, "tentative_species")
        self.assertTrue(is_species_level("Munidopsis cf, aries"))

    def test_sp_does_not_imply_genus_rank(self):
        """'X sp.' says nothing about X's rank.

        This batch writes 'Asteroidea sp.' (a class), 'Porifera sp.' (a phylum)
        and 'Scleralcyonacea sp.' (an order) in exactly that form, so assuming the
        first word is a genus mislabels 26 of 80 samples.
        """
        lineage = {"class": "Asteroidea"}
        self.assertEqual(parse_rank("Asteroidea sp.", lineage), ("class", "Asteroidea"))
        lineage = {"order": "Scleralcyonacea"}
        self.assertEqual(parse_rank("Scleralcyonacea sp.", lineage),
                         ("order", "Scleralcyonacea"))

    def test_resolved_lineage_beats_suffix_guess(self):
        """A resolved rank is a fact; a Latin ending is a guess.

        Hexactinellida (class) and Tetractinellida (order) share the -ida ending,
        so no suffix rule can separate them.
        """
        self.assertEqual(parse_rank("Hexactinellida", {"class": "Hexactinellida"}),
                         ("class", "Hexactinellida"))
        self.assertEqual(parse_rank("Tetractinellida", {"order": "Tetractinellida"}),
                         ("order", "Tetractinellida"))

    def test_unknown_higher_taxon_is_not_guessed_as_order(self):
        rank, _ = parse_rank("Hexactinellida")
        self.assertEqual(rank, "higher")

    def test_non_taxonomic_labels_are_unresolved(self):
        for label in ("Larval fish", "unknown", "unidentified", "", None, "?"):
            self.assertEqual(parse_rank(label)[0], "unresolved", repr(label))

    def test_trailing_parenthetical_note_is_stripped(self):
        self.assertEqual(parse_rank("Halieutaea sp (coccinea?)"), ("genus", "Halieutaea"))

    def test_nothing_above_species_is_species_level(self):
        for label in ("Zoanthus sp.", "Brisingidae", "Tetractinellida (Order)",
                      "Larval fish", "Asteroidea sp."):
            self.assertFalse(is_species_level(label), label)


if __name__ == "__main__":
    unittest.main()
