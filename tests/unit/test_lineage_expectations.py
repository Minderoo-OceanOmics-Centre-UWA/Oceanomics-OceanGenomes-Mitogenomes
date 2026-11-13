"""Genes expected outside the vertebrate set, keyed on taxonomic rank.

Two things this must get right, because getting either wrong fails a whole
lineage: mtMutS is an OCTOCORAL gene, not an anthozoan one (300/302 Octocorallia
reference records carry it, 0/335 Hexacorallia), and ATP9 in sponges is common
but not universal (81/87), so it can never be blocking.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from lineage_expectations import classify_missing, expected_genes, load_rules


class RuleAssetTests(unittest.TestCase):
    def test_asset_loads(self):
        self.assertTrue(load_rules(), "assets/lineage_gene_expectations.json is empty")

    def test_every_rule_carries_its_evidence(self):
        """A rule that fails a sample without a stated reason is not reviewable."""
        for rule in load_rules():
            self.assertTrue(str(rule.get("evidence", "")).strip(),
                            f"rule has no evidence: {rule}")


class OctocoralTests(unittest.TestCase):
    def test_octocoral_orders_require_mtmuts(self):
        for order in ("Malacalcyonacea", "Scleralcyonacea"):
            blocking, _ = classify_missing({"order": order, "class": "Anthozoa"}, ["CO1"])
            self.assertEqual(blocking, ["MTMUTS"], order)

    def test_hexacorals_do_not_require_mtmuts(self):
        """Keying this on class Anthozoa instead of the octocoral orders would
        wrongly fail every sea anemone and zoanthid in the batch."""
        blocking, advisory = classify_missing(
            {"order": "Zoantharia", "class": "Anthozoa"}, ["CO1"])
        self.assertEqual(blocking, [])
        self.assertEqual(advisory, [])

    def test_present_mtmuts_is_not_reported_missing(self):
        blocking, _ = classify_missing(
            {"order": "Malacalcyonacea"}, ["CO1", "MTMUTS"])
        self.assertEqual(blocking, [])


class PoriferaTests(unittest.TestCase):
    def test_atp9_is_advisory_never_blocking(self):
        for taxonomy in ({"class": "Demospongiae"}, {"phylum": "Porifera"},
                         {"class": "Hexactinellida"}):
            blocking, advisory = classify_missing(taxonomy, ["CO1"])
            self.assertEqual(blocking, [], taxonomy)
            self.assertEqual(advisory, ["ATP9"], taxonomy)


class RankResolutionTests(unittest.TestCase):
    def test_blank_ranks_match_nothing(self):
        """Most invertebrates in this batch carry no order label at all, so an
        order-keyed rule has to be inert rather than throwing."""
        blocking, advisory = classify_missing({"class": "", "order": None}, [])
        self.assertEqual((blocking, advisory), ([], []))

    def test_matching_is_case_insensitive(self):
        blocking, _ = classify_missing({"order": "malacalcyonacea"}, [])
        self.assertEqual(blocking, ["MTMUTS"])

    def test_most_specific_rank_wins(self):
        """A class rule must not dilute an order rule for the same gene."""
        genes = expected_genes({"order": "Malacalcyonacea", "class": "Anthozoa"})
        self.assertEqual(genes.get("MTMUTS"), "blocking")

    def test_unknown_taxon_matches_nothing(self):
        self.assertEqual(expected_genes({"class": "Actinopteri"}), {})


if __name__ == "__main__":
    unittest.main()
