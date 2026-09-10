"""Selection tests for bin/select_fallback_seed.py.

An empty GetOrganelle first pass has no sequence for SELECT_REFERENCE_DB to rank a
group database against, so it used to be reseeded from the WHOLE group. That is
harmless for a four-record database and bad for a 647-record one: INV12_ANOMURA was
handed all of Arthropoda while 20 Anomura references sat in the same manifest.
taxonomy_shortlist() is the replacement -- it bounds the group by the sample's own
taxonomy before reads rank what is left -- so these tests are the gate on it picking
the right tier, and on it never silently widening back to the whole group.

Deliberately assert against the REAL tracked manifests, not a fixture: the failure
being prevented is a mismatch between the cascade and the lineage strings those
manifests actually carry, which a hand-written fixture cannot reproduce. The gap
case is synthetic because, since the ctenophore rebuild, no tracked group still has
a taxon with zero records.

No Biopython import here, deliberately: the module keeps Bio lazy so this half stays
testable in a plain interpreter, and importing it under one is part of the contract.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# bin/ on sys.path so sibling imports resolve as they do in the task container.
sys.path.insert(0, str(ROOT / "bin"))

import refdb_record                                          # noqa: E402
from select_fallback_seed import (balanced_cap, norm,        # noqa: E402
                                  taxonomy_shortlist)

REFDB = ROOT / "assets" / "refdb"


def manifest(group):
    _fasta, manifest_path, _features = refdb_record.refdb_paths(REFDB / group, group)
    return refdb_record.read_manifest(manifest_path)


def organisms(rows, accessions):
    return sorted(rows[a]["organism"] for a in accessions)


class TaxonomyShortlist(unittest.TestCase):
    def test_anomura_is_bounded_to_its_order(self):
        """INV12_ANOMURA: 20 Anomura references, not all 647 arthropods."""
        rows = manifest("arthropoda")
        hits, tier, missed = taxonomy_shortlist(rows, order="Anomura", max_candidates=50)
        self.assertEqual(tier, "order")
        self.assertEqual(missed, [])
        self.assertLess(len(hits), len(rows) / 10)
        # every survivor is genuinely an anomuran, not merely a near-miss on the string
        for acc in hits:
            self.assertIn("anomura", [norm(t) for t in rows[acc]["lineage"]])

    def test_platyctenida_resolves_since_the_ctenophore_rebuild(self):
        """INV08_TJALFIELLA was called an unfixable reference gap.

        It was not: the group held four records only because the builder searched
        RefSeq alone. With refseq_only off for ctenophora the order resolves, and it
        carries the sample's OWN genus.
        """
        rows = manifest("ctenophora")
        hits, tier, missed = taxonomy_shortlist(rows, order="Platyctenida", max_candidates=50)
        self.assertEqual(tier, "order")
        self.assertEqual(missed, [])
        self.assertIn("Tjalfiella sp.", organisms(rows, hits))

    def test_nominal_taxon_beats_the_coarser_tiers(self):
        """A resolvable nominal name must not be skipped in favour of its order."""
        rows = manifest("ctenophora")
        hits, tier, _missed = taxonomy_shortlist(
            rows, taxon="Tjalfiella sp.", order="Platyctenida", max_candidates=50)
        self.assertEqual(tier, "nominal")
        self.assertEqual(set(organisms(rows, hits)), {"Tjalfiella sp."})

    def test_chiton_rescue_still_narrows(self):
        """INV10_CHITON was rescued by the old whole-mollusc seed; keep it bounded."""
        rows = manifest("mollusca")
        hits, tier, _missed = taxonomy_shortlist(
            rows, order="Chitonida", tax_class="Polyplacophora", max_candidates=50)
        self.assertEqual(tier, "order")
        self.assertLess(len(hits), len(rows) / 10)

    def test_class_tier_is_reached_when_order_misses(self):
        rows = manifest("mollusca")
        hits, tier, missed = taxonomy_shortlist(
            rows, order="Nonesuchida", tax_class="Polyplacophora", max_candidates=50)
        self.assertEqual(tier, "class")
        self.assertEqual(missed, ["order:Nonesuchida"])
        self.assertTrue(hits)

    def test_absent_taxon_reports_the_gap_and_falls_back_to_the_group(self):
        """The REFERENCE_GAP path. No tracked group still exercises it, so build one.

        INV05_FARREA is the live example: Farrea has no complete mitogenome in
        GenBank at all, so widening the search cannot help it and the honest answer
        is an audited gap rather than a silent whole-group seed.
        """
        rows = {
            "AA000001.1": {"organism": "Farrea sp.", "family": "Farreidae",
                           "lineage": ["Eukaryota", "Metazoa", "Porifera"], "length_bp": 10},
        }
        hits, tier, missed = taxonomy_shortlist(
            rows, taxon="Farrea occa", family="Farreidae2", order="Sceptrulophora",
            max_candidates=50)
        self.assertEqual(tier, "group")
        self.assertEqual(missed, ["nominal:Farrea occa", "family:Farreidae2",
                                  "order:Sceptrulophora"])
        self.assertEqual(hits, ["AA000001.1"])

    def test_unresolved_placeholders_are_skipped_not_matched(self):
        rows = manifest("ctenophora")
        for placeholder in ("", "unknown", "na", "none", "dropped", "  Unknown  "):
            _hits, tier, missed = taxonomy_shortlist(
                rows, order=placeholder, max_candidates=50)
            self.assertEqual(tier, "group", placeholder)
            # skipped entirely: an unresolved rank is not a reference gap
            self.assertEqual(missed, [], placeholder)


class BalancedCap(unittest.TestCase):
    def setUp(self):
        self.rows = manifest("mollusca")

    def test_never_exceeds_the_cap(self):
        for limit in (1, 5, 20, 50):
            self.assertEqual(len(balanced_cap(list(self.rows), self.rows, limit)), limit)

    def test_under_the_cap_is_returned_whole_and_sorted(self):
        subset = sorted(self.rows)[:5]
        self.assertEqual(balanced_cap(subset, self.rows, 50), sorted(subset))

    def test_is_deterministic_regardless_of_input_order(self):
        accs = list(self.rows)
        first = balanced_cap(accs, self.rows, 20)
        self.assertEqual(first, balanced_cap(list(reversed(accs)), self.rows, 20))

    def test_spreads_across_families_rather_than_taking_the_first_n(self):
        """The point of the cap: a bounded panel that is still taxonomically broad.

        Taking the first N accessions of a large group tends to collect siblings,
        which is what a top-n seed panel least needs.
        """
        accs = list(self.rows)
        picked = balanced_cap(accs, self.rows, 20)
        fam = lambda a: norm(self.rows[a].get("family")) or "_unknown"   # noqa: E731
        self.assertGreater(len({fam(a) for a in picked}),
                           len({fam(a) for a in sorted(accs)[:20]}))


class UnionCoverage(unittest.TestCase):
    """Moved to refdb_record so the fallback selector can use it without Biopython."""

    def test_merges_overlapping_and_abutting_intervals(self):
        self.assertEqual(refdb_record.union_coverage([]), 0)
        self.assertEqual(refdb_record.union_coverage([(1, 10)]), 10)
        self.assertEqual(refdb_record.union_coverage([(1, 10), (5, 20)]), 20)
        self.assertEqual(refdb_record.union_coverage([(1, 10), (11, 20)]), 20)
        self.assertEqual(refdb_record.union_coverage([(1, 10), (21, 30)]), 20)
        self.assertEqual(refdb_record.union_coverage([(21, 30), (1, 10)]), 20)


if __name__ == "__main__":
    unittest.main()
