"""Grading tests for bin/audit_reference_selection_diff.py.

The audit decides whether rebuilding a group database moved a sample's reference
CLOSER or FURTHER, and that verdict is the evidence a rebuild is committed on. So the
grading has to be right about the two things invertebrate taxonomy makes awkward:

  * Invertebrate orders have no shared suffix. reference_divergence_check grades order
    via order_from_lineage(), which only matches 'iformes' and therefore never fires
    for Anomura, Platyctenida or Scleractinia. This grader matches the sample's own
    order against the reference lineage instead, so the tier exists at all.
  * OceanOmics invert labels are frequently coarser than a species. INV01's
    nominal_species_id is 'Acanthogorgiidae' -- a family. Read as a genus it matches
    no reference genus, which would grade a perfectly good confamilial reference as
    DISTANT and turn a real improvement into a 'same' verdict.

Everything under test is pure, so it runs without biopython, blastn or a database.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import audit_reference_selection_diff as audit  # noqa: E402

# INV08, the sample the whole ctenophore rebuild was for.
TJALFIELLA = {"nominal_species_id": "Tjalfiella", "family": "Tjalfiellidae",
              "order": "Platyctenida", "class": "Tentaculata"}
TJALFIELLA_LINEAGE = ["Eukaryota", "Metazoa", "Ctenophora", "Tentaculata",
                      "Platyctenida", "Tjalfiellidae"]
# What INV08 was actually seeded from before the rebuild: another order entirely.
BEROE_LINEAGE = ["Eukaryota", "Metazoa", "Ctenophora", "Nuda", "Beroida", "Beroidae"]

# INV01, whose label is a family name rather than a genus.
ACANTHOGORGIA = {"nominal_species_id": "Acanthogorgiidae", "family": "Acanthogorgiidae",
                 "order": "Malacalcyonacea", "class": "Anthozoa"}
ACANTHOGORGIIDAE_LINEAGE = ["Eukaryota", "Metazoa", "Cnidaria", "Anthozoa",
                            "Octocorallia", "Malacalcyonacea", "Acanthogorgiidae"]


def result(acc, tier, state="SELECTED"):
    return {"acc": acc, "tier": tier, "score": audit.TIERS[tier], "state": state}


class Grading(unittest.TestCase):
    def test_congeneric_beats_every_other_tier(self):
        tier, score = audit.grade(TJALFIELLA, "Tjalfiella sp.", TJALFIELLA_LINEAGE)
        self.assertEqual(tier, "CONGENERIC")
        self.assertEqual(score, max(audit.TIERS.values()))

    def test_invertebrate_order_is_graded(self):
        """The tier reference_divergence_check cannot reach: Platyctenida has no
        'iformes' suffix, so order_from_lineage() returns '' for every ctenophore."""
        other_platyctenid = ["Eukaryota", "Metazoa", "Ctenophora", "Tentaculata",
                             "Platyctenida", "Lyroctenidae"]
        tier, _ = audit.grade(TJALFIELLA, "Lyrocteis imperatoris", other_platyctenid)
        self.assertEqual(tier, "SAME_ORDER")

    def test_the_pre_rebuild_inv08_pick_shares_only_the_phylum(self):
        """What INV08 was actually seeded from. Tjalfiella is class Tentaculata and
        Beroe is class Nuda, so the two share nothing below Ctenophora -- which is
        why DISTANT, inside a per-phylum database, is the floor rather than a
        statement that the reference is from another phylum."""
        tier, _ = audit.grade(TJALFIELLA, "Beroe forskalii", BEROE_LINEAGE)
        self.assertEqual(tier, "DISTANT")

    def test_same_class_sits_between_order_and_distant(self):
        another_tentaculate = ["Eukaryota", "Metazoa", "Ctenophora", "Tentaculata",
                               "Cydippida", "Pleurobrachiidae"]
        tier, _ = audit.grade(TJALFIELLA, "Pleurobrachia bachei", another_tentaculate)
        self.assertEqual(tier, "SAME_CLASS")
        self.assertLess(audit.TIERS["SAME_CLASS"], audit.TIERS["SAME_ORDER"])
        self.assertGreater(audit.TIERS["SAME_CLASS"], audit.TIERS["DISTANT"])

    def test_a_reference_with_no_organism_is_unknown_not_distant(self):
        """UNKNOWN must score below DISTANT: 'we could not grade it' is not evidence
        of an improvement, and it must never produce a 'closer' verdict."""
        tier, score = audit.grade(TJALFIELLA, "", [])
        self.assertEqual(tier, "UNKNOWN")
        self.assertLess(score, audit.TIERS["DISTANT"])

    def test_family_level_label_still_grades_confamilial(self):
        """INV01: a family-name label must not be read as a genus."""
        tier, _ = audit.grade(ACANTHOGORGIA, "Acanthogorgia sp.",
                              ACANTHOGORGIIDAE_LINEAGE)
        self.assertEqual(tier, "CONFAMILIAL")

    def test_sample_genus_rejects_family_shaped_labels(self):
        self.assertEqual(audit.sample_genus("Acanthogorgiidae", "Acanthogorgiidae"), "")
        self.assertEqual(audit.sample_genus("Umbellula", "Umbellulidae"), "Umbellula")
        self.assertEqual(audit.sample_genus("Tjalfiella sp.", "Tjalfiellidae"),
                         "Tjalfiella")
        self.assertEqual(audit.sample_genus("", "Tjalfiellidae"), "")

    def test_grading_is_case_insensitive(self):
        lower = [t.lower() for t in TJALFIELLA_LINEAGE]
        self.assertEqual(audit.grade(TJALFIELLA, "TJALFIELLA SP.", lower)[0],
                         "CONGENERIC")


class Verdicts(unittest.TestCase):
    def test_same_accession_is_equal_record_even_at_an_ungraded_tier(self):
        r = result("NC_038065.1", "UNKNOWN")
        self.assertEqual(audit.verdict_of(r, dict(r)), "equal_record")

    def test_the_inv08_rescue_reads_as_closer(self):
        self.assertEqual(
            audit.verdict_of(result("NC_038065.1", "SAME_CLASS"),
                             result("PP327218.1", "CONGENERIC")),
            "closer")

    def test_a_regression_reads_as_further(self):
        self.assertEqual(
            audit.verdict_of(result("PP327218.1", "CONGENERIC"),
                             result("NC_038065.1", "SAME_CLASS")),
            "further")

    def test_a_sibling_species_is_same_not_closer(self):
        self.assertEqual(
            audit.verdict_of(result("NC_038065.1", "CONFAMILIAL"),
                             result("MG655622.1", "CONFAMILIAL")),
            "same")

    def test_neither_side_selecting_is_not_an_unchanged_result(self):
        """INV08_TJALFIELLA aligns to no ctenophore record from either build. That is
        the absence of a measurement, not evidence the rebuild changed nothing."""
        none = result("", "UNKNOWN", state="NONE")
        self.assertEqual(audit.verdict_of(none, dict(none)), "no_selection")

    def test_a_rescue_from_nothing_still_reads_as_closer(self):
        self.assertEqual(
            audit.verdict_of(result("", "UNKNOWN", state="NONE"),
                             result("PP327218.1", "CONGENERIC")),
            "closer")

    def test_an_ungradeable_new_pick_never_reads_as_an_improvement(self):
        self.assertEqual(
            audit.verdict_of(result("NC_038065.1", "DISTANT"),
                             result("XX000000.1", "UNKNOWN")),
            "further")


class StatusParsing(unittest.TestCase):
    """The status line is select_reference_db.py's only machine-readable output."""

    def test_reference_mode_selected_line(self):
        line = "SELECTED\tPP327218.1 Tjalfiella sp. [Tjalfiellidae] cov=0.98 pid=97.4"
        path = Path(self.tmp) / "s.txt"
        path.write_text(line + "\n")
        state, detail = audit.parse_status(path)
        self.assertEqual(state, "SELECTED")
        self.assertEqual(detail.split()[0], "PP327218.1")
        self.assertEqual(audit.parse_metrics(detail), {"cov": "0.98", "pid": "97.4"})

    def test_low_confidence_still_parses_as_a_selection(self):
        detail = "NC_038065.1 Beroe forskalii [Beroidae] cov=0.31 pid=79.0"
        self.assertTrue("SELECTED_LOW_CONFIDENCE".startswith("SELECTED"))
        self.assertEqual(audit.parse_metrics(detail)["cov"], "0.31")

    def test_a_none_status_yields_no_metrics(self):
        self.assertEqual(audit.parse_metrics("empty assembly INV12.fasta"),
                         {"cov": "", "pid": ""})

    def test_missing_and_empty_status_files(self):
        missing = Path(self.tmp) / "absent.txt"
        self.assertEqual(audit.parse_status(missing), ("MISSING", ""))
        empty = Path(self.tmp) / "empty.txt"
        empty.write_text("")
        self.assertEqual(audit.parse_status(empty), ("MISSING", ""))

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)


class AssemblyDiscovery(unittest.TestCase):
    """A reseeded sample publishes its empty first pass beside the reseed that
    replaced it. Picking the wrong one audits an assembly the run discarded, and
    picking an empty one audits nothing at all."""

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def publish(self, sample, prefix, content):
        p = self.root / sample / prefix / "mtdna" / f"{prefix}.fasta"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
        return p

    def test_the_reseed_wins_over_an_empty_first_pass(self):
        self.publish("INV10_CHITON", "INV10_CHITON.ilmn.260902.getorg1770", "")
        reseed = self.publish("INV10_CHITON",
                              "INV10_CHITON.ilmn.260902.getorg1770reseed",
                              ">c\n" + "A" * 15000 + "\n")
        found, empty = audit.find_assemblies([self.root])
        self.assertEqual(found["INV10_CHITON"], reseed)
        self.assertNotIn("INV10_CHITON", empty)

    def test_a_sample_with_only_empty_passes_is_reported_not_dropped(self):
        self.publish("INV12_ANOMURA", "INV12_ANOMURA.ilmn.260902.getorg1770", "")
        self.publish("INV12_ANOMURA", "INV12_ANOMURA.ilmn.260902.getorg1770reseed", "")
        found, empty = audit.find_assemblies([self.root])
        self.assertNotIn("INV12_ANOMURA", found)
        self.assertIn("INV12_ANOMURA", empty)


if __name__ == "__main__":
    unittest.main()
