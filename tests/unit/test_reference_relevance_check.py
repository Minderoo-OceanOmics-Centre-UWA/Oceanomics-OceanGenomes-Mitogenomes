"""Unit tests for the reference grading in bin/reference_relevance_check.py.

Tests classify_relevance directly (pure function, no BLAST / biopython needed).
The numbers are the real refcov/pid values measured with dc-megablast over the
mitogenomes-missing-audit-5 cohort, which is what the thresholds were set from.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# bin/ on sys.path so the module's sibling import of species_name_utils resolves
# the same way it does in the task container (Nextflow puts bin/ on PATH there).
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "reference_relevance_check", ROOT / "bin" / "reference_relevance_check.py"
)
rrc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rrc)

VERT_PID = 82.0   # module default for vertebrates
CORAL_PID = 88.0  # module raises to this for invertebrates
MIN_COV = 0.70


class ClassifyRelevanceTests(unittest.TestCase):
    def state(self, ref_cov, mean_pid, congeneric=False, aln_bp=5000,
              min_cov=MIN_COV, min_pid=VERT_PID):
        return rrc.classify_relevance(ref_cov, mean_pid, aln_bp, congeneric,
                                      min_cov, min_pid)[0]

    # --- the false positives this rewrite exists to remove --------------------

    def test_congeneric_teleost_at_low_identity_is_not_a_mismatch(self):
        # OG54: Pseudolabrus biserialis vs P. eoethinus. Same genus, 87.5% identity
        # under the old coral-calibrated 88.0 floor -> was wrongly MISMATCH.
        self.assertEqual(self.state(0.95, 87.5, congeneric=True), "PASS")

    def test_congeneric_at_80_percent_is_divergent_not_mismatch(self):
        # OG834: Scorpaena sumptuosa vs S. agassizii, refcov 0.93 pid 80.0.
        self.assertEqual(self.state(0.93, 80.0, congeneric=True), "DIVERGENT")

    def test_well_covered_distant_relative_is_divergent(self):
        # OG65: 98% of the reference present at 76% identity. Right molecule,
        # distant relative -- the assembly is fine.
        self.assertEqual(self.state(0.98, 76.0), "DIVERGENT")

    def test_inflated_assembly_does_not_read_as_a_bad_reference(self):
        # OG852/OG853: control-region tandem repeat makes the assembly 19-22 kb, so
        # assembly-normalised coverage collapsed. Reference-normalised it is 0.95.
        self.assertEqual(self.state(0.95, 87.6, congeneric=True), "PASS")

    # --- the signal it must still catch ---------------------------------------

    def test_no_alignment_is_a_mismatch(self):
        self.assertEqual(self.state(0.0, 0.0, aln_bp=0), "MISMATCH")

    def test_wrong_family_reference_is_a_mismatch(self):
        # The coral case the check was built for: a wrong-family reference aligns
        # only in short low-identity patches.
        self.assertEqual(self.state(0.13, 81.5), "MISMATCH")

    def test_mismatch_requires_both_signals_to_fail(self):
        # Partly present but a good match -> a fragmented assembly, not a wrong
        # reference (OG829: a 2.4 kb assembly against a 16.8 kb reference at 86.9%).
        self.assertEqual(self.state(0.14, 86.9), "DIVERGENT")
        # Fully present but a poor match -> a distant relative, not a wrong one.
        self.assertEqual(self.state(1.00, 60.0), "DIVERGENT")

    # --- the congeneric veto ---------------------------------------------------

    def test_congeneric_veto_caps_at_divergent(self):
        # Same genus can never be the "wrong" reference, however bad the alignment.
        self.assertEqual(self.state(0.10, 60.0, congeneric=True), "DIVERGENT")

    def test_congeneric_veto_applies_when_nothing_aligns(self):
        self.assertEqual(self.state(0.0, 0.0, congeneric=True, aln_bp=0), "DIVERGENT")

    def test_non_congeneric_same_numbers_is_a_mismatch(self):
        self.assertEqual(self.state(0.10, 60.0, congeneric=False), "MISMATCH")

    # --- taxon-aware identity floor -------------------------------------------

    def test_coral_floor_flags_what_the_vertebrate_floor_passes(self):
        # 85% identity: normal for a teleost congener, too low for an anthozoan
        # (same-family corals sit at ~96.9%).
        self.assertEqual(self.state(0.95, 85.0, min_pid=VERT_PID), "PASS")
        self.assertEqual(self.state(0.95, 85.0, min_pid=CORAL_PID), "DIVERGENT")

    # --- boundaries ------------------------------------------------------------

    def test_thresholds_are_inclusive_lower_bounds(self):
        self.assertEqual(self.state(MIN_COV, VERT_PID), "PASS")
        self.assertEqual(self.state(MIN_COV - 0.01, VERT_PID), "DIVERGENT")
        self.assertEqual(self.state(MIN_COV, VERT_PID - 0.1), "DIVERGENT")


class SameGenusTests(unittest.TestCase):
    """The congeneric veto is only as good as the genus comparison behind it."""

    def test_matches_case_insensitively(self):
        self.assertTrue(rrc.same_genus("Scorpaena sumptuosa", "scorpaena agassizii"))

    def test_different_genus_does_not_match(self):
        self.assertFalse(rrc.same_genus("Vincentia punctata", "Phaeoptyx conklini"))

    def test_missing_sample_species_is_not_congeneric(self):
        # No label supplied -> no veto, so the BLAST evidence decides alone.
        self.assertFalse(rrc.same_genus("", "Scorpaena agassizii"))

    def test_open_nomenclature_name_still_yields_a_genus(self):
        self.assertTrue(rrc.same_genus("Hymenocephalus sp 4", "Hymenocephalus italicus"))


if __name__ == "__main__":
    unittest.main()
