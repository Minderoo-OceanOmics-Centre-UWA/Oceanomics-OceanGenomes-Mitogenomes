"""Start codons: what the table lists, what ENA accepts, and what is a defect.

Three tables used to disagree about what a valid initiation codon is. orf_utils
rejected code-9 ATT/ATC starts that process_files.py was already handling
correctly with a transl_except, so 8 complete echinoderm assemblies failed QC on
starts the submission layer had no problem with. And 9 decapods failed on CGA,
the documented crustacean COI initiator, which neither table listed.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from orf_utils import (PLAUSIBLE_MITO_STARTS, START_CODONS_BY_CODE,
                       canonical_gene, classify_cds, lineage_alt_starts)


class Code9StartTests(unittest.TestCase):
    def test_code9_lists_the_starts_that_actually_occur(self):
        """Observed in this batch's asteroids: ND3 ATT x7, ND4L ATC x5,
        ND1 GTG x6, plus TTG. All were rejected by a two-codon table."""
        for codon in ("ATT", "ATC", "TTG", "ATG", "GTG"):
            self.assertIn(codon, START_CODONS_BY_CODE[9], codon)

    def test_code5_and_code9_agree(self):
        """Both are invertebrate mitochondrial tables; a divergent initiator set
        between them is not defensible."""
        self.assertEqual(set(START_CODONS_BY_CODE[5]), set(START_CODONS_BY_CODE[9]))


class StartClassTests(unittest.TestCase):
    def test_table_listed_start_is_canonical(self):
        self.assertEqual(classify_cds("ATTGGGTAA", 9)["start_class"], "canonical")

    def test_plausible_but_unlisted_start_is_documented_not_broken(self):
        """TTA initiates in code 4 but not code 2. process_files.py declares it
        with transl_except aa:Met, so QC must not call it a defect."""
        info = classify_cds("TTAGGGTAA", 2)
        self.assertEqual(info["start_class"], "alternative_documented")
        self.assertTrue(info["start_ok"])
        self.assertIn("TTA", PLAUSIBLE_MITO_STARTS)

    def test_codon_in_no_table_is_implausible(self):
        info = classify_cds("TTTGGGTAA", 4)
        self.assertEqual(info["start_class"], "implausible")
        self.assertFalse(info["start_ok"])

    def test_start_codon_is_reported(self):
        self.assertEqual(classify_cds("TTTGGGTAA", 4)["start_codon"], "TTT")


class LineageAltStartTests(unittest.TestCase):
    def test_decapod_co1_on_cga_is_documented(self):
        info = classify_cds("CGAGGGTAA", 5, lineage="Malacostraca", gene="CO1")
        self.assertEqual(info["start_class"], "alternative_documented")

    def test_same_codon_without_the_lineage_is_still_a_defect(self):
        """The exception is scoped. Omitting the lineage can only make the
        verdict stricter, never looser."""
        self.assertEqual(classify_cds("CGAGGGTAA", 5)["start_class"], "implausible")

    def test_exception_does_not_leak_to_other_genes(self):
        info = classify_cds("CGAGGGTAA", 5, lineage="Malacostraca", gene="ND5")
        self.assertEqual(info["start_class"], "implausible")

    def test_mitos_and_emma_gene_spellings_both_resolve(self):
        self.assertEqual(canonical_gene("cox1"), "CO1")
        self.assertEqual(lineage_alt_starts("Malacostraca", "cox1"),
                         lineage_alt_starts("Malacostraca", "CO1"))


class StopTests(unittest.TestCase):
    def test_terminal_stop_accepted(self):
        self.assertTrue(classify_cds("ATGGGGTAA", 4)["stop_ok"])

    def test_in_frame_sense_codon_is_not_a_stop(self):
        """The residual Porifera defects all look like this: in frame (rem 0) but
        ending on a sense codon, which table2asn raises NoStop for and
        compute_transl_except_pos correctly declines to rescue."""
        self.assertFalse(classify_cds("ATGGGGAGA", 4)["stop_ok"])

    def test_polyadenylation_completed_stop_accepted(self):
        self.assertTrue(classify_cds("ATGGGGT", 4)["stop_ok"])


if __name__ == "__main__":
    unittest.main()
