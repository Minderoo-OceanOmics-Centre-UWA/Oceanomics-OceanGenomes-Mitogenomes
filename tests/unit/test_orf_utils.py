"""Unit tests for bin/orf_utils.py.

Pure stdlib -- no Biopython, no BLAST -- so these always run:
    python3 -m pytest tests/unit/test_orf_utils.py -v
"""

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "orf_utils", ROOT / "bin" / "orf_utils.py")
orf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(orf)

# The real OG2350 (Montipora grisea) ND1 region, EMMA-prefix genome coordinates.
# MITOS annotated ND1 as 3503-4450 (1-based); the spliced CDS below begins on AAA
# (no initiator under code 4) and ends on CCG (no stop), 948 nt / 316 codons, no
# internal stop -- exactly the mis-placed boundary that slipped past the old gate
# and blocked table2asn with SEQ_INST.BadProteinStart / SEQ_FEAT.NoStop.
OG2350_ND1_CDS = (
    "AAAATATTAATGGTCATAGTTCCATTGCTTATCACAGTAGCTTATTTAACTTTAGCCGAACGAAAGGTTTTAGGA"
    "TATATGCAGGCTAGAAAAGGGCCAAATGTAGTGGGGGTTAGTGGGCTGGCTCAGCCGTTTGCAGATGGCCTAAAA"
    "CTATTTACCAAAGAGATGGTGGTTCCGCATCAAACTAATTTGTTTATTTATATAGTGGCGCCGGTGTTTTCCTTT"
    "ATCTTGGCGTTAATCGTTTGAGGGGTGGTGCCTTATGAGAGAGGGGCTTTAATAAGTGATTTAAAAATAGGGGTT"
    "TTGTATATCTTGGCTGTTTCTTCGATAAGTGTCTATGCGATATTAATGTCTGGGTGGGCGAGTAATTCTAAATAT"
    "GCTTTTTTGGGGGCGATTCGGGCAGCTGCTCAAATGATTAGTTATGAAGTTTCGATTGGGCTGATCATCATTTCG"
    "GTTCTATTGTGTGTTGGCTCTTTGAGTGTTACTGAAATTGTTTTAGCGCAGAATAGTGGTGTTTGGTTTTTTTTT"
    "CCGTTATTTCCTGTTGCAATAATGTTTTTTGTTTCAGCTTTGGCGGAGACAAATCGAGCTCCCTTTGATTTAACA"
    "GAAGGAGAGTCGGAGCTAGTGTCGGGGTATAACGTAGAGTACGCCTCGATGTCTTTTGCTCTATTTTTCCTTGCT"
    "GAATATGCTCATATAATATTAATGAGTTGTTTAACAATTATCTTTTTTGGGGGGGGGTGACTTTCTCCAGTAAAA"
    "TATTTAAAAGGCGGGGCCGGGTGGTTTGGTCTAAAAGTTGTTTTAATAATTTTTCTTTTTATTTGGGTAAGGGCC"
    "TCTTTTCCTCGTATTCGCTACGACCAGCTTATGTCTTTGTTATGAAAGGCGTATTTACCTCTAAGTTTAGGGGTG"
    "GTGACTTTTGTGGCTAGTGTTCTTTTTGGGTTAAATGGCCCGCCTCCG"
)

# The same region with 33 nt of upstream flank and enough 3' flank to reach the
# true stop, as refine_orf sees it. The true ORF is ATG at +9 (relative to the
# MITOS start), running to TAA 6 nt past the MITOS end -> 314 aa.
OG2350_ND1_WINDOW = (
    "AAGTGTTCGTGGAAATAATTCATTTACTTTTT"
    + OG2350_ND1_CDS
    + "ATCTAAGACATTTTGTAATTTGTTGTTTGA"
)


class CodonTableTests(unittest.TestCase):
    def test_supported_codes_have_both_tables(self):
        for code in (2, 4, 5, 9, 13, 14):
            self.assertTrue(orf.start_codons(code))
            self.assertTrue(orf.stop_codons(code))

    def test_code4_stops_exclude_aga_agg(self):
        self.assertEqual(set(orf.stop_codons(4)), {"TAA", "TAG"})

    def test_code2_stops_include_aga_agg(self):
        self.assertIn("AGA", orf.stop_codons(2))

    def test_code14_only_stop_is_tag(self):
        self.assertEqual(set(orf.stop_codons(14)), {"TAG"})

    def test_translate_code4_reads_tga_as_trp(self):
        self.assertEqual(orf.translate("TGA", 4), "W")

    def test_translate_code9_reads_aaa_as_asn(self):
        self.assertEqual(orf.translate("AAA", 9), "N")


class ClassifyCdsTests(unittest.TestCase):
    def test_clean_atg_orf_passes(self):
        info = orf.classify_cds("ATGAAAGGGTAA", 4)
        self.assertTrue(info["start_ok"])
        self.assertTrue(info["stop_ok"])
        self.assertEqual(info["internal_stops"], 0)

    def test_gtg_alt_start_is_ok_under_code4(self):
        self.assertTrue(orf.classify_cds("GTGAAATAA", 4)["start_ok"])

    def test_non_start_first_codon_flagged(self):
        self.assertFalse(orf.classify_cds("AAAGGGTAA", 4)["start_ok"])

    def test_no_terminal_stop_flagged(self):
        self.assertFalse(orf.classify_cds("ATGAAAGGG", 4)["stop_ok"])

    def test_truncated_codon_counts_as_polyA_stop(self):
        # Trailing partial codon is the prefix of a stop (TAA/TAG under code 4),
        # i.e. a genuine poly-A-completed stop -> accepted.
        self.assertTrue(orf.classify_cds("ATGAAAT", 4)["stop_ok"])
        self.assertTrue(orf.classify_cds("ATGAAATA", 4)["stop_ok"])

    def test_truncated_codon_that_is_not_a_stop_prefix_is_flagged(self):
        # "GG" cannot become TAA/TAG by adding 3' A residues, so a length that is
        # not a multiple of 3 here is a frameshift or a mis-called boundary, not
        # the poly-A convention -- it must not be waved through.
        self.assertFalse(orf.classify_cds("ATGAAAGG", 4)["stop_ok"])
        # Under code 2 the same applies, and AGA/AGG being stops there means "AG"
        # IS a valid truncated-stop prefix.
        self.assertFalse(orf.classify_cds("ATGAAAGG", 2)["stop_ok"])
        self.assertTrue(orf.classify_cds("ATGAAAAG", 2)["stop_ok"])

    def test_internal_stop_counted(self):
        self.assertEqual(orf.classify_cds("ATGTAAGGGTAA", 4)["internal_stops"], 1)

    def test_og2350_nd1_is_all_fail(self):
        info = orf.classify_cds(OG2350_ND1_CDS, 4)
        self.assertFalse(info["start_ok"])
        self.assertFalse(info["stop_ok"])
        self.assertEqual(info["internal_stops"], 0)
        self.assertEqual(info["length_nt"], 948)


class RefineOrfTests(unittest.TestCase):
    def test_og2350_nd1_snaps_to_true_atg(self):
        # MITOS start sits at 1-based 33 in the window (32 nt of 5' flank).
        res = orf.refine_orf(OG2350_ND1_WINDOW, 33, len(OG2350_ND1_WINDOW), 4)
        self.assertIsNotNone(res)
        start, end, aa, start_codon, poly_a = res
        self.assertEqual(start_codon, "ATG")
        self.assertEqual(start, 42)          # +9 nt from the MITOS start
        self.assertEqual(len(aa), 314)
        self.assertTrue(aa.startswith("M"))
        self.assertFalse(poly_a)
        self.assertNotIn("*", aa)

    def test_code9_orf_uses_code9_start_set(self):
        # code 9 initiators are ATG/GTG only; ATT must not be taken as a start.
        window = "AAAATTAAAATGAAAGGGTAA"
        res = orf.refine_orf(window, 4, len(window), 9)
        self.assertIsNotNone(res)
        self.assertEqual(res[3], "ATG")


class UnknownCodeTests(unittest.TestCase):
    """An unrecognised table must raise, never fall back.

    The old fallback returned ("TAA","TAG") while its comment claimed the
    vertebrate code -- but code 2 also stops on AGA/AGG, so a mistyped table
    silently got the WRONG stop set and produced plausible-looking ORFs. Every
    caller now resolves the code from meta.genetic_code, so an unknown value is a
    wiring bug and has to surface as one.
    """

    def test_unsupported_code_raises(self):
        for code in (1, 7, 99):
            with self.subTest(code=code):
                self.assertRaises(ValueError, orf.stop_codons, code)
                self.assertRaises(ValueError, orf.start_codons, code)

    def test_non_integer_code_raises(self):
        self.assertRaises(ValueError, orf.stop_codons, None)
        self.assertRaises(ValueError, orf.start_codons, "vertebrate")

    def test_code_2_keeps_its_four_stops(self):
        self.assertEqual(set(orf.stop_codons(2)), {"TAA", "TAG", "AGA", "AGG"})


if __name__ == "__main__":
    unittest.main()
