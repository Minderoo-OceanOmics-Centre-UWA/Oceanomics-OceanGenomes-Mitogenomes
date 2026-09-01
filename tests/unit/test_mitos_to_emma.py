"""Unit tests for the .tbl writer in mitos_to_emma.py.

mitos_to_emma.py imports Biopython at module scope (it adapts MITOS2 output
using Bio.SeqIO/Seq), which isn't installed in every environment these unit
tests run in (e.g. a bare login-node python3 without the pipeline's
containers). Skip cleanly rather than failing when Biopython is absent.
"""

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

try:
    import Bio  # noqa: F401
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

if BIOPYTHON_AVAILABLE:
    SPEC = importlib.util.spec_from_file_location(
        "mitos_to_emma", ROOT / "bin" / "mitos_to_emma.py"
    )
    mitos_to_emma = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(mitos_to_emma)


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class MapGeneNameAnticodonTests(unittest.TestCase):
    def test_trna_anticodon_captured_as_rna(self):
        emma_name, ftype, frag, anticodon = mitos_to_emma.map_gene_name("trnW(tca)")
        self.assertEqual((emma_name, ftype, frag, anticodon), ("TW", "tRNA", None, "UCA"))

    def test_cds_has_no_anticodon(self):
        emma_name, ftype, frag, anticodon = mitos_to_emma.map_gene_name("cox1")
        self.assertEqual((emma_name, ftype, frag, anticodon), ("CO1", "CDS", None, None))

    def test_split_gene_fragment_still_parsed(self):
        emma_name, ftype, frag, anticodon = mitos_to_emma.map_gene_name("nad5-a")
        self.assertEqual((emma_name, ftype, frag), ("ND5", "CDS", "a"))


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class TrnaProductTests(unittest.TestCase):
    def test_known_suffix_with_anticodon(self):
        self.assertEqual(mitos_to_emma.trna_product("TW", "UCA"), "tRNA-Trp(UCA)")

    def test_arm_suffix(self):
        self.assertEqual(mitos_to_emma.trna_product("TS1", "GCU"), "tRNA-Ser(GCU)")

    def test_falls_back_when_unrecognised(self):
        self.assertEqual(mitos_to_emma.trna_product("TQQ", None), "tRNA-TQQ")


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class TblIntervalsTests(unittest.TestCase):
    def _exon(self, start, end, strand):
        return {"start": start, "end": end, "strand": strand}

    def test_plus_strand_simple(self):
        exons = [self._exon(71, 110, "+")]
        self.assertEqual(mitos_to_emma._tbl_intervals(exons, 126), [(71, 110)])

    def test_minus_strand_simple_is_high_low(self):
        # Mirrors real EMMA .tbl output: a minus-strand feature is written
        # high..low with no separate strand column.
        exons = [self._exon(13797, 14318, "-")]
        self.assertEqual(mitos_to_emma._tbl_intervals(exons, 20000), [(14318, 13797)])

    def test_plus_strand_origin_wrap_splits_in_genomic_order(self):
        exons = [self._exon(110, 9, "+")]
        self.assertEqual(mitos_to_emma._tbl_intervals(exons, 126), [(110, 126), (1, 9)])

    def test_minus_strand_origin_wrap_reverses_and_swaps(self):
        exons = [self._exon(111, 10, "-")]
        self.assertEqual(mitos_to_emma._tbl_intervals(exons, 126), [(10, 1), (126, 111)])

    def test_multi_exon_kept_in_transcript_order(self):
        exons = [self._exon(1, 50, "+"), self._exon(80, 120, "+")]
        self.assertEqual(mitos_to_emma._tbl_intervals(exons, 200), [(1, 50), (80, 120)])


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class WriteTblTests(unittest.TestCase):
    def _features(self):
        return {
            "CO1": [{"start": 71, "end": 110, "strand": "+", "ftype": "CDS", "frag": None}],
            "TW": [{"start": 20, "end": 29, "strand": "-", "ftype": "tRNA", "frag": None,
                    "anticodon": "UCA"}],
            "RNR1": [{"start": 1, "end": 10, "strand": "+", "ftype": "rRNA", "frag": None}],
        }

    def test_write_tbl_shape(self, tmp_path=None):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            tbl_path = Path(d) / "sample.tbl"
            mitos_to_emma.write_tbl(tbl_path, self._features(), "chr1", 126, "sample", 5)
            lines = tbl_path.read_text().splitlines()

        self.assertEqual(lines[0], ">Feature sample")
        joined = "\n".join(lines)
        # CDS emits gene -> mRNA -> CDS with product + transl_table, no protein_id
        # (process_files.py strips protein_id anyway; omitted here by design).
        self.assertIn("71\t110\tgene", joined)
        self.assertIn("71\t110\tmRNA", joined)
        self.assertIn("71\t110\tCDS", joined)
        self.assertIn("\t\t\ttransl_table\t5", joined)
        self.assertNotIn("protein_id", joined)
        # Minus-strand tRNA written high..low with its /product.
        self.assertIn("29\t20\ttRNA", joined)
        self.assertIn("\t\t\tproduct\ttRNA-Trp(UCA)", joined)
        # rRNA product comes from the shared PRODUCT map.
        self.assertIn("\t\t\tproduct\t12S ribosomal RNA", joined)


if __name__ == "__main__":
    unittest.main()
