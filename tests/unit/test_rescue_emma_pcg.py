"""Unit tests for bin/rescue_emma_pcg.py.

The ORF-refinement and coordinate logic is pure Python and always tested. The
end-to-end BLAST path needs Biopython + BLAST+ (the MITOS2 BioContainer) and is
skipped when either is missing.
"""

import importlib.util
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

try:
    import Bio  # noqa: F401
    BIOPYTHON = True
except ImportError:
    BIOPYTHON = False

TBLASTN = shutil.which("tblastn") is not None

def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# rescue_emma_pcg.py's ORF refinement now lives in bin/orf_utils.py, which is
# stdlib-only, so the ORF tests below run everywhere -- they no longer need the
# Biopython container just to exercise the boundary logic the rescue depends on.
sys.path.insert(0, str(ROOT / "bin"))
orf = _load("orf_utils")

rescue = None
if BIOPYTHON:
    rescue = _load("rescue_emma_pcg")


# A plausible vertebrate-mito ND4L: ATG start, 96 sense codons, TAA stop.
ND4L_BODY = "ATG" + ("ACC" * 96) + "TAA"          # 294 nt -> M + 96*T


class RefineOrfTests(unittest.TestCase):
    def test_clean_orf_from_internal_hsp(self):
        window = "GGGGGG" + ND4L_BODY + "ATGCTATG"   # 6 nt 5' pad, junk 3'
        # tblastn would land a few codons in; refine walks back to the ATG.
        res = orf.refine_orf(window, s_from=6 + 13, down_limit=len(window),
                                code=2)
        self.assertIsNotNone(res)
        cds_lo, cds_hi, aa, start_codon, poly_a = res
        self.assertEqual(start_codon, "ATG")
        self.assertFalse(poly_a)
        self.assertEqual(aa, "M" + "T" * 96)
        # CDS span (1-based inclusive) covers the ATG..TAA, stop included.
        self.assertEqual(window[cds_lo - 1:cds_hi], ND4L_BODY)

    def test_orf_ends_at_first_in_frame_stop(self):
        body = "ATG" + "ACC" * 10 + "TAA" + "ACC" * 20 + "TAA"
        window = "GG" + body
        res = orf.refine_orf(window, s_from=3, down_limit=len(window), code=2)
        self.assertIsNotNone(res)
        _lo, _hi, aa, _sc, poly_a = res
        self.assertFalse(poly_a)
        self.assertNotIn("*", aa)
        self.assertEqual(aa, "M" + "T" * 10)     # stops at the first TAA

    def test_polyadenylation_when_stop_runs_into_downstream_gene(self):
        # No stop codon before the downstream-gene boundary -> truncate + poly-A.
        body = "ATG" + "ACC" * 60           # 183 nt, no stop
        window = "GG" + body + "GGGGGG"
        down_limit = 2 + len(body)           # boundary right after the ORF bases
        res = orf.refine_orf(window, s_from=3, down_limit=down_limit, code=2)
        self.assertIsNotNone(res)
        _lo, _hi, aa, _sc, poly_a = res
        self.assertTrue(poly_a)
        self.assertEqual(aa, "M" + "T" * 60)

    def test_upstream_gtg_does_not_beat_the_true_atg(self):
        """Regression for OG811 (batch-11) ATP8.

        The genomic context there is ...TTG GTG AAA ATG CCT CAA...: an in-frame
        GTG sits two codons upstream of the real ATG. The old local refine_orf
        walked upstream keeping the FURTHEST initiator within 12 codons, so it
        started on the GTG and published a 57 aa VKMPQLNP... ATP8 instead of the
        canonical 55 aa MPQLNP... . orf_utils.refine_orf ranks a canonical ATG
        above an alternative initiator, so the true start wins from any in-frame
        HSP position.
        """
        body = "ATG" + "CCTCAATTAAACCCAACC" + "ACC" * 34 + "TAA"
        window = "CAACCACCCTTGGTGAAA" + body + "GGGGGG"
        atg_at = window.index(body) + 1        # 1-based
        for hsp_offset in (0, 3, 6, 9):        # wherever tblastn lands, in frame
            res = orf.refine_orf(window, s_from=atg_at + hsp_offset,
                                    down_limit=len(window), code=2)
            self.assertIsNotNone(res)
            cds_lo, _hi, aa, start_codon, _p = res
            self.assertEqual(start_codon, "ATG", f"hsp_offset={hsp_offset}")
            self.assertEqual(cds_lo, atg_at, f"hsp_offset={hsp_offset}")
            self.assertTrue(aa.startswith("MPQ"), aa[:6])

    def test_five_prime_partial_when_no_start(self):
        # HSP frame has no ATG upstream within reach -> flagged 'partial'.
        window = ("ACC" * 40) + "TAA" + "GGG"
        res = orf.refine_orf(window, s_from=61, down_limit=len(window), code=2)
        self.assertIsNotNone(res)
        _lo, _hi, _aa, start_codon, _p = res
        self.assertEqual(start_codon, "partial")


@unittest.skipUnless(BIOPYTHON, "Biopython not installed")
class EvalueFormatTests(unittest.TestCase):
    def test_emma_style_evalue(self):
        self.assertEqual(rescue.emma_evalue(6.5e-51), "6.5e-51")
        self.assertEqual(rescue.emma_evalue(1e-30), "1.0e-30")


@unittest.skipUnless(BIOPYTHON and TBLASTN, "needs Biopython + tblastn")
class EndToEndTests(unittest.TestCase):
    """Build a tiny genome with a real ND4L flanked by TR and ND4 markers,
    a matching GFF, and confirm rescue_gene() splices the feature in."""

    def setUp(self):
        import tempfile
        self.d = Path(tempfile.mkdtemp())
        ann = self.d / "annotation"
        (ann / "cds").mkdir(parents=True)
        (ann / "proteins").mkdir(parents=True)

        tr = "GG" * 35                       # 70 nt stand-in for tRNA-Arg
        nd4 = "ATG" + "CTT" * 60             # ND4 start marker
        spacer = "AA"
        genome = spacer + tr + ND4L_BODY + nd4 + spacer
        tr_start = len(spacer) + 1
        tr_end = tr_start + len(tr) - 1
        nd4l_start = tr_end + 1
        nd4_start = nd4l_start + len(ND4L_BODY)
        nd4_end = nd4_start + len(nd4) - 1

        (ann / "og.test.1.emma.gff").write_text(
            "##gff-version 3\n"
            f"##sequence-region og 1 {len(genome)}\n"
            f"og\tEmma\tgene\t{tr_start}\t{tr_end}\t.\t+\t.\tID=a;Name=MT-TR\n"
            f"og\tEmma\tgene\t{nd4_start}\t{nd4_end}\t.\t+\t.\tID=b;Name=MT-ND4\n")
        (ann / "og.test.1.emma.tbl").write_text(">Feature og\n")
        (ann / "og.test.1.emma.fa").write_text(f">og\n{genome}\n")

        # Reference: the exact ND4L protein (translate the body).
        from Bio.Seq import Seq
        prot = str(Seq(ND4L_BODY).translate(table=2)).rstrip("*")
        self.ref = self.d / "ref.faa"
        self.ref.write_text(f">ND4L_Test_species__NC_000000.1\n{prot}\n")
        self.ann = ann

    def test_nd4l_is_spliced_in(self):
        genes = rescue.parse_gff_genes(next(self.ann.glob("*.gff")))
        from Bio import SeqIO
        with open(next(self.ann.glob("*.fa"))) as fh:
            rec = next(SeqIO.parse(fh, "fasta"))
            seq, rid, slen = str(rec.seq).upper(), rec.id, len(rec.seq)
        state, gene, msg = rescue.rescue_gene(
            "ND4L", genes, seq, rid, slen,
            self.ann, "og.test.1.emma", self.ref, 2, 55.0, 0.75, 0.15)
        self.assertEqual((state, gene), ("RESCUED", "ND4L"), msg)
        self.assertIn("MT-ND4L", next(self.ann.glob("*.gff")).read_text())
        self.assertTrue((self.ann / "proteins" / "MT-ND4L.og.test.1.emma.fa").exists())


if __name__ == "__main__":
    unittest.main()
