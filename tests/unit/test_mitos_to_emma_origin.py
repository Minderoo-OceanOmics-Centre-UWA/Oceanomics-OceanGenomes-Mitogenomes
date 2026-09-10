"""Unit tests for the published-origin re-rotation in mitos_to_emma.py.

This code had NO tests before, which is how it shipped re-origining every
invertebrate to tRNA-Met while its docstring claimed that matched the NCBI coral
convention. Measured against the repo's own reference databases, trnM is the
deposited origin for 0-2% of most invertebrate phyla; it is the convention only
for Scleractinia (63.8%), which is where the rule came from.

The load-bearing test here is the round-trip: re-origining must move coordinates
without changing a single gene's sequence. That property covers the composition
of the reverse-complement mirror and the rotation, which is the part most likely
to be subtly wrong and least likely to be noticed.

Biopython is imported at module scope by mitos_to_emma.py and isn't installed in
every environment these tests run in, so skip cleanly when it is absent.
"""

import importlib.util
import sys
import unittest
import unittest.mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

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


def _exon(start, end, strand="+"):
    return {"start": start, "end": end, "strand": strand}


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class AnchorOriginTests(unittest.TestCase):
    def test_plus_strand_returns_start_and_no_rc(self):
        features = {"CO1": [_exon(101, 200)]}
        self.assertEqual(mitos_to_emma.anchor_origin(features, "CO1"), (101, False))

    def test_minus_strand_returns_end_and_requests_rc(self):
        # On the minus strand the 5' end is the HIGHER genomic coordinate.
        features = {"TM": [_exon(101, 170, "-")]}
        self.assertEqual(mitos_to_emma.anchor_origin(features, "TM"), (170, True))

    def test_multi_exon_anchors_on_the_first_exon_not_the_largest(self):
        # cox1 is intron-split in some sponge families and rrnL can be called in
        # pieces. The origin belongs at the start of the GENE, so exons[0] (the
        # 5'-most piece in MITOS transcript order) wins even though exon 2 is longer.
        features = {"CO1": [_exon(50, 100), _exon(400, 900)]}
        self.assertEqual(mitos_to_emma.anchor_origin(features, "CO1"), (50, False))

    def test_missing_anchor_returns_none(self):
        features = {"CO1": [_exon(1, 100)]}
        self.assertIsNone(mitos_to_emma.anchor_origin(features, "RNR2"))
        self.assertIsNone(mitos_to_emma.anchor_origin({}, "TM"))


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class ReoriginRoundTripTests(unittest.TestCase):
    """Re-origining must be a pure coordinate change: every gene keeps its sequence.

    This is the property the whole transform has to hold, and it is what makes the
    cds/ and proteins/ outputs rotation-invariant.
    """

    SEQ_LEN = 200

    def setUp(self):
        from Bio.Seq import Seq
        import random
        rng = random.Random(20260910)
        self.seq = Seq("".join(rng.choice("ACGT") for _ in range(self.SEQ_LEN)))
        self.features = {
            "CO1":  [_exon(10, 60)],
            "RNR2": [_exon(70, 120, "-")],
            "TF":   [_exon(130, 160)],
            # A feature that already wraps the origin, to prove the wrap survives.
            "TM":   [_exon(190, 200)],
        }

    def _gene_seqs(self, seq, features):
        """Spliced nucleotide sequence of every gene, exon by exon.

        gene_nt_seq takes a {chrom: SeqRecord} dict, which is how the script holds
        the genome, and already handles an origin-spanning feature (start > end).
        """
        from Bio.SeqRecord import SeqRecord
        genome = {"c": SeqRecord(seq, id="c")}
        out = {}
        for name, exons in features.items():
            pieces = [
                str(mitos_to_emma.gene_nt_seq(genome, "c", e["start"], e["end"], e["strand"]))
                for e in exons
            ]
            out[name] = "".join(pieces)
        return out

    def test_every_anchor_preserves_every_gene_sequence(self):
        before = self._gene_seqs(self.seq, self.features)
        for anchor in ("CO1", "RNR2", "TF", "TM"):
            with self.subTest(anchor=anchor):
                pos, rc = mitos_to_emma.anchor_origin(self.features, anchor)
                rotated = mitos_to_emma.reorigin_seq(self.seq, pos, rc, self.SEQ_LEN)
                shifted = mitos_to_emma.reorigin_features(
                    self.features, pos, rc, self.SEQ_LEN)
                self.assertEqual(len(rotated), self.SEQ_LEN)
                self.assertEqual(self._gene_seqs(rotated, shifted), before)

    def test_the_anchor_lands_at_position_one(self):
        for anchor in ("CO1", "RNR2", "TF"):
            with self.subTest(anchor=anchor):
                pos, rc = mitos_to_emma.anchor_origin(self.features, anchor)
                shifted = mitos_to_emma.reorigin_features(
                    self.features, pos, rc, self.SEQ_LEN)
                self.assertEqual(shifted[anchor][0]["start"], 1)

    def test_minus_strand_anchor_flips_every_strand(self):
        # RNR2 is on the minus strand, so the whole molecule is reverse-complemented
        # and every feature's strand flips with it.
        pos, rc = mitos_to_emma.anchor_origin(self.features, "RNR2")
        self.assertTrue(rc)
        shifted = mitos_to_emma.reorigin_features(self.features, pos, rc, self.SEQ_LEN)
        self.assertEqual(shifted["RNR2"][0]["strand"], "+")
        self.assertEqual(shifted["CO1"][0]["strand"], "-")

    def test_trnm_reproduces_the_pre_change_behaviour(self):
        """Regression pin: --origin-gene TM must be exactly what the old code did.

        The generalisation has to be a strict superset of the trnM-only version,
        because Scleractinia still depends on that exact behaviour.
        """
        pos, rc = mitos_to_emma.anchor_origin(self.features, "TM")
        # Old trnmet_origin() body, inlined.
        e = self.features["TM"][0]
        expected = (e["end"], True) if e["strand"] == "-" else (e["start"], False)
        self.assertEqual((pos, rc), expected)


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class AnchorLabelTests(unittest.TestCase):
    def test_labels_are_readable(self):
        self.assertEqual(mitos_to_emma.anchor_label("CO1"), "cox1")
        self.assertEqual(mitos_to_emma.anchor_label("RNR2"), "rrnL (16S)")
        self.assertIn("Met", mitos_to_emma.anchor_label("TM"))
        self.assertIn("Phe", mitos_to_emma.anchor_label("TF"))

    def test_unknown_key_falls_back_to_itself(self):
        self.assertEqual(mitos_to_emma.anchor_label("ND5"), "ND5")


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class OriginGeneIsRequiredTests(unittest.TestCase):
    """--origin-gene has no default, on purpose.

    There are two call sites (MITOS2 and CORAL_ANNOTATION_FIX, one per sample). A
    default would let one be updated without the other, silently publishing FIX and
    PASS anthozoans of the same species on different origins. Requiring the flag
    turns that into a loud failure.
    """

    def test_parser_rejects_a_missing_origin_gene(self):
        import contextlib
        import io

        # main() reads sys.argv directly, so drive it that way.
        argv = ["mitos_to_emma.py",
                "--bed", "b.bed", "--genome", "g.fa", "--prefix", "p",
                "--outdir", "o", "--code", "5"]
        stderr = io.StringIO()
        with self.assertRaises(SystemExit) as ctx, \
                contextlib.redirect_stderr(stderr), \
                unittest.mock.patch.object(sys, "argv", argv):
            mitos_to_emma.main()
        self.assertNotEqual(ctx.exception.code, 0)
        self.assertIn("origin-gene", stderr.getvalue())

    def test_parser_accepts_an_explicit_origin_gene(self):
        # Sanity check the opposite direction: with the flag present the parser gets
        # past argparse and fails later, on the missing BED, not on the arguments.
        import contextlib
        import io

        argv = ["mitos_to_emma.py",
                "--bed", "missing.bed", "--genome", "g.fa", "--prefix", "p",
                "--outdir", "o", "--code", "5", "--origin-gene", "TM"]
        stderr = io.StringIO()
        with self.assertRaises(SystemExit) as ctx, \
                contextlib.redirect_stderr(stderr), \
                unittest.mock.patch.object(sys, "argv", argv):
            mitos_to_emma.main()
        self.assertNotIn("origin-gene", stderr.getvalue())
        self.assertIn("BED not found", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
