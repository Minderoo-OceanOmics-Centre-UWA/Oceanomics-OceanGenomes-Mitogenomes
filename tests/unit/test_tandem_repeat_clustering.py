"""Unit tests for the tandem-repeat detector and over-length classification.

Needs blastn/makeblastdb on PATH; skipped otherwise.

Regression cover for NOVA_260909_LA. find_tandem_repeat used to reduce over EVERY
off-diagonal self-hit (min/max for the array bounds, the smallest offset as the
period), so two short unrelated repeats near opposite ends reported an "array"
spanning the whole molecule with an arbitrarily small period. Combined with a
concatemer rule that only looked at that array fraction and the length ratio, six of
seven over-length assemblies in that batch were handed a suggested_trim_region
overlapping annotated genes, and OG3004 (a 2.04x Munidopsis whose 13 CDS span all
30,972 bp) was called a concatemer against a 15,182 bp reference at 75% identity.
"""
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
SCRIPT = BIN / "check_getorganelle.py"
HAVE_BLAST = shutil.which("blastn") and shutil.which("makeblastdb")

sys.path.insert(0, str(BIN))


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_gb(path, seq):
    """Minimal GenBank; check_getorganelle reads the reference from ORIGIN only."""
    lines = ["LOCUS       REF   %d bp    DNA   circular   01-JAN-2026" % len(seq), "ORIGIN"]
    for i in range(0, len(seq), 60):
        lines.append("%9d %s" % (i + 1, seq[i:i + 60]))
    lines.append("//")
    path.write_text("\n".join(lines) + "\n")


@unittest.skipUnless(HAVE_BLAST, "blastn/makeblastdb not on PATH")
class FindTandemRepeatTests(unittest.TestCase):
    def find(self, seq, min_len=200, min_pident=70.0):
        from check_getorganelle import find_tandem_repeat
        with tempfile.TemporaryDirectory() as wd:
            return find_tandem_repeat(seq, wd, min_len, min_pident)

    def test_true_tandem_array_reports_its_period(self):
        unit = rand_seq(2000, seed=11)
        got = self.find(unit * 3)
        self.assertIsNotNone(got)
        self.assertAlmostEqual(got["unit_bp"], 2000, delta=100)
        self.assertGreaterEqual(got["array_len"], 5000)

    def test_scattered_repeats_do_not_inflate_the_array(self):
        """Two unrelated 300 bp repeats near opposite ends of a 30 kb molecule.

        The old reduction spanned the whole molecule; the clustered one must not.
        """
        rep_a, rep_b = rand_seq(300, seed=21), rand_seq(300, seed=22)
        body = rand_seq(30000, seed=23)
        seq = rep_a + body[:1000] + rep_b + body[1000:28000] + rep_a + body[28000:] + rep_b
        got = self.find(seq)
        if got is not None:
            self.assertLess(
                got["array_len"], 0.8 * len(seq),
                "array spans the molecule: the offset clustering is not holding",
            )


@unittest.skipUnless(HAVE_BLAST, "blastn/makeblastdb not on PATH")
class OverLengthClassificationTests(unittest.TestCase):
    def classify(self, seq, ref_seq):
        with tempfile.TemporaryDirectory() as wd:
            wd = Path(wd)
            fa, gb, ev = wd / "a.fa", wd / "ref.gb", wd / "ev.tsv"
            fa.write_text(">a\n%s\n" % seq)
            write_gb(gb, ref_seq)
            subprocess.run(
                [sys.executable, str(SCRIPT), "--fasta", str(fa), "--reference-gb", str(gb),
                 "--sample", "OGTEST", "--out-evidence", str(ev)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            rows = [l.split("\t") for l in ev.read_text().splitlines()]
            return dict(zip(rows[0], rows[1]))

    def test_real_dimer_is_still_a_concatemer_and_trims_the_second_copy(self):
        mono = rand_seq(16000, seed=31)
        out = self.classify(mono + mono, mono)
        self.assertEqual(out["anomaly_type"], "concatemer")
        start = int(out["suggested_trim_region"].split("-")[0])
        self.assertAlmostEqual(start, 16001, delta=1600,
                               msg="trim must start at the measured monomer period")

    def test_over_length_without_a_genome_period_is_not_a_concatemer(self):
        """2x the reference, but the excess is unique sequence, not a duplication.

        This is the OG3004 shape. It must not be called a concatemer, and must not
        propose any trim: at this stage there is no annotation to prove the region
        carries no genes.
        """
        ref = rand_seq(15000, seed=41)
        seq = ref + rand_seq(15000, seed=42)
        out = self.classify(seq, ref)
        self.assertNotEqual(out["anomaly_type"], "concatemer")
        self.assertEqual(out["suggested_trim_region"], "NA")

    def test_tandem_repeat_proposes_no_trim_without_gene_evidence(self):
        """A localised array in an over-length assembly is reported, never trimmed."""
        ref = rand_seq(16000, seed=51)
        vntr = rand_seq(200, seed=52) * 25
        out = self.classify(ref + vntr, ref)
        self.assertIn(out["anomaly_type"], ("tandem_repeat", "unresolved"))
        self.assertEqual(out["suggested_trim_region"], "NA")


if __name__ == "__main__":
    unittest.main()
