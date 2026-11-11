"""bin/check_getorganelle.py: a tiny GetOrganelle "circle" is a fragment, not circular.

OG3000's reseed published a 155 bp "circular genome" covering 2.1% of the
16,107 bp Glyptelasma reference. The check kept GetOrganelle's True verdict
unconditionally, so final_verdict_circular said True. BLAST and the repeat scan
are patched out: these tests pin the verdict logic, not the aligner.
"""

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("check_getorganelle", ROOT / "bin" / "check_getorganelle.py")
cg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cg)

REF_LEN = 16107


class CheckGetorganelleVerdictTests(unittest.TestCase):
    def verdict(self, length, getorg_circular, ref_cov, records=1):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = Path(tmp.name)
        fasta = d / "asm.fasta"
        fasta.write_text("".join(f">s{i}\n{'A' * length}\n" for i in range(records)))
        gb = d / "ref.gb"
        gb.write_text("LOCUS ref\nORIGIN\n        1 acgt\n//\n")
        out = d / "ev.tsv"
        argv = ["check_getorganelle.py", "--fasta", str(fasta), "--reference-gb", str(gb),
                "--sample", "OG3000.ilmn.260909.getorg1770reseed",
                "--getorg-circular", getorg_circular, "--out-evidence", str(out)]
        with mock.patch.object(cg, "gb_to_fasta", return_value=REF_LEN), \
             mock.patch.object(cg, "blast_reference_coverage", return_value=(ref_cov, 0)), \
             mock.patch.object(cg, "find_tandem_repeat", return_value=None), \
             mock.patch.object(sys, "argv", argv):
            self.assertEqual(cg.main(), 0)
        with open(out) as fh:
            return next(csv.DictReader(fh, delimiter="\t"))

    def test_og3000_tiny_circle_is_a_fragment(self):
        row = self.verdict(155, "true", 0.021)
        self.assertEqual(row["getorg_circular"], "True")
        self.assertEqual(row["final_verdict_circular"], "False")
        self.assertEqual(row["circular_corrected"], "yes")
        self.assertEqual(row["note"], "circular_fragment")

    def test_real_circle_stays_circular(self):
        row = self.verdict(16262, "true", 0.98)
        self.assertEqual(row["final_verdict_circular"], "True")
        self.assertEqual(row["circular_corrected"], "no")
        self.assertEqual(row["note"], "ok")

    def test_divergent_full_length_circle_stays_circular(self):
        # Low coverage alone is a distant reference, not a fragment.
        row = self.verdict(16000, "true", 0.30)
        self.assertEqual(row["final_verdict_circular"], "True")

    def test_short_circle_with_good_coverage_stays_circular(self):
        # Short alone is not enough either: both conditions must hold.
        row = self.verdict(7000, "true", 0.60)
        self.assertEqual(row["final_verdict_circular"], "True")

    def test_reference_confirmed_scaffold_still_corrected_to_circular(self):
        row = self.verdict(16100, "false", 0.99)
        self.assertEqual(row["final_verdict_circular"], "True")
        self.assertEqual(row["circular_corrected"], "yes")
        self.assertEqual(row["note"], "ok")


if __name__ == "__main__":
    unittest.main()
