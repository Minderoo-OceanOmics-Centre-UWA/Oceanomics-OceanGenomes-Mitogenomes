"""Unit tests for bin/collapse_concatemer.py.

Exercises the concatemer collapse end-to-end (needs blastn/makeblastdb on PATH;
skipped otherwise). Verifies a genuine head-to-tail dimer collapses to a monomer,
that a non-concatemer passes through untouched, and that an over-length chimera
mislabelled as a concatemer fails closed (is not corrupted).
"""
import csv
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "collapse_concatemer.py"
HAVE_BLAST = shutil.which("blastn") and shutil.which("makeblastdb")

HDR = "sample\tanomaly_type\tsuggested_trim_region\treference_length\tfinal_verdict_circular"


def monomer(n=16700, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@unittest.skipUnless(HAVE_BLAST, "blastn/makeblastdb not on PATH")
class CollapseConcatemerTests(unittest.TestCase):
    def run_collapse(self, seq, evidence_line, sample="OG"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        fasta = root / "in.fa"
        fasta.write_text(f">{sample}_ctg\n{seq}\n")
        ev = root / "ev.tsv"
        ev.write_text(HDR + "\n" + evidence_line + "\n")
        out = root / "out.fa"
        rep = root / "rep.tsv"
        subprocess.run(
            ["python3", str(SCRIPT), "--fasta", str(fasta), "--evidence", str(ev),
             "--sample", sample, "--out-fasta", str(out), "--out-report", str(rep)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        out_seq = "".join(l.strip() for l in out.read_text().splitlines() if not l.startswith(">"))
        with rep.open(newline="") as fh:
            report = list(csv.DictReader(fh, delimiter="\t"))[0]
        return out_seq, report

    def test_genuine_dimer_collapses(self):
        mono = monomer()
        seq, report = self.run_collapse(
            mono + mono,
            f"OG750\tconcatemer\t{len(mono)+1}-{2*len(mono)}\t16703\tTrue",
            sample="OG750",
        )
        self.assertEqual(report["action"], "collapsed")
        self.assertEqual(len(seq), len(mono))
        self.assertEqual(seq, mono)

    def test_non_concatemer_passthrough(self):
        mono = monomer()
        seq, report = self.run_collapse(mono + mono, "OGY\tnone\tNA\t16703\tTrue", sample="OGY")
        self.assertEqual(report["action"], "passthrough")
        self.assertEqual(len(seq), 2 * len(mono))

    def test_chimera_mislabelled_concatemer_fails_closed(self):
        mono = monomer()
        junk = "".join(random.Random(2).choice("ACGT") for _ in range(16000))
        seq, report = self.run_collapse(
            mono + junk,
            f"OGX\tconcatemer\t{len(mono)+1}-{len(mono)+len(junk)}\t16703\tTrue",
            sample="OGX",
        )
        self.assertEqual(report["action"], "passthrough")   # not corrupted
        self.assertEqual(len(seq), len(mono) + len(junk))


if __name__ == "__main__":
    unittest.main()
