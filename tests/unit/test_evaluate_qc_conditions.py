"""Unit tests for bin/evaluate_qc_conditions.py -- specifically the held-reason
output that feeds the run-level held_samples.tsv (Part E).

Pure stdlib; invokes the script as a subprocess.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "evaluate_qc_conditions.py"

BLAST_HEADER = "og_id\tx\tnom_species\ty\tfound_in_blast\n"


def run(tmp, blast_rows, passed, circular="null", circ_check=None, trna_advisory=None):
    tmp = Path(tmp)
    blast = tmp / "blast.tsv"
    blast.write_text(BLAST_HEADER + "".join(blast_rows))
    ann = tmp / "ann.csv"
    if trna_advisory is None:
        ann.write_text("passed\n" + passed + "\n")
    else:
        ann.write_text("passed,trna_advisory\n" + passed + "," + trna_advisory + "\n")
    args = [sys.executable, str(SCRIPT),
            "--blast-table", str(blast), "--annotation-csv", str(ann),
            "--circular", circular,
            "--output-species", str(tmp / "s.txt"),
            "--output-proceed", str(tmp / "p.txt"),
            "--output-circular", str(tmp / "c.txt"),
            "--output-reason", str(tmp / "r.txt"),
            "--output-versions", str(tmp / "v.yml")]
    if circ_check is not None:
        cc = tmp / "cc.tsv"
        cc.write_text(circ_check)
        args += ["--circularity-check", str(cc)]
    proc = subprocess.run(args, check=True, capture_output=True, text=True)
    run.last_stdout = proc.stdout
    return (tmp / "p.txt").read_text(), (tmp / "r.txt").read_text()


class HeldReasonTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.tmp = self.td.name

    def tearDown(self):
        self.td.cleanup()

    def test_clean_sample_proceeds_with_empty_reason(self):
        proceed, reason = run(self.tmp,
                              ["OG1\t.\tGadus morhua\t.\tyes\n"], "yes",
                              circular="true")
        self.assertEqual(proceed, "true")
        self.assertEqual(reason, "")

    def test_not_in_blast_reason(self):
        proceed, reason = run(self.tmp,
                              ["OG1\t.\tGadus morhua\t.\tno\n"], "yes",
                              circular="true")
        self.assertEqual(proceed, "false")
        self.assertEqual(reason, "species_not_in_blast")

    def test_trna_tolerated_pass_still_proceeds_and_is_noted(self):
        proceed, reason = run(self.tmp,
                              ["OG1\t.\tGadus morhua\t.\tyes\n"], "yes",
                              circular="true", trna_advisory="TP")
        self.assertEqual(proceed, "true")
        self.assertEqual(reason, "")
        self.assertIn("tolerated missing tRNAs: TP", run.last_stdout)

    def test_trna_advisory_no_is_not_noted(self):
        proceed, _ = run(self.tmp,
                         ["OG1\t.\tGadus morhua\t.\tyes\n"], "yes",
                         circular="true", trna_advisory="no")
        self.assertEqual(proceed, "true")
        self.assertNotIn("tolerated missing tRNAs", run.last_stdout)

    def test_multiple_causes_are_joined(self):
        proceed, reason = run(self.tmp,
                              ["OG1\t.\tGadus morhua\t.\tno\n"], "no",
                              circular="false")
        self.assertEqual(proceed, "false")
        self.assertEqual(
            reason.split(";"),
            ["species_not_in_blast", "annotation_failed", "not_circular"])

    def test_assembly_anomaly_reason_carries_type(self):
        cc = ("anomaly_type\tlength_anomaly\tfinal_verdict_circular\n"
              "concatemer\tno\tna\n")
        proceed, reason = run(self.tmp,
                              ["OG1\t.\tGadus morhua\t.\tyes\n"], "yes",
                              circular="true", circ_check=cc)
        self.assertEqual(proceed, "false")
        self.assertEqual(reason, "assembly_anomaly:concatemer")


if __name__ == "__main__":
    unittest.main()
