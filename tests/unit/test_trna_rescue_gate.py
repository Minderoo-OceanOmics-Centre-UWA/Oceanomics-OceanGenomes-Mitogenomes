"""Unit tests for bin/trna_rescue_gate.py (the tRNA rescue FIX/PASS gate).

Pure stdlib -- no tRNAscan-SE -- so these always run.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "trna_rescue_gate", ROOT / "bin" / "trna_rescue_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

REF = gate.REF_GENES


def gff_for(genes, start=100, step=100):
    lines = ["##gff-version 3", "##sequence-region chr 1 16500"]
    pos = start
    for g in genes:
        lines.append(
            f"chr\tEmma\tgene\t{pos}\t{pos + 50}\t.\t+\t.\tID=x{pos};Name=MT-{g}")
        pos += step
    return "\n".join(lines) + "\n"


class GateDecisionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _decide(self, genes):
        p = Path(self.tmp) / "a.gff"
        p.write_text(gff_for(genes))
        return gate.decide(p)

    def test_complete_passes(self):
        self.assertEqual(self._decide(REF), ("PASS", "-"))

    def test_single_trna_missing_is_fix(self):
        self.assertEqual(self._decide([g for g in REF if g != "TP"]), ("FIX", "TP"))

    def test_multiple_trna_missing_fix_in_ref_order(self):
        self.assertEqual(
            self._decide([g for g in REF if g not in ("TN", "TA", "TW")]),
            ("FIX", "TW,TA,TN"))

    def test_missing_pcg_is_pass(self):
        self.assertEqual(self._decide([g for g in REF if g != "ND4L"]), ("PASS", "-"))

    def test_missing_rrna_is_pass(self):
        self.assertEqual(self._decide([g for g in REF if g != "RNR1"]), ("PASS", "-"))

    def test_trna_missing_but_out_of_order_is_pass(self):
        genes = [g for g in REF if g != "TP"]
        genes[4], genes[7] = genes[7], genes[4]
        self.assertEqual(self._decide(genes), ("PASS", "-"))

    def test_malformed_gff_degrades_to_pass(self):
        p = Path(self.tmp) / "bad.gff"
        p.write_text("not\ta\tgff\n")
        self.assertEqual(gate.decide(p), ("PASS", "-"))


if __name__ == "__main__":
    unittest.main()
