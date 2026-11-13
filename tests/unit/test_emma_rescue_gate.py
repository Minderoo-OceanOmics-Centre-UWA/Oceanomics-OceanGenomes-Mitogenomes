"""Unit tests for bin/emma_rescue_gate.py (the ND4L/ATP8 rescue FIX/PASS gate).

Pure stdlib -- no Biopython, no BLAST -- so these always run.
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
    "emma_rescue_gate", ROOT / "bin" / "emma_rescue_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

REF = gate.REF_GENES


def gff_for(genes, start=100, step=100):
    """Minimal EMMA-style GFF with a `gene` line per name, in the given order,
    at ascending coordinates."""
    lines = ["##gff-version 3", "##sequence-region chr 1 20000"]
    pos = start
    for g in genes:
        lines.append(
            f"chr\tEmma\tgene\t{pos}\t{pos + 50}\t.\t+\t.\tID=x{pos};Name=MT-{g}")
        pos += step
    return "\n".join(lines) + "\n"


class GateDecisionTests(unittest.TestCase):
    def _decide(self, genes):
        p = Path(self.tmp) / "a.gff"
        p.write_text(gff_for(genes))
        return gate.decide(p)

    def setUp(self):
        import tempfile
        self.tmp = tempfile.mkdtemp()

    def test_complete_annotation_passes(self):
        self.assertEqual(self._decide(REF), ("PASS", "-"))

    def test_only_nd4l_missing_is_fix(self):
        self.assertEqual(self._decide([g for g in REF if g != "ND4L"]),
                         ("FIX", "ND4L"))

    def test_only_atp8_missing_is_fix(self):
        self.assertEqual(self._decide([g for g in REF if g != "ATP8"]),
                         ("FIX", "ATP8"))

    def test_both_missing_is_fix_in_ref_order(self):
        self.assertEqual(
            self._decide([g for g in REF if g not in ("ND4L", "ATP8")]),
            ("FIX", "ATP8,ND4L"))

    def test_other_gene_missing_is_pass(self):
        # ND1 gone -> not a rescuable case, even though ND4L is also dropped.
        self.assertEqual(
            self._decide([g for g in REF if g not in ("ND1", "ND4L")]),
            ("PASS", "-"))

    def test_missing_flank_is_pass(self):
        # ND4L and its downstream flank ND4 both gone -> window undefined.
        self.assertEqual(
            self._decide([g for g in REF if g not in ("ND4L", "ND4")]),
            ("PASS", "-"))

    def test_scrambled_order_is_pass(self):
        genes = [g for g in REF if g != "ND4L"]
        genes[5], genes[20] = genes[20], genes[5]  # disorder the present set
        self.assertEqual(self._decide(genes), ("PASS", "-"))

    def test_malformed_gff_degrades_to_pass(self):
        p = Path(self.tmp) / "bad.gff"
        p.write_text("not\ta\tgff\n")
        self.assertEqual(gate.decide(p), ("PASS", "-"))

    def test_tRNA_only_drop_is_pass(self):
        self.assertEqual(self._decide([g for g in REF if g != "TS1"]),
                         ("PASS", "-"))


if __name__ == "__main__":
    unittest.main()


class VariantOrderTests(unittest.TestCase):
    """The gate must judge order against the ACCEPTED order for the taxon.

    Left keyed on canonical only, a clade whose real gene order is non-canonical
    is judged out of order here and silently declined for rescue -- so an assembly
    missing ND4L would be held for the very defect this gate could have repaired.
    Preventing that drift between the QC step and the gates is why the order table
    and the GFF reader both live in bin/mito_gene_order.py.
    """

    def _imq(self, drop=None):
        genes = list(REF)
        i = genes.index("TQ")
        genes[i:i + 2] = ["TM", "TQ"]
        if drop:
            genes = [g for g in genes if g != drop]
        return genes

    def _decide(self, genes, taxon, tmp_path=None):
        import tempfile
        p = Path(tempfile.mkdtemp()) / "a.gff"
        p.write_text(gff_for(genes))
        return gate.decide(p, taxon)

    def test_a_variant_taxon_missing_nd4l_is_offered_for_rescue(self):
        state, targets = self._decide(self._imq(drop="ND4L"), {"genus": "Chlorurus"})
        self.assertEqual(state, "FIX")
        self.assertEqual(targets, "ND4L")

    def test_the_same_assembly_without_the_taxonomy_is_declined(self):
        # The regression: no taxon, so only the canonical order is accepted, the
        # IMQ order reads as out-of-order, and the rescue is suppressed.
        state, _ = self._decide(self._imq(drop="ND4L"), None)
        self.assertEqual(state, "PASS")

    def test_a_variant_taxon_with_a_canonical_assembly_still_works(self):
        genes = [g for g in REF if g != "ND4L"]
        state, targets = self._decide(genes, {"genus": "Chlorurus"})
        self.assertEqual(state, "FIX")
        self.assertEqual(targets, "ND4L")

    def test_an_unkeyed_genus_with_the_variant_order_is_still_declined(self):
        state, _ = self._decide(self._imq(drop="ND4L"), {"genus": "Epibulus"})
        self.assertEqual(state, "PASS")
