"""Unit tests for bin/annotation_qc_gate.py (the MITOS2 FIX/PASS gate).

Pure stdlib. Run:
    python3 -m pytest tests/unit/test_annotation_qc_gate.py -v
"""

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# annotation_qc_gate imports orf_utils as a sibling (Nextflow puts bin/ on PATH).
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "annotation_qc_gate", ROOT / "bin" / "annotation_qc_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

PREFIX = "OG1.ilmn.250101.getorg1770.mitos2110"

# A clean coral CDS per PCG: ATG ... TAA, no internal stop (code 4).
CLEAN_CDS = "ATG" + "AAA" * 20 + "TAA"
# ND1 mis-boundary flavours.
NO_START_CDS = "AAA" + "AAA" * 20 + "TAA"
NO_STOP_CDS = "ATG" + "AAA" * 20 + "AAA"
INTERNAL_STOP_CDS = "ATG" + "AAA" * 5 + "TAA" + "AAA" * 5 + "TAA"

PCGS = gate.PCGS


def write_annotation(tmp, genes=None, cds_overrides=None, nd5_aa=600,
                     drop_from_gff=()):
    """Lay out a MITOS2-style annotation dir. genes defaults to the full
    cnidarian core; cds_overrides maps gene -> CDS nt string."""
    tmp = Path(tmp)
    genes = list(gate.CNIDARIAN_CORE if genes is None else genes)
    cds_overrides = cds_overrides or {}

    gff = tmp / f"{PREFIX}.gff"
    lines = ["##gff-version 3", f"##sequence-region {PREFIX} 1 18000"]
    pos = 100
    for g in genes:
        if g in drop_from_gff:
            continue
        lines.append(
            f"{PREFIX}\tmitos\tgene\t{pos}\t{pos + 90}\t.\t+\t.\t"
            f"ID=gene-{g};Name=MT-{g};Product={g}")
        pos += 200
    gff.write_text("\n".join(lines) + "\n")

    cds_dir = tmp / "cds"
    prot_dir = tmp / "proteins"
    cds_dir.mkdir()
    prot_dir.mkdir()
    for g in PCGS:
        if g not in genes:
            continue
        nt = cds_overrides.get(g, CLEAN_CDS)
        (cds_dir / f"{g}.{PREFIX}.fa").write_text(f">{PREFIX}\n{nt}\n")
    # ND5 protein length drives the cnidarian truncation heuristic.
    if "ND5" in genes:
        (prot_dir / f"MT-ND5.{PREFIX}.fa").write_text(">x\nM" + "A" * (nd5_aa - 1) + "\n")
    return gff, prot_dir, cds_dir


class GateTests(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.tmp = self.td.name

    def tearDown(self):
        self.td.cleanup()

    def _eval(self, code=4, **kw):
        gff, prot, cds = write_annotation(self.tmp, **kw)
        return gate.evaluate(str(gff), str(prot), str(cds), code,
                             540, is_cnidarian=(code == 4))

    def test_clean_annotation_passes(self):
        self.assertEqual(self._eval()[0], "PASS")

    def test_missing_rnr2_is_fix(self):
        decision, reason = self._eval(
            genes=[g for g in gate.CNIDARIAN_CORE if g != "RNR2"])
        self.assertEqual(decision, "FIX")
        self.assertIn("RNR2", reason)

    def test_short_nd5_is_fix(self):
        decision, reason = self._eval(nd5_aa=400)
        self.assertEqual(decision, "FIX")
        self.assertIn("ND5_trunc", reason)

    def test_nd1_non_m_start_is_fix(self):
        decision, reason = self._eval(cds_overrides={"ND1": NO_START_CDS})
        self.assertEqual(decision, "FIX")
        self.assertIn("ND1_no_start", reason)

    def test_nd1_no_stop_is_fix(self):
        decision, reason = self._eval(cds_overrides={"ND1": NO_STOP_CDS})
        self.assertEqual(decision, "FIX")
        self.assertIn("ND1_no_stop", reason)

    def test_internal_stop_is_fix(self):
        decision, reason = self._eval(cds_overrides={"CO1": INTERNAL_STOP_CDS})
        self.assertEqual(decision, "FIX")
        self.assertIn("CO1_internal_stop", reason)

    def test_pcg_check_is_code_generic_non_cnidarian(self):
        # code 9 (echinoderm): the cnidarian core/ND5 heuristics are skipped,
        # but a broken PCG ORF still routes to FIX.
        decision, reason = self._eval(
            code=9, cds_overrides={"ND1": NO_START_CDS})
        self.assertEqual(decision, "FIX")
        self.assertIn("ND1_no_start", reason)

    def test_non_cnidarian_missing_core_gene_still_passes(self):
        # A non-cnidarian invert missing RNR2 must NOT be forced to FIX by the
        # Anthozoa-specific core list.
        decision, _ = self._eval(
            code=9, genes=[g for g in gate.CNIDARIAN_CORE if g != "RNR2"])
        self.assertEqual(decision, "PASS")

    def test_missing_gff_degrades_to_pass_and_exits_zero(self):
        # A parse error (here: no such GFF) must never drop a sample: main()
        # catches, writes PASS, and exits 0.
        out = Path(self.tmp) / "qc.txt"
        proc = subprocess.run(
            [sys.executable, str(ROOT / "bin" / "annotation_qc_gate.py"),
             "--gff", str(Path(self.tmp) / "nope.gff"),
             "--proteins", self.tmp, "--cds", self.tmp,
             "--genetic-code", "4", "--out", str(out)],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0)
        self.assertTrue(out.read_text().startswith("PASS"))


if __name__ == "__main__":
    unittest.main()
