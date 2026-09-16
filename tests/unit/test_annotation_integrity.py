"""The final annotation verdict, computed on what was actually published.

annotation_qc_gate.py cannot serve as this verdict: it runs on raw MITOS2 output
BEFORE any repair branch, covers only the fixed 13-PCG list, and always exits 0
with exceptions degrading to PASS. It is a router. This runs post-merge, over
every exported CDS, and is allowed to fail.
"""

import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

from annotation_stats import (check_cds_integrity, grade_gene_order,
                              read_conversion_audit)
from mito_gene_order import REF_GENES


def make_cds_dir(tmp, genes):
    d = Path(tmp) / "cds"
    d.mkdir()
    for gene, seq in genes.items():
        (d / f"{gene}.sample.fa").write_text(f">sample\n{seq}\n")
    return d


class IntegrityTests(unittest.TestCase):
    def test_clean_cds_is_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {"ND1": "ATG" + "GGG" * 5 + "TAA"})
            status, issues = check_cds_integrity(d, "sample", 4)
        self.assertEqual(status, "ok")
        self.assertEqual(issues, [])

    def test_implausible_start_is_broken(self):
        """The five identical Demospongiae ND4L defects look exactly like this."""
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {"ND4L": "TTT" + "GGG" * 5 + "TAA"})
            status, issues = check_cds_integrity(d, "sample", 4)
        self.assertEqual(status, "broken")
        self.assertTrue(any("implausible_start(TTT)" in i for i in issues))

    def test_in_frame_sense_codon_end_is_broken(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {"ATP8": "ATG" + "GGG" * 5 + "AGA"})
            status, issues = check_cds_integrity(d, "sample", 4)
        self.assertEqual(status, "broken")
        self.assertTrue(any("no_stop" in i for i in issues))

    def test_documented_alternative_start_is_advisory_not_broken(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {"CO1": "CGA" + "GGG" * 5 + "TAA"})
            status, issues = check_cds_integrity(d, "sample", 5,
                                                 lineage="Malacostraca")
        self.assertEqual(status, "advisory")
        self.assertTrue(any("alternative_start(CGA)" in i for i in issues))

    def test_every_exported_cds_is_inspected_not_a_fixed_list(self):
        """ATP9 and mtMutS are outside the vertebrate 13, so a fixed-list check
        would never look at them."""
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {
                "ND1": "ATG" + "GGG" * 5 + "TAA",
                "ATP9": "TTT" + "GGG" * 5 + "AGA",
            })
            status, issues = check_cds_integrity(d, "sample", 4)
        self.assertEqual(status, "broken")
        self.assertTrue(any(i.startswith("ATP9") for i in issues))

    def test_stop_excluded_convention_is_detected_not_assumed(self):
        """EMMA exports CDS without the terminator (an EMMA ATP6 of 681 nt is
        exactly 227 aa). Applying the MITOS convention reports a missing stop on
        every single vertebrate CDS: 13 false findings a sample."""
        genes = {g: "ATG" + "GGG" * 20 for g in
                 ("ND1", "ND2", "ND3", "ND4", "CO1", "CO2", "CO3", "CYTB")}
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, genes)
            status, issues = check_cds_integrity(d, "sample", 2)
        self.assertEqual(status, "advisory")
        self.assertIn("cds_convention:stop_codon_excluded", issues)
        self.assertFalse([i for i in issues if "no_stop" in i])

    def test_frame_is_still_checked_under_the_stop_excluded_convention(self):
        """Only the terminator check is suspended, not the reading frame."""
        genes = {g: "ATG" + "GGG" * 20 for g in
                 ("ND1", "ND2", "ND3", "ND4", "CO1", "CO2", "CO3")}
        genes["CYTB"] = "ATG" + "GGG" * 20 + "GG"
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, genes)
            status, issues = check_cds_integrity(d, "sample", 2)
        self.assertEqual(status, "broken")
        self.assertTrue(any("frame_remainder" in i for i in issues))

    def test_missing_cds_dir_is_not_evaluated_not_a_pass(self):
        status, issues = check_cds_integrity(None, "sample", 4)
        self.assertEqual(status, "not_evaluated")


class ConversionAuditTests(unittest.TestCase):
    def test_unmapped_and_auxiliary_are_read_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "audit.tsv"
            p.write_text(
                "raw_name\tinternal_name\tlocus\tfeature_type\tdisposition\treason\n"
                "cox1\tCO1\t1\tCDS\tretained\t\n"
                "lagli\tLAGLI\t1\tmisc_feature\tknown_auxiliary\t\n"
                "wibble\t\t\t\tunmapped\tno mapping\n")
            unmapped, auxiliary = read_conversion_audit(p)
        self.assertEqual(unmapped, ["wibble"])
        self.assertEqual(auxiliary, ["lagli"])

    def test_absent_audit_is_empty_not_an_error(self):
        self.assertEqual(read_conversion_audit(None), ([], []))


class OrderGradingTests(unittest.TestCase):
    def test_reference_order(self):
        ref = [g for g in REF_GENES]
        self.assertEqual(grade_gene_order(ref, ref), ("reference", []))

    def test_one_displaced_trna(self):
        found = list(REF_GENES)
        i = found.index("TV")
        found.insert(i + 3, found.pop(i))
        status, deviations = grade_gene_order(found, list(REF_GENES))
        self.assertEqual(status, "trna_displacement")
        self.assertIn("TV", deviations)

    def test_displacement_is_minimal_not_positional(self):
        """Moving one gene shifts every gene after it, so a position-by-position
        comparison would report a dozen deviations for a single translocation."""
        found = list(REF_GENES)
        i = found.index("TV")
        found.insert(i + 3, found.pop(i))
        _status, deviations = grade_gene_order(found, list(REF_GENES))
        self.assertLessEqual(len(deviations), 2, deviations)

    def test_pcg_block_move_is_a_rearrangement(self):
        found = list(REF_GENES)
        i, j = found.index("ND6"), found.index("CYTB")
        found[i], found[j] = found[j], found[i]
        status, _deviations = grade_gene_order(found, list(REF_GENES))
        self.assertIn(status, ("rearranged_block", "rearranged_major"))


if __name__ == "__main__":
    unittest.main()
