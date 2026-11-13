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
                              read_conversion_audit, score_cds_dir)
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


class ScoreCdsDirTests(unittest.TestCase):
    """score_cds_dir is the shared scorer bin/mito_run_snapshot.py reuses.

    The snapshot used to re-implement this loop without the terminator-convention
    detection, which flagged no_stop on all 13 PCGs of every EMMA sample (136
    false findings in NOVA_260909_LA). These tests pin the contract it relies on.
    """

    def test_per_gene_results_are_returned_in_filename_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_cds_dir(tmp, {
                "ND1": "ATG" + "GGG" * 5 + "TAA",
                "ATP6": "ATG" + "GGG" * 5 + "TAA",
            })
            per_gene, stop_excluded, _note = score_cds_dir(d, 4)
        self.assertEqual([g for g, _i, _a in per_gene], ["ATP6", "ND1"])
        self.assertFalse(stop_excluded)
        self.assertEqual([i for _g, i, _a in per_gene], [[], []])

    def test_stop_excluded_is_detected_and_suspends_only_the_terminator(self):
        """EMMA publishes CDS without the stop codon; that is a convention, not a defect."""
        with tempfile.TemporaryDirectory() as tmp:
            # Six in-frame genes, none carrying a terminator: the EMMA convention.
            d = make_cds_dir(tmp, {g: "ATG" + "GGG" * 5
                                   for g in ("ND1", "ND2", "ND3", "ND4", "ND5", "ND6")})
            per_gene, stop_excluded, _note = score_cds_dir(d, 4)
        self.assertTrue(stop_excluded)
        flat = [i for _g, issues, _a in per_gene for i in issues]
        self.assertEqual(flat, [], f"stop-excluded CDS must not be flagged: {flat}")
        # The convention note must NOT be attached to a gene: a per-gene consumer would
        # report it as a CDS defect on whichever gene sorted first.
        per_gene_advisories = [a for _g, _i, adv in per_gene for a in adv]
        self.assertEqual(per_gene_advisories, [])

    def test_frame_check_survives_stop_exclusion(self):
        """Only the terminator check is suspended, never the frame check."""
        with tempfile.TemporaryDirectory() as tmp:
            genes = {g: "ATG" + "GGG" * 5 for g in ("ND1", "ND2", "ND3", "ND4", "ND5")}
            genes["ND6"] = "ATG" + "GGG" * 5 + "G"   # remainder 1
            d = make_cds_dir(tmp, genes)
            per_gene, stop_excluded, _note = score_cds_dir(d, 4)
        self.assertTrue(stop_excluded)
        flat = [i for _g, issues, _a in per_gene for i in issues]
        self.assertEqual(flat, ["ND6:frame_remainder(1)"])

    def test_missing_dir_reports_not_evaluated(self):
        per_gene, stop_excluded, note = score_cds_dir(None, 4)
        self.assertIsNone(per_gene)
        self.assertFalse(stop_excluded)
        self.assertIsNone(note)
        self.assertEqual(check_cds_integrity(None, "sample", 4), ("not_evaluated", []))

def make_emma_pair(tmp, genes):
    """Build an EMMA-shaped annotation: a reoriented genome + GFF, plus stop-EXCLUDED CDS files.

    `genes` is {gene: full_nt_including_terminator}. The per-gene CDS files get that sequence
    with its last codon removed, which is what EMMA publishes and what makes stop_excluded
    fire; the GFF span covers the full sequence, terminator included.
    """
    d = Path(tmp)
    cds = d / "cds"
    cds.mkdir()
    genome, gff_rows, offset = [], [], 0
    for gene, full in genes.items():
        start = offset + 1
        end = offset + len(full)
        genome.append(full)
        offset = end
        gff_rows.append(
            f"asm\tEmma\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene};Name=MT-{gene}\n")
        # EMMA drops the terminator from the per-gene file: a whole codon for a real stop,
        # or the partial T/TA for a poly-A-completed one.
        trimmed = full[:-3] if len(full) % 3 == 0 else full[:len(full) - (len(full) % 3)]
        (cds / f"MT-{gene}.asm.emma102.fa").write_text(f">asm\n{trimmed}\n")
    (d / "asm.emma102.fa").write_text(">asm\n" + "".join(genome) + "\n")
    gff = d / "asm.emma102.gff"
    gff.write_text("##gff-version 3\n" + "".join(gff_rows))
    return cds, gff


class StopExcludedIsCheckedFromTheGffTests(unittest.TestCase):
    """The terminator check used to be SUSPENDED for a stop-excluded annotator.

    EMMA publishes annotation/cds/*.fa without the stop codon, so scoring those files reports
    no_stop on every CDS. The old code detected that and switched the terminator check off
    entirely -- which made a genuinely broken EMMA CDS invisible, a false negative on all 9
    fish in NOVA_260909_LA. The GFF span carries the terminator, so it can be checked instead.
    """

    # Six clean genes are enough for the >=3 / <=40% convention detection to fire.
    CLEAN = {f"ND{i}": "ATG" + "GGG" * 5 + "TAA" for i in range(1, 7)}

    def test_a_clean_emma_annotation_is_fully_checked_and_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, dict(self.CLEAN))
            status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        self.assertEqual(status, "advisory")
        # The suffix is the claim that matters: the terminator WAS checked.
        self.assertEqual(issues, ["cds_convention:stop_codon_excluded(checked_from_gff)"])

    def test_a_poly_a_completed_terminator_is_not_a_defect(self):
        """A span whose length is not a multiple of 3 is EMMA's other convention."""
        genes = dict(self.CLEAN)
        genes["ATP6"] = "ATG" + "GGG" * 5 + "T"      # 16 nt, truncated stop
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, genes)
            status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        self.assertEqual(status, "advisory")
        self.assertFalse([i for i in issues if "ATP6" in i], issues)

    def test_a_missing_terminator_is_now_caught(self):
        """This is the whole point: before, no_stop here was unreportable."""
        genes = dict(self.CLEAN)
        genes["CYTB"] = "ATG" + "GGG" * 5 + "CCC"    # in frame, ends on no stop
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, genes)
            status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        self.assertEqual(status, "broken")
        self.assertTrue([i for i in issues if i.startswith("CYTB:no_stop")], issues)

    def test_an_internal_stop_is_now_caught(self):
        genes = dict(self.CLEAN)
        genes["CO1"] = "ATG" + "GGG" + "TAA" + "GGG" + "TAA"
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, genes)
            status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        self.assertEqual(status, "broken")
        self.assertTrue([i for i in issues if i.startswith("CO1:internal_stops")], issues)

    def test_without_a_gff_the_old_suspension_still_applies(self):
        """The fallback has to stay: a gene the GFF cannot resolve must not go unchecked."""
        genes = dict(self.CLEAN)
        genes["CYTB"] = "ATG" + "GGG" * 5 + "CCC"
        with tempfile.TemporaryDirectory() as tmp:
            cds, _gff = make_emma_pair(tmp, genes)
            status, issues = check_cds_integrity(cds, "asm", 2)
        self.assertEqual(status, "advisory")
        self.assertEqual(issues, ["cds_convention:stop_codon_excluded"])

    def test_the_gff_is_read_against_the_reoriented_fasta_not_the_assembly(self):
        """A missing sibling .fa must fall back, never score against the wrong sequence."""
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, dict(self.CLEAN))
            Path(tmp, "asm.emma102.fa").unlink()
            status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        self.assertEqual(status, "advisory")
        self.assertEqual(issues, ["cds_convention:stop_codon_excluded"])

    def test_a_multi_exon_gene_falls_back_rather_than_guessing_exon_order(self):
        """An origin-spanning gene's pieces are high-coordinate first, not sorted by start."""
        with tempfile.TemporaryDirectory() as tmp:
            cds, gff = make_emma_pair(tmp, dict(self.CLEAN))
            gff.write_text(gff.read_text() +
                           "asm\tEmma\tCDS\t1\t9\t.\t+\t0\tID=ND1;Name=MT-ND1\n")
            _status, issues = check_cds_integrity(cds, "asm", 2, gff_path=gff)
        # ND1 now has two CDS rows, so it is not GFF-graded and the run-wide note drops the
        # (checked_from_gff) claim rather than overstating what was verified.
        self.assertEqual(issues, ["cds_convention:stop_codon_excluded"])

