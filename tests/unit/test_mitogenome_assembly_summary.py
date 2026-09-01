"""Unit tests for the mitogenome assembly-summary QC classification.

Covers the Phase 1 robustness changes: the blocking-vs-advisory split
(is_complete_core / blocking_reasons) and the rewritten ambiguous-graph
detector (getorganelle_graph_ambiguous). Rows are modelled on the real
manual-review samples from the mitogenomes-missing-audit-3 run.
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
    "mitogenome_assembly_summary", ROOT / "bin" / "mitogenome_assembly_summary.py"
)
mas = importlib.util.module_from_spec(SPEC)
# Register before exec so the module's @dataclass forward-ref resolution works.
sys.modules[SPEC.name] = mas
SPEC.loader.exec_module(mas)

THRESHOLDS = mas.Thresholds(
    min_mean_coverage=20,
    max_coverage_cv=1.0,
    min_length=10000,
    max_length=25000,
    expected_gene_count=37,
    expected_pcg_count=13,
)


def complete_row(**overrides):
    """A clean circular vertebrate mitogenome row (all 37 genes / 13 PCGs)."""
    row = {
        "final_length_bp": "16500",
        "circularised": "true",
        "num_candidate_contigs": "1",
        "num_final_contigs": "1",
        "num_genes": "37",
        "num_cds": "13",
        "missing_genes": "no",
        "frameshift_flag": "false",
        "mean_coverage": "300",
        "coverage_cv": "0.2",
        "numt_flag": "false",
        "reference_relevance": "",
        "reference_divergence": "",
        "anomaly_type": "none",
        "length_anomaly": "no",
        "manual_review_reason": "",
    }
    row.update(overrides)
    return row


def paths(*names):
    return [Path(n) for n in names]


class AmbiguousGraphTests(unittest.TestCase):
    def test_single_path_multi_segment_is_not_ambiguous(self):
        # The false-positive that dominated the audit run: one resolved circular
        # path whose GFA carries several S-segments (normal repeat structure).
        files = paths(
            "OG61.animal_mt.K115.complete.graph1.1.path_sequence.fasta",
            "OG61.animal_mt.K115.complete.graph1.selected_graph.gfa",
            "OG61.get_org.log.txt",
        )
        self.assertFalse(mas.getorganelle_graph_ambiguous(files))

    def test_multiple_path_sequences_is_ambiguous(self):
        files = paths(
            "OG100.complete.graph1.1.path_sequence.fasta",
            "OG100.complete.graph1.2.path_sequence.fasta",
            "OG100.selected_graph.gfa",
        )
        self.assertTrue(mas.getorganelle_graph_ambiguous(files))

    def test_multiple_selected_graphs_is_ambiguous(self):
        files = paths(
            "OG9.complete.graph1.1.path_sequence.fasta",
            "OG9.animal_mt.graph1.selected_graph.gfa",
            "OG9.animal_mt.graph2.selected_graph.gfa",
        )
        self.assertTrue(mas.getorganelle_graph_ambiguous(files))

    def test_duplicate_staged_names_not_double_counted(self):
        # Same basename staged from two publish locations must not read as two paths.
        files = paths(
            "a/OG61.complete.graph1.1.path_sequence.fasta",
            "b/OG61.complete.graph1.1.path_sequence.fasta",
        )
        self.assertFalse(mas.getorganelle_graph_ambiguous(files))


class BlockingAdvisoryTests(unittest.TestCase):
    def apply(self, row):
        mas.apply_qc(row, THRESHOLDS)
        return row.get("manual_review_reason", ""), row.get(mas.BLOCKING_KEY, "")

    def test_no_congeneric_on_complete_is_advisory(self):
        reason, blocking = self.apply(complete_row(reference_divergence="NON_CONGENERIC"))
        self.assertIn("no_congeneric_reference", reason)   # kept for transparency
        self.assertEqual(blocking, "")                     # but not blocking

    def test_reference_mismatch_on_complete_is_advisory(self):
        # A finished mitogenome is finished whatever reference built it. This is the
        # false positive that held 18 complete assemblies at manual_review in the
        # mitogenomes-missing-audit-5 run.
        reason, blocking = self.apply(complete_row(reference_relevance="MISMATCH"))
        self.assertIn("reference_mismatch", reason)        # still recorded
        self.assertEqual(blocking, "")                     # but not blocking

    def test_reference_divergent_on_complete_is_advisory(self):
        reason, blocking = self.apply(complete_row(reference_relevance="DIVERGENT"))
        self.assertIn("reference_divergent", reason)
        self.assertEqual(blocking, "")

    def test_reference_relevance_states_are_mutually_exclusive(self):
        # DIVERGENT must not also raise the mismatch reason, and vice versa.
        divergent, _ = self.apply(complete_row(reference_relevance="DIVERGENT"))
        self.assertNotIn("reference_mismatch", divergent)
        mismatch, _ = self.apply(complete_row(reference_relevance="MISMATCH"))
        self.assertNotIn("reference_divergent", mismatch)

    def test_reference_pass_raises_no_reason(self):
        reason, blocking = self.apply(complete_row(reference_relevance="PASS"))
        self.assertNotIn("reference", reason)
        self.assertEqual(blocking, "")

    def test_reference_mismatch_still_blocks_a_damaged_assembly(self):
        # OG2102-like: the reference really was too divergent and the assembly
        # collapsed. The damage blocks, even though the reference reason no longer
        # does on its own.
        _, blocking = self.apply(
            complete_row(reference_relevance="MISMATCH", num_genes="25", num_cds="7")
        )
        self.assertIn("missing_protein_coding_genes", blocking)

    def test_reference_mismatch_still_blocks_a_non_circular_assembly(self):
        # OG810-like: structural defect present, so the row stays in review.
        _, blocking = self.apply(
            complete_row(reference_relevance="MISMATCH", circularised="false")
        )
        self.assertIn("not_circularised", blocking)

    def test_low_coverage_on_complete_is_advisory(self):
        reason, blocking = self.apply(complete_row(mean_coverage="12"))
        self.assertIn("low_mean_coverage", reason)
        self.assertEqual(blocking, "")

    def test_single_trna_shortfall_is_advisory(self):
        reason, blocking = self.apply(complete_row(num_genes="36", missing_genes="TS1"))
        self.assertIn("missing_genes", reason)
        self.assertEqual(blocking, "")

    def test_missing_pcg_always_blocks(self):
        # Genuine collapse: only 7 CDS. Must stay blocking despite no-congeneric.
        _, blocking = self.apply(
            complete_row(num_genes="25", num_cds="7", missing_genes="CO2;CO3;ND3",
                         reference_divergence="NON_CONGENERIC")
        )
        self.assertIn("missing_protein_coding_genes", blocking)

    def test_low_coverage_when_not_circular_blocks(self):
        _, blocking = self.apply(complete_row(circularised="false", mean_coverage="12"))
        self.assertIn("low_mean_coverage", blocking)

    def test_length_anomaly_still_blocks_on_otherwise_complete(self):
        _, blocking = self.apply(
            complete_row(final_length_bp="32000", anomaly_type="concatemer")
        )
        self.assertIn("concatemer", blocking)

    def test_multiple_final_contigs_blocks(self):
        _, blocking = self.apply(complete_row(num_final_contigs="3", circularised="true"))
        self.assertIn("multiple_final_contigs", blocking)

    def test_low_coverage_fragmented_tagged_data_limited(self):
        # OG765-like: shallow HiC, non-circular, out-of-range length.
        reason, _ = self.apply(complete_row(
            circularised="false", final_length_bp="3407", num_final_contigs="2",
            mean_coverage="9", num_genes="", num_cds="", missing_genes=""))
        self.assertIn("data_limited", reason)

    def test_adequate_coverage_not_tagged_data_limited(self):
        # OG810-like: fragmented but coverage 26 (above 0.75*20=15) -> not data_limited.
        reason, _ = self.apply(complete_row(
            circularised="false", final_length_bp="9856", num_final_contigs="6",
            mean_coverage="26", num_genes="37", num_cds="13"))
        self.assertNotIn("data_limited", reason)

    def test_complete_low_coverage_not_data_limited(self):
        # A complete circular assembly at low coverage is advisory, never data_limited
        # (no fragmentation reason present).
        reason, _ = self.apply(complete_row(mean_coverage="12"))
        self.assertNotIn("data_limited", reason)


class CollapseProvenanceTests(unittest.TestCase):
    """A genuine collapse forks the molecule to <prefix>_collapsed, which carries its
    own identity end to end and therefore its own summary row. This run is then the
    superseded pre-collapse original, and must keep its OWN length rather than
    borrowing the monomer's."""

    def make_run(self, report_text):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        rep = root / "OG750.hifi.v323mitohifi.concatemer_collapse.tsv"
        rep.write_text(report_text)
        return mas.RunFiles(sample_id="OG750", prefix="OG750.hifi.v323mitohifi",
                            assembler="MitoHiFi", files=[rep])

    def test_collapsed_report_marks_superseded_and_keeps_own_length(self):
        run = self.make_run(
            "sample\taction\toriginal_length\tcollapsed_length\treference_length\ttail_identity\treason\n"
            "OG750\tcollapsed\t32672\t16466\t16703\t1.0\tcollapsed_2.0x_concatemer\n"
        )
        row = complete_row(final_length_bp="32672", anomaly_type="concatemer", length_anomaly="yes")
        mas.apply_collapse_provenance(row, run)
        # The monomer's 16466 belongs to the <prefix>_collapsed row, not this one.
        self.assertEqual(row["final_length_bp"], "32672")
        self.assertTrue(row[mas.SUPERSEDED_KEY])

    def test_superseded_original_leaves_the_review_queue(self):
        run = self.make_run(
            "sample\taction\toriginal_length\tcollapsed_length\treference_length\ttail_identity\treason\n"
            "OG750\tcollapsed\t32672\t16466\t16703\t1.0\tcollapsed_2.0x_concatemer\n"
        )
        row = complete_row(final_length_bp="32672", anomaly_type="concatemer", length_anomaly="yes")
        mas.apply_collapse_provenance(row, run)
        mas.apply_qc(row, THRESHOLDS)
        mas.finalise_status(row)
        # It is over-length and flagged, but triaging it would mean reviewing the
        # same molecule twice -- the monomer's row is the one to review.
        self.assertTrue(row.get(mas.BLOCKING_KEY))
        self.assertEqual(row["status"], "superseded")
        # apply_qc rebuilds manual_review_reason, so the marker has to survive it.
        self.assertIn("superseded_by_collapse", row["manual_review_reason"])

    def test_passthrough_report_leaves_row_untouched(self):
        run = self.make_run(
            "sample\taction\toriginal_length\tcollapsed_length\treference_length\ttail_identity\treason\n"
            "OGX\tpassthrough\t32000\t32000\t16703\t0.0\ttail_identity_below_threshold\n"
        )
        row = complete_row(final_length_bp="32000", anomaly_type="unresolved", length_anomaly="yes")
        mas.apply_collapse_provenance(row, run)
        self.assertEqual(row["final_length_bp"], "32000")
        self.assertEqual(row["anomaly_type"], "unresolved")
        self.assertFalse(row.get(mas.SUPERSEDED_KEY))


class CollapseChildEvidenceTests(unittest.TestCase):
    """The post-curation check is measured ON the monomer but written under the
    PRE-collapse prefix, so the monomer's own row has to reach across for it."""

    PARENT = "OG750.hifi.241004.v323mitohifi"
    HEADER = "sample\tfinal_verdict_circular\tlength_anomaly\tanomaly_type\n"

    def make(self, prefix, evidence_name):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        ev = root / evidence_name
        ev.write_text(self.HEADER + f"{self.PARENT}\tTrue\tno\tnone\n")
        run = mas.RunFiles(sample_id="OG750", prefix=prefix, assembler="MitoHiFi", files=[])
        return run, [ev]

    def test_collapsed_child_takes_the_post_curation_verdict(self):
        run, files = self.make(f"{self.PARENT}_collapsed", f"{self.PARENT}.post_curation_check.tsv")
        # anomaly_type "none" is dropped by first_value (it is a PLACEHOLDER_VALUE),
        # exactly as anomaly_for_run drops it; apply_qc treats "" and "none" alike.
        self.assertEqual(
            mas.collapse_child_evidence(run, files),
            {"circularised": "true", "length_anomaly": "no"},
        )

    def test_a_real_anomaly_on_the_monomer_is_carried_over(self):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        ev = root / f"{self.PARENT}.post_curation_check.tsv"
        ev.write_text(self.HEADER + f"{self.PARENT}\tFalse\tyes\tunresolved\n")
        run = mas.RunFiles(sample_id="OG750", prefix=f"{self.PARENT}_collapsed",
                           assembler="MitoHiFi", files=[])
        self.assertEqual(
            mas.collapse_child_evidence(run, [ev]),
            {"circularised": "false", "length_anomaly": "yes", "anomaly_type": "unresolved"},
        )

    def test_non_collapsed_run_reaches_for_nothing(self):
        run, files = self.make(self.PARENT, f"{self.PARENT}.post_curation_check.tsv")
        self.assertEqual(mas.collapse_child_evidence(run, files), {})

    def test_another_samples_evidence_is_not_borrowed(self):
        run, files = self.make(f"{self.PARENT}_collapsed", "OG751.hifi.241004.v323mitohifi.post_curation_check.tsv")
        self.assertEqual(mas.collapse_child_evidence(run, files), {})


class GetOrganelleEvidenceTests(unittest.TestCase):
    """GetOrganelle's log verdict must yield the same circularity evidence the
    other assemblers are judged on. It writes either "circular genome" or
    "N scaffold(s)"; only the former is a closed molecule."""

    def log(self, *lines):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "OG5.hic.getorg1770.get_org.log.txt"
        path.write_text("\n".join(lines) + "\n")
        return [path]

    def test_circular_genome_is_circularised(self):
        circ, _, failed = mas.getorganelle_evidence_from_log(
            self.log("2026-07-27 10:54:00 - INFO: Result status of animal_mt: circular genome"))
        self.assertEqual(circ, "true")
        self.assertFalse(failed)

    def test_scaffold_count_is_not_circularised(self):
        # The audit-5 OG5 / OG64 / OG869 case: previously left blank, which let the
        # row skip not_circularised entirely.
        for verdict in ("1 scaffold(s)", "2 scaffold(s)", "8 scaffold(s)"):
            with self.subTest(verdict=verdict):
                circ, _, _ = mas.getorganelle_evidence_from_log(
                    self.log(f"2026-07-27 10:54:00 - INFO: Result status of animal_mt: {verdict}"))
                self.assertEqual(circ, "false")

    def test_incomplete_is_not_read_as_complete(self):
        # Regression: "complete" is a substring of "incomplete", and the old parser
        # tested for it first, so an incomplete result resolved as complete.
        circ, _, _ = mas.getorganelle_evidence_from_log(
            self.log("2026-07-27 10:54:00 - INFO: Result status of animal_mt: incomplete genome"))
        self.assertEqual(circ, "false")

    def test_last_verdict_wins(self):
        circ, _, _ = mas.getorganelle_evidence_from_log(self.log(
            "INFO: Result status of animal_mt: 3 scaffold(s)",
            "INFO: Disentangling failed: retrying",
            "INFO: Result status of animal_mt: circular genome",
        ))
        self.assertEqual(circ, "true")

    def test_missing_verdict_leaves_circularity_unknown(self):
        # A truncated log means unknown topology, not "not circular" -- the same
        # state a MitoHiFi row with no contig-stats sidecar is left in.
        circ, _, failed = mas.getorganelle_evidence_from_log(
            self.log("2026-07-27 10:54:00 - INFO: Assembling reads"))
        self.assertEqual(circ, "")
        self.assertFalse(failed)

    def test_error_line_is_a_failure_but_disentangling_is_not(self):
        _, _, failed = mas.getorganelle_evidence_from_log(
            self.log("2026-07-27 10:54:00 - ERROR: No animal_mt seed reads found!"))
        self.assertTrue(failed)
        _, _, failed = mas.getorganelle_evidence_from_log(self.log(
            "INFO: Disentangling failed: cannot resolve, retrying",
            "INFO: Result status of animal_mt: circular genome",
        ))
        self.assertFalse(failed)

    def test_base_coverage_parsed(self):
        _, cov, _ = mas.getorganelle_evidence_from_log(self.log(
            "INFO: Average animal_mt base-coverage = 335.9",
            "INFO: Result status of animal_mt: circular genome",
        ))
        self.assertAlmostEqual(cov, 335.9)


class StatusVocabularyTests(unittest.TestCase):
    """`status` must mean the same thing whichever assembler produced the row."""

    VALUES = {"complete", "manual_review", "failed", "superseded"}

    def finalise(self, row, failed=False):
        mas.apply_qc(row, THRESHOLDS)
        mas.finalise_status(row, failed=failed)
        return row["status"]

    def test_clean_row_is_complete(self):
        self.assertEqual(self.finalise(complete_row()), "complete")

    def test_blocking_reason_is_manual_review(self):
        self.assertEqual(self.finalise(complete_row(num_final_contigs="3")), "manual_review")

    def test_no_final_assembly_is_failed(self):
        self.assertEqual(self.finalise(complete_row(final_length_bp="")), "failed")

    def test_assembler_failure_beats_other_evidence(self):
        row = complete_row()
        self.assertEqual(self.finalise(row, failed=True), "failed")
        self.assertIn("failed_run", row["manual_review_reason"])

    def test_advisory_flags_do_not_block(self):
        # A complete-core assembly at low coverage stays complete, with the flag
        # retained in manual_review_reason for transparency.
        row = complete_row(mean_coverage="12")
        self.assertEqual(self.finalise(row), "complete")
        self.assertIn("low_mean_coverage", row["manual_review_reason"])

    def test_vocabulary_is_closed_and_has_no_circular(self):
        rows = [
            complete_row(),
            complete_row(circularised="false"),
            complete_row(num_final_contigs="3"),
            complete_row(final_length_bp=""),
            complete_row(num_cds="7", missing_genes="CO2;CO3"),
        ]
        produced = {self.finalise(row) for row in rows}
        self.assertTrue(produced <= self.VALUES, f"unexpected status values: {produced - self.VALUES}")
        self.assertNotIn("circular", produced)


class CrossAssemblerStatusTests(unittest.TestCase):
    """The point of the change: equivalent assemblies get the same status whichever
    assembler produced them. Goes through the real per-assembler parsers, since that
    is where the two vocabularies and the two evidence rules used to diverge."""

    SEQ = "ACGT" * 4200   # 16800 bp, inside the expected range

    def root(self):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        return Path(tmp.name)

    def getorganelle_row(self, circular):
        root = self.root()
        prefix = "OG1.hic.getorg1770"
        fasta = root / f"{prefix}.fasta"
        fasta.write_text(f">OG1\n{self.SEQ}\n")
        log = root / f"{prefix}.get_org.log.txt"
        verdict = "circular genome" if circular else "1 scaffold(s)"
        log.write_text(f"INFO: Result status of animal_mt: {verdict}\n")
        files = [fasta, log]
        run = mas.RunFiles(sample_id="OG1", prefix=prefix, assembler="GetOrganelle", files=files)
        return mas.parse_getorganelle_run(run, THRESHOLDS, files)

    def mitohifi_row(self, circular):
        root = self.root()
        prefix = "OG1.hifi.v323mitohifi"
        fasta = root / f"{prefix}.fasta"
        fasta.write_text(f">OG1\n{self.SEQ}\n")
        stats = root / f"{prefix}.contigs_stats.tsv"
        stats.write_text(
            "contig_id\tlength_bp\twas_circular\tselected\n"
            f"OG1\t16800\t{'True' if circular else 'False'}\tTrue\n"
        )
        files = [fasta, stats]
        run = mas.RunFiles(sample_id="OG1", prefix=prefix, assembler="MitoHiFi", files=files)
        return mas.parse_mitohifi_run(run, THRESHOLDS, files)

    def oatk_row(self, circular):
        root = self.root()
        prefix = "OG1.hifi.oatk"
        fasta = root / f"{prefix}.fasta"
        fasta.write_text(f">OG1\n{self.SEQ}\n")
        gfa = root / f"{prefix}.gfa"
        # A self-link (from-segment == to-segment) is oatk's circularity marker.
        gfa.write_text(
            "S\tu1\t*\n" + ("L\tu1\t+\tu1\t+\t0M\n" if circular else "L\tu1\t+\tu2\t+\t0M\n")
        )
        files = [fasta, gfa]
        run = mas.RunFiles(sample_id="OG1", prefix=prefix, assembler="Oatk", files=files)
        return mas.parse_oatk_run(run, THRESHOLDS, files)

    def rows(self, circular):
        return {
            "GetOrganelle": self.getorganelle_row(circular),
            "MitoHiFi": self.mitohifi_row(circular),
            "Oatk": self.oatk_row(circular),
        }

    def test_circular_assemblies_agree(self):
        rows = self.rows(circular=True)
        self.assertEqual({name: row["circularised"] for name, row in rows.items()},
                         {"GetOrganelle": "true", "MitoHiFi": "true", "Oatk": "true"})
        self.assertEqual({row["status"] for row in rows.values()}, {"complete"})

    def test_non_circular_assemblies_agree(self):
        # Before the change GetOrganelle reported "complete" here (blank circularity
        # skipped not_circularised) while MitoHiFi reported manual_review.
        rows = self.rows(circular=False)
        self.assertEqual({name: row["circularised"] for name, row in rows.items()},
                         {"GetOrganelle": "false", "MitoHiFi": "false", "Oatk": "false"})
        self.assertEqual({row["status"] for row in rows.values()}, {"manual_review"})
        for name, row in rows.items():
            self.assertIn("not_circularised", row["manual_review_reason"], name)

    def test_no_assembler_emits_a_private_status_value(self):
        for circular in (True, False):
            for name, row in self.rows(circular).items():
                with self.subTest(assembler=name, circular=circular):
                    self.assertIn(row["status"], StatusVocabularyTests.VALUES)


class GetOrganelleRowStatusTests(unittest.TestCase):
    """End-to-end over parse_getorganelle_run, covering the audit-5 cases."""

    def make_run(self, verdict, check_row=None):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        prefix = "OG1161.hic.250624.getorg1770reseed"
        files = []

        fasta = root / f"{prefix}.fasta"
        fasta.write_text(">OG1161\n" + ("ACGT" * 4200) + "\n")   # 16800 bp, in range
        files.append(fasta)

        log = root / f"{prefix}.get_org.log.txt"
        log.write_text(
            "INFO: Average animal_mt base-coverage = 451.7\n"
            f"INFO: Result status of animal_mt: {verdict}\n"
        )
        files.append(log)

        if check_row is not None:
            check = root / f"{prefix}.getorg_check.tsv"
            check.write_text(
                "sample\tgetorg_circular\tfinal_verdict_circular\n"
                f"{prefix}\tFalse\t{check_row}\n"
            )
            files.append(check)

        run = mas.RunFiles(sample_id="OG1161", prefix=prefix, assembler="GetOrganelle", files=files)
        return mas.parse_getorganelle_run(run, THRESHOLDS, files)

    def test_circular_log_is_complete(self):
        row = self.make_run("circular genome")
        self.assertEqual(row["circularised"], "true")
        self.assertEqual(row["status"], "complete")

    def test_scaffold_log_without_check_is_manual_review(self):
        # OG5 / OG64 / OG869-like: never circularity-checked, so it must not sit at
        # the same status as a confirmed MitoHiFi assembly.
        row = self.make_run("1 scaffold(s)")
        self.assertEqual(row["circularised"], "false")
        self.assertEqual(row["status"], "manual_review")
        self.assertIn("not_circularised", row["manual_review_reason"])

    def test_check_sidecar_overrides_scaffold_verdict(self):
        # OG1161-like: the reference test confirms the scaffold is a linearised circle.
        row = self.make_run("1 scaffold(s)", check_row="True")
        self.assertEqual(row["circularised"], "true")
        self.assertEqual(row["status"], "complete")
        self.assertNotIn("not_circularised", row["manual_review_reason"])

    def test_no_row_reports_circular_status(self):
        for verdict in ("circular genome", "1 scaffold(s)"):
            with self.subTest(verdict=verdict):
                self.assertNotEqual(self.make_run(verdict)["status"], "circular")


DEPTH_HEADER = (
    "sample\ttarget_fasta\ttarget_length_bp\tn_contigs\tcircular_doubled\tpreset\t"
    "sequencing_type\tmean_depth\tmedian_depth\tsd_depth\tdepth_cv\tp10_depth\t"
    "min_depth\tbreadth_1x\tbreadth_10x\tbreadth_20x\tmito_mapped_reads\t"
    "reads_fail_identity\treads_fail_clip\tsupplementary_dropped\ttotal_reads\t"
    "mito_read_fraction\tmean_identity\tmin_identity\tmin_aligned_frac\tsubsampled\t"
    "subsample_fraction\tscale_factor\tdepth_method\tmean_coverage\tcoverage_cv"
)


def depth_tsv_text(prefix, mean_depth="412.5", cv="0.11"):
    return DEPTH_HEADER + "\n" + "\t".join([
        prefix, prefix + ".fasta", "16500", "1", "true", "sr", "ilmn",
        mean_depth, mean_depth, "45", cv, "380", "300", "1", "1", "1",
        "120000", "300", "12", "0", "900000", "0.133", "0.998", "0.95",
        "0.8", "false", "NA", "1", "remap_full_v1", mean_depth, cv,
    ]) + "\n"


class UniformDepthPrecedenceTests(unittest.TestCase):
    """The remap-based depth must win over every legacy coverage source."""

    def write(self, tmp, name, text):
        path = Path(tmp) / name
        path.write_text(text)
        return path

    def test_mito_depth_beats_getorganelle_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = "OG1.ilmn.240101.getorg1770"
            log = self.write(tmp, prefix + ".get_org.log.txt",
                             "INFO: Average animal_mt base-coverage = 335.9\n")
            depth = self.write(tmp, prefix + ".mito_depth.tsv", depth_tsv_text(prefix))
            mean, cv = mas.parse_coverage([log, depth])
            self.assertAlmostEqual(mean, 412.5)
            self.assertAlmostEqual(cv, 0.11)

    def test_mito_depth_beats_contigs_stats_with_coverage(self):
        """Regression: the old loose glob let with_coverage.tsv win by file order."""
        with tempfile.TemporaryDirectory() as tmp:
            prefix = "OG1.hifi.240101.v323mitohifi"
            stats = self.write(
                tmp, prefix + ".contigs_stats.with_coverage.tsv",
                "contig_id\twas_circular\tavg_coverage\tcoverage_cv\n"
                "final_mitogenome\tTrue\t88.2\t0.19\n",
            )
            depth = self.write(tmp, prefix + ".mito_depth.tsv", depth_tsv_text(prefix))
            # Both orderings must give the same answer.
            for files in ([stats, depth], [depth, stats]):
                mean, _cv = mas.parse_coverage(files)
                self.assertAlmostEqual(mean, 412.5)

    def test_with_coverage_reads_final_mitogenome_row_not_row_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = self.write(
                tmp, "OG1.hifi.240101.v323mitohifi.contigs_stats.with_coverage.tsv",
                "contig_id\twas_circular\tavg_coverage\tcoverage_cv\n"
                "ptg000001l\tTrue\t\t\n"
                "final_mitogenome\tTrue\t88.2\t0.19\n",
            )
            mean, cv = mas.parse_coverage([stats])
            self.assertAlmostEqual(mean, 88.2)
            self.assertAlmostEqual(cv, 0.19)

    def test_depth_suffix_stripped_so_no_phantom_run(self):
        """A reseed depth must join the reseed run, not spawn <prefix>.mito_depth."""
        self.assertEqual(
            mas.strip_known_suffix("OG1.ilmn.240101.getorg1770reseed.mito_depth.tsv"),
            "OG1.ilmn.240101.getorg1770reseed",
        )

    def test_depth_tsv_does_not_trip_numt_flag(self):
        """has_numt_signal greps .tsv text; the depth columns must be inert."""
        with tempfile.TemporaryDirectory() as tmp:
            prefix = "OG1.ilmn.240101.getorg1770"
            depth = self.write(tmp, prefix + ".mito_depth.tsv", depth_tsv_text(prefix))
            self.assertFalse(mas.has_numt_signal([depth]))


class CoverageVariabilityAdvisoryTests(unittest.TestCase):
    """high_coverage_variability must not newly fail finished mitogenomes.

    Before the uniform remap this reason could barely fire: coverage_cv was
    populated for MitoHiFi alone, and even there it was dominated by the
    triangular profile that mapping 15-20 kb reads to a LINEAR reference
    produces. Now that it is measured properly for all three assemblers it has to
    behave like low_mean_coverage: advisory on a finished mitogenome, blocking on
    anything less.
    """

    def test_high_cv_on_complete_assembly_is_advisory(self):
        row = complete_row(coverage_cv="2.5")
        mas.apply_qc(row, THRESHOLDS)
        mas.finalise_status(row)
        self.assertIn("high_coverage_variability", row["manual_review_reason"])
        self.assertEqual(row["status"], "complete")

    def test_high_cv_on_incomplete_assembly_still_blocks(self):
        row = complete_row(coverage_cv="2.5", circularised="false")
        mas.apply_qc(row, THRESHOLDS)
        mas.finalise_status(row)
        self.assertIn("high_coverage_variability", row["manual_review_reason"])
        self.assertEqual(row["status"], "manual_review")


class AnnotationJoinTests(unittest.TestCase):
    """Gene counts must come from the assembly's OWN re-annotation.

    Every case here is a real misattribution from the mitogenomes-missing-audit-6
    cohort, caused by the join being a substring test: og_id "OG5" is a substring
    of "OG58", "OG8" of "OG810" and "OG848", and code "getorg1770" of
    "getorg1770reseed".
    """

    HEADER = "og_id,tech,seq_date,code,annotation,missing_genes,num_cds,num_trna,num_rrna\n"

    def write(self, name, og_id, tech, seq_date, code, missing="no", cds=13, trna=22, rrna=2):
        import tempfile
        if not hasattr(self, "_root"):
            tmp = tempfile.TemporaryDirectory()
            self.addCleanup(tmp.cleanup)
            self._root = Path(tmp.name)
        path = self._root / name
        path.write_text(
            self.HEADER
            + f"{og_id},{tech},{seq_date},{code},emma102,{missing},{cds},{trna},{rrna}\n"
        )
        return path

    def test_shorter_sample_id_does_not_bleed_into_a_longer_one(self):
        og5 = self.write("OG5.hifi.230609.v323mitohifi.annotation_stats.csv",
                         "OG5", "hifi", "230609", "v323mitohifi")
        self.assertEqual(mas.parse_annotation_stats([og5], "OG58.hifi.250704.v323mitohifi"), {})

    def test_base_assembly_does_not_bleed_into_its_reseed(self):
        base = self.write("OG838.hic.250522.getorg1770.annotation_stats.csv",
                          "OG838", "hic", "250522", "getorg1770")
        self.assertEqual(mas.parse_annotation_stats([base], "OG838.hic.250522.getorg1770reseed"), {})

    def test_collapsed_variant_does_not_bleed_into_its_parent(self):
        child = self.write("OG750.hifi.241004.v323mitohifi_collapsed.annotation_stats.csv",
                           "OG750", "hifi", "241004", "v323mitohifi_collapsed",
                           missing="TS2;TD;CO2;TK;ATP8", cds=11, trna=19, rrna=2)
        self.assertEqual(mas.parse_annotation_stats([child], "OG750.hifi.241004.v323mitohifi"), {})
        self.assertEqual(
            mas.parse_annotation_stats([child], "OG750.hifi.241004.v323mitohifi_collapsed"),
            {
                "num_genes": "32",
                "num_cds": "11",
                "missing_genes": "TS2;TD;CO2;TK;ATP8",
                "frameshift_flag": "",
            },
        )

    def test_degenerate_row_matches_nothing(self):
        # bin/annotation_stats.py blanks og_id/code when the GFF stem is not five
        # dot-fields. A blank key must match nothing, not everything.
        blank = self.write("weird.annotation_stats.csv", "", "", "", "")
        self.assertEqual(mas.parse_annotation_stats([blank], "OG58.hifi.250704.v323mitohifi"), {})

    def test_own_annotation_is_found_by_filename(self):
        mine = self.write("OG58.hifi.250704.v10oatk.annotation_stats.csv",
                          "OG58", "hifi", "250704", "v10oatk")
        stats = mas.parse_annotation_stats([mine], "OG58.hifi.250704.v10oatk")
        self.assertEqual(stats["num_genes"], "37")
        self.assertEqual(stats["num_cds"], "13")

    def test_content_match_is_used_when_the_file_was_not_renamed(self):
        mine = self.write("emma_out.annotation_stats.csv", "OG58", "hifi", "250704", "v10oatk")
        stats = mas.parse_annotation_stats([mine], "OG58.hifi.250704.v10oatk")
        self.assertEqual(stats["num_cds"], "13")
        self.assertEqual(mas.parse_annotation_stats([mine], "OG58.hifi.250704.v323mitohifi"), {})


class RunDiscoveryTests(unittest.TestCase):
    """A per-run sidecar must not manufacture an assembly of its own. The two
    reference sidecars below produced 38 of the 250 rows in the
    mitogenomes-missing-audit-6 table, every one of them reported as `failed`."""

    def test_reference_sidecars_do_not_spawn_runs(self):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        prefix = "OG750.hifi.241004.v323mitohifi"
        (root / f"{prefix}.fasta").write_text(">c\n" + "ACGT" * 4000 + "\n")
        for sidecar in ("reference_ranking", "reference_candidates_status"):
            (root / f"{prefix}.{sidecar}.tsv").write_text("sample\tvalue\nOG750\t1\n")

        prefixes = [run.prefix for run in mas.discover_assembler_runs([root])]
        self.assertEqual(prefixes, [prefix])

    def test_unknown_sidecar_suffix_is_still_rejected(self):
        # The generic net: strip_known_suffix cannot list a suffix nobody has
        # written yet, so the shape of the prefix has to carry the decision.
        self.assertFalse(mas.is_assembly_run_prefix("OG750.hifi.241004.v323mitohifi.some_new_check"))

    def test_real_assembly_prefixes_survive(self):
        for prefix in (
            "OG750.hifi.241004.v323mitohifi",
            "OG750.hifi.241004.v323mitohifi_collapsed",
            "OG838.hic.250522.getorg1770",
            "OG838.hic.250522.getorg1770reseed",
            "OG778.hic.250624.getorg1770reseed_rgj",
            "OG58.hifi.250704.v10oatk",
            "OG750.hifi.v323mitohifi",
        ):
            self.assertTrue(mas.is_assembly_run_prefix(prefix), prefix)


class GeneCountProvenanceTests(unittest.TestCase):
    """num_genes and num_cds come from the same annotation or from neither."""

    def test_contigs_stats_gene_count_is_not_reported(self):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        stats = root / "OG750.hifi.241004.v323mitohifi.contigs_stats.tsv"
        stats.write_text(
            "contig_id\tframeshifts_found\tannotation_file\tlength(bp)\tnumber_of_genes\twas_circular\n"
            "final_mitogenome\tNo frameshift found\tfinal_mitogenome.gb\t32672\t56\tTrue\n"
        )
        parsed = mas.parse_mitohifi_stats([stats])
        # 56 is MitoHiFi's reference-guided count over the un-collapsed 2.14x
        # concatemer; the pipeline's gene count comes from EMMA / MITOS2 only.
        self.assertNotIn("num_genes", parsed)
        # The rest of the sidecar is still read.
        self.assertEqual(parsed["circularised"], "true")
        self.assertEqual(parsed["frameshift_flag"], "false")

    def test_unannotated_run_reports_neither_count(self):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        prefix = "OG750.hifi.241004.v323mitohifi"
        (root / f"{prefix}.fasta").write_text(">c\n" + "ACGT" * 4000 + "\n")
        (root / f"{prefix}.contigs_stats.tsv").write_text(
            "contig_id\tframeshifts_found\tannotation_file\tlength(bp)\tnumber_of_genes\twas_circular\n"
            "final_mitogenome\tNo frameshift found\tfinal_mitogenome.gb\t16000\t37\tTrue\n"
        )
        runs = mas.discover_assembler_runs([root])
        self.assertEqual(len(runs), 1)
        row = mas.parse_mitohifi_run(runs[0], THRESHOLDS, mas.collect_input_files([root]))
        self.assertEqual(row["num_genes"], "")
        self.assertEqual(row["num_cds"], "")


class PrefixBoundaryTests(unittest.TestCase):
    """The shared anchoring primitive: a prefix must end at a name/component
    boundary, so it never matches the curated variants built on top of it."""

    PARENT = "OG750.hifi.241004.v323mitohifi"

    def test_curated_variant_does_not_match_its_parent(self):
        for variant in ("_collapsed", "reseed", "reseed_rgj"):
            path = Path(f"/w/{self.PARENT}{variant}.annotation_stats.csv")
            self.assertFalse(mas.path_carries_prefix(path, self.PARENT), variant)

    def test_own_files_match(self):
        self.assertTrue(
            mas.path_carries_prefix(Path(f"/w/{self.PARENT}.contigs_stats.tsv"), self.PARENT))
        self.assertTrue(
            mas.path_carries_prefix(Path(f"/w/{self.PARENT}"), self.PARENT))
        self.assertTrue(
            mas.path_carries_prefix(Path(f"/w/{self.PARENT}/mtdna/final.fasta"), self.PARENT))

    def test_child_prefix_matches_only_its_own_files(self):
        child = f"{self.PARENT}_collapsed"
        self.assertTrue(
            mas.path_carries_prefix(Path(f"/w/{child}.annotation_stats.csv"), child))
        self.assertFalse(
            mas.path_carries_prefix(Path(f"/w/{self.PARENT}.contigs_stats.tsv"), child))


if __name__ == "__main__":
    unittest.main()
