"""Unit tests for the mitogenome assembly-summary QC classification.

Covers the Phase 1 robustness changes: the blocking-vs-advisory split
(is_complete_core / blocking_reasons) and the rewritten ambiguous-graph
detector (getorganelle_graph_ambiguous). Rows are modelled on the real
manual-review samples from the mitogenomes-missing-audit-3 run.
"""
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
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


class CollapseOverrideTests(unittest.TestCase):
    def make_run(self, report_text):
        import tempfile
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        rep = root / "OG750.hifi.v323mitohifi.concatemer_collapse.tsv"
        rep.write_text(report_text)
        return mas.RunFiles(sample_id="OG750", prefix="OG750.hifi.v323mitohifi",
                            assembler="MitoHiFi", files=[rep])

    def test_collapsed_report_overrides_length_and_clears_anomaly(self):
        run = self.make_run(
            "sample\taction\toriginal_length\tcollapsed_length\treference_length\ttail_identity\treason\n"
            "OG750\tcollapsed\t32672\t16466\t16703\t1.0\tcollapsed_2.0x_concatemer\n"
        )
        row = complete_row(final_length_bp="32672", anomaly_type="concatemer", length_anomaly="yes")
        mas.apply_collapse_override(row, run)
        self.assertEqual(row["final_length_bp"], "16466")
        self.assertEqual(row["anomaly_type"], "none")
        # And the resolved assembly must no longer be blocked.
        mas.apply_qc(row, THRESHOLDS)
        self.assertNotIn("concatemer", row.get(mas.BLOCKING_KEY, ""))
        self.assertNotIn("length_outside_expected_range", row.get(mas.BLOCKING_KEY, ""))

    def test_passthrough_report_leaves_row_untouched(self):
        run = self.make_run(
            "sample\taction\toriginal_length\tcollapsed_length\treference_length\ttail_identity\treason\n"
            "OGX\tpassthrough\t32000\t32000\t16703\t0.0\ttail_identity_below_threshold\n"
        )
        row = complete_row(final_length_bp="32000", anomaly_type="unresolved", length_anomaly="yes")
        mas.apply_collapse_override(row, run)
        self.assertEqual(row["final_length_bp"], "32000")
        self.assertEqual(row["anomaly_type"], "unresolved")


if __name__ == "__main__":
    unittest.main()
