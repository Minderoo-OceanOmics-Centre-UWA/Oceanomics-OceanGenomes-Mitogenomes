import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
# push_assembly_qc imports load_db_config from its sibling, exactly as it does at
# runtime with bin/ on PATH, so bin/ has to be importable here too.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location("push_assembly_qc", ROOT / "bin" / "push_assembly_qc.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

HEADER = ["assembly_prefix", "sample_id", "assembler", "status", "final_length_bp", "circularised",
          "num_candidate_contigs", "num_final_contigs", "num_genes", "num_cds", "missing_genes",
          "frameshift_flag", "mean_coverage", "coverage_cv", "reference_species",
          "reference_accession", "numt_flag", "manual_review_reason"]


def summary(tmp, *rows):
    path = Path(tmp) / "mitogenome_assembly_summary_mqc.tsv"
    lines = ["\t".join(HEADER)]
    for row in rows:
        lines.append("\t".join(row.get(column, "") for column in HEADER))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class FakeCursor:
    def __init__(self, missing):
        self.missing = missing
        self.rowcount = 0
        self.params = []

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, _query, params):
        self.params.append(params)
        self.rowcount = 0 if params["og_id"] in self.missing else 1


class FakeConnection:
    def __init__(self, missing=()):
        self.cursor_obj = FakeCursor(set(missing))
        self.committed = self.rolled_back = False

    def cursor(self):
        return self.cursor_obj

    def commit(self):
        self.committed = True

    def rollback(self):
        self.rolled_back = True

    def close(self):
        pass


class ReadVerdictTests(unittest.TestCase):
    def test_parses_the_verdict_and_keys_the_assembly(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, skipped = MODULE.read_verdicts(summary(tmp, {
                "assembly_prefix": "OG3078.ilmn.260909.getorg1770reseed", "status": "manual_review",
                "manual_review_reason": "numt_signal;low_coverage", "numt_flag": "true",
                "frameshift_flag": "false", "reference_species": "Pagrus auratus",
                "reference_accession": "NC_000001.1"}))
        self.assertEqual(skipped, [])
        record = records[0]
        self.assertEqual((record["og_id"], record["tech"], record["seq_date"], record["code"]),
                         ("OG3078", "ilmn", "260909", "getorg1770reseed"))
        self.assertEqual(record["qc_status"], "manual_review")
        self.assertIs(record["numt_flag"], True)
        self.assertIs(record["frameshift_flag"], False)

    def test_blank_flags_and_text_become_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, _ = MODULE.read_verdicts(summary(tmp, {
                "assembly_prefix": "OG1.hifi.250101.v323mitohifi", "status": "complete"}))
        record = records[0]
        for column in ("numt_flag", "frameshift_flag", "manual_review_reason",
                       "reference_species", "reference_accession"):
            self.assertIsNone(record[column], column)

    def test_prefix_without_four_fields_is_skipped_by_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, skipped = MODULE.read_verdicts(summary(tmp, {
                "assembly_prefix": "OG750.hifi.v323mitohifi", "status": "complete"}))
        self.assertEqual(records, [])
        self.assertEqual(skipped, ["OG750.hifi.v323mitohifi"])

    def test_unknown_status_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                MODULE.read_verdicts(summary(tmp, {
                    "assembly_prefix": "OG1.ilmn.250101.getorg1770", "status": "done"}))


class PushTests(unittest.TestCase):
    RECORD = {"og_id": "OG1", "tech": "ilmn", "seq_date": "250101", "code": "getorg1770",
              "qc_status": "complete", "manual_review_reason": None, "numt_flag": False,
              "frameshift_flag": None, "reference_species": None, "reference_accession": None}

    def test_updates_existing_rows_and_reports_missing_ones(self):
        connection = FakeConnection(missing={"OG2"})
        records = [self.RECORD, dict(self.RECORD, og_id="OG2")]
        outcome = MODULE.push_verdicts(records, {}, connect=lambda **_: connection)
        self.assertEqual(outcome["set"], ["OG1.ilmn.250101.getorg1770"])
        self.assertEqual(outcome["no_row"], ["OG2.ilmn.250101.getorg1770"])
        self.assertTrue(connection.committed)

    def test_dry_run_rolls_back(self):
        connection = FakeConnection()
        MODULE.push_verdicts([self.RECORD], {}, dry_run=True, connect=lambda **_: connection)
        self.assertTrue(connection.rolled_back)
        self.assertFalse(connection.committed)

    def test_query_never_inserts(self):
        self.assertNotIn("INSERT", MODULE.UPDATE_QUERY.upper())


if __name__ == "__main__":
    unittest.main()
