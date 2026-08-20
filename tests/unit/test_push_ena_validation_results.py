import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "push_ena_validation_results", ROOT / "bin" / "push_ena_validation_results.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    def __init__(self, row):
        self.row = row
        self.query = None
        self.params = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, params):
        self.query = query
        self.params = params

    def fetchone(self):
        return self.row


class FakeConnection:
    def __init__(self, row=(1, True), fail=False):
        self.cursor_instance = FakeCursor(row)
        self.fail = fail
        self.committed = self.rolled_back = self.closed = False

    def cursor(self):
        if self.fail:
            raise RuntimeError("database unavailable")
        return self.cursor_instance

    def commit(self):
        self.committed = True

    def rollback(self):
        self.rolled_back = True

    def close(self):
        self.closed = True


class PushEnaValidationTests(unittest.TestCase):
    def make_record(self, root):
        row = {column: "" for column in MODULE.INSERT_COLUMNS}
        row.update(
            full_seqid="OG1.hifi.260101.final.emma102", og_id="OG1", annotation="emma102",
            ena_study="PRJEB1",
            validation_mode="pipeline", validation_attempt="one", table2asn_status="PASS",
            conversion_status="PASS", preflight_status="NOT_APPLICABLE", webin_status="PASS",
            submission_ready="true", reject_count="0",
        )
        path = root / "record.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.INSERT_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        return path

    def test_parameter_mapping_and_insert(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = MODULE.read_record(self.make_record(Path(tmp)))
        self.assertTrue(record["submission_ready"])
        self.assertEqual(record["reject_count"], 0)
        connection = FakeConnection(row=(1, True))
        result = MODULE.upload_record(record, {"dbname": "test"}, connect=lambda **_kw: connection)
        self.assertEqual(result, "inserted")
        self.assertTrue(connection.committed)
        self.assertIn("ON CONFLICT", connection.cursor_instance.query)
        self.assertIn("DO UPDATE", connection.cursor_instance.query)
        self.assertNotIn("validation_attempt = EXCLUDED", connection.cursor_instance.query)
        self.assertEqual(connection.cursor_instance.params["full_seqid"], record["full_seqid"])
        self.assertIn("ON CONFLICT (full_seqid, ena_study, validation_attempt)", connection.cursor_instance.query)
        self.assertEqual(connection.cursor_instance.params["annotation"], "emma102")

    def test_conflict_on_earlier_attempt_is_updated(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = MODULE.read_record(self.make_record(Path(tmp)))
        connection = FakeConnection(row=(1, False))
        self.assertEqual(
            MODULE.upload_record(record, {}, connect=lambda **_kw: connection), "updated"
        )

    def test_upsert_does_not_consult_the_retired_selection_tables(self):
        """Submission state lives in a separate pipeline; nothing freezes a row here."""
        with tempfile.TemporaryDirectory() as tmp:
            record = MODULE.read_record(self.make_record(Path(tmp)))
        connection = FakeConnection(row=(1, False))
        MODULE.upload_record(record, {}, connect=lambda **_kw: connection)
        self.assertNotIn("ena_submission_selections", connection.cursor_instance.query)
        self.assertNotIn("WHERE NOT EXISTS", connection.cursor_instance.query)

    def test_database_failure_rolls_back_and_closes(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = MODULE.read_record(self.make_record(Path(tmp)))
        connection = FakeConnection(fail=True)
        with self.assertRaisesRegex(RuntimeError, "database unavailable"):
            MODULE.upload_record(record, {}, connect=lambda **_kw: connection)
        self.assertTrue(connection.rolled_back)
        self.assertTrue(connection.closed)


if __name__ == "__main__":
    unittest.main()
