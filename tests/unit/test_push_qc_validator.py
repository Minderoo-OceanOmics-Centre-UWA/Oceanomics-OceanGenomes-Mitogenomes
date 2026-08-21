import csv
import importlib.util
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
# push_qc_validator imports load_db_config from its sibling, exactly as it does
# at runtime with bin/ on PATH, so bin/ has to be importable here too.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "push_qc_validator", ROOT / "bin" / "push_qc_validator.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


RECORD_COLUMNS = MODULE.IDENTITY_COLUMNS + ["submission_ready"]


class FakeCursor:
    """Records every statement; the UPDATE's rowcount and the SELECT's row are injected."""

    def __init__(self, update_rowcount, existing_row):
        self.update_rowcount = update_rowcount
        self.existing_row = existing_row
        self.rowcount = 0
        self.queries = []
        self.params = []

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, params):
        self.queries.append(query)
        self.params.append(params)
        self.rowcount = self.update_rowcount if "UPDATE" in query else 1

    def fetchone(self):
        return self.existing_row


class FakeConnection:
    def __init__(self, update_rowcount=1, existing_row=None, fail=False):
        self.cursor_instance = FakeCursor(update_rowcount, existing_row)
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


def make_record(root, submission_ready="true"):
    row = {
        "full_seqid": "OG1.hifi.260101.v3mitohifi.emma102",
        "og_id": "OG1",
        "tech": "hifi",
        "seq_date": "260101",
        "code": "v3mitohifi",
        "annotation": "emma102",
        "submission_ready": submission_ready,
    }
    path = root / "record.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    return path


def read_record(submission_ready="true"):
    with tempfile.TemporaryDirectory() as tmp:
        return MODULE.read_record(make_record(Path(tmp), submission_ready))


def run_main(config_path, record_path):
    """Invoke main() with argv patched, capturing stdout."""
    argv = sys.argv
    sys.argv = ["push_qc_validator.py", str(config_path), str(record_path)]
    buffer = io.StringIO()
    try:
        with redirect_stdout(buffer):
            code = MODULE.main()
    finally:
        sys.argv = argv
    return code, buffer.getvalue()


class ReadRecordTests(unittest.TestCase):
    def test_identity_and_readiness_are_parsed(self):
        record = read_record()
        self.assertTrue(record["submission_ready"])
        self.assertEqual(record["og_id"], "OG1")
        self.assertEqual(record["annotation"], "emma102")

    def test_submission_ready_false_is_falsey(self):
        self.assertFalse(read_record("false")["submission_ready"])


class SetQcValidatorTests(unittest.TestCase):
    def test_blank_validator_2_is_filled(self):
        record = read_record()
        connection = FakeConnection(update_rowcount=1)
        status, detail = MODULE.set_qc_validator(
            record, {"dbname": "test"}, connect=lambda **_kw: connection
        )
        self.assertEqual((status, detail), ("set", "QCd-nf-core"))
        self.assertTrue(connection.committed)
        self.assertTrue(connection.closed)
        # One statement only: no SELECT is needed when the UPDATE lands.
        self.assertEqual(len(connection.cursor_instance.queries), 1)
        query = connection.cursor_instance.queries[0]
        self.assertIn("UPDATE lca_validation", query)
        # The guard is the whole point: the UPDATE must be unable to overwrite.
        self.assertIn("validator_2 IS NULL OR btrim(validator_2) = ''", query)
        params = connection.cursor_instance.params[0]
        self.assertEqual(params["validator_2"], "QCd-nf-core")
        self.assertEqual(params["annotation"], "emma102")

    def test_existing_validator_2_is_preserved(self):
        record = read_record()
        connection = FakeConnection(update_rowcount=0, existing_row=("nf-core", "TP"))
        status, detail = MODULE.set_qc_validator(
            record, {}, connect=lambda **_kw: connection
        )
        self.assertEqual((status, detail), ("preserved", "TP"))

    def test_missing_row_is_reported_not_inserted(self):
        record = read_record()
        connection = FakeConnection(update_rowcount=0, existing_row=None)
        status, detail = MODULE.set_qc_validator(
            record, {}, connect=lambda **_kw: connection
        )
        self.assertEqual((status, detail), ("no_row", None))
        self.assertNotIn(
            "INSERT", " ".join(connection.cursor_instance.queries)
        )

    def test_database_failure_rolls_back_and_closes(self):
        record = read_record()
        connection = FakeConnection(fail=True)
        with self.assertRaisesRegex(RuntimeError, "database unavailable"):
            MODULE.set_qc_validator(record, {}, connect=lambda **_kw: connection)
        self.assertTrue(connection.rolled_back)
        self.assertTrue(connection.closed)


class MainMarkerTests(unittest.TestCase):
    """The stdout markers are parsed by compile_upload_report.py, so pin them."""

    def test_not_submission_ready_issues_no_sql(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            record = make_record(root, submission_ready="false")

            def explode(**_kw):
                raise AssertionError("no database connection should be opened")

            original = MODULE.set_qc_validator
            MODULE.set_qc_validator = explode
            try:
                code, out = run_main(root / "missing.cfg", record)
            finally:
                MODULE.set_qc_validator = original
        self.assertEqual(code, 0)
        self.assertIn("not submission_ready", out)

    def test_database_error_marker_and_exit_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            record = make_record(root)
            # No config file exists, so load_db_config raises before any connect.
            code, out = run_main(root / "missing.cfg", record)
        self.assertEqual(code, 1)
        self.assertIn("❌ Database error", out)


if __name__ == "__main__":
    unittest.main()
