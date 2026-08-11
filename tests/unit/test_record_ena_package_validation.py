import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "record_ena_package_validation",
    ROOT / "bin" / "record_ena_package_validation.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    def __init__(self):
        self.calls = []
        self.fetchone_values = [("OG910.hifi.250101.v3mitohifi", "OG910", "PRJEB123419")]

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, values):
        self.calls.append((query, values))

    def fetchone(self):
        return self.fetchone_values.pop(0)


class FakeConnection:
    def __init__(self):
        self.cursor_instance = FakeCursor()

    def cursor(self):
        return self.cursor_instance


class RecordEnaPackageValidationTests(unittest.TestCase):
    def test_reads_structured_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "status.tsv"
            path.write_text(
                "full_seqid\tservice\tstatus\treason\twebin_exit\n"
                "OG910.hifi.250101.v3mitohifi.emma102\tproduction\tPASS\tvalidated\t0\n"
            )
            row = MODULE.read_status(path)
            self.assertEqual(row["service"], "production")

    def test_production_pass_updates_package_and_readiness(self):
        connection = FakeConnection()
        MODULE.persist(
            connection,
            {
                "full_seqid": "OG910.hifi.250101.v3mitohifi.emma102",
                "service": "production",
                "status": "PASS",
                "reason": "validated",
                "webin_exit": "0",
            },
        )
        self.assertEqual(len(connection.cursor_instance.calls), 2)
        self.assertIn(
            "webin_production_status",
            connection.cursor_instance.calls[0][0],
        )
        self.assertIn(
            "submission_ready",
            connection.cursor_instance.calls[1][0],
        )


if __name__ == "__main__":
    unittest.main()
