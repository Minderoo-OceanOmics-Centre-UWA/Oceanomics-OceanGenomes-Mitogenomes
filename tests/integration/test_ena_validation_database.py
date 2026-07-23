"""Opt-in PostgreSQL integration test.

Set ENA_TEST_DB_DSN to a disposable/test database connection string. The test
uses a temporary schema and rolls the migration and inserts back.
"""

import os
import sys
import unittest
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))


class _KeepOpenConnection:
    """Proxies a real connection but no-ops commit/rollback/close.

    push_ena_validation_results.upload_record() commits and closes the
    connection it's given, since it's normally a one-shot CLI invocation. This
    lets the test call it repeatedly on one connection and keep everything in
    a single transaction that the test itself rolls back at the end.
    """

    def __init__(self, real_connection):
        self._real = real_connection

    def cursor(self):
        return self._real.cursor()

    def commit(self):
        pass

    def rollback(self):
        pass

    def close(self):
        pass


@unittest.skipUnless(os.environ.get("ENA_TEST_DB_DSN"), "ENA_TEST_DB_DSN is not set")
class EnaValidationDatabaseIntegrationTests(unittest.TestCase):
    def test_migration_overwrites_until_submission_ready(self):
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        import push_ena_validation_results as uploader

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_validation_test_{uuid.uuid4().hex}"
        keep_open = _KeepOpenConnection(connection)
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                cursor.execute((ROOT / "sql" / "001_create_ena_validation_attempts.sql").read_text())
                cursor.execute(
                    (ROOT / "sql" / "002_ena_validation_attempts_single_row_per_attempt.sql").read_text()
                )

            record = {column: None for column in uploader.INSERT_COLUMNS}
            record.update(
                assembly_prefix="OG1.hifi.260101.final", og_id="OG1", ena_study="PRJEB1",
                validation_mode="pipeline", validation_attempt="integration",
                table2asn_status="PASS", conversion_status="FAIL",
                preflight_status="NOT_APPLICABLE", webin_status="NOT_RUN",
                submission_ready=False, result_digest="a" * 64,
            )

            # First failed attempt: inserted.
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "inserted"
            )

            # Rerun with a different digest but still not ready: overwrites in place.
            record["result_digest"] = "b" * 64
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )
            with connection.cursor() as cursor:
                cursor.execute("SELECT count(*), max(attempt_count) FROM ena_validation_attempts")
                self.assertEqual(cursor.fetchone(), (1, 2))

            # Now it passes: still overwrites (it was not ready before).
            record.update(
                conversion_status="PASS", webin_status="PASS", submission_ready=True, result_digest="c" * 64
            )
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )

            # A later rerun under the same key, even a regressed/failed one, is locked out.
            record.update(webin_status="NOT_RUN", submission_ready=False, result_digest="d" * 64)
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "locked"
            )
            with connection.cursor() as cursor:
                cursor.execute("SELECT count(*), bool_and(submission_ready) FROM ena_validation_attempts")
                self.assertEqual(cursor.fetchone(), (1, True))
        finally:
            connection.rollback()
            connection.close()


if __name__ == "__main__":
    unittest.main()
