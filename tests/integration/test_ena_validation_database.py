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


@unittest.skipUnless(os.environ.get("ENA_TEST_DB_DSN"), "ENA_TEST_DB_DSN is not set")
class EnaValidationDatabaseIntegrationTests(unittest.TestCase):
    def test_migration_deduplication_and_changed_history(self):
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        import push_ena_validation_results as uploader

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_validation_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET LOCAL search_path TO "{schema}"')
                cursor.execute((ROOT / "sql" / "001_create_ena_validation_attempts.sql").read_text())

            record = {column: None for column in uploader.INSERT_COLUMNS}
            record.update(
                assembly_prefix="OG1.hifi.260101.final", og_id="OG1", ena_study="PRJEB1",
                validation_mode="pipeline", validation_attempt="integration",
                table2asn_status="PASS", conversion_status="PASS",
                preflight_status="NOT_APPLICABLE", webin_status="PASS",
                submission_ready=True, result_digest="a" * 64,
            )

            # Keep one transaction open so the temporary schema and all inserts
            # are discarded by the final rollback.
            with connection.cursor() as cursor:
                columns = ", ".join(uploader.INSERT_COLUMNS)
                values = ", ".join(f"%({column})s" for column in uploader.INSERT_COLUMNS)
                query = (
                    f"INSERT INTO ena_validation_attempts ({columns}) VALUES ({values}) "
                    "ON CONFLICT (assembly_prefix, ena_study, validation_attempt, result_digest) "
                    "DO NOTHING RETURNING id"
                )
                cursor.execute(query, record)
                self.assertIsNotNone(cursor.fetchone())
                cursor.execute(query, record)
                self.assertIsNone(cursor.fetchone())
                record["result_digest"] = "b" * 64
                cursor.execute(query, record)
                self.assertIsNotNone(cursor.fetchone())
                cursor.execute("SELECT count(*) FROM ena_validation_attempts")
                self.assertEqual(cursor.fetchone()[0], 2)
        finally:
            connection.rollback()
            connection.close()


if __name__ == "__main__":
    unittest.main()
