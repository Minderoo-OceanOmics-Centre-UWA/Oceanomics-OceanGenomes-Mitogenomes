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


def _migration_sql(name: str) -> str:
    """Migration text with its own transaction control removed.

    The files each wrap themselves in BEGIN/COMMIT for deliberate manual
    application. Executed as-is here, that COMMIT would commit the test's
    transaction and leave the throwaway schema behind, so strip both and let the
    test's single rollback undo everything.
    """
    text = (ROOT / "sql" / name).read_text()
    return text.replace("BEGIN;", "", 1).replace("COMMIT;", "", 1)


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
    def test_reruns_always_overwrite_the_attempt_row(self):
        """Nothing freezes a validation row now that selection lives elsewhere.

        Migration 012 drops ena_submission_selections along with the upsert
        guard that consulted it, so the latest validation of an assembly is
        always the one on record.
        """
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
                for migration in (
                    "001_create_ena_validation_attempts.sql",
                    "002_ena_validation_attempts_single_row_per_attempt.sql",
                    "004_ena_candidate_packages.sql",
                    "005_ena_validation_attempts_genome_context.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    # 010 has to be in the list: it drops ena_candidate_loci,
                    # whose foreign key onto ena_candidate_packages would
                    # otherwise block 012's drop of that table.  The real chain
                    # always runs 010 first, so skipping it here tested an
                    # ordering that cannot occur.
                    "010_drop_ena_locus_tables.sql",
                    "011_drop_local_package_validation.sql",
                    "012_drop_ena_selection_layer.sql",
                    "014_ena_validation_attempts_full_seqid.sql",
                ):
                    cursor.execute(_migration_sql(migration))

                # 012 takes the selection layer with it.
                for relation in (
                    "ena_candidate_packages",
                    "ena_submission_selections",
                    "ena_submission_queue",
                ):
                    cursor.execute("SELECT to_regclass(%s)", (f"{schema}.{relation}",))
                    self.assertIsNone(cursor.fetchone()[0], relation)

            record = {column: None for column in uploader.INSERT_COLUMNS}
            record.update(
                full_seqid="OG1.hifi.260101.final.emma102", og_id="OG1",
                annotation="emma102", ena_study="PRJEB1",
                validation_mode="pipeline", validation_attempt="integration",
                table2asn_status="PASS", conversion_status="FAIL",
                preflight_status="NOT_APPLICABLE", webin_status="NOT_RUN",
                submission_ready=False,
            )

            # First failed attempt: inserted.
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "inserted"
            )

            # Rerun of the same failure: overwrites in place rather than adding history.
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )
            with connection.cursor() as cursor:
                cursor.execute("SELECT count(*), max(attempt_count) FROM ena_validation_attempts")
                self.assertEqual(cursor.fetchone(), (1, 2))

            # Now the flatfile clears every gate, so the row becomes ready.
            record.update(
                conversion_status="PASS", webin_status="PASS", submission_ready=True
            )
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )
            with connection.cursor() as cursor:
                cursor.execute("SELECT submission_ready FROM ena_validation_attempts")
                self.assertEqual(cursor.fetchone()[0], True)

            # A ready row is still editable: validation is not submission.
            record.update(webin_status="NOT_RUN", submission_ready=False)
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )

            # submission_ready cannot outrun the flatfile format check.
            with connection.cursor() as cursor:
                with self.assertRaises(psycopg2.errors.CheckViolation):
                    cursor.execute(
                        "UPDATE ena_validation_attempts"
                        " SET submission_ready = TRUE, webin_status = 'FAIL_WEBIN'"
                    )
        finally:
            connection.rollback()
            connection.close()

    def test_submission_ledger_records_what_happened_after_validation(self):
        """015: the ledger the downstream submitter closes a sequence out in.

        Deliberately a separate table from ena_validation_attempts, which is
        upserted with DO UPDATE SET across every non-key column: an accession
        stored there would be erased by the next validation rerun.
        """
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        seqid = "OG82.ilmn.240313.getorg1770.emma102"
        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_submissions_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                for migration in (
                    "001_create_ena_validation_attempts.sql",
                    "002_ena_validation_attempts_single_row_per_attempt.sql",
                    "004_ena_candidate_packages.sql",
                    "005_ena_validation_attempts_genome_context.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    "010_drop_ena_locus_tables.sql",
                    "011_drop_local_package_validation.sql",
                    "012_drop_ena_selection_layer.sql",
                    "014_ena_validation_attempts_full_seqid.sql",
                    "015_ena_submissions.sql",
                ):
                    cursor.execute(_migration_sql(migration))

                cursor.execute(
                    """
                    INSERT INTO ena_validation_attempts (
                        full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
                        validation_mode, validation_attempt, table2asn_status,
                        conversion_status, preflight_status, webin_status, submission_ready
                    ) VALUES (
                        %s, 'OG82', 'ilmn', '240313', 'getorg1770', 'emma102', 'PRJEB123421',
                        'pipeline', 'initial', 'PASS', 'PASS', 'NOT_APPLICABLE', 'PASS', TRUE
                    )
                    """,
                    (seqid,),
                )

                # Validated but never submitted reads NOT_SUBMITTED, not NULL.
                cursor.execute(
                    "SELECT submission_status, ena_analysis_accession"
                    " FROM ena_submission_status WHERE full_seqid = %s",
                    (seqid,),
                )
                self.assertEqual(cursor.fetchone(), ("NOT_SUBMITTED", None))

                cursor.execute(
                    """
                    INSERT INTO ena_submissions (
                        full_seqid, webin_mode, submission_status, submitted_at, submitted_by,
                        receipt_path, receipt_sha256, webin_cli_version, ena_study_accession,
                        ena_analysis_accession, biosample_accession, biosample_source,
                        locus_tag_prefix, run_accessions
                    ) VALUES (
                        %s, 'production', 'ACCESSION_ASSIGNED', CURRENT_TIMESTAMP, 'tester',
                        '/receipts/OG82/receipt.xml', %s, '9.0.3', 'PRJEB123421',
                        'ERZ26543210', 'SAMN40589646', 'sample.ncbi_biosample_id',
                        'OGMITO01', ARRAY['ERR13579246', 'ERR13579247']
                    )
                    """,
                    (seqid, "a" * 64),
                )

                # Identity is GENERATED ALWAYS: the submitter writes full_seqid only.
                cursor.execute(
                    "SELECT og_id, tech, seq_date, code, annotation FROM ena_submissions"
                    " WHERE full_seqid = %s",
                    (seqid,),
                )
                self.assertEqual(
                    cursor.fetchone(), ("OG82", "ilmn", "240313", "getorg1770", "emma102")
                )

                # A dotted annotation version stays whole rather than truncating at the dot.
                cursor.execute(
                    "INSERT INTO ena_submissions (full_seqid) VALUES (%s)",
                    ("OG5.hifi.260101.final.emma1.0.2",),
                )
                cursor.execute(
                    "SELECT code, annotation FROM ena_submissions WHERE og_id = 'OG5'"
                )
                self.assertEqual(cursor.fetchone(), ("final", "emma1.0.2"))

                # A test-service dry run cannot overwrite the production record.
                cursor.execute(
                    """
                    INSERT INTO ena_submissions (
                        full_seqid, webin_mode, submission_status, submitted_at,
                        ena_analysis_accession
                    ) VALUES (%s, 'test', 'SUBMITTED', CURRENT_TIMESTAMP, 'ERZ99999999')
                    """,
                    (seqid,),
                )
                cursor.execute(
                    "SELECT webin_mode, ena_analysis_accession FROM ena_submissions"
                    " WHERE full_seqid = %s ORDER BY webin_mode",
                    (seqid,),
                )
                self.assertEqual(
                    cursor.fetchall(),
                    [("production", "ERZ26543210"), ("test", "ERZ99999999")],
                )

                # The later accessions land on a second pass, in place.
                cursor.execute(
                    """
                    INSERT INTO ena_submissions (
                        full_seqid, webin_mode, submission_status, submitted_at,
                        ena_analysis_accession, ena_sequence_accession
                    ) VALUES (
                        %s, 'production', 'ACCESSION_ASSIGNED', CURRENT_TIMESTAMP,
                        'ERZ26543210', 'OU123456.1'
                    )
                    ON CONFLICT (full_seqid, webin_mode) DO UPDATE SET
                        ena_sequence_accession = EXCLUDED.ena_sequence_accession,
                        updated_at = CURRENT_TIMESTAMP
                    RETURNING (xmax <> 0) AS updated
                    """,
                    (seqid,),
                )
                self.assertTrue(cursor.fetchone()[0], "upsert should update in place")

                cursor.execute(
                    "SELECT submission_status, ena_analysis_accession, ena_sequence_accession,"
                    " biosample_accession, locus_tag_prefix"
                    " FROM ena_submission_status WHERE full_seqid = %s",
                    (seqid,),
                )
                self.assertEqual(
                    cursor.fetchone(),
                    (
                        "ACCESSION_ASSIGNED",
                        "ERZ26543210",
                        "OU123456.1",
                        "SAMN40589646",
                        "OGMITO01",
                    ),
                )

                # Every status has to carry the evidence it claims, and every
                # accession has to look like the archive that mints it.
                rejected = (
                    ("ACCESSION_ASSIGNED without an accession",
                     "INSERT INTO ena_submissions (full_seqid, submission_status, submitted_at)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'ACCESSION_ASSIGNED', CURRENT_TIMESTAMP)"),
                    ("SUBMITTED without submitted_at",
                     "INSERT INTO ena_submissions (full_seqid, submission_status)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'SUBMITTED')"),
                    ("FAILED without a reason",
                     "INSERT INTO ena_submissions (full_seqid, submission_status)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'FAILED')"),
                    ("status outside the vocabulary",
                     "INSERT INTO ena_submissions (full_seqid, submission_status)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'DONE')"),
                    ("webin_mode outside the vocabulary",
                     "INSERT INTO ena_submissions (full_seqid, webin_mode)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'staging')"),
                    ("a full_seqid that is not an OG id",
                     "INSERT INTO ena_submissions (full_seqid) VALUES ('not-an-og-id')"),
                    ("a GenBank accession in the ERZ column",
                     "INSERT INTO ena_submissions (full_seqid, ena_analysis_accession)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'OQ123456')"),
                    ("a malformed BioSample",
                     "INSERT INTO ena_submissions (full_seqid, biosample_accession)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', 'SAMX123')"),
                    ("an analysis accession among the runs",
                     "INSERT INTO ena_submissions (full_seqid, run_accessions)"
                     " VALUES ('OG9.ilmn.260101.final.emma102', ARRAY['ERZ1234'])"),
                )
                for label, statement in rejected:
                    with self.subTest(rejects=label):
                        cursor.execute("SAVEPOINT rejected_row")
                        with self.assertRaises(psycopg2.errors.CheckViolation):
                            cursor.execute(statement)
                        cursor.execute("ROLLBACK TO SAVEPOINT rejected_row")
        finally:
            connection.rollback()
            connection.close()

    def test_specimen_registry_accepts_every_insdc_biosample_prefix(self):
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_biosample_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                for migration in (
                    "004_ena_candidate_packages.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    "007_insdc_biosample_accessions.sql",
                    "011_drop_local_package_validation.sql",
                ):
                    cursor.execute(_migration_sql(migration))

                # The registry records what we hold, including the SAMN that
                # webin cannot yet reference; ena_candidate_packages is where
                # that unusability is expressed, as BLOCKED_NCBI_ONLY_BIOSAMPLE.
                for numeric, accession in enumerate(
                    ("SAMEA132129018", "SAMN40589646", "SAMD00000001"), start=1
                ):
                    with self.subTest(accession=accession):
                        cursor.execute(
                            """
                            INSERT INTO ena_specimen_accessions (
                                og_id, og_numeric, ena_biosample_accession
                            ) VALUES (%s, %s, %s)
                            """,
                            (f"OG{numeric}", numeric, accession),
                        )

                # A savepoint, not a rollback: a plain rollback here would also
                # undo the throwaway schema and the migrations under test.
                for accession in ("SAMX123", "SRS123456", "SAMN"):
                    with self.subTest(rejected=accession):
                        cursor.execute("SAVEPOINT rejected_accession")
                        with self.assertRaises(psycopg2.errors.CheckViolation):
                            cursor.execute(
                                """
                                INSERT INTO ena_specimen_accessions (
                                    og_id, og_numeric, ena_biosample_accession
                                ) VALUES ('OG99', 99, %s)
                                """,
                                (accession,),
                            )
                        cursor.execute("ROLLBACK TO SAVEPOINT rejected_accession")

                # The new blocked status must be storable, or a package built
                # for an NCBI-only specimen cannot be recorded at all.
                cursor.execute(
                    """
                    INSERT INTO ena_candidate_packages (
                        full_seqid, og_id, assembly_prefix, annotation_version,
                        ena_study_accession, package_path, sequence_sha256,
                        normalised_circular_sha256, package_status
                    ) VALUES (
                        'OG2.hifi.260101.final.emma102', 'OG2',
                        'OG2.hifi.260101.final', 'emma102', 'PRJEB123419',
                        '/tmp/package', %s, %s,
                        'BLOCKED_NCBI_ONLY_BIOSAMPLE'
                    )
                    """,
                    ("a" * 64, "b" * 64),
                )
        finally:
            connection.rollback()
            connection.close()

    def test_migration_010_archives_then_drops_the_locus_tables(self):
        """010 retires the registry now that tags are assigned downstream.

        The drop is irreversible and the (og_id, canonical_gene,
        gene_occurrence) -> gene_serial assignment behind already-published tags
        cannot be rebuilt from the flat files, so the archive is the part that
        actually matters here.
        """
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_locus_drop_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                for migration in (
                    "004_ena_candidate_packages.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    "007_insdc_biosample_accessions.sql",
                    "011_drop_local_package_validation.sql",
                ):
                    cursor.execute(_migration_sql(migration))

                cursor.execute(
                    "INSERT INTO ena_specimen_accessions (og_id, og_numeric)"
                    " VALUES ('OG5', 5)"
                )
                seqid = "OG5.hifi.260101.final.emma102"
                cursor.execute(
                    """
                    INSERT INTO ena_candidate_packages (
                        full_seqid, og_id, assembly_prefix, annotation_version,
                        ena_study_accession, package_path, sequence_sha256,
                        normalised_circular_sha256, package_status,
                        platform, mean_depth
                    ) VALUES (
                        %s, 'OG5', 'OG5.hifi.260101.final', 'emma102',
                        'PRJEB123419', '/published/OG5/ena/package', %s, %s,
                        'READY', 'PACBIO_SMRT', 120.5
                    )
                    """,
                    (seqid, "5" * 64, "f" * 64),
                )
                for serial, gene in ((1, "TF"), (2, "RNR1"), (3, "ND6")):
                    cursor.execute(
                        """
                        INSERT INTO ena_locus_registry (
                            og_id, gene_serial, canonical_gene, gene_occurrence,
                            feature_type, strand, coordinate_snapshot
                        ) VALUES ('OG5', %s, %s, 1, 'gene', '+', '1..70')
                        """,
                        (serial, gene),
                    )
                    cursor.execute(
                        """
                        INSERT INTO ena_candidate_loci (
                            full_seqid, og_id, gene_serial, locus_tag, feature_key
                        ) VALUES (%s, 'OG5', %s, %s, %s)
                        """,
                        (seqid, serial, f"OGMTHIFI_000005{serial:03d}", f"gene:{serial}"),
                    )

                cursor.execute(_migration_sql("010_drop_ena_locus_tables.sql"))

                cursor.execute("SELECT to_regclass('ena_locus_registry')")
                self.assertIsNone(cursor.fetchone()[0])
                cursor.execute("SELECT to_regclass('ena_candidate_loci')")
                self.assertIsNone(cursor.fetchone()[0])

                cursor.execute(
                    "SELECT canonical_gene, gene_serial"
                    " FROM ena_locus_registry_archive ORDER BY gene_serial"
                )
                self.assertEqual(
                    cursor.fetchall(), [("TF", 1), ("RNR1", 2), ("ND6", 3)]
                )
                cursor.execute(
                    "SELECT locus_tag FROM ena_candidate_loci_archive"
                    " ORDER BY gene_serial"
                )
                self.assertEqual(
                    [row[0] for row in cursor.fetchall()],
                    [
                        "OGMTHIFI_000005001",
                        "OGMTHIFI_000005002",
                        "OGMTHIFI_000005003",
                    ],
                )

                # Idempotent, and a re-run must not blank the archive now that
                # the source tables are gone.
                cursor.execute(_migration_sql("010_drop_ena_locus_tables.sql"))
                cursor.execute("SELECT count(*) FROM ena_locus_registry_archive")
                self.assertEqual(cursor.fetchone()[0], 3)
        finally:
            connection.rollback()
            connection.close()


if __name__ == "__main__":
    unittest.main()
