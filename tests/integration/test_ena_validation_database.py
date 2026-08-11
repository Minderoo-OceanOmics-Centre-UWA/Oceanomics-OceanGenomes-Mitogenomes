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
    def test_migration_overwrites_until_selection_is_archived(self):
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
                ):
                    cursor.execute(_migration_sql(migration))

            record = {column: None for column in uploader.INSERT_COLUMNS}
            record.update(
                assembly_prefix="OG1.hifi.260101.final", og_id="OG1", ena_study="PRJEB1",
                validation_mode="pipeline", validation_attempt="integration",
                table2asn_status="PASS", conversion_status="FAIL",
                preflight_status="NOT_APPLICABLE", webin_status="NOT_RUN",
                package_status="READY", local_package_status="PASS",
                webin_test_status="NOT_RUN", webin_production_status="NOT_RUN",
                overall_status="LOCAL_PACKAGE_READY",
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
                conversion_status="PASS", webin_status="PASS",
                webin_production_status="PASS", overall_status="PRODUCTION_VALIDATED",
                submission_ready=True, result_digest="c" * 64
            )
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )

            # A later rerun remains editable because validation is not submission.
            record.update(
                webin_status="NOT_RUN", webin_production_status="NOT_RUN",
                overall_status="LOCAL_PACKAGE_READY", submission_ready=False,
                result_digest="d" * 64
            )
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "updated"
            )
            with connection.cursor() as cursor:
                cursor.execute(
                    """
                    INSERT INTO ena_specimen_accessions (og_id, og_numeric)
                    VALUES ('OG1', 1);
                    INSERT INTO ena_candidate_packages (
                        full_seqid, og_id, assembly_prefix, annotation_version,
                        ena_study_accession, package_path, sequence_sha256,
                        normalised_circular_sha256, package_status,
                        local_validation_status
                    ) VALUES (
                        'OG1.hifi.260101.final.emma102', 'OG1',
                        'OG1.hifi.260101.final', 'emma102', 'PRJEB1', '/tmp/package',
                        %s, %s, 'READY', 'PASS'
                    );
                    INSERT INTO ena_submission_selections (
                        ena_study_accession, og_id, selected_full_seqid,
                        selection_status, selection_reason, selected_by,
                        archive_status
                    ) VALUES (
                        'PRJEB1', 'OG1', 'OG1.hifi.260101.final.emma102',
                        'SELECTED', 'reviewed', 'integration', 'SUBMITTED'
                    )
                    """,
                    ("e" * 64, "f" * 64),
                )
            record["result_digest"] = "e" * 64
            self.assertEqual(
                uploader.upload_record(record, {}, connect=lambda **_kw: keep_open), "locked"
            )
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
                        normalised_circular_sha256, package_status,
                        local_validation_status
                    ) VALUES (
                        'OG2.hifi.260101.final.emma102', 'OG2',
                        'OG2.hifi.260101.final', 'emma102', 'PRJEB123419',
                        '/tmp/package', %s, %s,
                        'BLOCKED_NCBI_ONLY_BIOSAMPLE', 'FAIL'
                    )
                    """,
                    ("a" * 64, "b" * 64),
                )
        finally:
            connection.rollback()
            connection.close()

    def test_submission_queue_holds_only_work_the_submitter_should_do(self):
        """The handoff surface for the separate ENA submission pipeline.

        A row means: this package, at this path, is the chosen candidate for
        this specimen in this technology, and nobody has submitted it yet.
        """
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_queue_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                for migration in (
                    "004_ena_candidate_packages.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    "007_insdc_biosample_accessions.sql",
                ):
                    cursor.execute(_migration_sql(migration))
                # The only two columns the view needs from sample. Stubbed so
                # the test never reads the real specimen table.
                cursor.execute(
                    "CREATE TABLE sample (og_id TEXT PRIMARY KEY, embargo_status TEXT)"
                )
                cursor.execute(_migration_sql("008_ena_submission_queue.sql"))

                for numeric, embargo in ((1, "Release"), (2, "Embargoed"), (3, "Release")):
                    cursor.execute(
                        "INSERT INTO sample (og_id, embargo_status) VALUES (%s, %s)",
                        (f"OG{numeric}", embargo),
                    )
                    cursor.execute(
                        """
                        INSERT INTO ena_specimen_accessions (og_id, og_numeric)
                        VALUES (%s, %s)
                        """,
                        (f"OG{numeric}", numeric),
                    )

                def add_package(numeric, package_status="READY", local="PASS"):
                    seqid = f"OG{numeric}.hifi.260101.final.emma102"
                    cursor.execute(
                        """
                        INSERT INTO ena_candidate_packages (
                            full_seqid, og_id, assembly_prefix, annotation_version,
                            ena_study_accession, package_path, sequence_sha256,
                            normalised_circular_sha256, package_status,
                            local_validation_status, platform, mean_depth
                        ) VALUES (
                            %s, %s, %s, 'emma102', 'PRJEB123419',
                            %s, %s, %s, %s, %s, 'PACBIO_SMRT', 120.5
                        )
                        """,
                        (
                            seqid,
                            f"OG{numeric}",
                            f"OG{numeric}.hifi.260101.final",
                            f"/published/OG{numeric}/ena/package",
                            # Both digest columns are CHECKed as lowercase hex.
                            str(numeric) * 64,
                            chr(ord("a") + numeric) * 64,
                            package_status,
                            local,
                        ),
                    )
                    return seqid

                def add_selection(numeric, seqid, selection_status, archive_status):
                    cursor.execute(
                        """
                        INSERT INTO ena_submission_selections (
                            ena_study_accession, og_id, selected_full_seqid,
                            selection_status, selection_reason, selected_by,
                            archive_status
                        ) VALUES (
                            'PRJEB123419', %s, %s, %s, 'test', 'integration', %s
                        )
                        """,
                        (f"OG{numeric}", seqid, selection_status, archive_status),
                    )

                add_selection(1, add_package(1), "SELECTED", "NOT_SUBMITTED")
                # Ambiguous candidates are the submitter's business only after a
                # human resolves them, so they must not surface here.
                add_selection(2, add_package(2), "MANUAL_REVIEW_REQUIRED", "NOT_SUBMITTED")
                add_selection(
                    3, add_package(3, "BLOCKED_METADATA", "FAIL"), "SELECTED", "NOT_SUBMITTED"
                )

                cursor.execute(
                    """
                    SELECT og_id, tech, full_seqid, package_path, platform,
                           mean_depth, embargo_status, archive_status
                    FROM ena_submission_queue
                    """
                )
                rows = cursor.fetchall()
                self.assertEqual(
                    rows,
                    [
                        (
                            "OG1", "hifi", "OG1.hifi.260101.final.emma102",
                            "/published/OG1/ena/package", "PACBIO_SMRT",
                            120.5, "Release", "NOT_SUBMITTED",
                        )
                    ],
                )

                # Embargo is exposed, never applied: an embargoed but otherwise
                # clean specimen still appears, and the submitter filters it.
                cursor.execute(
                    "UPDATE sample SET embargo_status = 'Embargoed' WHERE og_id = 'OG1'"
                )
                cursor.execute("SELECT embargo_status FROM ena_submission_queue")
                self.assertEqual(cursor.fetchall(), [("Embargoed",)])

                # Closing the row out is what removes the work item.
                cursor.execute(
                    """
                    UPDATE ena_submission_selections
                    SET archive_status = 'SUBMITTED'
                    WHERE og_id = 'OG1'
                    """
                )
                cursor.execute("SELECT count(*) FROM ena_submission_queue")
                self.assertEqual(cursor.fetchone()[0], 0)
        finally:
            connection.rollback()
            connection.close()

    def test_migration_renumbers_serials_into_coordinate_order(self):
        """009 rewrites serials allocated in Emma's string-sorted file order.

        The registry is keyed on (canonical_gene, gene_occurrence) and the
        allocator skips a gene that already holds a serial, so fixing the sort
        upstream only helps new specimens. This is what repairs the ones that
        were already allocated.
        """
        try:
            import psycopg2
        except ImportError:
            self.skipTest("psycopg2 is not installed")

        connection = psycopg2.connect(os.environ["ENA_TEST_DB_DSN"])
        schema = f"ena_locus_order_test_{uuid.uuid4().hex}"
        try:
            with connection.cursor() as cursor:
                cursor.execute(f'CREATE SCHEMA "{schema}"')
                cursor.execute(f'SET search_path TO "{schema}"')
                for migration in (
                    "004_ena_candidate_packages.sql",
                    "006_ena_tech_aware_locus_tags.sql",
                    "007_insdc_biosample_accessions.sql",
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
                        local_validation_status, platform, mean_depth
                    ) VALUES (
                        %s, 'OG5', 'OG5.hifi.260101.final', 'emma102',
                        'PRJEB123419', '/published/OG5/ena/package', %s, %s,
                        'READY', 'PASS', 'PACBIO_SMRT', 120.5
                    )
                    """,
                    (seqid, "5" * 64, "f" * 64),
                )

                # Serials as Emma's string sort produced them: TF at 1 first,
                # then ND4L at 10053, TV at 1027, and RNR1 at 71 last.
                allocated = [
                    (1, "TF", "1..70"),
                    (2, "ND4L", "10053..10349"),
                    (3, "TV", "1027..1098"),
                    (4, "ND6", "14279..13755"),
                    (5, "RNR1", "71..1026"),
                ]
                for serial, gene, snapshot in allocated:
                    cursor.execute(
                        """
                        INSERT INTO ena_locus_registry (
                            og_id, gene_serial, canonical_gene, gene_occurrence,
                            feature_type, strand, coordinate_snapshot
                        ) VALUES ('OG5', %s, %s, 1, 'gene', '+', %s)
                        """,
                        (serial, gene, snapshot),
                    )
                    cursor.execute(
                        """
                        INSERT INTO ena_candidate_loci (
                            full_seqid, og_id, gene_serial, locus_tag, feature_key
                        ) VALUES (%s, 'OG5', %s, %s, %s)
                        """,
                        (seqid, serial, f"OGMTHIFI_000005{serial:03d}", f"gene:{serial}"),
                    )

                cursor.execute(_migration_sql("009_ena_locus_registry_canonical_order.sql"))

                cursor.execute(
                    "SELECT canonical_gene, gene_serial FROM ena_locus_registry"
                    " WHERE og_id = 'OG5' ORDER BY gene_serial"
                )
                self.assertEqual(
                    cursor.fetchall(),
                    # ND6 is stored 14279..13755, so it sorts on 13755.
                    [("TF", 1), ("RNR1", 2), ("TV", 3), ("ND4L", 4), ("ND6", 5)],
                )

                # The rendered tag carries the serial as its last three digits,
                # so it has to travel with it.
                cursor.execute(
                    "SELECT canonical_gene, loci.locus_tag"
                    " FROM ena_candidate_loci loci"
                    " JOIN ena_locus_registry USING (og_id, gene_serial)"
                    " ORDER BY gene_serial"
                )
                self.assertEqual(
                    cursor.fetchall(),
                    [
                        ("TF", "OGMTHIFI_000005001"),
                        ("RNR1", "OGMTHIFI_000005002"),
                        ("TV", "OGMTHIFI_000005003"),
                        ("ND4L", "OGMTHIFI_000005004"),
                        ("ND6", "OGMTHIFI_000005005"),
                    ],
                )

                # Idempotent: the second run is a no-op, not another rotation.
                cursor.execute(_migration_sql("009_ena_locus_registry_canonical_order.sql"))
                cursor.execute(
                    "SELECT canonical_gene, gene_serial FROM ena_locus_registry"
                    " WHERE og_id = 'OG5' ORDER BY gene_serial"
                )
                self.assertEqual(
                    cursor.fetchall(),
                    [("TF", 1), ("RNR1", 2), ("TV", 3), ("ND4L", 4), ("ND6", 5)],
                )

                # A published tag is frozen, so the migration must refuse.
                cursor.execute(
                    """
                    INSERT INTO ena_submission_selections (
                        ena_study_accession, og_id, selected_full_seqid,
                        selection_status, selection_reason, selected_by,
                        archive_status
                    ) VALUES (
                        'PRJEB123419', 'OG5', %s, 'SELECTED', 'test',
                        'integration', 'SUBMITTED'
                    )
                    """,
                    (seqid,),
                )
                with self.assertRaises(psycopg2.errors.RaiseException):
                    cursor.execute(
                        _migration_sql("009_ena_locus_registry_canonical_order.sql")
                    )
        finally:
            connection.rollback()
            connection.close()


if __name__ == "__main__":
    unittest.main()
