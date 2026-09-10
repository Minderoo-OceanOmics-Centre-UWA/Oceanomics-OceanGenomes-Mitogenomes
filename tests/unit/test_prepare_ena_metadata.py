import importlib.util
import io
import os
import sys
import tempfile
import unittest
from argparse import Namespace
from contextlib import redirect_stderr
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "prepare_ena_metadata", ROOT / "bin" / "prepare_ena_metadata.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    """Dispatch on the queried table so query order is not baked into the test."""

    def __init__(self):
        self.queries = []
        self.result = None
        self.empty = False

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, values):
        self.queries.append(" ".join(query.split()))
        if "FROM mitogenome_data" in query:
            self.result = [] if self.empty else [(123.5,)]
        else:
            raise AssertionError(f"unexpected query: {query}")

    def fetchone(self):
        return self.result[0] if self.result else None

    def fetchall(self):
        return list(self.result or [])


class FakeConnection:
    def __init__(self):
        self.cursor_instance = FakeCursor()

    def cursor(self):
        return self.cursor_instance


def fetch(depth_tsv=None):
    connection = FakeConnection()
    args = Namespace(
        depth_tsv=depth_tsv,
        og_id="OG910",
        assembly_prefix="OG910.hifi.241127.v3mitohifi",
        annotation_version="emma102",
        full_seqid="OG910.hifi.241127.v3mitohifi.emma102",
        tech="hifi",
        seq_date="241127",
        code="v3mitohifi",
        study="PRJEB123419",
        scientific_name="Choerodon rubescens",
    )
    return MODULE.fetch_metadata(connection, args), connection.cursor_instance


class PrepareEnaMetadataTests(unittest.TestCase):
    def test_the_submitter_owned_accessions_are_not_recorded(self):
        """The BioSample, the study target and the runs are registered downstream.

        Carrying a BioSample here produced a SAMPLE key webin could not resolve
        and a package that read as blocked on something this pipeline could fix.
        """
        metadata, cursor = fetch()
        for absent in (
            "biosample_accession",
            "biosample_source",
            "run_accessions",
            "study",
        ):
            self.assertNotIn(absent, metadata)
        for table in ("FROM sample", "ena_specimen_accessions", "ena_candidate_runs"):
            self.assertFalse(any(table in query for query in cursor.queries))

    def test_the_study_is_recorded_as_the_validation_study(self):
        """--ena_study names the study validation ran against, not a target."""
        metadata, _ = fetch()
        self.assertEqual(metadata["validation_study"], "PRJEB123419")
        self.assertEqual(metadata["schema_version"], 3)

    def test_other_metadata_still_collected(self):
        metadata, _ = fetch()
        self.assertEqual(metadata["mean_depth"], 123.5)
        self.assertEqual(metadata["mean_depth_source"], "database")
        self.assertEqual(metadata["program"], "MitoHiFi 3")
        self.assertEqual(metadata["platform"], "PACBIO_SMRT")
        self.assertEqual(metadata["scientific_name"], "Choerodon rubescens")

    def test_program_mapping(self):
        self.assertEqual(MODULE.assembly_program("v323mitohifi"), "MitoHiFi 3.2.3")
        self.assertEqual(MODULE.assembly_program("getorg1771"), "GetOrganelle 1.7.7.1")
        self.assertEqual(MODULE.assembly_program("v10oatk"), "Oatk 1.0")

    def test_platform_mapping(self):
        self.assertEqual(MODULE.platform_for_tech("hifi"), "PACBIO_SMRT")
        self.assertEqual(MODULE.platform_for_tech("ilmn"), "ILLUMINA")
        self.assertEqual(MODULE.platform_for_tech("hic"), "ILLUMINA")


class MeanDepthSourceTests(unittest.TestCase):
    """Coverage comes from this run's measurement, and the database is the fallback.

    The precedence is the whole reason submission prep can run with
    --skip_upload_results: reading the file instead of the row this run writes is
    what removed the ordering dependency on the upload.
    """

    def _depth_file(self, text):
        handle = tempfile.NamedTemporaryFile(
            "w", suffix=".mito_depth.tsv", delete=False
        )
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_the_pipeline_measurement_wins_over_the_stored_row(self):
        """A re-assembled molecule must not report the previous run's depth."""
        path = self._depth_file("sample\tmean_depth\nOG910\t812.5\n")
        metadata, cursor = fetch(depth_tsv=path)
        self.assertEqual(metadata["mean_depth"], 812.5)
        self.assertEqual(metadata["mean_depth_source"], "pipeline")
        # The stored 123.5 was never even queried.
        self.assertEqual(cursor.queries, [])

    def test_the_database_is_used_when_there_is_no_measurement(self):
        """--skip_mitogenome_depth and precomputed assemblies land here."""
        path = self._depth_file("sample\tmean_depth\n")
        metadata, cursor = fetch(depth_tsv=path)
        self.assertEqual(metadata["mean_depth"], 123.5)
        self.assertEqual(metadata["mean_depth_source"], "database")
        self.assertTrue(any("FROM mitogenome_data" in q for q in cursor.queries))

    def test_neither_source_is_reported_as_none_not_zero(self):
        """A false COVERAGE 0 is worse than no COVERAGE at all."""
        connection = FakeConnection()
        connection.cursor_instance.empty = True
        args = Namespace(
            depth_tsv=None,
            og_id="OG910",
            assembly_prefix="OG910.hifi.241127.v3mitohifi",
            annotation_version="emma102",
            full_seqid="OG910.hifi.241127.v3mitohifi.emma102",
            tech="hifi",
            seq_date="241127",
            code="v3mitohifi",
            study="PRJEB123419",
            scientific_name="Choerodon rubescens",
        )
        with redirect_stderr(io.StringIO()) as err:
            metadata = MODULE.fetch_metadata(connection, args)
        self.assertIsNone(metadata["mean_depth"])
        self.assertEqual(metadata["mean_depth_source"], "none")
        self.assertIn("no mean_depth", err.getvalue())


if __name__ == "__main__":
    unittest.main()
