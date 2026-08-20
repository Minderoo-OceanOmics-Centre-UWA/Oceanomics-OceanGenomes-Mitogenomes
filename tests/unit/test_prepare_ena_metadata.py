import importlib.util
import unittest
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "prepare_ena_metadata", ROOT / "bin" / "prepare_ena_metadata.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    """Dispatch on the queried table so query order is not baked into the test."""

    def __init__(self, biosample):
        self.biosample = biosample
        self.queries = []
        self.result = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, values):
        self.queries.append(" ".join(query.split()))
        if "FROM sample" in query:
            self.result = None if self.biosample is MISSING else [(self.biosample,)]
        elif "FROM mitogenome_data" in query:
            self.result = [(123.5,)]
        else:
            raise AssertionError(f"unexpected query: {query}")

    def fetchone(self):
        return self.result[0] if self.result else None

    def fetchall(self):
        return list(self.result or [])


class FakeConnection:
    def __init__(self, biosample):
        self.cursor_instance = FakeCursor(biosample)

    def cursor(self):
        return self.cursor_instance


MISSING = object()


def fetch(biosample):
    connection = FakeConnection(biosample)
    args = Namespace(
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
    def test_reads_biosample_from_sample_table(self):
        # sample.ncbi_biosample_id is the source of truth; the derived
        # ena_specimen_accessions cache must not be consulted.
        metadata, cursor = fetch("SAMN40589646")
        self.assertEqual(metadata["biosample_accession"], "SAMN40589646")
        self.assertEqual(metadata["biosample_source"], "sample.ncbi_biosample_id")
        self.assertTrue(any("FROM sample" in query for query in cursor.queries))
        self.assertFalse(
            any("ena_specimen_accessions" in query for query in cursor.queries)
        )

    def test_unregistered_specimen_normalises_to_none(self):
        # Unregistered specimens hold '' or whitespace rather than NULL. These
        # must read as absent so the package is WAITING_FOR_BIOSAMPLE rather
        # than BLOCKED_METADATA.
        for raw in ("", "   ", None):
            with self.subTest(raw=raw):
                metadata, _ = fetch(raw)
                self.assertIsNone(metadata["biosample_accession"])

    def test_missing_sample_row_is_not_an_error(self):
        metadata, _ = fetch(MISSING)
        self.assertIsNone(metadata["biosample_accession"])

    def test_surrounding_whitespace_is_stripped(self):
        metadata, _ = fetch("  SAMN40589646 ")
        self.assertEqual(metadata["biosample_accession"], "SAMN40589646")

    def test_other_metadata_still_collected(self):
        metadata, _ = fetch("SAMN40589646")
        self.assertEqual(metadata["mean_depth"], 123.5)
        self.assertEqual(metadata["run_accessions"], [])
        self.assertEqual(metadata["program"], "MitoHiFi 3")
        self.assertEqual(metadata["platform"], "PACBIO_SMRT")


    def test_run_accessions_are_not_queried(self):
        # Run accessions belong to the downstream submitter; ena_candidate_runs
        # was retired with the rest of the selection layer and must not be read.
        metadata, cursor = fetch("SAMN40589646")
        self.assertEqual(metadata["run_accessions"], [])
        self.assertFalse(
            any("ena_candidate_runs" in query for query in cursor.queries)
        )

    def test_program_mapping(self):
        self.assertEqual(MODULE.assembly_program("v323mitohifi"), "MitoHiFi 3.2.3")
        self.assertEqual(MODULE.assembly_program("getorg1771"), "GetOrganelle 1.7.7.1")
        self.assertEqual(MODULE.assembly_program("v10oatk"), "Oatk 1.0")

    def test_platform_mapping(self):
        self.assertEqual(MODULE.platform_for_tech("hifi"), "PACBIO_SMRT")
        self.assertEqual(MODULE.platform_for_tech("ilmn"), "ILLUMINA")
        self.assertEqual(MODULE.platform_for_tech("hic"), "ILLUMINA")


if __name__ == "__main__":
    unittest.main()
