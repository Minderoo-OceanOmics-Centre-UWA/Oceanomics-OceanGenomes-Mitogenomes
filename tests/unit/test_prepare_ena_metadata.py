import importlib.util
import sys
import unittest
from argparse import Namespace
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

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, values):
        self.queries.append(" ".join(query.split()))
        if "FROM mitogenome_data" in query:
            self.result = [(123.5,)]
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


def fetch():
    connection = FakeConnection()
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
        self.assertEqual(metadata["schema_version"], 2)

    def test_other_metadata_still_collected(self):
        metadata, _ = fetch()
        self.assertEqual(metadata["mean_depth"], 123.5)
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


if __name__ == "__main__":
    unittest.main()
