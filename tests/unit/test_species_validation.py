import configparser
import csv
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "species_validation", ROOT / "bin" / "species_validation.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    """Routes fetchone() by the most recent query so one fake DB can back both
    get_species_for_ogid (SELECT ... FROM sample) and upsert_lca_validation
    (SELECT ... FROM lca_validation / INSERT INTO lca_validation)."""

    def __init__(self, sample_species, existing_lca_row, log):
        self.sample_species = sample_species
        self.existing_lca_row = existing_lca_row
        self.log = log
        self._last_query = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, params=None):
        self._last_query = query
        self.log.append((query, params))

    def fetchone(self):
        if "FROM sample" in self._last_query:
            return (self.sample_species,) if self.sample_species is not None else None
        if "FROM lca_validation" in self._last_query:
            return self.existing_lca_row
        return None


class FakeConnection:
    def __init__(self, sample_species, existing_lca_row, log):
        self.cursor_instance = FakeCursor(sample_species, existing_lca_row, log)

    def cursor(self):
        return self.cursor_instance

    def commit(self):
        pass

    def rollback(self):
        pass

    def close(self):
        pass


def make_config(root):
    path = root / "oceanomics.cfg"
    config = configparser.ConfigParser()
    config["postgres"] = {
        "dbname": "test", "user": "u", "password": "p", "host": "h", "port": "5432",
    }
    with path.open("w") as fh:
        config.write(fh)
    return path


class NoNominalSpeciesTests(unittest.TestCase):
    """No nominal_species_id on the sample row: still write the summary file
    (with N/A comparison columns) and still upsert an lca_validation row with
    validated_species_name=None, so PUSH_LCA_BLAST_RESULTS still runs and the
    sample isn't silently absent from lca_validation."""

    def _run(self, force):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = "OG470.ilmn.230607.getorg1770"

            lca_path = root / "lca.12s.tsv"
            with lca_path.open("w", newline="") as fh:
                fh.write("seq_id\tspecies_in_LCA\n")
                fh.write(f"{prefix}.emma102.001\tGenus species\n")

            blast_path = root / "blast.12s.tsv"
            with blast_path.open("w") as fh:
                fh.write(f"{prefix}.emma102\tsubject1\tsome\tblast\tfields\n")

            config_path = make_config(root)

            log = []

            def fake_connect(**_kwargs):
                return FakeConnection(sample_species=None, existing_lca_row=None, log=log)

            class FakePsycopg2:
                connect = staticmethod(fake_connect)

            original_psycopg2 = MODULE.psycopg2
            original_cwd = Path.cwd()
            MODULE.psycopg2 = FakePsycopg2
            os.chdir(root)
            try:
                buffer = io.StringIO()
                with redirect_stdout(buffer):
                    MODULE.compare_lca_and_blast(
                        str(config_path),
                        "OG470",
                        [str(lca_path)],
                        [str(blast_path)],
                        f"lca_results.{prefix}.tsv",
                        assembly_prefix=prefix,
                        force=force,
                    )
            finally:
                MODULE.psycopg2 = original_psycopg2
                os.chdir(original_cwd)

            output = buffer.getvalue()
            with (root / f"lca_results.{prefix}.tsv").open() as fh:
                summary_rows = list(csv.reader(fh, delimiter="\t"))
            return output, summary_rows, log

    def test_writes_na_summary_columns(self):
        output, rows, _log = self._run(force=False)
        self.assertIn("has no nominal_species_id", output)
        self.assertEqual(
            rows[0],
            ["og_id", "LCA_result", "nom_species_id", "Match_YN", "Found_in_blast_YN"],
        )
        self.assertEqual(rows[1][2:], ["N/A", "N/A", "N/A"])

    def test_upserts_lca_validation_row_with_null_species(self):
        _output, _rows, log = self._run(force=False)
        insert_calls = [
            (query, params)
            for query, params in log
            if query and "INSERT INTO lca_validation" in query
        ]
        self.assertEqual(len(insert_calls), 1)
        _query, params = insert_calls[0]
        self.assertIsNone(params["validated_species_name"])
        self.assertEqual(params["og_id"], "OG470")
        self.assertEqual(params["tech"], "ilmn")
        self.assertEqual(params["annotation"], "emma102")
        self.assertEqual(params["validator"], "nf-core")


if __name__ == "__main__":
    unittest.main()
