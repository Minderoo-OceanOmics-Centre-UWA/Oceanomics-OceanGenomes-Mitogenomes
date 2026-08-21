"""The species table is not the only source of a sample's lineage.

A miss there used to leave class='unknown', which every downstream chooser reads
as "vertebrate": genetic code 2, EMMA instead of MITOS2, the curated fish BLAST
DB, the vertebrate 22-tRNA completeness expectation. These tests pin the NCBI
taxdump fallback that fills the gap, and the rule that the curated database still
wins whenever it has an answer.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"

sys.path.insert(0, str(BIN))
from taxdump_lineage import TaxdumpLineage  # noqa: E402


def load_create_samplesheet():
    """Import create_samplesheet.py without requiring psycopg2 to be installed."""
    spec = importlib.util.spec_from_file_location(
        "create_samplesheet_taxdump_under_test", BIN / "create_samplesheet.py"
    )
    module = importlib.util.module_from_spec(spec)
    patches = {}
    try:
        import psycopg2  # noqa: F401
    except ImportError:
        patches["psycopg2"] = types.ModuleType("psycopg2")
    sys.path.insert(0, str(BIN))
    try:
        with mock.patch.dict(sys.modules, patches):
            spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(BIN))
    return module


cs = load_create_samplesheet()


# A miniature taxdump: Anthozoa > Scleractinia > Acroporidae > Acropora >
# Acropora tenuis, plus a cross-kingdom homonym at genus rank ('Morus').
NODES = [
    (1, 1, "no rank"),
    (2, 1, "class"),
    (3, 2, "order"),
    (4, 3, "family"),
    (5, 4, "genus"),
    (6, 5, "species"),
    (7, 5, "species"),
    (10, 1, "class"),
    (11, 10, "genus"),
    (20, 1, "class"),
    (21, 20, "genus"),
]

NAMES = [
    (1, "root", "scientific name"),
    (2, "Anthozoa", "scientific name"),
    (3, "Scleractinia", "scientific name"),
    (4, "Acroporidae", "scientific name"),
    (5, "Acropora", "scientific name"),
    (6, "Acropora tenuis", "scientific name"),
    (6, "Madrepora tenuis", "synonym"),
    (7, "Acropora sp.", "scientific name"),
    (10, "Aves", "scientific name"),
    (11, "Morus", "scientific name"),
    (20, "Magnoliopsida", "scientific name"),
    (21, "Morus", "scientific name"),
]


def write_taxdump(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "nodes.dmp", "w") as handle:
        for taxid, parent, rank in NODES:
            handle.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\n")
    with open(directory / "names.dmp", "w") as handle:
        for taxid, name, name_class in NAMES:
            handle.write(f"{taxid}\t|\t{name}\t|\t\t|\t{name_class}\t|\n")
    return directory


class TaxdumpLineageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import tempfile

        cls._tmp = tempfile.TemporaryDirectory()
        cls.taxdump = write_taxdump(Path(cls._tmp.name) / "taxonkit_dbs")
        cls.resolver = TaxdumpLineage(str(cls.taxdump))

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_binomial_resolves_the_full_lineage(self):
        lineage = self.resolver.lineage_for_name("Acropora tenuis")
        self.assertEqual(lineage["class"], "Anthozoa")
        self.assertEqual(lineage["order"], "Scleractinia")
        self.assertEqual(lineage["family"], "Acroporidae")
        self.assertEqual(lineage["matched_rank"], "species")

    def test_unknown_species_falls_back_to_its_genus(self):
        # The species table misses far more often than NCBI misses a genus.
        lineage = self.resolver.lineage_for_name("Acropora notarealspecies")
        self.assertEqual(lineage["class"], "Anthozoa")
        self.assertEqual(lineage["family"], "Acroporidae")
        self.assertEqual(lineage["matched_rank"], "genus")

    def test_open_nomenclature_resolves_at_genus_not_the_placeholder_node(self):
        # NCBI's 'Acropora sp.' node stands for one submitter's unidentified
        # organism, so matching it would claim our sample is that record.
        lineage = self.resolver.lineage_for_name("Acropora sp.")
        self.assertEqual(lineage["matched_name"], "Acropora")
        self.assertEqual(lineage["matched_rank"], "genus")

    def test_family_name_alone_still_pins_a_class(self):
        lineage = self.resolver.lineage_for_name("Acroporidae")
        self.assertEqual(lineage["class"], "Anthozoa")
        self.assertEqual(lineage["family"], "Acroporidae")

    def test_cross_kingdom_homonym_resolves_to_nothing(self):
        # 'Morus' is both a gannet and a mulberry. No lineage beats the wrong one.
        self.assertEqual(self.resolver.lineage_for_name("Morus"), {})

    def test_absent_name_resolves_to_nothing(self):
        self.assertEqual(self.resolver.lineage_for_name("Nothing here"), {})

    def test_synonyms_are_not_indexed(self):
        # Only scientific names are indexed, so a synonym is a miss rather than a
        # silently different taxon.
        self.assertEqual(self.resolver.lineage_for_name("Madrepora tenuis"), {})

    def test_available_is_false_without_the_dmp_files(self):
        self.assertFalse(TaxdumpLineage("/nonexistent/taxdump").available)


class ResolveSpeciesInfoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import tempfile

        cls._tmp = tempfile.TemporaryDirectory()
        cls.taxdump = write_taxdump(Path(cls._tmp.name) / "taxonkit_dbs")
        cls.resolver = TaxdumpLineage(str(cls.taxdump))

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def resolve(self, db_result, resolver=None):
        with mock.patch.object(cs, "query_species_info", return_value=db_result):
            return cs.resolve_species_info(None, "OG1", resolver)

    def test_curated_database_answer_is_left_alone(self):
        # The species table is authoritative; the taxdump is only a backstop.
        result = self.resolve(
            ("Acropora tenuis", "Anthozoa", "Acroporidae", "Scleractinia", "Acropora tenuis"),
            self.resolver,
        )
        self.assertEqual(result[1:5],
                         ("Anthozoa", "Acroporidae", "Scleractinia", "Acropora tenuis"))
        self.assertEqual(result[5], "db")

    def test_unknown_class_is_filled_from_the_taxdump(self):
        result = self.resolve(
            ("Acropora tenuis", "unknown", "", "", "Acropora tenuis"), self.resolver
        )
        self.assertEqual(result[1], "Anthozoa")
        self.assertEqual(result[2], "Acroporidae")
        self.assertEqual(result[3], "Scleractinia")
        self.assertEqual(result[5], "taxdump")

    def test_blank_family_and_order_are_filled_without_touching_the_class(self):
        # family/order drive the CROSS_ORDER tier in reference_divergence_check;
        # blank ones collapse every reference to NON_CONGENERIC.
        result = self.resolve(
            ("Acropora tenuis", "Anthozoa", "", "", "Acropora tenuis"), self.resolver
        )
        self.assertEqual(result[1], "Anthozoa")
        self.assertEqual(result[2], "Acroporidae")
        self.assertEqual(result[3], "Scleractinia")
        self.assertEqual(result[5], "db+taxdump")

    def test_blank_reference_species_id_gains_the_ncbi_name(self):
        result = self.resolve(
            ("Acropora tenuis", "unknown", "", "", ""), self.resolver
        )
        self.assertEqual(result[4], "Acropora tenuis")

    def test_unresolvable_sample_stays_unknown_and_is_flagged(self):
        result = self.resolve(("Nothing here", "unknown", "", "", ""), self.resolver)
        self.assertEqual(result[1], "unknown")
        self.assertEqual(result[5], "unresolved")

    def test_no_taxdump_still_returns_the_database_answer(self):
        result = self.resolve(
            ("Acropora tenuis", "unknown", "", "", "Acropora tenuis"), None
        )
        self.assertEqual(result[1], "unknown")
        self.assertEqual(result[5], "unresolved")

    def test_unknown_class_is_never_treated_as_an_invertebrate(self):
        # The bug this whole guard exists for: 'unknown' is not 'vertebrate', but
        # is_invertebrate() cannot say so, which is why the run must abort.
        self.assertEqual(cs.is_invertebrate("unknown"), "false")
        self.assertEqual(cs.is_invertebrate("Anthozoa"), "true")


if __name__ == "__main__":
    unittest.main()
