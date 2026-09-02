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


# A miniature taxdump: Metazoa > Cnidaria > Anthozoa > Scleractinia >
# Acroporidae > Acropora > Acropora tenuis, a vertebrate branch under
# Vertebrata, a plant branch outside Metazoa, and two homonyms -- 'Morus'
# (gannet vs mulberry, one animal candidate) and 'Vertebrata' itself (the clade
# vs a plant genus), which is why the ancestry anchors are picked by lineage
# rather than by taking the first match.
NODES = [
    (1, 1, "no rank"),
    (100, 1, "kingdom"),      # Metazoa
    (101, 100, "clade"),      # Vertebrata
    (30, 100, "phylum"),      # Cnidaria
    (2, 30, "class"),         # Anthozoa
    (3, 2, "order"),
    (4, 3, "family"),
    (5, 4, "genus"),
    (6, 5, "species"),
    (7, 5, "species"),
    (40, 100, "phylum"),      # Porifera -- a phylum with no class below it here
    (10, 101, "class"),       # Aves
    (11, 10, "genus"),        # Morus, the gannet
    (20, 1, "class"),         # Magnoliopsida, outside Metazoa
    (21, 20, "genus"),        # Morus, the mulberry
    (22, 20, "genus"),        # Vertebrata, the red alga
]

NAMES = [
    (1, "root", "scientific name"),
    (100, "Metazoa", "scientific name"),
    (101, "Vertebrata", "scientific name"),
    (30, "Cnidaria", "scientific name"),
    (2, "Anthozoa", "scientific name"),
    (3, "Scleractinia", "scientific name"),
    (4, "Acroporidae", "scientific name"),
    (5, "Acropora", "scientific name"),
    (6, "Acropora tenuis", "scientific name"),
    (6, "Madrepora tenuis", "synonym"),
    (7, "Acropora sp.", "scientific name"),
    (40, "Porifera", "scientific name"),
    (10, "Aves", "scientific name"),
    (11, "Morus", "scientific name"),
    (20, "Magnoliopsida", "scientific name"),
    (21, "Morus", "scientific name"),
    (22, "Vertebrata", "scientific name"),
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

    def test_cross_kingdom_homonym_resolves_to_the_animal(self):
        # 'Morus' is both a gannet and a mulberry. This pipeline sequences animals
        # and never plants, so the animal candidate is the answer -- dropping the
        # name outright cost real samples their lineage (the sponge genus
        # Acanthella and the barnacle genus Calantica are both plant homonyms).
        lineage = self.resolver.lineage_for_name("Morus")
        self.assertEqual(lineage["class"], "Aves")
        self.assertTrue(lineage["is_animal"])
        self.assertTrue(lineage["is_vertebrate"])

    def test_the_vertebrata_anchor_is_the_animal_one(self):
        # 'Vertebrata' is itself a homonym (a red algal genus), so an anchor taken
        # by first match would classify every vertebrate as an invertebrate.
        self.assertTrue(self.resolver.lineage_for_name("Acropora tenuis")["is_animal"])
        self.assertFalse(self.resolver.lineage_for_name("Acropora tenuis")["is_vertebrate"])
        self.assertTrue(self.resolver.lineage_for_name("Morus")["is_vertebrate"])

    def test_a_phylum_only_name_still_resolves(self):
        # 'Porifera' pins no class, but the phylum is enough to choose the genetic
        # code and the annotation route, and is what the operator has when the
        # sponge has not been identified further.
        lineage = self.resolver.lineage_for_name("Porifera")
        self.assertEqual(lineage["phylum"], "Porifera")
        self.assertNotIn("class", lineage)
        self.assertTrue(lineage["is_animal"])
        self.assertFalse(lineage["is_vertebrate"])

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

    def test_ancestry_beats_the_class_allow_list(self):
        # A class nobody has added to INVERT_CLASSES is still an invertebrate when
        # the taxdump says it is an animal outside Vertebrata. Without this, a
        # missing class silently takes the vertebrate path -- how barnacles
        # (Thecostraca) were being annotated as fish.
        self.assertNotIn("Nothingoidea", cs.INVERT_CLASSES)
        self.assertEqual(
            cs.is_invertebrate("Nothingoidea",
                               {"is_animal": True, "is_vertebrate": False}),
            "true")
        self.assertEqual(
            cs.is_invertebrate("Anthozoa",
                               {"is_animal": True, "is_vertebrate": True}),
            "false")

    def test_the_class_list_is_still_the_fallback_without_a_lineage(self):
        # A class straight from the species table has no taxdump lineage behind it.
        self.assertEqual(cs.is_invertebrate("Anthozoa", {}), "true")
        self.assertEqual(cs.is_invertebrate("Actinopteri", None), "false")

    def test_a_phylum_only_sample_gets_the_phylum_as_its_class(self):
        result = self.resolve(("Porifera", "unknown", "", "", ""), self.resolver)
        self.assertEqual(result[1], "Porifera")
        self.assertEqual(result[5], "taxdump-phylum")
        self.assertEqual(cs.is_invertebrate(result[1], result[6]), "true")


if __name__ == "__main__":
    unittest.main()
