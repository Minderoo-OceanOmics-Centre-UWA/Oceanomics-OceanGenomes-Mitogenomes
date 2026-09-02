"""assets/mito_genetic_codes.json is the one class -> genetic code map.

It is read by both bin/create_samplesheet.py (which fills the samplesheet's
genetic_code column) and lib/InvertTaxonGroups.groovy (which prepare_samplesheet
uses to resolve meta.genetic_code). These tests pin the asset's internal
consistency and, most importantly, that it covers every class
create_samplesheet.py is willing to call an invertebrate -- a class that is in
INVERT_CLASSES but in neither the code map nor the documented `ambiguous` block
aborts every run that contains one, with no way for the operator to see it
coming from the sheet.
"""

import importlib.util
import json
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
ASSET = ROOT / "assets" / "mito_genetic_codes.json"

# The genetic codes the pipeline will actually annotate under; mirrors
# SUPPORTED_GENETIC_CODES in subworkflows/local/mitogenome_annotation_lca.
SUPPORTED_CODES = {2, 4, 5, 9, 13, 14, 21, 24, 33}


def load_create_samplesheet():
    """Import create_samplesheet.py without requiring psycopg2 to be installed."""
    spec = importlib.util.spec_from_file_location(
        "create_samplesheet_codes_under_test", BIN / "create_samplesheet.py"
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


class AssetShapeTests(unittest.TestCase):
    def setUp(self):
        self.codes, self.ambiguous = cs.load_genetic_codes(ASSET)

    def test_asset_parses_and_is_not_empty(self):
        self.assertTrue(self.codes)
        self.assertTrue(self.ambiguous)

    def test_every_code_is_one_the_pipeline_supports(self):
        # A code outside this set passes the samplesheet but is rejected by the
        # assertion in mitogenome_annotation_lca, i.e. mid-run instead of here.
        self.assertTrue(set(self.codes.values()) <= SUPPORTED_CODES,
                        f"unsupported code(s): {set(self.codes.values()) - SUPPORTED_CODES}")

    def test_class_names_are_normalised(self):
        for key in list(self.codes) + list(self.ambiguous):
            self.assertEqual(key, key.strip().lower())

    def test_every_entry_records_its_basis(self):
        payload = json.loads(ASSET.read_text())
        for entry in payload["codes"]:
            self.assertIn(entry.get("confidence"), {"ncbi", "convention"})
            self.assertTrue(entry.get("basis"))
        for reason in payload["ambiguous"].values():
            self.assertTrue(reason.strip())


class InvertClassCoverageTests(unittest.TestCase):
    """Every invertebrate class must resolve to a code or be a documented abort."""

    def setUp(self):
        self.codes, self.ambiguous = cs.load_genetic_codes(ASSET)
        source = (BIN / "create_samplesheet.py").read_text()
        block = re.search(r"INVERT_CLASSES = frozenset\(\{(.*?)\}\)", source, re.S).group(1)
        # Strip comments first: the Craniata note quotes the name it warns about.
        body = "\n".join(line.split("#")[0] for line in block.splitlines())
        self.invert_classes = {c.lower() for c in re.findall(r"'([^']+)'", body)}

    def test_invert_classes_is_not_empty(self):
        self.assertGreater(len(self.invert_classes), 50)

    def test_every_invert_class_is_mapped_or_documented_ambiguous(self):
        uncovered = sorted(self.invert_classes - set(self.codes) - set(self.ambiguous))
        self.assertEqual(uncovered, [], f"classes with no code and no reason: {uncovered}")

    def test_craniata_is_the_brachiopod_class(self):
        # Craniata is a homonym: NCBI has it as both the brachiopod class (115366)
        # and the vertebrate subphylum (89593). taxdump_lineage indexes only
        # species/genus/family/order/class, so the subphylum is unreachable and a
        # resolved class of 'Craniata' can only mean the brachiopod. It must stay
        # an invertebrate class -- dropping it makes craniid brachiopods silently
        # vertebrate, the worse of the two failures.
        self.assertIn("craniata", self.invert_classes)
        self.assertEqual(self.codes.get("craniata"), 5)

    def test_current_crustacean_classes_are_covered(self):
        # Barnacles resolve to Thecostraca and copepods to Copepoda; Maxillopoda
        # and Hexanauplia are retired names the taxdump no longer carries. A
        # missing class here does not abort -- it reads as invertebrates=false and
        # the sample is annotated as a vertebrate.
        for tax_class in ["thecostraca", "copepoda", "malacostraca", "ostracoda"]:
            self.assertIn(tax_class, self.invert_classes, tax_class)
            self.assertEqual(self.codes.get(tax_class), 5, tax_class)

    def test_known_codes(self):
        for tax_class, expected in [
            ("anthozoa", 4), ("porifera", 4), ("ctenophora", 4),
            ("bivalvia", 5), ("gastropoda", 5), ("polychaeta", 5), ("nematoda", 5),
            ("asteroidea", 9), ("rhabditophora", 9),
            ("ascidiacea", 13), ("thaliacea", 13),
            ("appendicularia", 5), ("tunicata", 5),
            ("trematoda", 9), ("platyhelminthes", 9), ("placozoa", 4),
        ]:
            self.assertEqual(self.codes.get(tax_class), expected, tax_class)

    def test_classes_spanning_several_codes_stay_unmapped(self):
        # Only the pterobranchs genuinely vary below class rank: nodes.dmp gives
        # Rhabdopleuridae 5 and Cephalodiscidae 33. NCBI is unambiguous on the
        # flatworms and Placozoa, so those are mapped rather than left to abort.
        for tax_class in ["pterobranchia", "hemichordata"]:
            self.assertNotIn(tax_class, self.codes)
            self.assertIn(tax_class, self.ambiguous)

    def test_the_map_agrees_with_ncbis_own_per_taxon_assignment(self):
        """Every mapped class must match nodes.dmp field 8 where NCBI names it.

        The map is only the fallback for a class with no taxdump lineage behind
        it, so it must not contradict the per-taxon codes the resolver returns
        for everything else. This test is what caught Appendicularia and Tunicata
        being listed as code 13 when NCBI assigns both 5.
        """
        taxdump = Path("/scratch/pawsey1348/tpeirce/taxonkit_dbs")
        if not (taxdump / "nodes.dmp").is_file():
            self.skipTest("no taxdump available")
        sys.path.insert(0, str(ROOT / "bin"))
        try:
            from taxdump_lineage import TaxdumpLineage
        finally:
            sys.path.remove(str(ROOT / "bin"))
        resolver = TaxdumpLineage(str(taxdump))
        resolver.load()
        mismatches = []
        for tax_class, code in self.codes.items():
            taxid = resolver._resolve_taxid(tax_class)
            if taxid is None:
                continue                      # retired or ambiguous name
            ncbi = resolver._mito_code.get(taxid)
            if ncbi is not None and ncbi != code:
                mismatches.append((tax_class, code, ncbi))
        self.assertEqual(mismatches, [], f"map disagrees with NCBI: {mismatches}")


class LoaderContradictionTests(unittest.TestCase):
    def _load(self, payload):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
            json.dump(payload, handle)
            path = handle.name
        return cs.load_genetic_codes(path)

    def test_a_class_under_two_codes_is_rejected(self):
        with self.assertRaises(ValueError):
            self._load({"codes": [{"code": 4, "classes": ["anthozoa"]},
                                  {"code": 5, "classes": ["anthozoa"]}]})

    def test_the_same_code_twice_is_fine(self):
        codes, _ = self._load({"codes": [{"code": 4, "classes": ["anthozoa"]},
                                         {"code": 4, "classes": ["anthozoa", "hydrozoa"]}]})
        self.assertEqual(codes, {"anthozoa": 4, "hydrozoa": 4})

    def test_mapped_and_ambiguous_is_rejected(self):
        with self.assertRaises(ValueError):
            self._load({"codes": [{"code": 5, "classes": ["trematoda"]}],
                        "ambiguous": {"trematoda": "table 21"}})

    def test_no_asset_yields_empty_maps(self):
        self.assertEqual(cs.load_genetic_codes(None), ({}, {}))


class ResolveGeneticCodeTests(unittest.TestCase):
    def setUp(self):
        self.codes, _ = cs.load_genetic_codes(ASSET)

    def test_mapped_class_resolves_to_a_string_code(self):
        self.assertEqual(cs.resolve_genetic_code("Anthozoa", self.codes), "4")

    def test_case_and_whitespace_insensitive(self):
        self.assertEqual(cs.resolve_genetic_code("  BIVALVIA ", self.codes), "5")

    def test_unmapped_and_unknown_classes_resolve_blank(self):
        # Blank defers to prepare_samplesheet: the vertebrate default for a
        # vertebrate, an abort for an invertebrate.
        for tax_class in ["Pterobranchia", "Actinopteri", "unknown", "", None]:
            self.assertEqual(cs.resolve_genetic_code(tax_class, self.codes), "")


class ReportUnresolvedGeneticCodeTests(unittest.TestCase):
    def _report(self, rows, ambiguous=None):
        import io
        from contextlib import redirect_stderr
        buffer = io.StringIO()
        with redirect_stderr(buffer):
            cs.report_unresolved_genetic_code(rows, ambiguous or {})
        return buffer.getvalue()

    def test_unresolved_invertebrate_is_named_with_its_reason(self):
        out = self._report(
            [{"sample": "OGX", "invertebrates": "true", "class": "Pterobranchia",
              "genetic_code": ""}],
            {"pterobranchia": "spans codes 24 and 33"})
        self.assertIn("OGX", out)
        self.assertIn("Pterobranchia", out)
        self.assertIn("spans codes 24 and 33", out)

    def test_resolved_and_vertebrate_rows_are_silent(self):
        self.assertEqual(self._report([
            {"sample": "OGA", "invertebrates": "true", "class": "Anthozoa", "genetic_code": "4"},
            {"sample": "OGB", "invertebrates": "false", "class": "Actinopteri", "genetic_code": ""},
        ]), "")

    def test_a_sample_is_reported_once(self):
        rows = [{"sample": "OGX", "invertebrates": "true", "class": "Placozoa",
                 "genetic_code": ""}] * 3
        self.assertEqual(self._report(rows).count("OGX"), 1)


if __name__ == "__main__":
    unittest.main()
