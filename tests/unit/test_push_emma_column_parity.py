"""Every annotation column the uploader INSERTs must also be refreshed on conflict, and
every value it uploads must land in exactly one table.

A column present in the INSERT list but missing from the ON CONFLICT ... DO UPDATE
list is written once and then never again: the row keeps whatever the first run put
there, forever, while every other column moves on. That is precisely how a stale
`annotation_integrity_issues` string ends up sitting next to an annotation that a
later run repaired.

The only columns legitimately absent from the DO UPDATE list are the conflict key
itself, which cannot be updated because it is what the rows are matched on. Since Phase R
that key is (assembly_id, annotation): one row per annotation, so a different annotator adds
a row instead of overwriting the previous annotation (docs/assembly_relations_plan.md §3.3).

Parsed from the SQL text rather than asserted against a hand-maintained list, so
adding a column to one place and forgetting the other fails here rather than in
production six months later.
"""

import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "push_emma_annotation_results.py"


def _load_module():
    """Import the uploader without pandas/psycopg2 installed (see test_push_mtdna_assm_decide_write)."""
    stubbed = [dep for dep in ("pandas", "numpy", "psycopg2") if dep not in sys.modules]
    for dep in stubbed:
        sys.modules[dep] = types.ModuleType(dep)
    try:
        spec = importlib.util.spec_from_file_location("push_emma_annotation_results", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        for dep in stubbed:
            del sys.modules[dep]


MODULE = _load_module()
CONFLICT_KEY = {"assembly_id", "annotation"}
ASSEMBLY_KEY = {"og_id", "tech", "seq_date", "code"}


def _sql():
    return SCRIPT.read_text()


class ColumnParityTests(unittest.TestCase):
    def setUp(self):
        text = _sql()
        insert = re.search(r"INSERT INTO mito_annotation \((.*?)\)\s*VALUES", text, re.S)
        update = re.search(r"DO UPDATE SET(.*?)RETURNING", text, re.S)
        self.assertIsNotNone(insert, "could not find the INSERT column list")
        self.assertIsNotNone(update, "could not find the DO UPDATE list")
        self.insert_cols = [
            c.strip() for c in insert.group(1).replace("\n", " ").split(",") if c.strip()
        ]
        self.update_cols = [
            m.group(1) for m in re.finditer(r"(\w+)\s*=\s*EXCLUDED\.", update.group(1))
        ]

    def test_every_inserted_column_is_refreshed_except_the_conflict_key(self):
        stale = set(self.insert_cols) - set(self.update_cols)
        self.assertEqual(
            stale,
            CONFLICT_KEY,
            "these columns are INSERTed but never refreshed on conflict, so a rerun "
            f"cannot correct them: {sorted(stale - CONFLICT_KEY)}",
        )

    def test_no_column_is_refreshed_without_being_inserted(self):
        orphan = set(self.update_cols) - set(self.insert_cols)
        self.assertEqual(orphan, set(), f"refreshed but never inserted: {sorted(orphan)}")

    def test_values_placeholder_count_matches_insert_columns(self):
        text = _sql()
        values = re.search(r"VALUES \((.*?)\)\s*\n\s*ON CONFLICT", text, re.S)
        self.assertIsNotNone(values, "could not find the VALUES list")
        placeholders = [v.strip() for v in values.group(1).split(",") if v.strip()]
        self.assertEqual(
            len(placeholders),
            len(self.insert_cols),
            "INSERT column count and VALUES placeholder count disagree: "
            f"{len(self.insert_cols)} columns vs {len(placeholders)} placeholders",
        )

    def test_integrity_columns_are_refreshed(self):
        """The columns most damaging to leave stale are explicitly covered."""
        for col in ("annotation_integrity_status", "annotation_integrity_issues",
                    "order_status", "missing_expected_lineage_genes",
                    "mitos_unmapped_features", "annotation_updated_at"):
            self.assertIn(col, self.update_cols, f"{col} is not refreshed on conflict")


    def test_conflict_target_is_the_annotation_not_the_assembly(self):
        """The overwrite bug: keyed on the assembly, a second annotator replaced the first."""
        self.assertRegex(_sql(), r"ON CONFLICT \(assembly_id, annotation\)")
        self.assertNotRegex(_sql(), r"ON CONFLICT \(og_id, tech, seq_date, code\)\s*DO UPDATE")

    def test_every_uploaded_value_lands_in_exactly_one_place(self):
        """Each params key is an assembly key, a summary column, or a gene length -- once."""
        text = _sql()
        keys = set(re.findall(r'^\s*"(\w+)":\s*(?:row_dict|seq_date|ANNOTATION_STATS_VERSION)', text, re.M))
        genes = set(MODULE.GENE_COLUMNS)
        translated = {f"{g}_trans" for g in MODULE.TRANSLATED_GENES}
        summary = set(self.insert_cols) - CONFLICT_KEY
        homes = [ASSEMBLY_KEY, summary, genes, translated]
        for key in sorted(keys | {"annotation"}):
            hits = sum(key in home for home in homes) + (key in CONFLICT_KEY)
            self.assertEqual(hits, 1, f"{key} lands in {hits} places")
        for col in summary:
            self.assertIn(col, keys | {"annotation_updated_at"}, f"{col} is inserted but never uploaded")

    def test_selection_flag_is_never_written(self):
        """is_selected is a person's choice (owner decision 3); the uploader must not touch it."""
        code = "\n".join(line for line in _sql().splitlines() if not line.lstrip().startswith("#"))
        self.assertNotIn("is_selected", code)
        self.assertNotIn("selected_by", code)


class GeneRowTests(unittest.TestCase):
    def test_only_genes_with_a_value_get_a_row(self):
        params = {"cox1": 1551, "cox1_trans": 516, "rrna12s": None, "trna_phe": 69,
                  "atp9_trans": 254}
        rows = MODULE.gene_rows(7, params)
        self.assertEqual(sorted(rows), sorted([
            (7, "cox1", 1551, 516), (7, "trna_phe", 69, None), (7, "atp9", None, 254)]))

    def test_rna_genes_never_carry_a_translation(self):
        rows = MODULE.gene_rows(1, {"trna_phe": 69, "trna_phe_trans": 99})
        self.assertEqual(rows, [(1, "trna_phe", 69, None)])


if __name__ == "__main__":
    unittest.main()
