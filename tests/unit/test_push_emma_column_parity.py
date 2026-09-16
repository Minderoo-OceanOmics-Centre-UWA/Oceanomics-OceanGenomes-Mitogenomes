"""Every annotation column the uploader INSERTs must also be refreshed on conflict.

A column present in the INSERT list but missing from the ON CONFLICT ... DO UPDATE
list is written once and then never again: the row keeps whatever the first run put
there, forever, while every other column moves on. That is precisely how a stale
`annotation_integrity_issues` string ends up sitting next to an annotation that a
later run repaired.

The only columns legitimately absent from the DO UPDATE list are the conflict key
itself, which cannot be updated because it is what the rows are matched on.

Parsed from the SQL text rather than asserted against a hand-maintained list, so
adding a column to one place and forgetting the other fails here rather than in
production six months later.
"""

import re
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "push_emma_annotation_results.py"
CONFLICT_KEY = {"og_id", "tech", "seq_date", "code"}


def _sql():
    return SCRIPT.read_text()


class ColumnParityTests(unittest.TestCase):
    def setUp(self):
        text = _sql()
        insert = re.search(r"INSERT INTO mitogenome_data \((.*?)\)\s*VALUES", text, re.S)
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


if __name__ == "__main__":
    unittest.main()
