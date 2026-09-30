"""The run snapshot is a comparison tool, so a false regression in it is expensive.

Both defects these tests pin were live at once in NOVA_260909_LA, and together they
manufactured two findings that sent a remediation after bugs that did not exist:
3 mtMutS and 61 OH annotations reported dropped when all were published, and 136
flagged CDS reported across 9 EMMA samples when the annotator had merely omitted
the terminators. Neither was a pipeline fault.
"""

import sys
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

from mito_run_snapshot import LINEAGE_FEATURES, bare, published_name


class PublishedNameTests(unittest.TestCase):
    def test_renamed_genes_resolve_to_their_published_name(self):
        """MITOS emits msh1; mitos_to_emma.py publishes it as MTMUTS."""
        self.assertEqual(published_name("msh1"), "MTMUTS")
        self.assertEqual(published_name("mutS"), "MTMUTS")
        self.assertEqual(published_name("mtMutS"), "MTMUTS")

    def test_non_gene_features_resolve_to_their_published_name(self):
        self.assertEqual(published_name("OH"), "REP_ORIGIN_H")
        self.assertEqual(published_name("OL"), "REP_ORIGIN_L")
        self.assertEqual(published_name("gpI"), "INTRON_GPI")
        self.assertEqual(published_name("gpII"), "INTRON_GPII")

    def test_names_the_rename_preserves_are_unchanged(self):
        """These two are why the bug hid: they matched case-insensitively anyway."""
        self.assertEqual(published_name("atp9"), "ATP9")
        self.assertEqual(published_name("lagli"), "LAGLI")

    def test_unknown_names_pass_through(self):
        self.assertEqual(published_name("cox1"), "cox1")
        self.assertEqual(published_name("nad5"), "nad5")

    def test_every_tracked_lineage_feature_resolves(self):
        """A feature tracked as lineage-specific but not mapped would read as lost."""
        for raw in LINEAGE_FEATURES:
            resolved = published_name(bare(raw))
            self.assertTrue(resolved, f"{raw} resolved to nothing")
            self.assertEqual(resolved, published_name(bare(raw).lower()),
                             f"{raw} must resolve case-insensitively")


if __name__ == "__main__":
    unittest.main()
