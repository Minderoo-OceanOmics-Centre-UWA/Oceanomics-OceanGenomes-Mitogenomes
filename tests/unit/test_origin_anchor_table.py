"""Validation of the generated asset assets/taxonomy/mito_origin_anchors.json.

Modelled on test_mito_genetic_codes.py, which cross-checks the Groovy sets from
Python by regex-parsing lib/InvertTaxonGroups.groovy. Same reasoning applies
here: the table is read by BOTH a Groovy loader and a Python one, and the two
must agree about what a valid table is.

The asset is GENERATED, so the strongest check available is that regenerating it
reproduces the tracked bytes. That is what makes "regenerate after a database
rebuild" enforceable rather than aspirational.
"""

import json
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from mito_gene_order import is_valid_anchor, load_origin_anchors, resolve_origin_anchor  # noqa: E402

ASSET = ROOT / "assets" / "taxonomy" / "mito_origin_anchors.json"
GROOVY = ROOT / "lib" / "InvertTaxonGroups.groovy"
TAXDUMP = Path("/scratch/pawsey1348/tpeirce/lca_cache/taxdump")


@unittest.skipUnless(ASSET.exists(), "anchor table not generated")
class AssetShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = json.loads(ASSET.read_text())
        cls.orders, cls.classes, cls.groups, cls.policy = load_origin_anchors(ASSET)

    def test_every_anchor_is_an_addressable_mitos_key(self):
        for level in ("orders", "classes", "groups"):
            for taxon, entry in self.raw[level].items():
                with self.subTest(level=level, taxon=taxon):
                    self.assertTrue(
                        is_valid_anchor(entry["anchor"]),
                        f"{level}[{taxon}] anchor {entry['anchor']!r} is not a key "
                        f"mitos_to_emma can look up")

    def test_no_ambiguous_leu_or_ser_anchor(self):
        # The reference DBs record bare "tRNA-Leu"/"tRNA-Ser" without a copy number,
        # so TL/TS can win a tally. They are not addressable in a MITOS annotation
        # and the generator must have fallen back a level instead of shipping one.
        for level in ("orders", "classes", "groups"):
            for taxon, entry in self.raw[level].items():
                self.assertNotIn(entry["anchor"], ("TL", "TS"),
                                 f"{level}[{taxon}] shipped an ambiguous tRNA anchor")

    def test_order_level_is_present(self):
        """Guards the silent taxdump degradation.

        A table generated without --taxdump-dir has no order level, which resolves
        Scleractinia from the Anthozoa class aggregate (no majority) and rotates
        every submitted stony coral off tRNA-Met.
        """
        self.assertEqual(self.policy["levels"][0], "order")
        self.assertGreater(len(self.orders), 10)

    def test_scleractinia_is_trnM(self):
        """The row that must never regress.

        Stony corals deposit from tRNA-Met (63.8% of 58 records) and everything
        already submitted was published that way. Getting this wrong changes those
        sequences' checksums and ENA treats them as different records.
        """
        self.assertEqual(self.orders["scleractinia"], "TM")

    def test_anthozoan_orders_do_not_all_agree(self):
        """The reason the table is keyed on order at all."""
        anthozoan = {o: self.orders[o] for o in
                     ("scleractinia", "malacalcyonacea", "zoantharia") if o in self.orders}
        self.assertEqual(len(anthozoan), 3)
        self.assertGreater(len(set(anthozoan.values())), 1,
                           "if these agreed, a class-level table would have sufficed")

    def test_rows_are_internally_consistent(self):
        # `level` names which level supplied the ANSWER, which is not always the
        # level the row lives at: a class row that inherits its phylum's anchor is
        # level="group" while carrying its own (small) n. Only a row resolved at its
        # OWN level is claiming to have cleared the bars.
        own_level = {"orders": "order", "classes": "class", "groups": "group"}
        for level in ("orders", "classes", "groups"):
            for taxon, entry in self.raw[level].items():
                with self.subTest(level=level, taxon=taxon):
                    self.assertIn("level", entry)
                    self.assertIn("n", entry)
                    if entry["level"] == own_level[level]:
                        self.assertGreaterEqual(entry["n"], self.policy["min_records"])
                        self.assertGreaterEqual(entry["fraction"],
                                                self.policy["min_fraction"])
                    else:
                        # Inherited or defaulted rows must say why.
                        self.assertTrue(entry.get("reason"),
                                        f"{level}[{taxon}] level={entry['level']} "
                                        f"carries no reason")
                    counts = entry.get("counts")
                    if counts:
                        self.assertLessEqual(sum(counts.values()), entry["n"])

    def test_inherited_rows_actually_match_what_they_inherit_from(self):
        """A row that says it took the group anchor must carry the group's anchor."""
        for taxon, entry in self.raw["classes"].items():
            if entry["level"] != "group":
                continue
            with self.subTest(taxon=taxon):
                group = entry.get("group")
                self.assertIsNotNone(group, f"classes[{taxon}] inherits but names no group")
                self.assertEqual(entry["anchor"], self.raw["groups"][group]["anchor"])

    def test_default_anchor_is_cox1(self):
        self.assertEqual(self.policy["default_anchor"], "CO1")


@unittest.skipUnless(ASSET.exists() and GROOVY.exists(), "asset or lib missing")
class GroovyAgreementTests(unittest.TestCase):
    """Python and Groovy must recognise the same class vocabulary."""

    @staticmethod
    def _groovy_classes():
        text = GROOVY.read_text()
        names = set()
        for set_name in ("CNIDARIA", "PORIFERA", "MOLLUSCA", "ARTHROPODA",
                         "ECHINODERMATA", "CTENOPHORA", "TUNICATA", "ANNELIDA"):
            m = re.search(
                rf"static\s+final\s+Set<String>\s+{set_name}_CLASSES\s*=\s*\[(.*?)\]\s*as\s+Set",
                text, re.S)
            assert m, f"could not parse {set_name}_CLASSES"
            names.update(n.strip().lower() for n in re.findall(r"'([^']+)'", m.group(1)))
        return names

    def test_every_invertebrate_class_resolves_to_some_anchor(self):
        """The 'never aborts' contract, unlike the genetic-code lookup."""
        orders, classes, groups, policy = load_origin_anchors(ASSET)
        groovy_classes = self._groovy_classes()
        self.assertGreater(len(groovy_classes), 30)
        # class -> group, mirroring InvertTaxonGroups.seedDbGroup().
        text = GROOVY.read_text()
        class_to_group = {}
        for set_name, group in (("CNIDARIA", "anthozoa"), ("PORIFERA", "porifera"),
                                ("MOLLUSCA", "mollusca"), ("ARTHROPODA", "arthropoda"),
                                ("ECHINODERMATA", "echinodermata"),
                                ("CTENOPHORA", "ctenophora"), ("TUNICATA", "tunicata"),
                                ("ANNELIDA", "annelida")):
            m = re.search(
                rf"static\s+final\s+Set<String>\s+{set_name}_CLASSES\s*=\s*\[(.*?)\]\s*as\s+Set",
                text, re.S)
            for n in re.findall(r"'([^']+)'", m.group(1)):
                class_to_group[n.strip().lower()] = group

        for taxon_class in sorted(groovy_classes):
            with self.subTest(taxon_class=taxon_class):
                anchor = resolve_origin_anchor("", taxon_class, orders, classes,
                                               groups, policy, class_to_group)
                self.assertTrue(is_valid_anchor(anchor))

    def test_a_class_with_no_records_still_inherits_its_phylum(self):
        """Calcarea has zero records in the porifera DB, so it has no class row.

        Without the group step it would take the cox1 default instead of the rrnL
        its phylum uses 90% of the time.
        """
        orders, classes, groups, policy = load_origin_anchors(ASSET)
        self.assertNotIn("calcarea", classes)
        self.assertEqual(
            resolve_origin_anchor("", "Calcarea", orders, classes, groups, policy,
                                  {"calcarea": "porifera"}),
            "RNR2")


@unittest.skipUnless(ASSET.exists() and TAXDUMP.is_dir(), "taxdump not available")
class AssetIsUpToDateTests(unittest.TestCase):
    def test_regenerating_reproduces_the_tracked_file(self):
        """--check is what makes 'regenerate on rebuild' enforceable."""
        result = subprocess.run(
            [sys.executable, str(ROOT / "bin" / "build_origin_anchor_table.py"),
             "--taxdump-dir", str(TAXDUMP), "--check"],
            capture_output=True, text=True, timeout=900)
        self.assertEqual(result.returncode, 0,
                         f"anchor table is stale:\n{result.stdout}\n{result.stderr}")


if __name__ == "__main__":
    unittest.main()
