import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "allocate_ena_locus_tags", ROOT / "bin" / "allocate_ena_locus_tags.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


TBL = """>Feature internal
1\t100\tgene
\t\t\tgene\tCOX1
1\t100\tCDS
\t\t\tgene\tCOX1
\t\t\tproduct\tcytochrome c oxidase subunit I
150\t220\tgene
\t\t\tgene\tTRNF
150\t220\ttRNA
\t\t\tgene\tTRNF
230\t300\tcontrol_region
\t\t\tnote\tD-loop
"""


# Emma orders its feature table by the start coordinate as a string, so a
# molecule arrives as 1, 10053, 1027 ... rather than in coordinate order. The
# minus-strand gene is written start > end the way a feature table encodes it,
# and the joined CDS carries a bare second interval line.
UNSORTED_TBL = """>Feature internal
1\t70\tgene
\t\t\tgene\tTF
1\t70\ttRNA
\t\t\tgene\tTF
10053\t10349\tgene
\t\t\tgene\tND4L
10053\t10349\tCDS
\t\t\tgene\tND4L
1027\t1098\tgene
\t\t\tgene\tTV
1027\t1098\ttRNA
\t\t\tgene\tTV
14279\t13755\tgene
\t\t\tgene\tND6
14279\t13755\tCDS
\t\t\tgene\tND6
71\t1026\tgene
\t\t\tgene\tRNR1
71\t900\trRNA
901\t1026
\t\t\tgene\tRNR1
"""


def _serials(registry):
    return {gene: int(row["gene_serial"]) for (gene, _occurrence), row in registry.items()}


class AllocateEnaLocusTagsTests(unittest.TestCase):
    def test_tags_gene_and_children_but_not_control_region(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(TBL)
            preamble, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            registry = MODULE.allocate_file_registry(
                roots, "OG910", root / "registry.tsv"
            )
            MODULE.render_tags(registry, "OG910", "OGMTHIFI")
            MODULE.inject_tags(preamble, features, registry, root / "tagged.tbl")
            text = (root / "tagged.tbl").read_text()
            self.assertEqual(text.count("OGMTHIFI_000910001"), 2)
            self.assertEqual(text.count("OGMTHIFI_000910002"), 2)
            control = text.split("230\t300\tcontrol_region", 1)[1]
            self.assertNotIn("locus_tag", control)

    def test_registry_stores_serials_not_rendered_tags(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(TBL)
            _, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            registry_path = root / "registry.tsv"
            MODULE.allocate_file_registry(roots, "OG910", registry_path)
            with registry_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertNotIn("locus_tag", rows[0])
            self.assertEqual([row["gene_serial"] for row in rows], ["1", "2"])

    def test_same_serial_renders_per_technology_prefix(self):
        # One shared specimen serial, one tag per technology study: ENA forbids
        # two published records sharing a locus tag.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(TBL)
            _, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            registry_path = root / "registry.tsv"
            hifi = MODULE.render_tags(
                MODULE.allocate_file_registry(roots, "OG910", registry_path),
                "OG910",
                "OGMTHIFI",
            )[("COX1", 1)]["locus_tag"]
            hic = MODULE.render_tags(
                MODULE.allocate_file_registry(roots, "OG910", registry_path),
                "OG910",
                "OGMTHIC",
            )[("COX1", 1)]["locus_tag"]
            self.assertEqual(hifi, "OGMTHIFI_000910001")
            self.assertEqual(hic, "OGMTHIC_000910001")

    def test_registry_is_reused_across_versions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(TBL)
            _, first_features = MODULE.parse_tbl(source)
            first_roots = MODULE.assign_locus_identities(first_features)
            registry_path = root / "registry.tsv"
            first = MODULE.allocate_file_registry(first_roots, "OG910", registry_path)
            first_serial = first[("COX1", 1)]["gene_serial"]
            changed = TBL.replace("1\t100", "5\t104").replace("150\t220", "154\t224")
            source.write_text(changed)
            _, second_features = MODULE.parse_tbl(source)
            second_roots = MODULE.assign_locus_identities(second_features)
            second = MODULE.allocate_file_registry(second_roots, "OG910", registry_path)
            self.assertEqual(first_serial, second[("COX1", 1)]["gene_serial"])
            with registry_path.open(newline="") as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 2)

    def test_serials_follow_coordinates_not_file_order(self):
        # Emma's string sort put ND4L at 10053 ahead of TV at 1027 and RNR1 at
        # 71 last. Serials have to ascend with position on the molecule, not
        # with the order the annotator happened to write.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(UNSORTED_TBL)
            _, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            serials = _serials(
                MODULE.allocate_file_registry(roots, "OG910", root / "registry.tsv")
            )
            self.assertEqual(
                serials, {"TF": 1, "RNR1": 2, "TV": 3, "ND4L": 4, "ND6": 5}
            )

    def test_minus_strand_gene_sorts_on_its_lower_coordinate(self):
        # ND6 is written 14279..13755. Sorting on the first column alone would
        # be harmless here, but a gene starting below a neighbour's end would
        # land in the wrong place, so the key is the lower coordinate.
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.tbl"
            source.write_text(UNSORTED_TBL)
            _, features = MODULE.parse_tbl(source)
            nd6 = next(f for f in features if f.key == "gene" and f.start == 14279)
            self.assertEqual(nd6.coordinate_key, (13755, 14279))
            self.assertEqual(nd6.strand, "-")

    def test_tagged_output_and_mapping_are_coordinate_ordered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(UNSORTED_TBL)
            preamble, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            registry = MODULE.render_tags(
                MODULE.allocate_file_registry(roots, "OG910", root / "registry.tsv"),
                "OG910",
                "OGMTHIFI",
            )
            for feature in features:
                if feature.canonical_gene:
                    feature.locus_tag = registry[
                        (feature.canonical_gene, feature.occurrence)
                    ]["locus_tag"]
            MODULE.inject_tags(preamble, features, registry, root / "tagged.tbl")
            MODULE.write_mapping(root / "mapping.tsv", "OG910.hifi.x", features)

            tagged = (root / "tagged.tbl").read_text().splitlines()
            self.assertEqual(tagged[0], ">Feature internal")
            lows = [
                min(int(cols[0]), int(cols[1]))
                for cols in (line.split("\t") for line in tagged)
                if len(cols) >= 3 and cols[0][:1].isdigit()
            ]
            self.assertEqual(lows, sorted(lows))

            with (root / "mapping.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            tags = [row["locus_tag"] for row in rows]
            self.assertEqual(tags, sorted(tags))
            self.assertEqual(rows[0]["canonical_gene"], "TF")
            self.assertEqual(rows[-1]["canonical_gene"], "ND6")

    def test_joined_interval_stays_with_its_feature(self):
        # The bare "901 1026" continuation line has no feature key, so it must
        # ride along with the rRNA block rather than being sorted on its own.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tbl"
            source.write_text(UNSORTED_TBL)
            preamble, features = MODULE.parse_tbl(source)
            roots = MODULE.assign_locus_identities(features)
            registry = MODULE.render_tags(
                MODULE.allocate_file_registry(roots, "OG910", root / "registry.tsv"),
                "OG910",
                "OGMTHIFI",
            )
            MODULE.inject_tags(preamble, features, registry, root / "tagged.tbl")
            lines = (root / "tagged.tbl").read_text().splitlines()
            self.assertEqual(lines[lines.index("901\t1026") - 1], "71\t900\trRNA")

    def test_orphan_child_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.tbl"
            path.write_text(">Feature x\n1\t100\tCDS\n\t\t\tgene\tCOX1\n")
            _, features = MODULE.parse_tbl(path)
            with self.assertRaisesRegex(ValueError, "No gene feature"):
                MODULE.assign_locus_identities(features)


if __name__ == "__main__":
    unittest.main()
