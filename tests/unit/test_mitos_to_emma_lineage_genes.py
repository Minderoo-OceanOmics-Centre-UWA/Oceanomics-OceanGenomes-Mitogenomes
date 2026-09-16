"""Non-vertebrate features must survive the MITOS -> EMMA conversion.

Across NOVA_260909_LA's 58 MITOS annotations, every one of these was silently
destroyed by map_gene_name returning None and parse_bed calling `continue`:
22 atp9, 3 msh1, 27 OH, 4 OL, 4 gpI/gpII and 2 lagli. The only signal was an
aggregate warning that fires when ZERO features survive, so 24 samples lost real
genes without a word in the logs.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import mitos_to_emma as m
from gene_names import submission_gene_name


def write_bed(tmp, rows):
    p = Path(tmp) / "result.bed"
    p.write_text("\n".join("\t".join(str(c) for c in r) for r in rows) + "\n")
    return p


class VocabularyTests(unittest.TestCase):
    def test_atp9_maps_to_a_cds(self):
        self.assertEqual(m.map_gene_name("atp9")[:2], ("ATP9", "CDS"))

    def test_every_mtmuts_alias_maps(self):
        for alias in ("msh1", "MSH1", "mutS", "mtMutS", "MTMUTS"):
            self.assertEqual(m.map_gene_name(alias)[:2], ("MTMUTS", "CDS"), alias)

    def test_replication_origins_are_rep_origin_not_genes(self):
        self.assertEqual(m.map_gene_name("OH")[:2], ("REP_ORIGIN_H", "rep_origin"))
        self.assertEqual(m.map_gene_name("OL")[:2], ("REP_ORIGIN_L", "rep_origin"))

    def test_introns_and_auxiliary_orf(self):
        self.assertEqual(m.map_gene_name("gpI")[:2], ("INTRON_GPI", "intron"))
        self.assertEqual(m.map_gene_name("gpII")[:2], ("INTRON_GPII", "intron"))
        self.assertEqual(m.map_gene_name("lagli")[:2], ("LAGLI", "misc_feature"))

    def test_genuinely_unknown_name_still_returns_none(self):
        self.assertEqual(m.map_gene_name("wibble")[0], None)


class ProductTests(unittest.TestCase):
    def test_products_defined_for_the_new_genes(self):
        self.assertEqual(m.product_for("ATP9", "CDS"), "ATP synthase F0 subunit 9")
        self.assertIn("mismatch repair", m.product_for("MTMUTS", "CDS"))

    def test_missing_product_raises_instead_of_shipping_the_gene_name(self):
        """PRODUCT.get(name, name) used to fall back to the bare gene name, which
        sends /product=ATP9 -- a gene name, not a product -- to table2asn."""
        with self.assertRaises(KeyError):
            m.product_for("NOTAGENE", "CDS")

    def test_non_genes_have_no_product(self):
        self.assertEqual(m.product_for("REP_ORIGIN_H", "rep_origin"), "")


class SubmissionNameTests(unittest.TestCase):
    def test_mtmuts_is_restored_to_insdc_case(self):
        """The internal key is uppercase so it can index features and filenames;
        INSDC spells it mtMutS and will not accept MTMUTS."""
        self.assertEqual(submission_gene_name("MTMUTS"), "mtMutS")
        self.assertEqual(submission_gene_name("MT-MTMUTS"), "mtMutS")

    def test_atp9_passes_through_unchanged(self):
        self.assertEqual(submission_gene_name("ATP9"), "ATP9")
        self.assertEqual(submission_gene_name("MT-CO1"), "CO1")


class LocusClusteringTests(unittest.TestCase):
    def _rows(self, spans):
        return [{"start": s, "end": e, "strand": st, "frag": f}
                for s, e, st, f in spans]

    def test_single_row_is_one_locus(self):
        self.assertEqual(len(m.cluster_loci(self._rows([(1, 60, "+", None)]))), 1)

    def test_explicit_part_suffixes_stay_one_locus(self):
        """MITOS writing nad5-a/nad5-b is the annotator stating they are exons of
        one gene, however far apart they sit."""
        rows = self._rows([(1, 500, "+", "a"), (9000, 9500, "+", "b")])
        self.assertEqual(len(m.cluster_loci(rows)), 1)

    def test_distant_unsuffixed_copies_split_into_separate_loci(self):
        """Two separate copies of a gene used to merge into one impossible feature
        spanning the gap between them."""
        rows = self._rows([(1, 70, "+", None), (15000, 15070, "+", None)])
        self.assertEqual(len(m.cluster_loci(rows, max_intron=10000)), 2)

    def test_opposite_strands_never_merge(self):
        rows = self._rows([(1, 70, "+", None), (200, 270, "-", None)])
        self.assertEqual(len(m.cluster_loci(rows, max_intron=10000)), 2)

    def test_nearby_same_strand_rows_are_one_intron_split_locus(self):
        rows = self._rows([(1, 500, "+", None), (900, 1400, "+", None)])
        self.assertEqual(len(m.cluster_loci(rows, max_intron=10000)), 1)


class AuditTests(unittest.TestCase):
    def test_audit_records_every_row_including_unmapped(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            bed = write_bed(tmp, [
                ["c", 0, 60, "cox1", 0, "+"],
                ["c", 100, 160, "atp9", 0, "+"],
                ["c", 200, 260, "OH", 0, "+"],
                ["c", 300, 360, "lagli", 0, "+"],
                ["c", 400, 460, "wibble", 0, "+"],
            ])
            _features, audit = m.parse_bed(bed)
        self.assertEqual(len(audit), 5, "every raw row must appear in the audit")
        by_disposition = {}
        for row in audit:
            by_disposition.setdefault(row["disposition"], []).append(row["raw_name"])
        self.assertIn("atp9", by_disposition["retained"])
        self.assertIn("OH", by_disposition["retained_non_gene"])
        self.assertIn("lagli", by_disposition["known_auxiliary"])
        self.assertEqual(by_disposition["unmapped"], ["wibble"])


if __name__ == "__main__":
    unittest.main()
