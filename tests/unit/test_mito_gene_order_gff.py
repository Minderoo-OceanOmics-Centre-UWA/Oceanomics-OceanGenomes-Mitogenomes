"""Unit tests for the shared GFF reader in bin/mito_gene_order.py.

The reference order was already shared to stop the QC step and the rescue gates
drifting on what "present and in order" means. Agreeing on the list was not
enough: they also have to agree on how a GFF becomes an ordered gene list.

parse_gff_attributes had three copies and genes_by_coord had two byte-identical
ones, while annotation_stats.py kept a third, DIFFERENT implementation inline that
deduped a repeated gene by first line seen rather than by lowest start. On a gene
split across the origin those two readings disagree, so one GFF could be judged
in-order by a rescue gate and out-of-order by the QC step -- and a failed order
check both fails QC and silently suppresses the rescue.

Pure stdlib.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mgo = _load("mito_gene_order")
stats = _load("annotation_stats")
emma_gate = _load("emma_rescue_gate")
trna_gate = _load("trna_rescue_gate")


def gene_line(name, start, end, strand="+"):
    return (f"chr\tEmma\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID=gene-{name}-{start};Name=MT-{name}")


class ParseGffAttributesTests(unittest.TestCase):
    def test_parses_key_value_pairs(self):
        self.assertEqual(
            mgo.parse_gff_attributes("ID=gene-CO1;Name=MT-CO1"),
            {"ID": "gene-CO1", "Name": "MT-CO1"})

    def test_ignores_entries_with_no_equals(self):
        self.assertEqual(mgo.parse_gff_attributes("ID=x;junk;Name=y"),
                         {"ID": "x", "Name": "y"})

    def test_empty_attribute_column(self):
        self.assertEqual(mgo.parse_gff_attributes(""), {})


class GeneEntriesByCoordTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _write(self, lines):
        p = Path(self.tmp) / "t.gff"
        p.write_text("##gff-version 3\n##sequence-region chr 1 16500\n"
                     + "\n".join(lines) + "\n")
        return str(p)

    def test_strips_the_MT_prefix_and_sorts_by_start(self):
        p = self._write([gene_line("CO2", 300, 400), gene_line("CO1", 100, 200)])
        self.assertEqual(mgo.genes_by_coord(p), ["CO1", "CO2"])

    def test_split_gene_is_kept_once_at_its_lowest_start(self):
        # An origin-spanning ND5 written as two gene lines, with the HIGH half
        # listed first -- the file order that made first-seen dedup wrong.
        p = self._write([
            gene_line("ND5", 16000, 16400),
            gene_line("CO1", 500, 600),
            gene_line("ND5", 100, 300),
        ])
        entries = mgo.gene_entries_by_coord(p)
        names = [e[0] for e in entries]
        self.assertEqual(names, ["ND5", "CO1"])
        self.assertEqual(names.count("ND5"), 1)
        # The retained entry is the LOWER fragment, which is what makes the
        # ordering meaningful.
        self.assertEqual(entries[0][1], 100)
        self.assertEqual(entries[0][2], 300)

    def test_gene_lines_without_a_name_are_skipped(self):
        p = self._write([
            "chr\tEmma\tgene\t100\t200\t.\t+\t.\tID=gene-nameless",
            gene_line("CO1", 300, 400),
        ])
        self.assertEqual(mgo.genes_by_coord(p), ["CO1"])

    def test_non_gene_feature_lines_are_ignored(self):
        p = self._write([
            gene_line("CO1", 100, 200),
            "chr\tEmma\tCDS\t100\t200\t.\t+\t.\tID=feat-CO1;Name=MT-CO1",
            "chr\tEmma\ttRNA\t300\t400\t.\t+\t.\tID=feat-TF;Name=MT-TF",
        ])
        # TF has a tRNA line but no gene line, so it is not a gene entry.
        self.assertEqual(mgo.genes_by_coord(p), ["CO1"])

    def test_strand_is_carried_through(self):
        p = self._write([gene_line("ND6", 100, 200, strand="-")])
        self.assertEqual(mgo.gene_entries_by_coord(p)[0][3], "-")


class ConsumersAgreeTests(unittest.TestCase):
    """The regression this file exists for.

    All three consumers must read one GFF into one gene list. Before the reader
    was shared, annotation_stats.py disagreed with both gates on a split gene.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _canonical_gff_with_split_origin_gene(self):
        # A complete, canonically ordered mitogenome whose FIRST gene (TF) spans
        # the origin and so is written as two gene lines, the high half first.
        lines = []
        span = 16500
        lines.append(gene_line("TF", span - 40, span))
        pos = 100
        for g in mgo.REF_GENES:
            end = pos + 50
            if g == "TF":
                lines.append(gene_line("TF", pos, end))
            else:
                lines.append(gene_line(g, pos, end))
            ftype = ("rRNA" if g.startswith("RNR")
                     else "tRNA" if g.startswith("T") else "CDS")
            lines.append(f"chr\tEmma\t{ftype}\t{pos}\t{end}\t.\t+\t.\t"
                         f"ID=feat-{g};Parent=gene-{g}-{pos};Name=MT-{g}")
            pos += 100
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text("##gff-version 3\n##sequence-region chr 1 %d\n" % span
                     + "\n".join(lines) + "\n")
        return p

    def test_all_three_consumers_read_the_same_gene_list(self):
        p = self._canonical_gff_with_split_origin_gene()
        shared = mgo.genes_by_coord(str(p))
        self.assertEqual(emma_gate.genes_by_coord(str(p)), shared)
        self.assertEqual(trna_gate.genes_by_coord(str(p)), shared)
        # annotation_stats reaches the same list through process_gff.
        summary = stats.process_gff(str(p), p.stem, genetic_code=2)
        self.assertEqual(shared, list(mgo.REF_GENES))
        # The whole point: the split origin gene does NOT make this look
        # out of order to the QC step while the gates call it fine.
        self.assertEqual(summary["order_correct"], "yes")
        self.assertEqual(summary["passed"], "yes")


if __name__ == "__main__":
    unittest.main()
