"""Unit tests for filename-derived metadata in process_files.py."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "process_files", ROOT / "bin" / "process_files.py"
)
process_files = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(process_files)


class AssemblyMethodTests(unittest.TestCase):
    def test_getorganelle(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1.ilmn.260101.v177getorg.emma102"),
            "GetOrganelle v.1.7.7",
        )

    def test_mitohifi(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1.hifi.260101.v323mitohifi.emma102"),
            "MitoHifi v.3.2.3",
        )

    def test_oatk(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1422.hifi.260227.v10oatk.emma102"),
            "Oatk v.1.0",
        )

    def test_unknown_assembler_is_rejected(self):
        with self.assertRaisesRegex(SystemExit, "Unknown assembler code"):
            process_files.derive_assembly_method("OG1.hifi.260101.v10unknown.emma102")


class ProcessTblFileTests(unittest.TestCase):
    TBL = (
        ">Feature OG1.hic.260101.v177getorg\n"
        "2847\t3821\tgene\n"
        "\t\t\tgene\tMT-ND1\n"
        "2847\t3821\tCDS\n"
        "\t\t\tproduct\tNADH dehydrogenase subunit 1\n"
        "\t\t\ttransl_table\t2\n"
        "\t\t\tprotein_id\tgnl|Emma|69ef424b-9f7c-5045-9ab9-9c90d1db603a\n"
    )

    def test_protein_id_placeholder_is_stripped(self):
        with tempfile.TemporaryDirectory() as tmp:
            tbl_in = Path(tmp) / "in.tbl"
            tbl_out = Path(tmp) / "out.tbl"
            tbl_in.write_text(self.TBL)

            process_files.process_tbl_gb_file(tbl_in, tbl_out, "OG1.hic.260101.v177getorg")

            out_text = tbl_out.read_text()
            self.assertNotIn("protein_id", out_text)
            self.assertIn("NADH dehydrogenase subunit 1", out_text)


class ProcessGffFileTests(unittest.TestCase):
    GFF = (
        "##gff-version 3\n"
        "##sequence-region\tOG1.hifi.260101.v323mitohifi\t1\t100\n"
        "# retain this comment\n"
        "\n"
        "OG1.hifi.260101.v323mitohifi\tEmma\tregion\t1\t100\t.\t+\t0\tIs_circular=true\n"
        "OG1.hifi.260101.v323mitohifi\tEmma\tgene\t1\t68\t.\t+\t.\t"
        "ID=gene-1;Name=12srna;Note=putative mitochondrial gene\n"
    )

    def test_all_sequence_ids_use_annotation_prefix(self):
        annotation_prefix = "OG1.hifi.260101.v323mitohifi.emma102"

        with tempfile.TemporaryDirectory() as tmp:
            gff_in = Path(tmp) / f"{annotation_prefix}.gff"
            gff_out = Path(tmp) / "processed.gff"
            gff_in.write_text(self.GFF)

            process_files.process_gff_file(gff_in, gff_out, annotation_prefix)

            lines = gff_out.read_text().splitlines()
            self.assertEqual(
                lines[1],
                f"##sequence-region\t{annotation_prefix}\t1\t100",
            )
            self.assertEqual(lines[2], "# retain this comment")
            self.assertEqual(lines[3], "")

            feature_lines = [
                line.split("\t")
                for line in lines
                if line and not line.startswith("#")
            ]
            self.assertEqual({fields[0] for fields in feature_lines}, {annotation_prefix})
            self.assertEqual(feature_lines[0][2], "region")
            self.assertIn("Name=RNR1", feature_lines[1][8])
            self.assertNotIn("putative", feature_lines[1][8])


# Emma sorts feature starts as strings, so 10053 lands ahead of 1027 and the
# 12S gene at 71 lands last. ND6 is minus strand, written start > end.
UNSORTED_TBL = """>Feature emma
1\t70\tgene
\t\t\tgene\tMT-TF
1\t70\ttRNA
\t\t\tproduct\ttRNA-Phe(GAA)
10053\t10349\tgene
\t\t\tgene\tMT-ND4L
10053\t10349\tCDS
\t\t\tproduct\tputative NADH dehydrogenase subunit 4L
\t\t\tprotein_id\tgnl|Emma|0f0e5b4c-0000-0000-0000-000000000000
1027\t1098\tgene
\t\t\tgene\tMT-TV
1027\t1098\ttRNA
\t\t\tproduct\ttRNA-Val(UAC)
14279\t13755\tgene
\t\t\tgene\tMT-ND6
14279\t13755\tCDS
\t\t\tproduct\tNADH dehydrogenase subunit 6
71\t1026\tgene
\t\t\tgene\tMT-RNR1
71\t1026\trRNA
\t\t\tproduct\t12S rRNA
"""

UNSORTED_GFF = """##gff-version 3
##sequence-region\temma\t1\t16745
emma\tEmma\tregion\t1\t16745\t.\t+\t0\tIs_circular=true
emma\tEmma\tgene\t1\t70\t.\t+\t.\tID=g1;Name=MT-TF
emma\tEmma\ttRNA\t1\t70\t.\t+\t.\tID=r1;Parent=g1;Name=MT-TF
emma\tEmma\tgene\t10053\t10349\t.\t+\t.\tID=g2;Name=MT-ND4L
emma\tEmma\tmRNA\t10053\t10349\t.\t+\t.\tID=m2;Parent=g2;Name=MT-ND4L
emma\tEmma\tCDS\t10053\t10349\t.\t+\t0\tID=c2;Parent=m2;Name=MT-ND4L
emma\tEmma\tgene\t1027\t1098\t.\t+\t.\tID=g3;Name=trnV-UAC
emma\tEmma\ttRNA\t1027\t1098\t.\t+\t.\tID=r3;Parent=g3;Name=trnV-UAC
emma\tEmma\tgene\t71\t1026\t.\t+\t.\tID=g4;Name=12srna
emma\tEmma\trRNA\t71\t1026\t.\t+\t.\tID=r4;Parent=g4;Name=12srna
"""


def _feature_lows(text):
    lows = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) >= 3 and cols[0][:1].isdigit():
            lows.append(min(int(cols[0]), int(cols[1])))
    return lows


class ProcessTblOrderingTests(unittest.TestCase):
    """Emma's string sort has to be undone before anything downstream reads it.

    The submission pipeline allocates locus tags by walking the feature table,
    so whatever order leaves here is the order the published tags carry.
    """

    def _process(self, tbl=UNSORTED_TBL, seq=None, genetic_code=2):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "in.tbl"
            source.write_text(tbl)
            output = root / "out.tbl"
            process_files.process_tbl_gb_file(
                source, output, "emma", seq=seq, genetic_code=genetic_code
            )
            return output.read_text()

    def test_features_are_written_in_coordinate_order(self):
        # Each locus contributes a gene line and one child line.
        text = self._process()
        self.assertEqual(
            _feature_lows(text),
            [1, 1, 71, 71, 1027, 1027, 10053, 10053, 13755, 13755],
        )

    def test_feature_header_stays_first(self):
        self.assertTrue(self._process().startswith(">Feature emma\n"))

    def test_qualifiers_stay_with_their_feature(self):
        # RNR1 moves from last to second; its product qualifier has to move too.
        lines = self._process().splitlines()
        self.assertEqual(
            lines[lines.index("71\t1026\trRNA") + 1], "\t\t\tproduct\t12S rRNA"
        )

    def test_existing_normalisation_survives_the_sort(self):
        text = self._process()
        self.assertNotIn("MT-", text)
        self.assertNotIn("putative ", text)
        self.assertNotIn("protein_id", text)

    def test_transl_except_lands_on_the_right_cds(self):
        # A polyA note with no transl_except triggers the stop-codon fix-up.
        # ND4L is extended to 10053..10350 so its length leaves one base over,
        # and that base is a T, i.e. the truncated TAA the poly(A) tail
        # completes. After sorting, the qualifier has to still sit under that
        # CDS, which the sort moved from second to fourth.
        seq = list("A" * 16745)
        seq[10350 - 1] = "T"
        tbl = UNSORTED_TBL.replace("10053\t10349\tCDS", "10053\t10350\tCDS").replace(
            "\t\t\tproduct\tputative NADH dehydrogenase subunit 4L\n",
            "\t\t\tproduct\tputative NADH dehydrogenase subunit 4L\n"
            "\t\t\tnote\tTAA stop codon is completed by the addition of "
            f"{process_files.POLYA_NOTE_MARK} to the mRNA\n",
        )
        lines = self._process(tbl=tbl, seq="".join(seq)).splitlines()
        cds = lines.index("10053\t10350\tCDS")
        block = lines[cds:cds + 4]
        self.assertIn("\t\t\ttransl_except\t(pos:10350,aa:TERM)", block)
        # and it is the ND4L block, not whichever CDS happened to sort there
        self.assertTrue(any("subunit 4L" in line for line in block))

    def _with_start_note(self, tbl, cds_line, product_line):
        """Add EMMA's non-standard-start note to the CDS whose product line matches."""
        return tbl.replace(product_line, product_line + (
            f"\t\t\tnote\t{process_files.START_NOTE_MARK} TTG\n"
        ))

    def test_start_transl_except_on_a_plus_strand_cds(self):
        # ND4L runs 10053..10349 on the plus strand. Put a TTG at its first codon
        # and mark it the way EMMA does; the initiation codon has to be declared
        # as aa:Met so table2asn stops emitting a gap symbol for residue 1.
        seq = list("A" * 16745)
        seq[10053 - 1:10053 + 2] = list("TTG")
        tbl = self._with_start_note(
            UNSORTED_TBL, "10053\t10349\tCDS",
            "\t\t\tproduct\tputative NADH dehydrogenase subunit 4L\n",
        )
        lines = self._process(tbl=tbl, seq="".join(seq)).splitlines()
        cds = lines.index("10053\t10349\tCDS")
        block = lines[cds:cds + 4]
        self.assertIn("\t\t\ttransl_except\t(pos:10053..10055,aa:Met)", block)
        self.assertTrue(any("subunit 4L" in line for line in block))

    def test_start_transl_except_on_a_minus_strand_cds(self):
        # ND6 is written 14279..13755, so its 5' end is the HIGH coordinate and
        # the initiation codon is read off the reverse complement.
        seq = list("A" * 16745)
        seq[14277 - 1:14277 + 2] = list("CAA")   # revcomp -> TTG
        tbl = self._with_start_note(
            UNSORTED_TBL, "14279\t13755\tCDS",
            "\t\t\tproduct\tNADH dehydrogenase subunit 6\n",
        )
        lines = self._process(tbl=tbl, seq="".join(seq)).splitlines()
        cds = lines.index("14279\t13755\tCDS")
        block = lines[cds:cds + 4]
        self.assertIn("\t\t\ttransl_except\t(pos:complement(14277..14279),aa:Met)", block)
        self.assertTrue(any("subunit 6" in line for line in block))

    def test_canonical_start_codon_is_left_alone(self):
        # A note can be stale or the codon already legal; only a codon outside the
        # table's start set earns a transl_except.
        seq = list("A" * 16745)
        seq[10053 - 1:10053 + 2] = list("ATG")
        tbl = self._with_start_note(
            UNSORTED_TBL, "10053\t10349\tCDS",
            "\t\t\tproduct\tputative NADH dehydrogenase subunit 4L\n",
        )
        self.assertNotIn("aa:Met", self._process(tbl=tbl, seq="".join(seq)))

    def test_start_transl_except_without_the_note(self):
        # EMMA only sometimes notes a non-standard start (it stayed silent on
        # OG663's CTG-initiated ATP6), so the sequence has to be what triggers the
        # declaration. The note the pipeline writes in EMMA's place names the codon.
        seq = list("A" * 16745)
        seq[10053 - 1:10053 + 2] = list("CTG")
        lines = self._process(seq="".join(seq)).splitlines()
        cds = lines.index("10053\t10349\tCDS")
        block = lines[cds:cds + 4]
        self.assertIn("\t\t\ttransl_except\t(pos:10053..10055,aa:Met)", block)
        self.assertIn(f"\t\t\tnote\t{process_files.START_NOTE_MARK} CTG", block)

    def test_implausible_start_codon_is_left_to_fail(self):
        # AAA is an initiator in no mitochondrial code, so it reads as a broken CDS
        # boundary or a bad base call rather than alternative initiation. Declaring
        # aa:Met there would hide the problem, so it is left for table2asn to flag.
        self.assertNotIn("aa:Met", self._process(seq="A" * 16745))

    def test_partial_start_cds_is_left_alone(self):
        # A CDS already marked 5'-partial has no initiation codon to declare.
        seq = list("A" * 16745)
        seq[10053 - 1:10053 + 2] = list("TTG")
        tbl = self._with_start_note(
            UNSORTED_TBL, "10053\t10349\tCDS",
            "\t\t\tproduct\tputative NADH dehydrogenase subunit 4L\n",
        ).replace("10053\t10349\tCDS", "<10053\t10349\tCDS")
        self.assertNotIn("aa:Met", self._process(tbl=tbl, seq="".join(seq)))


class ProcessGffOrderingTests(unittest.TestCase):
    def _process(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "in.gff"
            source.write_text(UNSORTED_GFF)
            output = root / "out.gff"
            process_files.process_gff_file(source, output, "emma")
            return output.read_text()

    def _records(self):
        return [
            line.split("\t")
            for line in self._process().splitlines()
            if not line.startswith("#") and len(line.split("\t")) == 9
        ]

    def test_records_are_written_in_coordinate_order(self):
        starts = [int(fields[3]) for fields in self._records()]
        # Leading 1 is the whole-molecule region record.
        self.assertEqual(starts, [1, 1, 1, 71, 71, 1027, 1027, 10053, 10053, 10053])

    def test_directives_and_region_stay_ahead_of_the_features(self):
        lines = self._process().splitlines()
        self.assertEqual(lines[0], "##gff-version 3")
        self.assertEqual(lines[2].split("\t")[2], "region")

    def test_children_stay_with_their_parent(self):
        self.assertEqual(
            [fields[2] for fields in self._records()],
            ["region", "gene", "tRNA", "gene", "rRNA", "gene", "tRNA",
             "gene", "mRNA", "CDS"],
        )

    def test_name_map_still_applies(self):
        text = self._process()
        self.assertIn("Name=RNR1", text)
        self.assertIn("Name=TV", text)


if __name__ == "__main__":
    unittest.main()
