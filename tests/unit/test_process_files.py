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


if __name__ == "__main__":
    unittest.main()
