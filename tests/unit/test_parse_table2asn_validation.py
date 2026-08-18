import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PARSER = ROOT / "bin" / "parse_table2asn_validation.py"


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ParseTable2AsnValidationTests(unittest.TestCase):
    def run_parser(self, val_text="", dr_text=""):
        tmp = tempfile.TemporaryDirectory()
        root = Path(tmp.name)
        val = root / "sample.val"
        dr = root / "sample.dr"
        val.write_text(val_text)
        dr.write_text(dr_text)
        findings = root / "findings.tsv"
        status = root / "status.tsv"
        flags = root / "flags.tsv"
        subprocess.run(
            [
                str(PARSER),
                "--sample", "sample",
                "--circular", "true",
                "--val", str(val),
                "--dr", str(dr),
                "--findings", str(findings),
                "--status", str(status),
                "--qc-flags", str(flags),
            ],
            check=True,
        )
        return tmp, read_rows(findings), read_rows(status)[0], read_rows(flags)[0]

    def test_clean_reports_pass(self):
        tmp, findings, status, flags = self.run_parser()
        self.addCleanup(tmp.cleanup)
        self.assertEqual(findings, [])
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["error_count"], "0")
        self.assertEqual(flags["nostop_count"], "0")

    def test_warning_is_reported_but_passes(self):
        tmp, findings, status, _ = self.run_parser(
            "WARNING: valid [SEQ_INST.ShortSeq] Protein is shorter than 10 aa\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["warning_count"], "1")
        self.assertEqual(status["warning_codes"], "SEQ_INST.ShortSeq")
        self.assertEqual(findings[0]["severity"], "WARNING")

    def test_error_reject_and_nostop_quarantine_sample(self):
        tmp, _, status, flags = self.run_parser(
            "ERROR: valid [SEQ_FEAT.NoStop] Missing stop codon FEATURE: CDS: COX1 [lcl|x]\n"
            "REJECT: valid [SEQ_DESCR.NoOrgFound] No organism name\n"
            "ERROR: valid [SEQ_FEAT.NoStop] Missing stop codon FEATURE: CDS: NAD5 [lcl|x]\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "FAIL_TABLE2ASN")
        self.assertEqual(status["error_count"], "2")
        self.assertEqual(status["reject_count"], "1")
        self.assertEqual(status["blocking_codes"], "SEQ_FEAT.NoStop,SEQ_DESCR.NoOrgFound")
        self.assertEqual(flags["nostop_count"], "2")
        self.assertEqual(flags["nostop_features"], "COX1;NAD5")

    def test_fatal_discrepancy_quarantines_sample(self):
        tmp, findings, status, _ = self.run_parser(
            dr_text="DISC_TEST: FATAL: MISSING_GENES: expected mitochondrial genes absent\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "FAIL_TABLE2ASN")
        self.assertEqual(status["fatal_discrepancy_count"], "1")
        self.assertEqual(findings[0]["severity"], "FATAL")

    def test_no_locus_tags_discrepancy_is_advisory(self):
        tmp, findings, status, _ = self.run_parser(
            dr_text="FATAL: NO_LOCUS_TAGS: None of 37 genes has locus tag.\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["fatal_discrepancy_count"], "0")
        self.assertEqual(status["warning_count"], "1")
        self.assertEqual(status["warning_codes"], "NO_LOCUS_TAGS")
        self.assertEqual(findings[0]["severity"], "WARNING")

    def test_missing_protein_id_discrepancy_is_advisory(self):
        tmp, findings, status, _ = self.run_parser(
            dr_text="FATAL: MISSING_PROTEIN_ID: 13 proteins have invalid IDs.\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["fatal_discrepancy_count"], "0")
        self.assertEqual(status["warning_count"], "1")
        self.assertEqual(status["warning_codes"], "MISSING_PROTEIN_ID")
        self.assertEqual(findings[0]["severity"], "WARNING")

    def test_mixed_case_severity_is_recognised(self):
        # table2asn writes "Error:", not "ERROR:". Matching only upper case sent
        # every real finding into the UNPARSED fallback and reported PASS.
        tmp, findings, status, _ = self.run_parser(
            "Error: valid [SEQ_FEAT.StartCodon] Illegal start codon used. "
            "Wrong genetic code [2] or protein should be partial FEATURE: CDS: ATP6 [lcl|x]\n"
            "Error: valid [SEQ_INST.BadProteinStart] gap symbol at start of protein "
            "sequence (ATP6) BIOSEQ: lcl|x\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "FAIL_TABLE2ASN")
        self.assertEqual(status["error_count"], "2")
        self.assertEqual(status["info_count"], "0")
        self.assertEqual(
            status["blocking_codes"], "SEQ_FEAT.StartCodon,SEQ_INST.BadProteinStart"
        )
        self.assertTrue(all(f["severity"] == "ERROR" for f in findings))

    def test_ncbi_only_validator_codes_are_advisory(self):
        # These fire against NCBI submission rules that the ENA route does not
        # apply, so they must stay visible without quarantining the assembly.
        tmp, findings, status, _ = self.run_parser(
            "Info: valid [SEQ_DESCR.LatLonWater] Lat_lon '14.0 S 121.8 E' is closest to "
            "'Australia' at distance 221 km, but in water 'Indian Ocean'\n"
            "Error: valid [SEQ_DESCR.OrganismIsUndefinedSpecies] Organism 'Chaunax sp.' "
            "is undefined species and does not have a specific identifier.\n"
            "Error: valid [SEQ_FEAT.GeneXrefWithoutLocus] Feature has Gene Xref with "
            "locus_tag but no locus FEATURE: tRNA: Phe [lcl|x]\n"
        )
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["error_count"], "0")
        self.assertEqual(status["warning_count"], "3")
        self.assertTrue(all(f["severity"] == "WARNING" for f in findings))

    def test_unparsed_line_is_a_warning_not_info(self):
        # An unrecognised line means the output format moved, not that the record
        # is clean; it has to stay visible.
        tmp, findings, status, _ = self.run_parser("something entirely unexpected\n")
        self.addCleanup(tmp.cleanup)
        self.assertEqual(status["status"], "PASS")
        self.assertEqual(status["warning_count"], "1")
        self.assertEqual(findings[0]["code"], "UNPARSED")
        self.assertEqual(findings[0]["severity"], "WARNING")


if __name__ == "__main__":
    unittest.main()
