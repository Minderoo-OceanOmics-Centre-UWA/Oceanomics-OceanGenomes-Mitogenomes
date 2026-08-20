import gzip
import importlib.util
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ena_package", ROOT / "bin" / "ena_package.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EnaPackageTests(unittest.TestCase):
    def test_circular_and_reverse_complement_equivalence(self):
        sequence = "AACCGT"
        rotated = "CGTAAC"
        reverse = MODULE.reverse_complement(sequence)
        self.assertEqual(
            MODULE.sequence_digests(sequence)[1],
            MODULE.sequence_digests(rotated)[1],
        )
        self.assertEqual(
            MODULE.sequence_digests(sequence)[1],
            MODULE.sequence_digests(reverse)[1],
        )
        self.assertNotEqual(
            MODULE.sequence_digests(sequence)[0],
            MODULE.sequence_digests(rotated)[0],
        )

    def test_manifest_uses_existing_full_seqid_and_mean_depth(self):
        text = MODULE.manifest_text(
            study="PRJEB123419",
            biosample="SAMEA123",
            full_seqid="OG910.hifi.241127.v3mitohifi.emma102",
            coverage=123.5,
            program="MitoHiFi 3.2.3",
            platform="PACBIO_SMRT",
            flatfile_name="record.embl.gz",
            chromosome_list_name="record.chromosome_list.tsv.gz",
            scientific_name="Choerodon rubescens",
        )
        self.assertIn("ASSEMBLYNAME\tOG910.hifi.241127.v3mitohifi.emma102\n", text)
        self.assertIn("COVERAGE\t123.5\n", text)
        self.assertIn("ASSEMBLY_TYPE\tclone or isolate\n", text)

    def test_build_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">internal\nAACCGT\n")
            seqid = "OG910.hifi.241127.v3mitohifi.emma102"
            embl = root / "input.embl"
            embl.write_text(
                f"ID   {seqid}; SV 1; circular; genomic DNA; STD; UNC; 6 BP.\n"
                f"AC * _{seqid}\n"
                "XX\nSQ   Sequence 6 BP;\n     aaccgt 6\n//\n"
            )
            package = root / "package"
            args = Namespace(
                og_id="OG910",
                assembly_prefix="OG910.hifi.241127.v3mitohifi",
                annotation_version="emma102",
                full_seqid=seqid,
                fasta=str(fasta),
                embl=str(embl),
                study="PRJEB123419",
                biosample="SAMEA123",
                coverage=100.0,
                program="MitoHiFi 3.2.3",
                platform="PACBIO_SMRT",
                scientific_name="Choerodon rubescens",
                run_accession=[],
                outdir=str(package),
            )
            self.assertEqual(MODULE.build_package(args), 0)
            with gzip.open(package / f"{seqid}.chromosome_list.tsv.gz", "rt") as handle:
                self.assertEqual(
                    handle.read(),
                    f"{seqid}\tMT\tCircular-Chromosome\tMitochondrion\n",
                )
            metadata = json.loads(
                (package / f"{seqid}.package_metadata.json").read_text()
            )
            # Status fields the package deliberately does not store: readiness
            # is derived from biosample_accession (it changes without the file
            # changing), and the durable path is the caller's to know.
            for absent in (
                "local_validation_status",
                "package_status",
                "metadata_blocker",
                "published_package_path",
            ):
                self.assertNotIn(absent, metadata)
            self.assertEqual(list(package.glob("*.local_validation.tsv")), [])
            # No --flatfile-status supplied, so no verdict is invented.
            self.assertEqual(
                metadata["flatfile_validation"]["status"], "NOT_REQUESTED"
            )

    def test_flatfile_verdict_is_recorded_in_the_package(self):
        """The pipeline's last gate travels with the package it validated."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">x\nAACCGT\n")
            seqid = "OG910.hifi.241127.v3mitohifi.emma102"
            embl = root / "input.embl"
            embl.write_text(
                f"ID   {seqid}; SV 1; circular; genomic DNA; STD; UNC; 6 BP.\n"
                f"AC * _{seqid}\n"
                "XX\nSQ   Sequence 6 BP;\n     aaccgt 6\n//\n"
            )
            status = root / "webin_status.tsv"
            status.write_text(
                "sample\tstatus\treason\twebin_exit\tvalidation_attempt"
                "\terror_count\twarning_count\twebin_cli_version\n"
                "OG910.hifi.241127.v3mitohifi\tPASS\tvalidated\t0\tinitial"
                "\t0\t4\t9.0.3\n"
            )
            package = root / "package"
            args = Namespace(
                og_id="OG910",
                assembly_prefix="OG910.hifi.241127.v3mitohifi",
                annotation_version="emma102",
                full_seqid=seqid,
                fasta=str(fasta),
                embl=str(embl),
                study="PRJEB123419",
                biosample="SAMEA123",
                coverage=100.0,
                program="MitoHiFi 3.2.3",
                platform="PACBIO_SMRT",
                scientific_name="Choerodon rubescens",
                run_accession=[],
                outdir=str(package),
                flatfile_status=str(status),
            )
            self.assertEqual(MODULE.build_package(args), 0)
            metadata_path = package / f"{seqid}.package_metadata.json"
            verdict = json.loads(metadata_path.read_text())["flatfile_validation"]
            self.assertEqual(verdict["status"], "PASS")
            self.assertEqual(verdict["reason"], "validated")
            self.assertEqual(verdict["error_count"], 0)
            self.assertEqual(verdict["warning_count"], 4)
            self.assertEqual(verdict["webin_cli_version"], "9.0.3")
            # A later manifest refresh must not drop it.
            refreshed = MODULE.refresh_package(package, {"mean_depth": 42.0})
            self.assertEqual(refreshed["flatfile_validation"]["status"], "PASS")

    def test_a_failing_flatfile_still_produces_a_package_that_says_so(self):
        self.assertEqual(
            MODULE.flatfile_validation(None),
            {
                "status": "NOT_REQUESTED",
                "reason": "not_requested",
                "error_count": None,
                "warning_count": None,
                "webin_cli_version": None,
            },
        )

    def test_missing_biosample_builds_blocked_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">x\nAACCGT\n")
            seqid = "OG5.ilmn.260101.getorg1770.emma102"
            embl = root / "input.embl"
            embl.write_text(f"ID   {seqid};\nAC * _{seqid}\n//\n")
            package = root / "package"
            args = Namespace(
                og_id="OG5",
                assembly_prefix="OG5.ilmn.260101.getorg1770",
                annotation_version="emma102",
                full_seqid=seqid,
                fasta=str(fasta),
                embl=str(embl),
                study="PRJEB123419",
                biosample="",
                coverage=5.0,
                program="GetOrganelle 1.7.7.1",
                platform="ILLUMINA",
                scientific_name="Testus organismus",
                run_accession=[],
                outdir=str(package),
            )
            MODULE.build_package(args)
            metadata = json.loads(
                (package / f"{seqid}.package_metadata.json").read_text()
            )
            self.assertIsNone(metadata["biosample_accession"])
            sequence_digest = metadata["sequence_sha256"]
            refreshed = MODULE.refresh_package(
                package, {"biosample_accession": "SAMEA123"}
            )
            self.assertEqual(refreshed["biosample_accession"], "SAMEA123")
            self.assertEqual(refreshed["sequence_sha256"], sequence_digest)
            self.assertIn(
                "SAMPLE\tSAMEA123",
                (package / f"{seqid}.manifest.txt").read_text(),
            )

    def manifest_with_biosample(self, biosample):
        return MODULE.manifest_text(
            study="PRJEB123419",
            biosample=biosample,
            full_seqid="OG910.hifi.241127.v3mitohifi.emma102",
            coverage=123.5,
            program="MitoHiFi 3.2.3",
            platform="PACBIO_SMRT",
            flatfile_name="record.embl.gz",
            chromosome_list_name="record.chromosome_list.tsv.gz",
            scientific_name="Choerodon rubescens",
        )

    def test_manifest_requires_an_ena_registered_biosample(self):
        # webin resolves SAMPLE against ENA's submission sample service, which
        # only knows Webin-registered samples. Verified against
        # ena-webin-cli 9.0.3 -context genome -validate -test: SAMEA132129018
        # validates, SAMN40589646 fails with "sample is null".
        self.assertIn(
            "SAMPLE\tSAMEA132129018\n",
            self.manifest_with_biosample("SAMEA132129018"),
        )

    def test_manifest_rejects_biosample_registered_outside_ena(self):
        for accession in ("SAMN40589646", "SAMD00000001"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError) as caught:
                    self.manifest_with_biosample(accession)
                self.assertIn("registered outside ENA", str(caught.exception))

    def test_manifest_rejects_non_biosample_values(self):
        for accession in ("", "SAMX123", "SAMN", "SAMN123x", "SRS123456", "123456"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError) as caught:
                    self.manifest_with_biosample(accession)
                self.assertIn("Invalid or missing", str(caught.exception))

    def test_ncbi_only_biosample_blocks_with_its_own_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">x\nAACCGT\n")
            seqid = "OG5.ilmn.260101.getorg1770.emma102"
            embl = root / "input.embl"
            embl.write_text(f"ID   {seqid};\nAC * _{seqid}\n//\n")
            package = root / "package"
            args = Namespace(
                og_id="OG5",
                assembly_prefix="OG5.ilmn.260101.getorg1770",
                annotation_version="emma102",
                full_seqid=seqid,
                fasta=str(fasta),
                embl=str(embl),
                study="PRJEB123419",
                biosample="SAMN40589646",
                coverage=5.0,
                program="GetOrganelle 1.7.7.1",
                platform="ILLUMINA",
                scientific_name="Testus organismus",
                run_accession=[],
                outdir=str(package),
            )
            MODULE.build_package(args)
            metadata = json.loads(
                (package / f"{seqid}.package_metadata.json").read_text()
            )
            self.assertEqual(metadata["biosample_accession"], "SAMN40589646")
            # A blocked manifest must never carry a SAMPLE line webin would reject.
            manifest = (package / f"{seqid}.manifest.txt").read_text()
            self.assertTrue(manifest.startswith("# BLOCKED:"))
            self.assertNotIn("SAMPLE\t", manifest)

            # Brokering the specimen into ENA unblocks it in place.
            refreshed = MODULE.refresh_package(
                package, {"biosample_accession": "SAMEA132129018"}
            )
            self.assertEqual(refreshed["biosample_accession"], "SAMEA132129018")
            self.assertIn(
                "SAMPLE\tSAMEA132129018",
                (package / f"{seqid}.manifest.txt").read_text(),
            )


class CollaboratorArtefactTests(unittest.TestCase):
    seqid = "OG910.hifi.241127.v3mitohifi.emma102"

    def build_args(self, root, package, **overrides):
        fasta = root / "input.fa"
        fasta.write_text(f">{self.seqid}\nAACCGT\n")
        embl = root / "input.embl"
        embl.write_text(
            f"ID   {self.seqid}; SV 1; circular; genomic DNA; STD; UNC; 6 BP.\n"
            f"AC * _{self.seqid}\n"
            "XX\nSQ   Sequence 6 BP;\n     aaccgt 6\n//\n"
        )
        args = Namespace(
            og_id="OG910",
            assembly_prefix="OG910.hifi.241127.v3mitohifi",
            annotation_version="emma102",
            full_seqid=self.seqid,
            fasta=str(fasta),
            embl=str(embl),
            study="PRJEB123419",
            biosample="SAMEA123",
            coverage=100.0,
            program="MitoHiFi 3.2.3",
            platform="PACBIO_SMRT",
            scientific_name="Choerodon rubescens",
            run_accession=[],
            outdir=str(package),
        )
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_package_carries_the_fasta_and_an_untouched_gff(self):
        """Locus tags are assigned downstream, so the GFF is passed through as-is.

        Byte-identical is the assertion that matters: the packaging step used to
        rewrite ID=/Parent= from locus tags, and nothing here may reintroduce a
        rewrite that the downstream pipeline would then have to undo.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff_text = (
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tgene\t1\t68\t.\t+\t.\tID=uuid-gene-tf;Name=TF\n"
            )
            gff = root / "in.gff"
            gff.write_text(gff_text)
            package = root / "package"
            args = self.build_args(root, package, gff=str(gff))
            self.assertEqual(MODULE.build_package(args), 0)

            self.assertEqual(
                (package / f"{self.seqid}.fa").read_text(), f">{self.seqid}\nAACCGT\n"
            )
            packaged_gff = (package / f"{self.seqid}.gff").read_text()
            self.assertEqual(packaged_gff, gff_text)
            self.assertNotIn("locus_tag", packaged_gff)

            metadata = json.loads((package / f"{self.seqid}.package_metadata.json").read_text())
            self.assertNotIn("gff_locus_tag_coverage", metadata)
            # Both collaborator files are covered by the package checksums.
            checksums = (package / "checksums.sha256").read_text()
            self.assertIn(f"{self.seqid}.fa\n", checksums)
            self.assertIn(f"{self.seqid}.gff\n", checksums)

    def test_package_carries_the_genes_fasta_renamed_onto_full_seqid(self):
        """The extractor names it on the assembly prefix; the package uses full_seqid.

        Every other file in the directory shares that stem, and a collaborator
        should not have to know that this one file was keyed differently.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes_text = f">{self.seqid}|1-6|+|COX1 [organism=x]\nAACCGT\n"
            genes = root / "OG910.hifi.241127.v3mitohifi.genes.fa"
            genes.write_text(genes_text)
            package = root / "package"
            args = self.build_args(root, package, genes=str(genes))
            self.assertEqual(MODULE.build_package(args), 0)

            self.assertEqual(
                (package / f"{self.seqid}.genes.fa").read_text(), genes_text
            )
            self.assertIn(
                f"{self.seqid}.genes.fa\n", (package / "checksums.sha256").read_text()
            )

    def test_package_holds_no_locus_tag_mapping(self):
        """The mapping TSV is gone: nothing in this pipeline allocates tags."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            package = root / "package"
            self.assertEqual(MODULE.build_package(self.build_args(root, package)), 0)
            self.assertEqual(
                list(package.glob("*.locus_tag_mapping.tsv")), []
            )

    def test_collaborator_files_do_not_move_the_package_digest(self):
        """package_digest means 'the ENA submission content changed'.

        The downstream submission pipeline uses it to tell a rebuilt package
        apart from an unchanged one, so a convenience artefact must not shift it.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tgene\t1\t68\t.\t+\t.\tID=uuid-gene-tf;Name=TF\n"
            )
            genes = root / "OG910.hifi.241127.v3mitohifi.genes.fa"
            genes.write_text(f">{self.seqid}|1-6|+|COX1\nAACCGT\n")
            digests = []
            for name, extra in (
                ("without", {}),
                ("with", {"gff": str(gff), "genes": str(genes)}),
            ):
                package = root / name
                MODULE.build_package(self.build_args(root, package, **extra))
                digests.append(
                    json.loads(
                        (package / f"{self.seqid}.package_metadata.json").read_text()
                    )["package_digest"]
                )
            self.assertEqual(digests[0], digests[1])


if __name__ == "__main__":
    unittest.main()
