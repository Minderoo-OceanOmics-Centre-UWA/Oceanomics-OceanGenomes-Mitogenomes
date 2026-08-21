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

    def manifest(self, **overrides):
        fields = dict(
            full_seqid="OG910.hifi.241127.v3mitohifi.emma102",
            coverage=123.5,
            program="MitoHiFi 3.2.3",
            platform="PACBIO_SMRT",
            flatfile_name="record.embl.gz",
            chromosome_list_name="record.chromosome_list.tsv.gz",
            scientific_name="Choerodon rubescens",
        )
        fields.update(overrides)
        return MODULE.manifest_fields(**fields)

    def test_manifest_uses_existing_full_seqid_and_mean_depth(self):
        fields = self.manifest()
        self.assertEqual(
            fields["ASSEMBLYNAME"], "OG910.hifi.241127.v3mitohifi.emma102"
        )
        self.assertEqual(fields["COVERAGE"], "123.5")
        self.assertEqual(fields["ASSEMBLY_TYPE"], "clone or isolate")
        self.assertEqual(
            fields["DESCRIPTION"], "Choerodon rubescens mitochondrial genome"
        )

    def test_manifest_omits_the_keys_the_submitter_owns(self):
        """STUDY, SAMPLE and RUN_REF are registered downstream, so we do not guess.

        Emitting them empty or as a placeholder would read as a value this
        pipeline had and could not fill, which is a different claim.
        """
        for key in ("STUDY", "SAMPLE", "RUN_REF"):
            self.assertNotIn(key, self.manifest())

    def test_an_unusable_value_drops_only_its_own_key(self):
        """A bad field must not cost the eight good ones next to it."""
        self.assertNotIn("PLATFORM", self.manifest(platform="OXFORD_NANOPORE"))
        self.assertNotIn("COVERAGE", self.manifest(coverage=None))
        self.assertNotIn("COVERAGE", self.manifest(coverage=float("nan")))
        self.assertNotIn("COVERAGE", self.manifest(coverage=-1.0))
        self.assertNotIn("PROGRAM", self.manifest(program="  "))
        self.assertNotIn("DESCRIPTION", self.manifest(scientific_name=""))
        for fields in (
            self.manifest(platform="OXFORD_NANOPORE"),
            self.manifest(coverage=None),
        ):
            self.assertEqual(
                fields["ASSEMBLYNAME"], "OG910.hifi.241127.v3mitohifi.emma102"
            )
            self.assertEqual(fields["FLATFILE"], "record.embl.gz")
            self.assertEqual(fields["MOLECULETYPE"], "genomic DNA")

    def test_specimen_facts_come_from_the_packaged_flatfile(self):
        """Read back off the flatfile so the two can never disagree."""
        with tempfile.TemporaryDirectory() as tmp:
            embl = Path(tmp) / "record.embl"
            embl.write_text(
                "FH   Key             Location/Qualifiers\n"
                "FH\n"
                "FT   source          1..16631\n"
                'FT                   /organism="Ophthalmolepis lineolata"\n'
                'FT                   /organelle="mitochondrion"\n'
                'FT                   /mol_type="genomic DNA"\n'
                'FT                   /isolate="OG51"\n'
                'FT                   /tissue_type="Gills"\n'
                'FT                   /geo_loc_name="Australia: WA, New Year Island"\n'
                'FT                   /collection_date="29-Mar-2023"\n'
                "FT   gene            1..68\n"
                'FT                   /gene="TF"\n'
            )
            self.assertEqual(
                MODULE.source_qualifiers(embl),
                {
                    "organism": "Ophthalmolepis lineolata",
                    "organelle": "mitochondrion",
                    "mol_type": "genomic DNA",
                    "isolate": "OG51",
                    "tissue_type": "Gills",
                    "geo_loc_name": "Australia: WA, New Year Island",
                    "collection_date": "29-Mar-2023",
                },
            )

    def test_a_wrapped_qualifier_is_rejoined(self):
        """EMBL wraps a long value across lines; half of one is not the value."""
        with tempfile.TemporaryDirectory() as tmp:
            embl = Path(tmp) / "record.embl"
            embl.write_text(
                "FT   source          1..6\n"
                'FT                   /geo_loc_name="Australia: New South Wales,\n'
                'FT                   Changte Shoal, east of Coffs Harbour"\n'
                'FT                   /isolate="OG193"\n'
            )
            self.assertEqual(
                MODULE.source_qualifiers(embl)["geo_loc_name"],
                "Australia: New South Wales, Changte Shoal, east of Coffs Harbour",
            )

    def test_an_absent_qualifier_is_left_out_not_blanked(self):
        """Most specimens have no lat_lon; an empty one would read as recorded."""
        with tempfile.TemporaryDirectory() as tmp:
            embl = Path(tmp) / "record.embl"
            embl.write_text(
                "FT   source          1..6\n"
                'FT                   /organism="Testus organismus"\n'
                'FT                   /isolate="OG5"\n'
            )
            qualifiers = MODULE.source_qualifiers(embl)
            self.assertNotIn("lat_lon", qualifiers)
            self.assertEqual(qualifiers["isolate"], "OG5")

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
                coverage=100.0,
                program="MitoHiFi 3.2.3",
                platform="PACBIO_SMRT",
                scientific_name="Choerodon rubescens",
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
            # Status fields the package deliberately does not store: the
            # durable path is the caller's to know, and the BioSample, study and
            # run accessions belong to the submission pipeline.
            for absent in (
                "local_validation_status",
                "package_status",
                "metadata_blocker",
                "published_package_path",
                "biosample_accession",
                "biosample_source",
                "run_accessions",
                "study",
            ):
                self.assertNotIn(absent, metadata)
            self.assertEqual(metadata["schema_version"], 3)
            self.assertEqual(metadata["validation_study"], "PRJEB123419")
            self.assertEqual(list(package.glob("*.manifest.txt")), [])
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
                coverage=100.0,
                program="MitoHiFi 3.2.3",
                platform="PACBIO_SMRT",
                scientific_name="Choerodon rubescens",
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

    def test_refresh_rewrites_the_manifest_and_clears_a_legacy_file(self):
        """A package built before the manifest moved into the metadata."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">x\nAACCGT\n")
            seqid = "OG5.ilmn.260101.getorg1770.emma102"
            embl = root / "input.embl"
            embl.write_text(
                f"ID   {seqid};\nAC * _{seqid}\n"
                "FT   source          1..6\n"
                'FT                   /organism="Testus organismus"\n'
                "//\n"
            )
            package = root / "package"
            args = Namespace(
                og_id="OG5",
                assembly_prefix="OG5.ilmn.260101.getorg1770",
                annotation_version="emma102",
                full_seqid=seqid,
                fasta=str(fasta),
                embl=str(embl),
                study="PRJEB123419",
                coverage=5.0,
                program="GetOrganelle 1.7.7.1",
                platform="ILLUMINA",
                scientific_name="Testus organismus",
                outdir=str(package),
            )
            MODULE.build_package(args)
            self.assertEqual(
                json.loads((package / f"{seqid}.package_metadata.json").read_text())[
                    "specimen"
                ],
                {"organism": "Testus organismus"},
            )
            legacy = package / f"{seqid}.manifest.txt"
            legacy.write_text("STUDY\tPRJEB110568\n")
            sequence_digest = json.loads(
                (package / f"{seqid}.package_metadata.json").read_text()
            )["sequence_sha256"]

            refreshed = MODULE.refresh_package(package, {"mean_depth": 42.0})
            self.assertEqual(refreshed["manifest"]["COVERAGE"], "42")
            self.assertEqual(refreshed["sequence_sha256"], sequence_digest)
            self.assertFalse(legacy.exists())
            self.assertNotIn(
                f"{seqid}.manifest.txt", (package / "checksums.sha256").read_text()
            )

    def test_refresh_carries_an_older_package_forward(self):
        """An on-disk schema 2 package gains the new shape rather than a hybrid."""
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
                coverage=5.0,
                program="GetOrganelle 1.7.7.1",
                platform="ILLUMINA",
                scientific_name="Testus organismus",
                outdir=str(package),
            )
            MODULE.build_package(args)
            metadata_path = package / f"{seqid}.package_metadata.json"
            stale = json.loads(metadata_path.read_text())
            stale["schema_version"] = 2
            stale["study"] = stale.pop("validation_study")
            stale["biosample_accession"] = "SAMN40589646"
            stale["biosample_source"] = "sample.ncbi_biosample_id"
            stale["run_accessions"] = []
            del stale["manifest"]
            metadata_path.write_text(json.dumps(stale, indent=2, sort_keys=True))

            refreshed = MODULE.refresh_package(package)
            self.assertEqual(refreshed["schema_version"], 3)
            self.assertEqual(refreshed["validation_study"], "PRJEB123419")
            for retired in (
                "study",
                "biosample_accession",
                "biosample_source",
                "run_accessions",
            ):
                self.assertNotIn(retired, refreshed)
            self.assertEqual(refreshed["manifest"]["PLATFORM"], "ILLUMINA")


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
            coverage=100.0,
            program="MitoHiFi 3.2.3",
            platform="PACBIO_SMRT",
            scientific_name="Choerodon rubescens",
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
