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
    def test_og_derived_locus_tags(self):
        self.assertEqual(MODULE.locus_tag("OG910", 1, "OGMTHIFI"), "OGMTHIFI_000910001")
        self.assertEqual(MODULE.locus_tag("OG5", 12, "OGMTILMN"), "OGMTILMN_000005012")
        # Same specimen serial, one tag per technology study.
        self.assertNotEqual(
            MODULE.locus_tag("OG910", 1, "OGMTHIFI"),
            MODULE.locus_tag("OG910", 1, "OGMTHIC"),
        )
        with self.assertRaises(ValueError):
            MODULE.locus_tag("OG1000000", 1, "OGMTHIFI")
        with self.assertRaises(ValueError):
            MODULE.locus_tag("OG910", 1000, "OGMTHIFI")
        with self.assertRaises(TypeError):
            MODULE.locus_tag("OG910", 1)

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

    def test_build_and_validate_package(self):
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
            self.assertEqual(MODULE.validate_package(package), [])
            with gzip.open(package / f"{seqid}.chromosome_list.tsv.gz", "rt") as handle:
                self.assertEqual(
                    handle.read(),
                    f"{seqid}\tMT\tCircular-Chromosome\tMitochondrion\n",
                )
            metadata = json.loads(
                (package / f"{seqid}.package_metadata.json").read_text()
            )
            self.assertEqual(metadata["package_status"], "READY")
            self.assertEqual(metadata["local_validation_status"], "PASS")
            # No --publish-root supplied, so nothing is claimed about where the
            # package will end up.
            self.assertIsNone(metadata["published_package_path"])

    def test_publish_root_is_recorded_and_survives_refresh(self):
        """The durable path must outlive the work directory the build ran in."""
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
            package = root / "package"
            published = "/published/mitogenomes/OG910/OG910.hifi.241127.v3mitohifi/ena/package"
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
                publish_root=f"{published}/",
            )
            self.assertEqual(MODULE.build_package(args), 0)
            metadata_path = package / f"{seqid}.package_metadata.json"
            metadata = json.loads(metadata_path.read_text())
            # Trailing separator normalised away so string comparison against a
            # stored path cannot fail on cosmetics.
            self.assertEqual(metadata["published_package_path"], published)
            refreshed = MODULE.refresh_package(package, {"mean_depth": 42.0})
            self.assertEqual(refreshed["published_package_path"], published)
            self.assertEqual(
                json.loads(metadata_path.read_text())["published_package_path"],
                published,
            )

    def test_publish_root_must_be_absolute(self):
        self.assertIsNone(MODULE.normalise_publish_root(None))
        self.assertIsNone(MODULE.normalise_publish_root("  "))
        self.assertEqual(MODULE.normalise_publish_root("/a/b/"), "/a/b")
        # A relative outdir resolves against whoever launched the run, which is
        # exactly the path that looks fine now and is unopenable later.
        with self.assertRaises(ValueError):
            MODULE.normalise_publish_root("results/ena/package")

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
            self.assertEqual(metadata["package_status"], "WAITING_FOR_BIOSAMPLE")
            sequence_digest = metadata["sequence_sha256"]
            refreshed = MODULE.refresh_package(
                package, {"biosample_accession": "SAMEA123"}
            )
            self.assertEqual(refreshed["package_status"], "READY")
            self.assertEqual(refreshed["local_validation_status"], "PASS")
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

    def test_biosample_block_status_separates_the_three_cases(self):
        self.assertEqual(
            MODULE.biosample_block_status(None), "WAITING_FOR_BIOSAMPLE"
        )
        self.assertEqual(MODULE.biosample_block_status("  "), "WAITING_FOR_BIOSAMPLE")
        self.assertEqual(
            MODULE.biosample_block_status("SAMN40589646"),
            "BLOCKED_NCBI_ONLY_BIOSAMPLE",
        )
        self.assertEqual(
            MODULE.biosample_block_status("SRS123456"), "BLOCKED_METADATA"
        )

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
            self.assertEqual(
                metadata["package_status"], "BLOCKED_NCBI_ONLY_BIOSAMPLE"
            )
            # A blocked manifest must never carry a SAMPLE line webin would reject.
            manifest = (package / f"{seqid}.manifest.txt").read_text()
            self.assertTrue(manifest.startswith("# BLOCKED:"))
            self.assertNotIn("SAMPLE\t", manifest)

            # Brokering the specimen into ENA unblocks it in place.
            refreshed = MODULE.refresh_package(
                package, {"biosample_accession": "SAMEA132129018"}
            )
            self.assertEqual(refreshed["package_status"], "READY")
            self.assertIn(
                "SAMPLE\tSAMEA132129018",
                (package / f"{seqid}.manifest.txt").read_text(),
            )


MAPPING_HEADER = (
    "full_seqid\tfeature_key\tfeature_type\tcanonical_gene\tgene_occurrence"
    "\tstart\tend\tstrand\tlocus_tag\n"
)


def write_mapping(path, seqid, rows):
    """Build a locus-tag mapping TSV in the shape allocate_ena_locus_tags emits."""
    lines = [MAPPING_HEADER]
    for index, (feature_type, gene, occurrence, start, end, strand, tag) in enumerate(rows, 1):
        lines.append(
            f"{seqid}\t{feature_type}:{index}\t{feature_type}\t{gene}\t{occurrence}"
            f"\t{start}\t{end}\t{strand}\t{tag}\n"
        )
    path.write_text("".join(lines))


def attributes_of(line):
    return dict(
        item.split("=", 1) for item in line.split("\t")[8].split(";") if "=" in item
    )


class TagGffTests(unittest.TestCase):
    seqid = "OG910.hifi.241127.v3mitohifi.emma102"

    def test_coordinates_join_and_ids_are_rebuilt_from_locus_tags(self):
        """A collaborator's GFF must key on the tags the record was published under."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"##sequence-region\t{self.seqid}\t1\t16777\n"
                "##organism Choerodon rubescens\n"
                "##isolate OG910\n"
                f"{self.seqid}\tEmma\tregion\t1\t16777\t.\t+\t0\tIs_circular=true\n"
                f"{self.seqid}\tEmma\tgene\t1\t68\t.\t+\t.\tID=uuid-gene-tf;Name=TF\n"
                f"{self.seqid}\tEmma\ttRNA\t1\t68\t3.0e-20\t+\t.\t"
                "ID=uuid-trna-tf;Parent=uuid-gene-tf;Name=TF;Product=tRNA-Phe(GAA)\n"
                f"{self.seqid}\tEmma\tgene\t3920\t3990\t.\t-\t.\tID=uuid-gene-tq;Name=TQ\n"
                f"{self.seqid}\tEmma\ttRNA\t3920\t3990\t1.4e-24\t-\t.\t"
                "ID=uuid-trna-tq;Parent=uuid-gene-tq;Name=TQ;Product=tRNA-Gln(UUG)\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping,
                self.seqid,
                [
                    ("gene", "TF", 1, 1, 68, "+", "OGMTHIFI_000910001"),
                    ("tRNA", "TF", 1, 1, 68, "+", "OGMTHIFI_000910001"),
                    # Minus strand arrives from the .tbl with start > end.
                    ("gene", "TQ", 1, 3990, 3920, "-", "OGMTHIFI_000910008"),
                    ("tRNA", "TQ", 1, 3990, 3920, "-", "OGMTHIFI_000910008"),
                ],
            )
            out = root / "out.gff"
            self.assertEqual(MODULE.tag_gff(gff, mapping, out), (4, 0))
            lines = out.read_text().splitlines()

            # Directives and the whole-molecule region record survive untouched.
            self.assertEqual(lines[:5], gff.read_text().splitlines()[:5])
            self.assertNotIn("uuid-", out.read_text())

            gene_tf, trna_tf, gene_tq, trna_tq = lines[5:9]
            self.assertEqual(
                attributes_of(gene_tf),
                {"ID": "gene-OGMTHIFI_000910001", "Name": "TF", "locus_tag": "OGMTHIFI_000910001"},
            )
            self.assertEqual(attributes_of(trna_tf)["ID"], "trna-OGMTHIFI_000910001")
            self.assertEqual(attributes_of(trna_tf)["Parent"], "gene-OGMTHIFI_000910001")
            # Name and Product keep their positions; locus_tag is appended.
            self.assertTrue(trna_tf.endswith(";Product=tRNA-Phe(GAA);locus_tag=OGMTHIFI_000910001"))
            self.assertEqual(attributes_of(gene_tq)["locus_tag"], "OGMTHIFI_000910008")
            self.assertEqual(attributes_of(trna_tq)["Parent"], "gene-OGMTHIFI_000910008")

    def test_cds_chains_through_its_mrna_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tgene\t2873\t3847\t.\t+\t.\tID=g;Name=ND1\n"
                f"{self.seqid}\tEmma\tmRNA\t2873\t3847\t.\t+\t.\tID=m;Parent=g;Name=ND1\n"
                f"{self.seqid}\tEmma\tCDS\t2873\t3847\t1.6e-191\t+\t0\tID=c;Parent=m;Name=ND1\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping,
                self.seqid,
                [
                    ("gene", "ND1", 1, 2873, 3847, "+", "OGMTHIFI_000910006"),
                    ("mRNA", "ND1", 1, 2873, 3847, "+", "OGMTHIFI_000910006"),
                    ("CDS", "ND1", 1, 2873, 3847, "+", "OGMTHIFI_000910006"),
                ],
            )
            out = root / "out.gff"
            self.assertEqual(MODULE.tag_gff(gff, mapping, out), (3, 0))
            cds = out.read_text().splitlines()[3]
            self.assertEqual(attributes_of(cds)["ID"], "cds-OGMTHIFI_000910006")
            self.assertEqual(attributes_of(cds)["Parent"], "mrna-OGMTHIFI_000910006")

    def test_multi_exon_records_stay_unique_under_one_tag(self):
        """The .tbl records one CDS block per gene; the GFF records one per exon."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tmitos\tgene\t100\t400\t.\t+\t.\tID=gene-ND5;Name=MT-ND5\n"
                f"{self.seqid}\tmitos\tCDS\t100\t200\t.\t+\t0\tID=ND5.1;Parent=gene-ND5;Name=MT-ND5\n"
                f"{self.seqid}\tmitos\tCDS\t300\t400\t.\t+\t0\tID=ND5.2;Parent=gene-ND5;Name=MT-ND5\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping,
                self.seqid,
                [
                    ("gene", "ND5", 1, 100, 400, "+", "OGMTHIFI_000910020"),
                    # Only the first interval reaches the mapping from the .tbl.
                    ("CDS", "ND5", 1, 100, 200, "+", "OGMTHIFI_000910020"),
                ],
            )
            out = root / "out.gff"
            self.assertEqual(MODULE.tag_gff(gff, mapping, out), (3, 0))
            first, second = out.read_text().splitlines()[2:4]
            self.assertEqual(attributes_of(first)["ID"], "cds-OGMTHIFI_000910020")
            self.assertEqual(attributes_of(second)["ID"], "cds-OGMTHIFI_000910020.2")
            self.assertEqual(attributes_of(second)["Parent"], "gene-OGMTHIFI_000910020")

    def test_gene_identity_is_the_fallback_when_coordinates_disagree(self):
        """MITOS clamps an origin-spanning feature's GFF end to the sequence length.

        bin/mitos_to_emma.py::_gff_span reports the wrapped end in a Note instead,
        so its coordinates cannot equal the .tbl the mapping came from.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tmitos\tgene\t16500\t16777\t.\t+\t.\t"
                "ID=gene-COX1;Name=MT-COX1;Note=origin-spanning;wrapped_end=120\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping,
                self.seqid,
                [("gene", "COX1", 1, 16500, 120, "+", "OGMTHIFI_000910030")],
            )
            out = root / "out.gff"
            self.assertEqual(MODULE.tag_gff(gff, mapping, out), (1, 0))
            gene = out.read_text().splitlines()[1]
            self.assertEqual(attributes_of(gene)["ID"], "gene-OGMTHIFI_000910030")
            self.assertEqual(attributes_of(gene)["locus_tag"], "OGMTHIFI_000910030")
            # The wrapped-end Note is the only record of the real coordinate.
            self.assertIn("wrapped_end=120", gene)

    def test_untagged_features_are_counted_not_raised_on(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tregion\t1\t16777\t.\t+\t0\tIs_circular=true\n"
                f"{self.seqid}\tEmma\tcontrol_region\t15000\t16777\t.\t+\t.\tID=uuid-cr\n"
            )
            mapping = root / "map.tsv"
            write_mapping(mapping, self.seqid, [])
            out = root / "out.gff"
            # region is structural and never counts; control_region is genuinely untagged.
            self.assertEqual(MODULE.tag_gff(gff, mapping, out), (0, 1))
            self.assertIn("ID=uuid-cr", out.read_text())


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

    def test_package_carries_the_fasta_and_a_locus_tagged_gff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tgene\t1\t68\t.\t+\t.\tID=uuid-gene-tf;Name=TF\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping, self.seqid, [("gene", "TF", 1, 1, 68, "+", "OGMTHIFI_000910001")]
            )
            package = root / "package"
            args = self.build_args(root, package, gff=str(gff), locus_map=str(mapping))
            self.assertEqual(MODULE.build_package(args), 0)
            self.assertEqual(MODULE.validate_package(package), [])

            self.assertEqual(
                (package / f"{self.seqid}.fa").read_text(), f">{self.seqid}\nAACCGT\n"
            )
            packaged_gff = (package / f"{self.seqid}.gff").read_text()
            self.assertIn("locus_tag=OGMTHIFI_000910001", packaged_gff)
            self.assertNotIn("uuid-", packaged_gff)

            metadata = json.loads((package / f"{self.seqid}.package_metadata.json").read_text())
            self.assertEqual(
                metadata["gff_locus_tag_coverage"], {"tagged": 1, "untagged": 0}
            )
            # Both new files are covered by the package checksums.
            checksums = (package / "checksums.sha256").read_text()
            self.assertIn(f"{self.seqid}.fa\n", checksums)
            self.assertIn(f"{self.seqid}.gff\n", checksums)

    def test_collaborator_files_do_not_move_the_package_digest(self):
        """package_digest means 'the ENA submission content changed'.

        It is persisted to ena_candidate_packages and drives reselection in
        bin/select_ena_submission.py, so a convenience artefact must not shift it.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "in.gff"
            gff.write_text(
                "##gff-version 3\n"
                f"{self.seqid}\tEmma\tgene\t1\t68\t.\t+\t.\tID=uuid-gene-tf;Name=TF\n"
            )
            mapping = root / "map.tsv"
            write_mapping(
                mapping, self.seqid, [("gene", "TF", 1, 1, 68, "+", "OGMTHIFI_000910001")]
            )
            digests = []
            for name, extra in (
                ("without", {"locus_map": str(mapping)}),
                ("with", {"locus_map": str(mapping), "gff": str(gff)}),
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
