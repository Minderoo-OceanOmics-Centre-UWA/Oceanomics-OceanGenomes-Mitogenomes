#!/usr/bin/env python3
"""Build and validate self-contained ENA mitogenome candidate packages."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path


FULL_SEQID_PATTERN = re.compile(r"^(OG[0-9]+)[.][A-Za-z0-9._-]+$")
# BioSample accessions are INSDC-wide: SAMEA (EBI), SAMN (NCBI), SAMD (DDBJ).
# OceanOmics registers specimens at NCBI, so every accession the database holds
# is SAMN. INSDC sharing mirrors those records into EBI BioSamples -- the ENA
# browser resolves them -- but webin-cli resolves SAMPLE against ENA's own
# submission sample service, which only knows samples registered through Webin.
# Verified against ena-webin-cli 9.0.3 -context genome -validate -test:
#   SAMEA132129018 -> "Submission(s) validated successfully."
#   SAMN40589646   -> "Failed to initialise validator ... sample is null"
# So a SAMN is a real registration that ENA cannot reference yet, which is a
# different situation from a specimen with no BioSample and from a malformed
# value. Keep the three apart so the blocked report says which action is needed.
INSDC_BIOSAMPLE_PATTERN = re.compile(r"SAM(?:EA|N|D)[0-9]+")
ENA_BIOSAMPLE_PATTERN = re.compile(r"SAMEA[0-9]+")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
VALID_BASES = frozenset("ACGTRYSWKMBDHVN")
PLATFORMS = frozenset({"PACBIO_SMRT", "ILLUMINA"})


def full_seqid_og_id(full_seqid: str) -> str:
    match = FULL_SEQID_PATTERN.fullmatch(full_seqid.strip())
    if not match:
        raise ValueError(f"Invalid full OceanOmics SeqID: {full_seqid!r}")
    return match.group(1)


def flatfile_validation(status_path: str | None) -> dict[str, object]:
    """The sequence-context Webin verdict, as recorded in the package.

    This is the pipeline's last ENA gate, so the package carries its own result
    rather than making a reader join the status TSV back on by file name.  A
    missing file means validation was switched off, not that it failed.
    """
    absent = {
        "status": "NOT_REQUESTED",
        "reason": "not_requested",
        "error_count": None,
        "warning_count": None,
        "webin_cli_version": None,
    }
    if not status_path:
        return absent
    try:
        with Path(status_path).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return absent
    if not rows:
        return absent
    row = {str(k): (v or "").strip() for k, v in rows[0].items() if k}

    def count(name: str) -> int | None:
        try:
            return int(row.get(name, ""))
        except ValueError:
            return None

    return {
        "status": row.get("status") or "MALFORMED_STATUS",
        "reason": row.get("reason") or "unspecified",
        "error_count": count("error_count"),
        "warning_count": count("warning_count"),
        "webin_cli_version": row.get("webin_cli_version") or None,
    }


def read_single_fasta(path: Path) -> tuple[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    bases: list[str] = []
    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(bases).upper()))
                name = line[1:].split()[0]
                bases = []
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence appears before header in {path}")
                bases.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(bases).upper()))
    if len(records) != 1:
        raise ValueError(f"Expected one mitochondrial FASTA record in {path}; found {len(records)}")
    invalid = sorted(set(records[0][1]) - VALID_BASES)
    if invalid:
        raise ValueError(f"Invalid IUPAC bases in {path}: {''.join(invalid)}")
    if not records[0][1]:
        raise ValueError(f"Empty mitochondrial sequence in {path}")
    return records[0]


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return sequence.translate(table)[::-1]


def minimal_rotation(sequence: str) -> str:
    """Return a lexicographically minimal rotation using Booth's algorithm."""
    if not sequence:
        return sequence
    doubled = sequence + sequence
    length = len(sequence)
    i, j, offset = 0, 1, 0
    while i < length and j < length and offset < length:
        left, right = doubled[i + offset], doubled[j + offset]
        if left == right:
            offset += 1
            continue
        if left > right:
            i = i + offset + 1
            if i <= j:
                i = j + 1
        else:
            j = j + offset + 1
            if j <= i:
                j = i + 1
        offset = 0
    start = min(i, j)
    return doubled[start : start + length]


def sequence_digests(sequence: str) -> tuple[str, str]:
    ordinary = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    forward = minimal_rotation(sequence)
    reverse = minimal_rotation(reverse_complement(sequence))
    circular = hashlib.sha256(min(forward, reverse).encode("ascii")).hexdigest()
    return ordinary, circular


def write_gzip_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as reader, destination.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as writer:
                shutil.copyfileobj(reader, writer)
    else:
        with source.open("rb") as reader, destination.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as writer:
                shutil.copyfileobj(reader, writer)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chromosome_list_text(full_seqid: str) -> str:
    full_seqid_og_id(full_seqid)
    return f"{full_seqid}\tMT\tCircular-Chromosome\tMitochondrion\n"


def manifest_text(
    *,
    study: str,
    biosample: str,
    full_seqid: str,
    coverage: float | None,
    program: str,
    platform: str,
    flatfile_name: str,
    chromosome_list_name: str,
    scientific_name: str,
    run_accessions: list[str] | None = None,
) -> str:
    full_seqid_og_id(full_seqid)
    if not study.strip():
        raise ValueError("STUDY is required")
    accession = biosample.strip()
    if not ENA_BIOSAMPLE_PATTERN.fullmatch(accession):
        if INSDC_BIOSAMPLE_PATTERN.fullmatch(accession):
            raise ValueError(
                f"BioSample {accession} is registered outside ENA and cannot be "
                "referenced by a Webin submission; register the specimen in ENA "
                "or ask ENA to broker the existing accession"
            )
        raise ValueError(f"Invalid or missing ENA BioSample accession: {biosample!r}")
    if coverage is None or not math.isfinite(float(coverage)) or float(coverage) < 0:
        raise ValueError(f"COVERAGE must be a finite non-negative mean depth: {coverage}")
    if not program.strip():
        raise ValueError("PROGRAM is required")
    platforms = [value.strip() for value in platform.split(",") if value.strip()]
    invalid_platforms = [value for value in platforms if value not in PLATFORMS]
    if not platforms or invalid_platforms:
        raise ValueError(f"Invalid PLATFORM value(s): {platform!r}")
    fields = [
        ("STUDY", study.strip()),
        ("SAMPLE", accession),
        ("ASSEMBLYNAME", full_seqid),
        ("ASSEMBLY_TYPE", "clone or isolate"),
        ("COVERAGE", f"{float(coverage):g}"),
        ("PROGRAM", program.strip()),
        ("PLATFORM", ",".join(platforms)),
        ("MOLECULETYPE", "genomic DNA"),
        ("DESCRIPTION", f"{scientific_name.strip()} mitochondrial genome"),
    ]
    runs = sorted(set(run_accessions or []))
    if runs:
        fields.append(("RUN_REF", ",".join(runs)))
    fields.extend(
        [
            ("FLATFILE", flatfile_name),
            ("CHROMOSOME_LIST", chromosome_list_name),
        ]
    )
    return "".join(f"{key}\t{value}\n" for key, value in fields)


def refresh_package(
    package_dir: Path, metadata_updates: dict[str, object] | None = None
) -> dict[str, object]:
    """Refresh mutable manifest metadata without rebuilding sequence/annotation."""
    metadata_paths = sorted(package_dir.glob("*.package_metadata.json"))
    if len(metadata_paths) != 1:
        raise ValueError(
            f"Expected one package metadata JSON in {package_dir}; found {len(metadata_paths)}"
        )
    metadata_path = metadata_paths[0]
    metadata = json.loads(metadata_path.read_text())
    metadata.update(
        {
            key: value
            for key, value in (metadata_updates or {}).items()
            if value is not None
        }
    )
    full_seqid = str(metadata["full_seqid"])
    embl = package_dir / f"{full_seqid}.embl.gz"
    chromosomes = package_dir / f"{full_seqid}.chromosome_list.tsv.gz"
    manifest = package_dir / f"{full_seqid}.manifest.txt"
    blocker = ""
    try:
        text = manifest_text(
            study=str(metadata.get("study") or ""),
            biosample=str(metadata.get("biosample_accession") or ""),
            full_seqid=full_seqid,
            coverage=metadata.get("mean_depth"),
            program=str(metadata.get("program") or ""),
            platform=str(metadata.get("platform") or ""),
            flatfile_name=embl.name,
            chromosome_list_name=chromosomes.name,
            scientific_name=str(metadata.get("scientific_name") or ""),
            run_accessions=list(metadata.get("run_accessions") or []),
        )
    except ValueError as error:
        blocker = str(error)
        text = (
            f"# BLOCKED: {blocker}\n"
            f"STUDY\t{metadata.get('study') or ''}\n"
            f"ASSEMBLYNAME\t{full_seqid}\n"
            f"FLATFILE\t{embl.name}\n"
            f"CHROMOSOME_LIST\t{chromosomes.name}\n"
        )
    manifest.write_text(text)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    primary_artifacts = sorted(
        path
        for path in package_dir.iterdir()
        if path.is_file()
        and (
            path.name.endswith(".embl.gz")
            or path.name.endswith(".chromosome_list.tsv.gz")
            or path.name.endswith(".manifest.txt")
            or path.name.endswith(".tbl")
        )
    )
    metadata["package_digest"] = hashlib.sha256(
        "".join(
            f"{sha256_file(path)}  {path.name}\n" for path in primary_artifacts
        ).encode()
    ).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    checksummed = sorted(
        path
        for path in package_dir.iterdir()
        if path.is_file() and path.name != "checksums.sha256"
    )
    (package_dir / "checksums.sha256").write_text(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in checksummed)
    )
    return metadata


def build_package(args: argparse.Namespace) -> int:
    supplied_metadata: dict[str, object] = {}
    if getattr(args, "metadata_input", None):
        supplied_metadata = json.loads(Path(args.metadata_input).read_text())
        for argument, key in (
            ("og_id", "og_id"),
            ("assembly_prefix", "assembly_prefix"),
            ("annotation_version", "annotation_version"),
            ("full_seqid", "full_seqid"),
            ("study", "study"),
            ("biosample", "biosample_accession"),
            ("coverage", "mean_depth"),
            ("program", "program"),
            ("platform", "platform"),
            ("scientific_name", "scientific_name"),
        ):
            current = getattr(args, argument, None)
            if current in (None, "") and supplied_metadata.get(key) is not None:
                setattr(args, argument, supplied_metadata[key])
        if not args.run_accession:
            args.run_accession = list(supplied_metadata.get("run_accessions") or [])
    required_arguments = (
        "og_id",
        "assembly_prefix",
        "annotation_version",
        "full_seqid",
        "study",
        "program",
        "platform",
        "scientific_name",
    )
    missing = [name for name in required_arguments if getattr(args, name, None) in (None, "")]
    if missing:
        raise ValueError(f"Missing package metadata: {', '.join(missing)}")
    full_seqid = args.full_seqid.strip()
    og_id = full_seqid_og_id(full_seqid)
    if args.og_id != og_id:
        raise ValueError(f"--og-id {args.og_id} does not match --full-seqid {full_seqid}")
    _, sequence = read_single_fasta(Path(args.fasta))
    sequence_sha, circular_sha = sequence_digests(sequence)
    package_dir = Path(args.outdir)
    package_dir.mkdir(parents=True, exist_ok=True)
    embl_name = f"{full_seqid}.embl.gz"
    chromosome_name = f"{full_seqid}.chromosome_list.tsv.gz"
    manifest_name = f"{full_seqid}.manifest.txt"
    write_gzip_copy(Path(args.embl), package_dir / embl_name)
    with (package_dir / chromosome_name).open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(chromosome_list_text(full_seqid).encode("ascii"))
    manifest_error = ""
    try:
        manifest = manifest_text(
            study=args.study,
            biosample=args.biosample or "",
            full_seqid=full_seqid,
            coverage=args.coverage,
            program=args.program,
            platform=args.platform,
            flatfile_name=embl_name,
            chromosome_list_name=chromosome_name,
            scientific_name=args.scientific_name,
            run_accessions=args.run_accession,
        )
    except ValueError as error:
        manifest_error = str(error)
        manifest = (
            f"# BLOCKED: {manifest_error}\n"
            f"STUDY\t{args.study}\n"
            f"ASSEMBLYNAME\t{full_seqid}\n"
            f"FLATFILE\t{embl_name}\n"
            f"CHROMOSOME_LIST\t{chromosome_name}\n"
        )
    (package_dir / manifest_name).write_text(manifest)
    metadata = {
        # 2: dropped package_status, metadata_blocker and published_package_path;
        # added flatfile_validation.  Readiness is derived from
        # biosample_accession and the durable path is the caller's to know.
        "schema_version": 2,
        "full_seqid": full_seqid,
        "og_id": og_id,
        "assembly_prefix": args.assembly_prefix,
        "annotation_version": args.annotation_version,
        "study": args.study,
        "biosample_accession": args.biosample or None,
        "mean_depth": args.coverage,
        "program": args.program,
        "platform": args.platform,
        "scientific_name": args.scientific_name,
        "run_accessions": sorted(set(args.run_accession or [])),
        "sequence_length": len(sequence),
        "sequence_sha256": sequence_sha,
        "normalised_circular_sha256": circular_sha,
        "flatfile_validation": flatfile_validation(getattr(args, "flatfile_status", None)),
    }
    metadata_path = package_dir / f"{full_seqid}.package_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if getattr(args, "tbl", None):
        shutil.copy2(args.tbl, package_dir / f"{full_seqid}.tbl")
    # Collaborators get sequence plus annotation, not an ENA submission bundle.
    # Both live here so the package directory is the single thing to hand over.
    # The GFF is copied through as the annotator produced it: locus tags are
    # assigned downstream, so there is nothing here to key it on.
    packaged_fasta = package_dir / f"{full_seqid}.fa"
    if Path(args.fasta).resolve() != packaged_fasta.resolve():
        shutil.copy2(args.fasta, packaged_fasta)
    if getattr(args, "gff", None):
        packaged_gff = package_dir / f"{full_seqid}.gff"
        if Path(args.gff).resolve() != packaged_gff.resolve():
            shutil.copy2(args.gff, packaged_gff)
    # Extracted gene sequences arrive named on the assembly prefix; renaming them
    # onto full_seqid keeps every file in the package sharing one stem.
    if getattr(args, "genes", None):
        packaged_genes = package_dir / f"{full_seqid}.genes.fa"
        if Path(args.genes).resolve() != packaged_genes.resolve():
            shutil.copy2(args.genes, packaged_genes)
    primary_artifacts = sorted(
        path for path in package_dir.iterdir()
        if path.is_file()
        and (
            path.name.endswith(".embl.gz")
            or path.name.endswith(".chromosome_list.tsv.gz")
            or path.name.endswith(".manifest.txt")
            or path.name.endswith(".tbl")
        )
    )
    package_digest = hashlib.sha256(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in primary_artifacts).encode()
    ).hexdigest()
    metadata["package_digest"] = package_digest
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    checksummed = sorted(
        path for path in package_dir.iterdir()
        if path.is_file() and path.name != "checksums.sha256"
    )
    checksum_path = package_dir / "checksums.sha256"
    checksum_path.write_text(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in checksummed)
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build", help="Build one ENA candidate package")
    build.add_argument("--metadata-input")
    build.add_argument("--og-id")
    build.add_argument("--assembly-prefix")
    build.add_argument("--annotation-version")
    build.add_argument("--full-seqid")
    build.add_argument("--fasta", required=True)
    build.add_argument("--embl", required=True)
    build.add_argument("--study")
    build.add_argument("--biosample", default="")
    build.add_argument("--coverage", type=float)
    build.add_argument("--program")
    build.add_argument("--platform")
    build.add_argument("--scientific-name")
    build.add_argument("--run-accession", action="append", default=[])
    build.add_argument("--tbl")
    build.add_argument(
        "--gff",
        help="Processed annotation GFF; packaged verbatim for collaborators.",
    )
    build.add_argument(
        "--genes",
        help="Concatenated gene FASTA; packaged for collaborators as <full_seqid>.genes.fa.",
    )
    build.add_argument("--outdir", required=True)
    build.add_argument(
        "--flatfile-status",
        help="WEBIN_VALIDATE status TSV; its verdict is recorded in the package metadata.",
    )
    build.set_defaults(func=build_package)
    args = parser.parse_args()
    try:
        return args.func(args)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
