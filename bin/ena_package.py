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
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
VALID_BASES = frozenset("ACGTRYSWKMBDHVN")
PLATFORMS = frozenset({"PACBIO_SMRT", "ILLUMINA"})
# Qualifiers of the EMBL source feature that a submitter would otherwise have to
# re-derive from the database. Order is the order they appear in the flatfile.
SOURCE_QUALIFIERS = (
    "organism",
    "organelle",
    "mol_type",
    "isolate",
    "tissue_type",
    "geo_loc_name",
    "collection_date",
    "lat_lon",
)


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


def source_qualifiers(embl_path: Path) -> dict[str, str]:
    """Specimen facts read back out of the packaged flatfile's source feature.

    Read from the flatfile rather than re-queried from the database so the
    package cannot disagree with itself: what is recorded here is exactly what
    was submitted.  A qualifier the flatfile does not carry is left out rather
    than emitted empty, so absent is distinguishable from blank.
    """
    opener = gzip.open if embl_path.suffix == ".gz" else open
    entries: list[str] = []
    in_source = False
    with opener(embl_path, "rt") as handle:
        for raw in handle:
            if not raw.startswith("FT"):
                continue
            body = raw[5:].rstrip("\n")
            if raw[5:6] != " ":
                # A feature key starts in column 6; source is the first feature,
                # so the next key after it ends the block.
                if in_source:
                    break
                in_source = bool(body.split()) and body.split()[0] == "source"
                continue
            if not in_source:
                continue
            stripped = body.strip()
            if stripped.startswith("/"):
                entries.append(stripped)
            elif entries:
                # A value long enough to wrap continues on the following line.
                entries[-1] = f"{entries[-1]} {stripped}"
    wanted = set(SOURCE_QUALIFIERS)
    found: dict[str, str] = {}
    for entry in entries:
        match = re.fullmatch(r'/([A-Za-z_]+)="?(.*?)"?', entry)
        if match and match.group(1) in wanted:
            found[match.group(1)] = match.group(2)
    return {name: found[name] for name in SOURCE_QUALIFIERS if name in found}


def manifest_fields(
    *,
    full_seqid: str,
    coverage: float | None,
    program: str,
    platform: str,
    flatfile_name: str,
    chromosome_list_name: str,
    scientific_name: str,
) -> dict[str, str]:
    """The Webin genome-context manifest keys this pipeline can actually fill.

    STUDY, SAMPLE and RUN_REF are deliberately absent: the study and the sample
    are registered by the submission pipeline and the runs belong to the raw-read
    submissions, so none of the three is known here.  A key whose value does not
    validate is omitted rather than emitted wrong, and never fails the build:
    the package is still worth handing over without it.
    """
    full_seqid_og_id(full_seqid)
    fields: dict[str, str] = {
        "ASSEMBLYNAME": full_seqid,
        "ASSEMBLY_TYPE": "clone or isolate",
    }
    if coverage is not None and math.isfinite(float(coverage)) and float(coverage) >= 0:
        fields["COVERAGE"] = f"{float(coverage):g}"
    if program.strip():
        fields["PROGRAM"] = program.strip()
    platforms = [value.strip() for value in platform.split(",") if value.strip()]
    if platforms and not [value for value in platforms if value not in PLATFORMS]:
        fields["PLATFORM"] = ",".join(platforms)
    fields["MOLECULETYPE"] = "genomic DNA"
    if scientific_name.strip():
        fields["DESCRIPTION"] = f"{scientific_name.strip()} mitochondrial genome"
    fields["FLATFILE"] = flatfile_name
    fields["CHROMOSOME_LIST"] = chromosome_list_name
    return fields


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
    # A package built under an earlier schema is carried forward rather than
    # left half-refreshed: the retired keys go, the study is renamed to what it
    # always was, and the specimen block is read back off the flatfile.
    if int(metadata.get("schema_version") or 0) < 3:
        if "study" in metadata:
            metadata["validation_study"] = metadata.pop("study")
        for retired in ("biosample_accession", "biosample_source", "run_accessions"):
            metadata.pop(retired, None)
        metadata["schema_version"] = 3
    if int(metadata.get("schema_version") or 0) < 4:
        # A package built before coverage had a provenance cannot have one
        # reconstructed: it predates the distinction.
        metadata.setdefault("mean_depth_source", "unknown")
        metadata["schema_version"] = 4
    if embl.exists():
        metadata["specimen"] = source_qualifiers(embl)
    metadata["manifest"] = manifest_fields(
        full_seqid=full_seqid,
        coverage=metadata.get("mean_depth"),
        program=str(metadata.get("program") or ""),
        platform=str(metadata.get("platform") or ""),
        flatfile_name=embl.name,
        chromosome_list_name=chromosomes.name,
        scientific_name=str(metadata.get("scientific_name") or ""),
    )
    # Packages built before the manifest moved into the metadata still carry a
    # standalone manifest file; leaving it would keep a stale copy in the
    # checksums of a package that has just been refreshed.
    (package_dir / f"{full_seqid}.manifest.txt").unlink(missing_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    primary_artifacts = sorted(
        path
        for path in package_dir.iterdir()
        if path.is_file()
        and (
            path.name.endswith(".embl.gz")
            or path.name.endswith(".chromosome_list.tsv.gz")
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
            ("study", "validation_study"),
            ("coverage", "mean_depth"),
            ("program", "program"),
            ("platform", "platform"),
            ("scientific_name", "scientific_name"),
        ):
            current = getattr(args, argument, None)
            if current in (None, "") and supplied_metadata.get(key) is not None:
                setattr(args, argument, supplied_metadata[key])
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
    write_gzip_copy(Path(args.embl), package_dir / embl_name)
    with (package_dir / chromosome_name).open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(chromosome_list_text(full_seqid).encode("ascii"))
    metadata = {
        # 2: dropped package_status, metadata_blocker and published_package_path;
        # added flatfile_validation.  The durable path is the caller's to know.
        # 3: the standalone <full_seqid>.manifest.txt is gone and this file is the
        # whole handoff.  "manifest" renders the Webin genome-context keys this
        # pipeline can fill, so a submitter reads them under the names Webin uses;
        # "specimen" carries the source-feature facts so the flatfile need not be
        # parsed for them.  STUDY, SAMPLE and RUN_REF are absent because the
        # study, the sample and the read submissions are registered downstream,
        # so biosample_accession and run_accessions went with them.  "study" is
        # "validation_study": it is the study sequence-context validation ran
        # against, never a submission target.
        # 4: added mean_depth_source, so a package records whether its COVERAGE
        # came from the run that built it or from a previously stored value.
        "schema_version": 4,
        "full_seqid": full_seqid,
        "og_id": og_id,
        "assembly_prefix": args.assembly_prefix,
        "annotation_version": args.annotation_version,
        "validation_study": args.study,
        "mean_depth": args.coverage,
        "mean_depth_source": str(supplied_metadata.get("mean_depth_source") or "unknown"),
        "program": args.program,
        "platform": args.platform,
        "scientific_name": args.scientific_name,
        "manifest": manifest_fields(
            full_seqid=full_seqid,
            coverage=args.coverage,
            program=args.program,
            platform=args.platform,
            flatfile_name=embl_name,
            chromosome_list_name=chromosome_name,
            scientific_name=args.scientific_name,
        ),
        "specimen": source_qualifiers(package_dir / embl_name),
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
    build.add_argument(
        "--study",
        help="Study sequence-context validation ran against; recorded as validation_study.",
    )
    build.add_argument("--coverage", type=float)
    build.add_argument("--program")
    build.add_argument("--platform")
    build.add_argument("--scientific-name")
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
