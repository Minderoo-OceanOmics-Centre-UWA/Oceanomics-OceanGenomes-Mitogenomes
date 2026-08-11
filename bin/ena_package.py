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


OG_PATTERN = re.compile(r"^OG([0-9]+)$")
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


def og_numeric(og_id: str) -> int:
    match = OG_PATTERN.fullmatch(og_id.strip())
    if not match:
        raise ValueError(f"Invalid OceanOmics specimen identifier: {og_id!r}")
    value = int(match.group(1))
    if value > 999999:
        raise ValueError(f"OG number exceeds six-digit locus-tag capacity: {og_id}")
    return value


def locus_tag(og_id: str, gene_serial: int, prefix: str) -> str:
    """Render a tag under the prefix registered to this candidate's ENA study.

    There is no default prefix: one prefix belongs to one study, so a fallback
    would tag a candidate with another study's namespace.
    """
    if not re.fullmatch(r"[A-Z][A-Z0-9]{2,11}", prefix):
        raise ValueError(f"Invalid INSDC locus-tag prefix: {prefix!r}")
    if not 1 <= int(gene_serial) <= 999:
        raise ValueError(f"Gene serial must be between 1 and 999: {gene_serial}")
    return f"{prefix}_{og_numeric(og_id):06d}{int(gene_serial):03d}"


def normalise_gene_name(value: str) -> str:
    """Canonicalise a gene name the way allocate_ena_locus_tags does.

    Kept byte-compatible with normalize_gene() there: the mapping TSV's
    canonical_gene column is produced by that function, so the GFF's Name= has
    to be folded the same way for the name fallback below to line up.
    """
    return re.sub(r"\s+", "", value.strip().replace("MT-", "")).upper()


def gff_attributes(field: str) -> list[tuple[str, str | None]]:
    """Split a GFF3 attribute column into ordered key/value pairs.

    Order is preserved and rebuilt verbatim so a record that gains nothing from
    tagging comes back out byte-identical.  A value of None means the attribute
    carried no '=' at all, which is not the same as an empty value and must not
    grow one on the way back out.
    """
    pairs: list[tuple[str, str | None]] = []
    for item in field.split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        pairs.append((key, value) if sep else (item, None))
    return pairs


def render_gff_attributes(pairs: list[tuple[str, str | None]]) -> str:
    return ";".join(key if value is None else f"{key}={value}" for key, value in pairs)


def read_locus_map(path: Path) -> tuple[dict, dict]:
    """Index a locus-tag mapping TSV for joining onto a GFF.

    Two indices, because the two files agree on coordinates in every case seen
    on disk but cannot be relied on to: the MITOS2 fallback clamps an
    origin-spanning feature's GFF end to the sequence length (see _gff_span in
    bin/mitos_to_emma.py) while the .tbl this mapping came from keeps the real
    wrapped end.  Coordinates are the primary key; gene identity is the backstop.

    The mapping is .tbl-derived, so a minus-strand feature arrives as
    start > end.  Both indices store the span low-to-high to match GFF columns
    4 and 5, which are always ascending with the strand in column 7.
    """
    by_coordinate: dict[tuple[str, int, int, str], str] = {}
    by_name: dict[tuple[str, str, int], str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tag = (row.get("locus_tag") or "").strip()
            if not tag:
                continue
            feature_type = row["feature_type"]
            start, end = int(row["start"]), int(row["end"])
            span = (min(start, end), max(start, end))
            by_coordinate[(feature_type, span[0], span[1], row["strand"])] = tag
            name_key = (
                feature_type,
                normalise_gene_name(row.get("canonical_gene") or ""),
                int(row.get("gene_occurrence") or 1),
            )
            by_name.setdefault(name_key, tag)
    return by_coordinate, by_name


def tag_gff(gff_path: Path, locus_map_path: Path, out_path: Path) -> tuple[int, int]:
    """Write a collaborator-facing GFF keyed on locus tags.

    The processed GFF identifies features by the annotator's own IDs -- Emma
    emits UUIDs -- which are meaningless next to a published INSDC record.  This
    rewrites ID=/Parent= from the locus tags the record was submitted under, so
    a collaborator can line the GFF up against the entry at ENA, and adds an
    explicit locus_tag= attribute for tools that read it directly.

    Records that carry no locus tag (the whole-molecule region, control_region,
    D-loop) pass through untouched, ID included: they have no tag to be renamed
    after, and rewriting them would break Parent references that point at them.

    Returns (tagged, untagged) record counts.  An untagged taggable record is
    reported rather than raised on: this file is a convenience artefact and must
    never be able to fail an ENA candidate build.
    """
    by_coordinate, by_name = read_locus_map(locus_map_path)
    lines = gff_path.read_text().splitlines()

    # Pass 1: resolve a tag and a new ID for every record, and remember the
    # old-to-new ID mapping so Parent= can be repointed in pass 2.
    resolved: dict[int, tuple[str, str]] = {}
    id_rewrites: dict[str, str] = {}
    name_occurrence: dict[tuple[str, str], int] = {}
    id_counts: dict[str, int] = {}
    for index, line in enumerate(lines):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] == "region":
            continue
        feature_type, start, end, strand = fields[2], fields[3], fields[4], fields[6]
        attributes = dict(gff_attributes(fields[8]))
        gene_name = normalise_gene_name(attributes.get("Name") or "")
        try:
            span = (int(start), int(end))
        except ValueError:
            continue
        tag = by_coordinate.get((feature_type, min(span), max(span), strand))
        if tag is None and gene_name:
            occurrence_key = (feature_type, gene_name)
            name_occurrence[occurrence_key] = name_occurrence.get(occurrence_key, 0) + 1
            tag = by_name.get((feature_type, gene_name, name_occurrence[occurrence_key]))
        if tag is None:
            continue
        stem = f"gene-{tag}" if feature_type == "gene" else f"{feature_type.lower()}-{tag}"
        id_counts[stem] = id_counts.get(stem, 0) + 1
        # Multi-exon genes emit several records of one type under one tag; suffix
        # from the second so every ID in the file stays unique.
        new_id = stem if id_counts[stem] == 1 else f"{stem}.{id_counts[stem]}"
        resolved[index] = (tag, new_id)
        old_id = attributes.get("ID") or ""
        if old_id:
            id_rewrites[old_id] = new_id

    # Pass 2: rewrite in place.
    tagged = untagged = 0
    output: list[str] = []
    for index, line in enumerate(lines):
        fields = line.split("\t")
        if line.startswith("#") or len(fields) != 9:
            output.append(line)
            continue
        if index not in resolved:
            if fields[2] != "region":
                untagged += 1
            output.append(line)
            continue
        tag, new_id = resolved[index]
        pairs = []
        for key, value in gff_attributes(fields[8]):
            if key == "ID":
                value = new_id
            elif key == "Parent":
                value = id_rewrites.get(value, value)
            elif key == "locus_tag":
                continue
            pairs.append((key, value))
        pairs.append(("locus_tag", tag))
        fields[8] = render_gff_attributes(pairs)
        output.append("\t".join(fields))
        tagged += 1
    out_path.write_text("\n".join(output) + "\n")
    return tagged, untagged


def full_seqid_og_id(full_seqid: str) -> str:
    match = FULL_SEQID_PATTERN.fullmatch(full_seqid.strip())
    if not match:
        raise ValueError(f"Invalid full OceanOmics SeqID: {full_seqid!r}")
    return match.group(1)


def normalise_publish_root(publish_root: str | None) -> str | None:
    """Where this package will live once Nextflow has published it.

    Nextflow stages package files into a task work directory, so nothing the
    build sees on disk survives a cleanup.  The publish target is known up
    front from publishDir, so record it here and let downstream consumers
    (selection, and ultimately the ENA submitter) resolve a durable path.

    Requiring it to be absolute is deliberate: a relative outdir would resolve
    against whatever launched the run, which is exactly the kind of path that
    looks fine in the database and cannot be opened weeks later.
    """
    root = (publish_root or "").strip().rstrip("/")
    if not root:
        return None
    if not root.startswith("/"):
        raise ValueError(f"--publish-root must be an absolute path; got {publish_root!r}")
    return root


def biosample_block_status(biosample: str | None) -> str:
    """Classify why a BioSample cannot be put on a manifest.

    The three cases need different follow-up: nothing registered anywhere,
    registered at another INSDC archive and awaiting an ENA-referenceable
    accession, or a value that is not an accession at all.
    """
    accession = (biosample or "").strip()
    if not accession:
        return "WAITING_FOR_BIOSAMPLE"
    if INSDC_BIOSAMPLE_PATTERN.fullmatch(accession):
        return "BLOCKED_NCBI_ONLY_BIOSAMPLE"
    return "BLOCKED_METADATA"


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


def embl_entry_name(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    found: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("AC * "):
                value = line[len("AC * ") :].strip()
                if value.startswith("_"):
                    value = value[1:]
                found.append(value.rstrip(";"))
    if len(found) != 1 or not found[0]:
        raise ValueError(f"Expected exactly one populated 'AC * _<entry>' line in {path}")
    return found[0]


def parse_manifest(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open() as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(f"Malformed manifest line {number}: {line!r}")
            if parts[0] in result:
                raise ValueError(f"Duplicate manifest field: {parts[0]}")
            result[parts[0]] = parts[1]
    return result


def validate_package(package_dir: Path) -> list[str]:
    metadata_path = next(iter(sorted(package_dir.glob("*.package_metadata.json"))), None)
    if metadata_path is None:
        return ["missing_package_metadata"]
    metadata = json.loads(metadata_path.read_text())
    full_seqid = metadata.get("full_seqid", "")
    errors: list[str] = []
    try:
        expected_og = full_seqid_og_id(full_seqid)
        if expected_og != metadata.get("og_id"):
            errors.append("og_id_full_seqid_mismatch")
    except ValueError:
        errors.append("invalid_full_seqid")
    embl = package_dir / f"{full_seqid}.embl.gz"
    chromosomes = package_dir / f"{full_seqid}.chromosome_list.tsv.gz"
    manifest = package_dir / f"{full_seqid}.manifest.txt"
    for path, code in (
        (embl, "missing_embl"),
        (chromosomes, "missing_chromosome_list"),
        (manifest, "missing_manifest"),
    ):
        if not path.is_file():
            errors.append(code)
    if errors:
        return errors
    try:
        if embl_entry_name(embl) != full_seqid:
            errors.append("embl_entry_name_mismatch")
    except (OSError, ValueError):
        errors.append("invalid_embl_entry_name")
    try:
        with gzip.open(chromosomes, "rt", encoding="ascii") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if rows != [[full_seqid, "MT", "Circular-Chromosome", "Mitochondrion"]]:
            errors.append("invalid_chromosome_list")
    except (OSError, UnicodeError):
        errors.append("unreadable_chromosome_list")
    try:
        fields = parse_manifest(manifest)
        required = {
            "STUDY",
            "SAMPLE",
            "ASSEMBLYNAME",
            "ASSEMBLY_TYPE",
            "COVERAGE",
            "PROGRAM",
            "PLATFORM",
            "FLATFILE",
            "CHROMOSOME_LIST",
        }
        if required - fields.keys():
            errors.append("missing_manifest_fields")
        if fields.get("ASSEMBLYNAME") != full_seqid:
            errors.append("manifest_assembly_name_mismatch")
        if fields.get("FLATFILE") != embl.name:
            errors.append("manifest_flatfile_mismatch")
        if fields.get("CHROMOSOME_LIST") != chromosomes.name:
            errors.append("manifest_chromosome_list_mismatch")
    except ValueError:
        errors.append("invalid_manifest")
    return errors


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
    package_status = "READY"
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
        package_status = biosample_block_status(metadata.get("biosample_accession"))
        text = (
            f"# BLOCKED: {blocker}\n"
            f"STUDY\t{metadata.get('study') or ''}\n"
            f"ASSEMBLYNAME\t{full_seqid}\n"
            f"FLATFILE\t{embl.name}\n"
            f"CHROMOSOME_LIST\t{chromosomes.name}\n"
        )
    manifest.write_text(text)
    metadata["package_status"] = package_status
    metadata["metadata_blocker"] = blocker or None
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    errors = validate_package(package_dir) if package_status == "READY" else [package_status]
    local_status = "PASS" if not errors else "FAIL"
    with (package_dir / f"{full_seqid}.local_validation.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["full_seqid", "status", "errors"])
        writer.writerow([full_seqid, local_status, ",".join(errors)])
    primary_artifacts = sorted(
        path
        for path in package_dir.iterdir()
        if path.is_file()
        and (
            path.name.endswith(".embl.gz")
            or path.name.endswith(".chromosome_list.tsv.gz")
            or path.name.endswith(".manifest.txt")
            or path.name.endswith(".locus_tag_mapping.tsv")
            or path.name.endswith(".tbl")
        )
    )
    metadata["local_validation_status"] = local_status
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
    published_package_path = normalise_publish_root(getattr(args, "publish_root", None))
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
    package_status = "READY"
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
        package_status = biosample_block_status(args.biosample)
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
        "schema_version": 1,
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
        "package_status": package_status,
        "metadata_blocker": manifest_error or None,
        "published_package_path": published_package_path,
    }
    metadata_path = package_dir / f"{full_seqid}.package_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if getattr(args, "locus_map", None):
        shutil.copy2(args.locus_map, package_dir / f"{full_seqid}.locus_tag_mapping.tsv")
    if getattr(args, "tagged_tbl", None):
        shutil.copy2(args.tagged_tbl, package_dir / f"{full_seqid}.tbl")
    # Collaborators get sequence plus annotation, not an ENA submission bundle.
    # Both live here so the package directory is the single thing to hand over,
    # and the GFF is keyed on the same locus tags as the submitted record.
    packaged_fasta = package_dir / f"{full_seqid}.fa"
    if Path(args.fasta).resolve() != packaged_fasta.resolve():
        shutil.copy2(args.fasta, packaged_fasta)
    if getattr(args, "gff", None) and getattr(args, "locus_map", None):
        gff_tagged, gff_untagged = tag_gff(
            Path(args.gff), Path(args.locus_map), package_dir / f"{full_seqid}.gff"
        )
        metadata["gff_locus_tag_coverage"] = {
            "tagged": gff_tagged,
            "untagged": gff_untagged,
        }
        if gff_untagged:
            print(
                f"WARNING: {gff_untagged} GFF feature(s) in {full_seqid} carry no locus tag",
                file=sys.stderr,
            )
    errors = validate_package(package_dir) if package_status == "READY" else [package_status]
    local_status = "PASS" if not errors else "FAIL"
    status_path = package_dir / f"{full_seqid}.local_validation.tsv"
    with status_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["full_seqid", "status", "errors"])
        writer.writerow([full_seqid, local_status, ",".join(errors)])
    primary_artifacts = sorted(
        path for path in package_dir.iterdir()
        if path.is_file()
        and (
            path.name.endswith(".embl.gz")
            or path.name.endswith(".chromosome_list.tsv.gz")
            or path.name.endswith(".manifest.txt")
            or path.name.endswith(".locus_tag_mapping.tsv")
            or path.name.endswith(".tbl")
        )
    )
    package_digest = hashlib.sha256(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in primary_artifacts).encode()
    ).hexdigest()
    metadata["local_validation_status"] = local_status
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


def validate_command(args: argparse.Namespace) -> int:
    errors = validate_package(Path(args.package_dir))
    if errors:
        print("\n".join(errors))
        return 1
    print("PASS")
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
    build.add_argument("--locus-map")
    build.add_argument("--tagged-tbl")
    build.add_argument(
        "--gff",
        help="Processed annotation GFF; packaged with locus tags for collaborators.",
    )
    build.add_argument("--outdir", required=True)
    build.add_argument(
        "--publish-root",
        help="Absolute directory this package will occupy after publishDir copies it.",
    )
    build.set_defaults(func=build_package)
    validate = subparsers.add_parser("validate", help="Validate one package directory")
    validate.add_argument("package_dir")
    validate.set_defaults(func=validate_command)
    args = parser.parse_args()
    try:
        return args.func(args)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
