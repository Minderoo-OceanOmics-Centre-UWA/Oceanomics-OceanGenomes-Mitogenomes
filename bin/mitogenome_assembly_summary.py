#!/usr/bin/env python3
"""Build a MultiQC custom-content table for mitogenome assembly QC."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


COLUMNS = [
    "assembly_prefix",
    "sample_id",
    "assembler",
    "status",
    "final_length_bp",
    "circularised",
    "num_candidate_contigs",
    "num_final_contigs",
    "num_genes",
    "num_cds",
    "missing_genes",
    "frameshift_flag",
    "mean_coverage",
    "coverage_cv",
    "reference_species",
    "reference_accession",
    "numt_flag",
    "manual_review_reason",
]

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")
MISSING = ""
PLACEHOLDER_VALUES = {"", ".", "NA", "na", "N/A", "n/a", "none", "None", "null", "unknown", "Unknown"}


@dataclass
class Thresholds:
    min_mean_coverage: float | None = None
    max_coverage_cv: float | None = None
    min_length: int | None = None
    max_length: int | None = None
    expected_gene_count: int | None = None
    expected_pcg_count: int | None = None
    trna_tolerance: int = 2


@dataclass
class RunFiles:
    sample_id: str
    prefix: str
    assembler: str
    files: list[Path] = field(default_factory=list)

    def add(self, path: Path) -> None:
        if path not in self.files:
            self.files.append(path)


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def normalise_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def parse_bool(value: object) -> str:
    text = str(value or "").strip().lower()
    if not text or text in {"na", "nan", "none", "null"}:
        return MISSING
    if text in {"true", "yes", "y", "1", "circular", "circularized", "circularised"}:
        return "true"
    if "no frameshift" in text or "no frame shift" in text:
        return "false"
    if "frameshift" in text or "frame shift" in text:
        return "true"
    if text in {"false", "no", "n", "0", "linear", "not_circular", "not circular"}:
        return "false"
    return MISSING


def first_value(row: dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        value = row.get(name)
        if value is None:
            continue
        text = str(value).strip()
        if text not in PLACEHOLDER_VALUES:
            return text
    return MISSING


def is_missing(value: object) -> bool:
    return str(value or "").strip() in PLACEHOLDER_VALUES


def first_numeric(row: dict[str, str], names: Iterable[str]) -> float | None:
    for name in names:
        value = row.get(name)
        if value in (None, "", "NA", "na", "null"):
            continue
        number = parse_number(value)
        if number is not None:
            return number
    return None


def parse_number(value: object) -> float | None:
    text = str(value or "").strip().replace(",", "")
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def format_number(value: float | int | None) -> str:
    if value is None:
        return MISSING
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def parse_fasta(path: Path) -> tuple[int | None, int | None]:
    seq_count = 0
    total_length = 0
    saw_sequence = False
    try:
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    seq_count += 1
                    continue
                saw_sequence = True
                total_length += len(re.sub(r"\s+", "", line))
    except OSError as error:
        warn(f"Could not read FASTA {path}: {error}")
        return None, None
    if not saw_sequence and seq_count == 0:
        return None, None
    return total_length, seq_count


def read_table(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as handle:
            lines = [line for line in handle if not line.lstrip().startswith("#")]
            sample = "".join(lines[:20])
            delimiter = "\t" if "\t" in sample else ","
            reader = csv.DictReader(lines, delimiter=delimiter)
            if not reader.fieldnames:
                return []
            reader.fieldnames = [normalise_header(name) for name in reader.fieldnames]
            parsed_rows = []
            for row in reader:
                parsed_rows.append(
                    {
                        normalise_header(str(key)): (value or "").strip()
                        for key, value in row.items()
                        if key is not None
                    }
                )
            return parsed_rows
    except OSError as error:
        warn(f"Could not read table {path}: {error}")
        return []


def read_text(path: Path, max_chars: int = 1_000_000) -> str:
    try:
        return path.read_text(errors="replace")[:max_chars]
    except OSError as error:
        warn(f"Could not read text file {path}: {error}")
        return ""


def collect_input_files(inputs: Iterable[Path]) -> list[Path]:
    files: list[Path] = []
    seen_dirs: set[Path] = set()

    def add_from_directory(directory: Path) -> None:
        try:
            resolved = directory.resolve()
        except OSError:
            resolved = directory
        if resolved in seen_dirs:
            return
        seen_dirs.add(resolved)

        for path in directory.rglob("*"):
            if path.is_file():
                files.append(path)
            elif path.is_dir() and path.is_symlink():
                add_from_directory(path)

    for root in inputs:
        if root.is_file():
            files.append(root)
        elif root.exists():
            add_from_directory(root)
        else:
            warn(f"Input path does not exist: {root}")

    return files


def sample_from_prefix(prefix: str) -> str:
    return prefix.split(".")[0] if "." in prefix else prefix


def strip_known_suffix(name: str) -> str:
    suffixes = [
        ".contigs_stats.with_coverage.tsv",
        ".contigs_stats.tsv",
        ".circularity_check.tsv",
        ".getorg_check.tsv",
        ".concatemer_collapse.tsv",
        ".post_curation_check.tsv",
        ".oatk_status.tsv",
        ".assembly_status.tsv",
        ".findmitoreference_status.tsv",
        ".reference_candidates_status.tsv",
        ".reference_ranking.tsv",
        ".scaffold_join.tsv",
        ".mito_depth.tsv",
        ".coverage.tsv",
        ".reference_relevance.txt",
        ".reference_divergence.txt",
        ".get_org.log.txt",
        ".hifiasm.log",
        ".log",
        ".gb",
        ".gbk",
        ".gbf",
        ".fasta",
        ".fna",
        ".fa",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def infer_prefix(path: Path, assembler: str) -> str:
    name = path.name
    if assembler == "GetOrganelle":
        for marker in (".animal_mt.", ".embplant_mt.", ".embplant_pt.", ".fungus_mt."):
            if marker in name:
                return name.split(marker, 1)[0]
        if ".extended_K" in name:
            return name.split(".extended_K", 1)[0]
    if assembler == "Oatk":
        # Oatk outputs are <prefix>.fasta / <prefix>.gfa / <prefix>.mito.ctg.fasta /
        # <prefix>.oatk.log plus the shared <prefix>.annotation_stats.csv; the run
        # prefix is the part before the first oatk-specific suffix.
        for marker in (".mito.ctg.", ".mito.", ".oatk.", ".annot_mito."):
            if marker in name:
                return name.split(marker, 1)[0]
    return strip_known_suffix(name)


def classify_file(path: Path) -> str | None:
    name = path.name.lower()
    parts = {part.lower() for part in path.parts}
    ignored_parts = {"annotation", "emma", "lca", "genbank", "cds", "proteins", "genes"}
    ignored_prefixes = (
        "mt-",
        "blast.",
        "filtered_summary.",
        "lca.",
        "lca_raw.",
        "lca_short.",
    )
    if name.endswith(".annotation_stats.csv"):
        return None
    # Reference GenBanks (the relabelled <prefix>.reference.gb and the published
    # MitoReference/reference_seed records) are inputs to the summary, not assembly
    # runs of their own; never let them spawn a spurious run row.
    if is_reference_genbank(path):
        return None
    if parts.intersection(ignored_parts) or name.startswith(ignored_prefixes):
        return None
    # Oatk (reference-free HiFi fallback) emits <prefix>.fasta / <prefix>.gfa /
    # <prefix>.mito.ctg.fasta / <prefix>.oatk.log, and its prefix carries the "oatk"
    # assembler tag. Check before MitoHiFi so an oatk run is never misfiled (its
    # prefix never contains "mitohifi").
    if "oatk" in name or ".mito.ctg." in name:
        return "Oatk"
    if "mitohifi" in name or "hifiasm" in name or "contigs_stats" in name:
        return "MitoHiFi"
    if (
        "getorg" in name
        or "get_org.log.txt" in name
        or "path_sequence" in name
        or "selected_graph" in name
        or "assembly_graph.fastg" in name
        or "getorganelle" in parts
    ):
        return "GetOrganelle"
    return None


def discover_sample_dirs(inputs: Iterable[Path]) -> list[Path]:
    sample_dirs: set[Path] = set()
    for root in inputs:
        if root.is_file():
            root = root.parent
        for path in root.rglob("*") if root.exists() else []:
            if path.is_dir() and (path / "mtdna").is_dir():
                sample_dirs.add(path)
    return sorted(sample_dirs)


# Tokens that mark the assembler-code field of an assembly prefix. Same tokens
# classify_file keys on, because they are the same naming convention: the
# mt_assembly_prefix's last dot-field is always <version><assembler>, e.g.
# v323mitohifi, v323mitohifi_collapsed, getorg1770, getorg1770reseed_rgj, v10oatk.
ASSEMBLER_PREFIX_TOKENS = ("mitohifi", "hifiasm", "getorg", "oatk")


def has_assembler_token(field: str) -> bool:
    lowered = field.lower()
    return any(token in lowered for token in ASSEMBLER_PREFIX_TOKENS)


def is_assembly_run_prefix(prefix: str) -> bool:
    """Reject a prefix that is a real assembly prefix with a sidecar suffix glued on.

    strip_known_suffix falls back to Path(name).stem for any suffix it does not
    know, so a new per-run sidecar silently becomes an assembly of its own:
    <prefix>.reference_ranking.tsv used to manufacture a whole extra row, reported
    as `failed` because it has no FASTA, and inheriting its parent's annotation via
    the substring join. That was 38 of 250 rows in the mitogenomes-missing-audit-6
    cohort. Listing each new suffix in strip_known_suffix only fixes the ones we
    already know about, so also reject structurally: if an EARLIER dot-field
    carries the assembler token and the LAST one does not, the tail is a suffix,
    not an assembly.

    Deliberately conservative -- a prefix with no assembler token anywhere is left
    alone rather than dropped, because losing a real assembly row is worse than
    keeping a spurious one.
    """
    fields = prefix.split(".")
    if len(fields) < 2 or has_assembler_token(fields[-1]):
        return True
    return not any(has_assembler_token(field) for field in fields[:-1])


def discover_assembler_runs(inputs: Iterable[Path]) -> list[RunFiles]:
    runs: dict[tuple[str, str], RunFiles] = {}
    files = collect_input_files(inputs)

    for path in files:
        assembler = classify_file(path)
        if not assembler:
            continue
        prefix = infer_prefix(path, assembler)
        # The file stays in the flat list collect_input_files returned, so
        # files_for_run / reference_for_run can still see it; it just no longer
        # spawns an assembly row of its own.
        if not is_assembly_run_prefix(prefix):
            continue
        key = (assembler, prefix)
        runs.setdefault(key, RunFiles(sample_from_prefix(prefix), prefix, assembler)).add(path)

    return sorted(runs.values(), key=lambda run: (run.sample_id, run.assembler, run.prefix))


def choose_final_fasta(files: Iterable[Path], assembler: str, prefix: str) -> Path | None:
    candidates = [
        path
        for path in files
        if path.suffix.lower() in FASTA_EXTENSIONS
        and "emma" not in path.name.lower()
        and "all_potential_contigs" not in path.name.lower()
        and "assembly_graph" not in path.name.lower()
        and "selected_graph" not in path.name.lower()
    ]
    if not candidates:
        return None

    def score(path: Path) -> tuple[int, int, str]:
        name = path.name.lower()
        value = 0
        if assembler == "GetOrganelle":
            if "path_sequence" in name:
                value += 50
            if "complete" in name:
                value += 25
            if path.name == f"{prefix}.fasta":
                value += 40
        elif assembler == "Oatk":
            # The module republishes the contig as <prefix>.fasta (annotated); the
            # native <prefix>.mito.ctg.fasta is the fallback.
            if path.name == f"{prefix}.fasta":
                value += 60
            elif name.endswith(".mito.ctg.fasta"):
                value += 40
        else:
            if path.name == f"{prefix}.fasta":
                value += 50
            if "mitohifi" in name:
                value += 20
        return (-value, len(path.name), path.name)

    return sorted(candidates, key=score)[0]


def annotation_row_prefix(row: dict[str, str]) -> str:
    """Rebuild the assembly prefix an annotation_stats.csv row describes.

    bin/annotation_stats.py splits the GFF stem into og_id/tech/seq_date/code/
    annotation, so those four fields reassemble the mt_assembly_prefix exactly.
    Returns "" when any component is blank -- which is what that script writes
    when the stem is not 5 dot-fields -- so a degenerate row matches nothing
    rather than matching everything.
    """
    og_id = first_value(row, ["sample", "sample_id", "og_id"])
    tech = first_value(row, ["tech", "sequencing_type"])
    seq_date = first_value(row, ["seq_date", "date"])
    code = first_value(row, ["code", "assembly", "assembler", "mt_assembly_prefix"])
    if not all((og_id, tech, seq_date, code)):
        return ""
    return f"{og_id}.{tech}.{seq_date}.{code}"


def parse_annotation_stats(files: Iterable[Path], prefix: str) -> dict[str, str]:
    """Gene counts for ONE assembly, from that assembly's own re-annotation.

    Matching is exact, never a substring test. The old `og_id not in prefix` /
    `code not in prefix` guards leaked badly in both directions: "OG5" is a
    substring of "OG58", "OG8" of "OG810" and "OG848", "OG10" of "OG107", and
    "getorg1770" of "getorg1770reseed" -- so one sample's annotation was reported
    against another's assembly (45 rows of the mitogenomes-missing-audit-6 cohort
    carried counts they never earned). Because the loop returned the FIRST match
    over an unsorted rglob, which sample won was filesystem-order dependent too.
    """
    # EMMA renames its output to exactly <mt_assembly_prefix>.annotation_stats.csv
    # (modules/local/upload_results/emma/main.nf), so the filename IS the join key.
    # Sorted so a duplicate stage-in can never make the result order dependent.
    candidates = sorted(
        (path for path in files if path.name == f"{prefix}.annotation_stats.csv"),
        key=lambda path: str(path),
    )
    rows = []
    for path in candidates:
        rows.extend(read_table(path))

    if not rows:
        # Fall back to the row's own contents for inputs that were not renamed,
        # still requiring whole-prefix equality.
        for path in sorted(
            (path for path in files if path.name.endswith(".annotation_stats.csv")),
            key=lambda path: str(path),
        ):
            rows.extend(row for row in read_table(path) if annotation_row_prefix(row) == prefix)

    for row in rows:
        num_cds = first_numeric(row, ["num_cds"])
        num_trna = first_numeric(row, ["num_trna"])
        num_rrna = first_numeric(row, ["num_rrna"])
        num_genes = first_numeric(row, ["num_genes", "gene_count"])
        if num_genes is None and None not in (num_cds, num_trna, num_rrna):
            num_genes = (num_cds or 0) + (num_trna or 0) + (num_rrna or 0)

        return {
            "num_genes": format_number(num_genes),
            "num_cds": format_number(num_cds),
            "missing_genes": first_value(row, ["missing_genes", "num_missing"]),
            "frameshift_flag": parse_bool(first_value(row, ["frameshift", "frameshift_flag", "frameshifts"])),
            # Which completeness profile annotation_stats.py judged this assembly
            # under. Absent on rows written before that column existed, which
            # read_expected_gene_count() treats as the vertebrate profile.
            "completeness_profile": first_value(row, ["completeness_profile"]),
        }
    return {}


def _coverage_from_row(row: dict[str, str]) -> tuple[float | None, float | None]:
    mean = first_numeric(row, ["mean_depth", "mean_coverage", "avg_coverage", "average_coverage", "coverage"])
    cv = first_numeric(row, ["depth_cv", "coverage_cv", "cv"])
    return mean, cv


def _raw_depth_stream(path: Path) -> tuple[float | None, float | None]:
    """Last-ditch: treat the file as raw `samtools depth` output."""
    depths = []
    try:
        with path.open() as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                parts = re.split(r"[\t, ]+", line.strip())
                if len(parts) >= 3 and parts[1].isdigit():
                    depth = parse_number(parts[2])
                    if depth is not None:
                        depths.append(depth)
    except OSError:
        return None, None
    if depths:
        mean = statistics.fmean(depths)
        return mean, (statistics.pstdev(depths) / mean if mean else None)
    return None, None


def parse_coverage(files: Iterable[Path]) -> tuple[float | None, float | None]:
    """Mean depth and its CV for a run, most trustworthy source first.

    Ordered explicitly rather than by a substring match. The previous predicate
    (`endswith(".coverage.tsv") or "coverage" in name and suffix == ".tsv"`) also
    matched `.contigs_stats.with_coverage.tsv` through its second clause -- note
    that file does NOT end in ".coverage.tsv", the character before "coverage" is
    an underscore -- so whichever of the two happened to come first in `files`
    won, non-deterministically.

    Precedence:
      1. <prefix>.mito_depth.tsv  -- the uniform remap-based measurement. Same
         definition for every assembler and platform, so it is authoritative
         wherever it exists.
      2. <prefix>.coverage.tsv    -- MitoHiFi's depth over reference-recruited
         reads only. Legacy; kept so reruns of old output still report something.
      3. <prefix>.contigs_stats.with_coverage.tsv -- the same number joined into
         the contig table, read from the final_mitogenome row rather than row 0.
    """
    files = list(files)

    for path in files:
        if path.name.endswith(".mito_depth.tsv"):
            rows = read_table(path)
            if rows:
                mean, cv = _coverage_from_row(rows[0])
                if mean is not None or cv is not None:
                    return mean, cv

    for path in files:
        if path.name.endswith(".coverage.tsv"):
            rows = read_table(path)
            if rows:
                mean, cv = _coverage_from_row(rows[0])
                if mean is not None or cv is not None:
                    return mean, cv
            mean, cv = _raw_depth_stream(path)
            if mean is not None:
                return mean, cv

    for path in files:
        if path.name.endswith(".contigs_stats.with_coverage.tsv"):
            rows = read_table(path)
            # Prefer the final_mitogenome row; row 0 is only a coincidence.
            target = next(
                (row for row in rows if str(row.get("contig_id", "")).strip() == "final_mitogenome"),
                rows[0] if rows else None,
            )
            if target:
                mean, cv = _coverage_from_row(target)
                if mean is not None or cv is not None:
                    return mean, cv

    return None, None


def parse_genbank_reference(path: Path) -> tuple[str, str]:
    species = MISSING
    accession = MISSING
    for line in read_text(path).splitlines():
        if line.startswith("VERSION"):
            parts = line.split()
            if len(parts) > 1:
                accession = parts[1]
        elif line.strip().startswith("ORGANISM"):
            species = line.strip().replace("ORGANISM", "", 1).strip()
        if species and accession:
            break
    return species, accession


GENBANK_SUFFIXES = {".gb", ".gbk", ".gbf"}


def is_reference_genbank(path: Path) -> bool:
    """True for a GenBank file that holds a *reference* mitogenome (not the assembly's own).

    Recognises the relabelled file the pipeline stages into the summary inputs
    (``<assembly_prefix>.reference.gb``), the per-sample published layouts
    (``mtdna/MitoReference/NC_*.gb`` for MitoHiFi, ``mtdna/reference_seed/NC_*.gb``
    for the GetOrganelle reseed), and bare RefSeq-accession-named records. The
    assembly's own annotation (``<assembly_prefix>.gb``) is deliberately excluded.
    """
    if path.suffix.lower() not in GENBANK_SUFFIXES:
        return False
    name = path.name.lower()
    parts = {part.lower() for part in path.parts}
    if name.endswith(tuple(f".reference{suffix}" for suffix in GENBANK_SUFFIXES)):
        return True
    if "mitoreference" in parts or "reference_seed" in parts:
        return True
    if name.startswith(("nc_", "nw_", "nz_")):
        return True
    return False


def path_carries_prefix(path: Path, prefix: str) -> bool:
    """True when `path` is named for the assembly `prefix`, not merely for one
    whose name starts with it.

    A bare `prefix in str(path)` also matches every curated variant built from it
    -- <prefix>_collapsed, <prefix>reseed, <prefix>reseed_rgj -- which pulls the
    child's files into the parent's evidence pool. Require the prefix to end at a
    filename or path-component boundary instead.
    """
    if path.name == prefix or path.name.startswith(f"{prefix}."):
        return True
    return any(part == prefix for part in path.parts)


def files_for_run(files: Iterable[Path], run: "RunFiles") -> list[Path]:
    """Restrict a flat file list to those belonging to a single assembly run.

    The summary module stages every sample into one flat directory, so the
    per-sample directory structure is gone in production; association then falls
    back to the sample id that prefixes every relabelled reference filename. The
    per-sample published layout (where files still sit under ``<sample>/``) is
    matched by the directory checks.
    """
    sample_id = run.sample_id
    matches = []
    for path in files:
        if (
            path_carries_prefix(path, run.prefix)
            or f"/{sample_id}/" in str(path)
            or path.name.startswith(f"{sample_id}.")
        ):
            matches.append(path)
    return matches


def find_reference_genbank(files: Iterable[Path]) -> tuple[str, str]:
    genbanks = [path for path in files if is_reference_genbank(path)]

    def score(path: Path) -> tuple[int, str]:
        lowered = [part.lower() for part in path.parts]
        name = path.name.lower()
        value = 0
        if name.endswith(tuple(f".reference{suffix}" for suffix in GENBANK_SUFFIXES)):
            value -= 150
        if "mitoreference" in lowered or "reference_seed" in lowered:
            value -= 100
        if name.startswith(("nc_", "nw_", "nz_")):
            value -= 50
        return value, str(path)

    for path in sorted(genbanks, key=score):
        species, accession = parse_genbank_reference(path)
        if not is_missing(species) and not is_missing(accession):
            return species, accession
    return MISSING, MISSING


def reference_for_run(run: "RunFiles", all_files: Iterable[Path]) -> tuple[str, str]:
    """Find the reference species/accession for one run, scoped to that run's files.

    The lookup is strictly scoped to files belonging to this run: a run with no
    reference of its own must report nothing rather than borrow another sample's
    reference (every sample is staged into one flat directory in production).
    """
    scoped = [path for path in files_for_run(all_files, run) if is_reference_genbank(path)]
    if not scoped:
        return MISSING, MISSING
    return find_reference_genbank(scoped)


def parse_mitohifi_stats(files: Iterable[Path]) -> dict[str, str]:
    stats: dict[str, str] = {}
    rows = []
    for path in files:
        if "contigs_stats" in path.name.lower():
            rows.extend(read_table(path))
    if not rows:
        return stats
    unique_rows = []
    seen_rows = set()
    for row in rows:
        row_key = (
            first_value(row, ["contig_id", "contig", "name"]),
            first_value(row, ["annotation_file", "file"]),
            first_value(row, ["length_bp", "length", "len", "size", "bp", "contig_length"]),
        )
        if row_key in seen_rows:
            continue
        seen_rows.add(row_key)
        unique_rows.append(row)
    rows = unique_rows

    stats["num_candidate_contigs"] = str(len(rows))
    final_rows = [
        row
        for row in rows
        if "final" in " ".join(row.values()).lower()
        or parse_bool(first_value(row, ["selected", "final", "chosen"])) == "true"
    ]
    final_row = final_rows[0] if final_rows else rows[0]

    length = first_numeric(final_row, ["final_length_bp", "length", "len", "size", "bp", "contig_length"])
    circular = parse_bool(first_value(final_row, ["circularised", "circularized", "circular", "is_circular", "was_circular"]))
    frameshift = parse_bool(first_value(final_row, ["frameshift_flag", "frameshift", "frameshifts", "frameshifts_found"]))
    mean_cov = first_numeric(final_row, ["mean_coverage", "avg_coverage", "average_coverage"])

    if length is not None:
        stats["final_length_bp"] = format_number(length)
    # num_genes is deliberately NOT taken from contigs_stats' number_of_genes.
    # That is MitoHiFi's own reference-guided annotation, which every assembly is
    # re-annotated over by EMMA / MITOS2 downstream -- see the note at
    # subworkflows/local/mitogenome_assembly/mitohifi/main.nf:64-68. Reading it
    # here made num_genes and num_cds come from two different annotations of two
    # different molecules: OG750's pre-collapse concatemer reported 56 genes
    # counted over the 2.14x multimer, next to a num_cds from the collapsed
    # monomer. Both counts now come from the same *.annotation_stats.csv row, and
    # an assembly that never reached annotation reports neither.
    if circular:
        stats["circularised"] = circular
    if frameshift:
        stats["frameshift_flag"] = frameshift
    if mean_cov is not None:
        stats["mean_coverage"] = format_number(mean_cov)
    reference_species = first_value(final_row, ["reference_species", "ref_species", "species"])
    reference_accession = first_value(final_row, ["reference_accession", "ref_accession", "accession"])
    if reference_species:
        stats["reference_species"] = reference_species
    if reference_accession:
        stats["reference_accession"] = reference_accession
    return {key: value for key, value in stats.items() if value}


def has_numt_signal(files: Iterable[Path]) -> bool:
    # Note: do NOT match the bare token "numt" or the "contigs_filtering"
    # directory name -- both appear in every MitoHiFi run as part of its routine
    # output (the standard "filter ... to avoid NUMTS" log line and the
    # contigs_filtering/ working dir), so matching them flags ~every HiFi
    # assembly as a false positive. Keep only patterns that indicate an actual
    # nuclear-mitochondrial / low-confidence call.
    patterns = [
        "nuclear mitochondrial",
        "rejected",
        "low-confidence",
        "low confidence",
        "low_coverage",
        "divergent",
    ]
    for path in files:
        text = str(path).lower()
        if path.suffix.lower() in {".tsv", ".txt", ".log"}:
            text += "\n" + read_text(path, max_chars=200_000).lower()
        if any(pattern in text for pattern in patterns):
            return True
    return False


def anomaly_for_run(run: "RunFiles") -> dict[str, str]:
    """Read the circularity-check sidecar for this run and return its
    length/repeat anomaly fields. The MitoHiFi check writes
    <prefix>.circularity_check.tsv and the GetOrganelle check writes
    <prefix>.getorg_check.tsv; both carry anomaly_type / length_anomaly. The check
    flags over-length assemblies as a concatemer (tandem genome duplication), a
    control_region_repeat (D-loop VNTR) or unresolved; apply_qc turns these into
    manual-review reasons."""
    for path in run.files:
        if path.name.endswith((".circularity_check.tsv", ".getorg_check.tsv")):
            rows = read_table(path)
            if rows:
                row = rows[0]
                return {
                    "anomaly_type": first_value(row, ["anomaly_type"]),
                    "length_anomaly": first_value(row, ["length_anomaly"]),
                }
    return {}


COLLAPSED_SUFFIX = "_collapsed"


def collapse_child_evidence(run: "RunFiles", all_files: Iterable[Path]) -> dict[str, str]:
    """Circularity / anomaly for a collapsed monomer, from the post-curation check.

    COLLAPSE_CONCATEMER re-runs the circularity and length checks on the molecule
    it produced, but writes the result under the PRE-collapse prefix
    (<parent>.post_curation_check.tsv). So the monomer's own row -- keyed
    <parent>_collapsed -- has no circularity evidence of its own at all, and would
    be reported as non-circular despite the check having confirmed it is
    (OG750: final_verdict_circular = True on a 15293 bp monomer).

    Reach across to the parent for exactly this one file and nothing else.
    Everything else stays strictly scoped to the run's own prefix; inheriting a
    parent's evidence wholesale is the bug this is carved out of, not a pattern to
    extend. It is sound here only because the check was measured ON the child.
    """
    if not run.prefix.endswith(COLLAPSED_SUFFIX):
        return {}
    parent = run.prefix[: -len(COLLAPSED_SUFFIX)]
    for path in sorted(all_files, key=lambda item: str(item)):
        if path.name != f"{parent}.post_curation_check.tsv":
            continue
        rows = read_table(path)
        if not rows:
            continue
        row = rows[0]
        evidence = {
            "circularised": parse_bool(first_value(row, ["final_verdict_circular"])),
            "anomaly_type": first_value(row, ["anomaly_type"]),
            "length_anomaly": first_value(row, ["length_anomaly"]),
        }
        return {key: value for key, value in evidence.items() if value}
    return {}


def collapse_for_run(run: "RunFiles") -> dict[str, str] | None:
    """Read the concatemer-collapse report (<prefix>.concatemer_collapse.tsv) for
    this run, if the auto-curation step ran. The report is written under the
    PRE-collapse prefix, so finding one on a run means this run is the molecule
    the collapse acted on -- never the monomer it produced."""
    for path in run.files:
        if path.name.endswith(".concatemer_collapse.tsv"):
            rows = read_table(path)
            if rows:
                return rows[0]
    return None


def apply_collapse_provenance(row: dict[str, str], run: "RunFiles") -> None:
    """Mark a pre-collapse concatemer as superseded by the monomer it produced.

    COLLAPSE_CONCATEMER renames a genuinely collapsed assembly to
    <prefix>_collapsed.fasta, and every downstream stage forks on that basename,
    so the monomer now carries its own identity end to end: its own publish dir,
    its own mitogenome_data row (see the provenance rows in
    workflows/oceangenomesmitogenomes.nf) and its own remap depth. It therefore
    gets its own row here too, and this run is the superseded original.

    This function used to write the monomer's length onto THIS row instead --
    correct back when the curated monomer was filed under the original's name,
    but since the fork it put the collapsed length (OG750: 15293) next to the
    pre-collapse assembly's contig counts and reference stats, while the row that
    was actually named for the monomer sat empty. So the length is left alone;
    the row keeps its own real 32672 bp and is flagged superseded so it drops out
    of the manual-review queue rather than being triaged twice.

    A passthrough report means nothing was rewritten, no second row exists, and
    the row is left untouched.
    """
    collapse = collapse_for_run(run)
    if not collapse or (collapse.get("action") or "").strip().lower() != "collapsed":
        return
    # Only the marker here. apply_qc rebuilds manual_review_reason from scratch and
    # runs after this, so the reason itself is recorded in finalise_status.
    row[SUPERSEDED_KEY] = "true"


def getorg_circular_override(run: "RunFiles") -> str:
    """Corrected circular verdict from the GetOrganelle check sidecar
    (<prefix>.getorg_check.tsv -> final_verdict_circular). 'true'/'false' or ''.
    Lets a single scaffold the reference test confirms circular be recorded as
    circularised even though GetOrganelle's log said 'N scaffold(s)'."""
    for path in run.files:
        if path.name.endswith(".getorg_check.tsv"):
            rows = read_table(path)
            if rows:
                return parse_bool(first_value(rows[0], ["final_verdict_circular"])) or ""
    return ""


def reference_relevance_for_run(run: "RunFiles") -> str:
    """Read the REFERENCE_RELEVANCE flag for this run, or ''.

    PASS | DIVERGENT | MISMATCH | UNKNOWN. The per-sample
    <prefix>.reference_relevance.txt records how well the resolved reference
    corresponds to the assembly: MISMATCH means the reference neither covers nor
    matches it, so the species label most likely pointed findMitoReference at a
    wrong-family reference; DIVERGENT means it is the right molecule but a distant
    relative, which degrades seeding and annotation transfer without making the
    assembly wrong.
    """
    for path in run.files:
        if path.name.endswith(".reference_relevance.txt"):
            text = read_text(path).strip()
            if text:
                return text.split("\t", 1)[0].strip().upper()
    return ""


def reference_divergence_for_run(run: "RunFiles") -> str:
    """Read the REFERENCE_DIVERGENCE tier for this run, or ''.

    The per-sample <prefix>.reference_divergence.txt records, pre-assembly, how
    closely the resolved reference is related to the sample by taxonomy
    (CONGENERIC | CONFAMILIAL | DIFFERENT_FAMILY | NON_CONGENERIC | UNKNOWN).
    Anything other than CONGENERIC means findMitoReference fell back to a
    non-congeneric reference, a known cause of gene-incomplete MitoHiFi collapses
    for deep-sea / poorly-sampled taxa.
    """
    for path in run.files:
        if path.name.endswith(".reference_divergence.txt"):
            text = read_text(path).strip()
            if text:
                return text.split("\t", 1)[0].strip().upper()
    return ""


def apply_qc(row: dict[str, str], thresholds: Thresholds) -> None:
    reasons = []

    length = parse_number(row.get("final_length_bp"))
    mean_cov = parse_number(row.get("mean_coverage"))
    cov_cv = parse_number(row.get("coverage_cv"))
    num_candidate = parse_number(row.get("num_candidate_contigs"))
    num_final = parse_number(row.get("num_final_contigs"))
    num_genes = parse_number(row.get("num_genes"))
    num_cds = parse_number(row.get("num_cds"))

    if not row.get("final_length_bp"):
        reasons.append("missing_final_fasta")
    if row.get("circularised") == "false":
        reasons.append("not_circularised")
    # hifiasm routinely emits many candidate contigs and MitoHiFi picks the best
    # one, so num_candidate > 1 is the norm, not a defect. Only flag for review
    # when the *final* assembly is not clean (multiple final contigs or not
    # circular); a single circular final contig is fine regardless of candidates.
    clean_final = num_final == 1 and row.get("circularised") == "true"
    if num_candidate is not None and num_candidate > 1 and not clean_final:
        reasons.append("multiple_candidate_contigs")
    if num_final is not None and num_final > 1:
        reasons.append("multiple_final_contigs")
    if row.get("missing_genes") and row["missing_genes"].lower() not in {"no", "0", "none", "na"}:
        reasons.append("missing_genes")
    if row.get("frameshift_flag") == "true":
        reasons.append("frameshift_detected")
    if thresholds.min_mean_coverage is not None and mean_cov is not None and mean_cov < thresholds.min_mean_coverage:
        reasons.append("low_mean_coverage")
    if thresholds.max_coverage_cv is not None and cov_cv is not None and cov_cv > thresholds.max_coverage_cv:
        reasons.append("high_coverage_variability")
    if row.get("numt_flag") == "true":
        reasons.append("possible_numt")
    if row.get("reference_relevance") == "MISMATCH":
        reasons.append("reference_mismatch")
    elif row.get("reference_relevance") == "DIVERGENT":
        reasons.append("reference_divergent")
    # Pre-assembly taxonomy guard: a non-congeneric reference is the leading cause
    # of gene-incomplete MitoHiFi collapses for taxa with no close NCBI relative.
    # UNKNOWN (unparseable reference / missing taxonomy) is not treated as a defect.
    if row.get("reference_divergence") in {"CONFAMILIAL", "DIFFERENT_FAMILY", "NON_CONGENERIC", "CROSS_ORDER"}:
        reasons.append("no_congeneric_reference")
    # Length / tandem-repeat anomaly from the circularity check. The specific type
    # (concatemer / control_region_repeat / unresolved) is the curation action, so
    # surface it verbatim; fall back to a generic flag if only length_anomaly is set.
    anomaly_type = (row.get("anomaly_type") or "").strip().lower()
    if anomaly_type and anomaly_type not in PLACEHOLDER_VALUES and anomaly_type != "none":
        reasons.append(anomaly_type)
    elif (row.get("length_anomaly") or "").strip().lower() == "yes":
        reasons.append("length_anomaly")
    if (
        thresholds.min_length is not None
        and length is not None
        and length < thresholds.min_length
    ) or (
        thresholds.max_length is not None
        and length is not None
        and length > thresholds.max_length
    ):
        reasons.append("length_outside_expected_range")
    expected_genes = expected_gene_count_for(row, thresholds)
    if expected_genes is not None and num_genes is not None and num_genes < expected_genes:
        reasons.append("missing_genes")
    # Protein-coding-gene check. More robust than the total gene count: a collapse
    # can drop several PCGs while tRNAs keep num_genes near the expected total, so
    # gate the CDS count directly (e.g. < 13 PCGs for a vertebrate mitogenome).
    if thresholds.expected_pcg_count is not None and num_cds is not None and num_cds < thresholds.expected_pcg_count:
        reasons.append("missing_protein_coding_genes")

    # Data-limited tag: a fragmented / non-circular assembly at coverage well below
    # the minimum is a sequencing-depth problem (low mito content, e.g. HiC or
    # low-yield libraries), not something the pipeline can assemble its way out of.
    # Tagging it lets triage separate data-limited samples from ones the pipeline
    # could still improve. It never changes status on its own.
    if (
        thresholds.min_mean_coverage is not None and mean_cov is not None
        and mean_cov < thresholds.min_mean_coverage * DATA_LIMITED_COVERAGE_FRACTION
        and any(reason in FRAGMENTATION_REASONS for reason in reasons)
    ):
        reasons.append("data_limited")

    deduped = []
    for reason in reasons:
        if reason not in deduped:
            deduped.append(reason)
    row["manual_review_reason"] = ";".join(deduped)
    row[BLOCKING_KEY] = ";".join(blocking_reasons(row, deduped, thresholds))


# Reasons that describe a *complete* mitogenome rather than a defect: they are
# retained in manual_review_reason for transparency but, on an assembly that
# passes the complete-core guard, they no longer force a manual_review status.
#
# The reference reasons are here because they describe the *reference*, not the
# assembly. A circular molecule of the expected length carrying all 13 PCGs and
# both rRNAs is a finished mitogenome whatever reference was used to build it, so
# the reference verdict is provenance metadata for curation (and a signal that a
# closer reference would help on a rerun), not evidence of an assembly defect.
# Anything that does describe the assembly -- missing_protein_coding_genes and
# every structural flag -- still blocks, so a genuinely bad reference that damaged
# the assembly is still caught, by the damage rather than by the reference.
#
# high_coverage_variability belongs here for the same reason as low_mean_coverage:
# depth variability describes the library, not the assembly. It was previously
# blocking only because it could not fire in practice -- coverage_cv was populated
# for MitoHiFi alone, and even there it was measured on a LINEAR reference, so
# 15-20 kb HiFi reads produced a triangular depth profile whose CV was mostly an
# artefact of linearisation. The uniform remap folds circular molecules, which
# removes that artefact and switches the metric on for GetOrganelle and Oatk as
# well. Leaving it blocking would mean measuring coverage properly caused
# previously-passing finished mitogenomes to start failing, which inverts the
# intent. It remains visible in manual_review_reason, and still blocks on anything
# that is not a complete-core assembly.
ADVISORY_WHEN_COMPLETE = {
    "no_congeneric_reference", "low_mean_coverage", "high_coverage_variability",
    "reference_mismatch", "reference_divergent",
}
# Number of tRNAs a complete-core assembly may be missing (annotation limitation,
# not an assembly defect) while all 13 PCGs + 2 rRNAs are still present.
# Default only. The effective value comes from --trna-tolerance (wired to
# params.annotation_trna_tolerance), the SAME knob annotation_stats.py uses for
# its pass/hold decision -- a hardcoded constant here would silently disagree
# with the gate the moment that param is changed.
TRNA_TOLERANCE = 2
# A fragmented / non-circular assembly whose mean coverage is below this fraction
# of the minimum coverage threshold is treated as sequencing-depth limited rather
# than a pipeline defect, and tagged data_limited for triage.
DATA_LIMITED_COVERAGE_FRACTION = 0.75
# Reasons that are symptoms of shallow / low-mito-content data rather than a
# pipeline shortcoming.
FRAGMENTATION_REASONS = {
    "not_circularised", "multiple_final_contigs", "multiple_candidate_contigs",
    "length_outside_expected_range",
}
# Row key holding the reasons that actually force manual_review. Not a summary
# column, so it is dropped by the DictWriter (extrasaction="ignore") on write.
BLOCKING_KEY = "_blocking_reason"
# Row key marking an assembly that a curation step replaced with a molecule of its
# own (currently only the concatemer collapse). Same deal: not a summary column.
SUPERSEDED_KEY = "_superseded"


def expected_gene_count_for(row: dict[str, str], thresholds: Thresholds) -> int | None:
    """The total-gene expectation that applies to THIS assembly, or None.

    --expected-gene-count is a run-level knob, but the 37-gene total is a
    vertebrate figure: cnidarians carry ~15 genes because most of their tRNAs are
    nuclear-encoded, and no invertebrate follows the vertebrate complement. Applying
    37 to them reported missing_genes on finished mitogenomes and contradicted
    annotation_stats.py, which had already passed the same assembly on the
    PCG+rRNA core. Rows judged under the core profile therefore get no total-gene
    expectation; the protein-coding-gene check (13 PCGs, true for these lineages
    too) still applies and still blocks.
    """
    if (row.get("completeness_profile") or "").strip().lower() == "core":
        return None
    return thresholds.expected_gene_count


def is_complete_core(row: dict[str, str], thresholds: Thresholds) -> bool:
    """A finished mitogenome: circular, all protein-coding genes present, length
    in the expected range, and at most TRNA_TOLERANCE tRNAs short. Soft flags on
    such an assembly are advisory, not blocking."""
    if row.get("circularised") != "true":
        return False
    length = parse_number(row.get("final_length_bp"))
    num_cds = parse_number(row.get("num_cds"))
    num_genes = parse_number(row.get("num_genes"))
    if thresholds.expected_pcg_count is not None:
        if num_cds is None or num_cds < thresholds.expected_pcg_count:
            return False
    if thresholds.min_length is not None and (length is None or length < thresholds.min_length):
        return False
    if thresholds.max_length is not None and (length is None or length > thresholds.max_length):
        return False
    expected_genes = expected_gene_count_for(row, thresholds)
    if expected_genes is not None and num_genes is not None:
        if num_genes < expected_genes - thresholds.trna_tolerance:
            return False
    return True


def blocking_reasons(row: dict[str, str], reasons: list[str], thresholds: Thresholds) -> list[str]:
    """Filter the review reasons down to those that should force manual_review.

    On a complete-core assembly, coverage / no-congeneric flags and a tRNA-only
    gene shortfall are downgraded to advisory. missing_protein_coding_genes and
    every structural flag (length, contigs, circularity, repeats) always block.
    """
    if not is_complete_core(row, thresholds):
        return list(reasons)
    num_cds = parse_number(row.get("num_cds"))
    num_genes = parse_number(row.get("num_genes"))
    expected_genes = expected_gene_count_for(row, thresholds)
    trna_only_shortfall = (
        thresholds.expected_pcg_count is not None and num_cds is not None
        and num_cds >= thresholds.expected_pcg_count
        and expected_genes is not None and num_genes is not None
        and 0 < (expected_genes - num_genes) <= thresholds.trna_tolerance
    )
    advisory = set(ADVISORY_WHEN_COMPLETE)
    if trna_only_shortfall:
        advisory.add("missing_genes")
    return [reason for reason in reasons if reason not in advisory]


def finalise_status(row: dict[str, str], failed: bool = False) -> None:
    """Terminal QC verdict. Identical for every assembler so the column can be
    sorted / filtered / counted across a mixed cohort:

      superseded     - a curation step replaced this assembly with one that has
                       its own row, so this one is provenance, not a deliverable
      failed         - the assembler errored, or produced no final assembly
      manual_review  - a final assembly exists but a blocking reason stands
      complete       - a final assembly exists and nothing blocking survived QC

    Topology is deliberately NOT encoded here (an earlier GetOrganelle-only
    "circular" status made the column incomparable); read `circularised` for it.

    `superseded` is checked first and beats every other verdict: a pre-collapse
    concatemer is over-length and non-circular by definition, so triaging it
    would mean reviewing the same molecule twice -- once as the raw concatemer
    and once as the monomer that replaced it.

    Must run after apply_qc, which is what populates BLOCKING_KEY.
    """
    if row.get(SUPERSEDED_KEY):
        row["status"] = "superseded"
        row["manual_review_reason"] = add_reason(
            row["manual_review_reason"], "superseded_by_collapse")
    elif failed or not row["final_length_bp"]:
        row["status"] = "failed"
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "failed_run")
    elif row.get(BLOCKING_KEY):
        row["status"] = "manual_review"
    else:
        row["status"] = "complete"


def parse_mitohifi_run(run: RunFiles, thresholds: Thresholds, all_files: Iterable[Path]) -> dict[str, str]:
    row = {column: MISSING for column in COLUMNS}
    row.update({"sample_id": run.sample_id, "assembly_prefix": run.prefix, "assembler": "MitoHiFi"})
    row.update(parse_mitohifi_stats(run.files))

    final_fasta = choose_final_fasta(run.files, "MitoHiFi", run.prefix)
    if final_fasta:
        length, count = parse_fasta(final_fasta)
        row["final_length_bp"] = row["final_length_bp"] or format_number(length)
        row["num_final_contigs"] = format_number(count)
    else:
        warn(f"No MitoHiFi final FASTA detected for {run.prefix}")

    potential = [path for path in run.files if path.name == "all_potential_contigs.fa"]
    if potential:
        _, count = parse_fasta(potential[0])
        row["num_candidate_contigs"] = format_number(count)

    mean_cov, cov_cv = parse_coverage(run.files)
    if mean_cov is not None:
        row["mean_coverage"] = format_number(mean_cov)
    if cov_cv is not None:
        row["coverage_cv"] = format_number(cov_cv)

    species, accession = reference_for_run(run, all_files)
    if species:
        row["reference_species"] = species
    if accession:
        row["reference_accession"] = accession

    row.update({key: value for key, value in parse_annotation_stats(all_files, run.prefix).items() if value})
    row["numt_flag"] = "true" if has_numt_signal(run.files) else "false"
    row["reference_relevance"] = reference_relevance_for_run(run)
    row["reference_divergence"] = reference_divergence_for_run(run)
    row.update(anomaly_for_run(run))
    row.update(collapse_child_evidence(run, all_files))
    apply_collapse_provenance(row, run)

    apply_qc(row, thresholds)
    finalise_status(row)
    return row


def add_reason(existing: str, reason: str) -> str:
    reasons = [item for item in existing.split(";") if item]
    if reason not in reasons:
        reasons.append(reason)
    return ";".join(reasons)


def getorganelle_evidence_from_log(files: Iterable[Path]) -> tuple[str, float | None, bool]:
    """Circularity, coverage and failure evidence from the GetOrganelle log.

    Returns (circularised, mean_coverage, failed). GetOrganelle states its verdict
    on a "Result status of animal_mt: ..." line, and in practice writes one of two
    forms: "circular genome" or "N scaffold(s)". Only the circular form is
    evidence of a closed molecule; every other parseable verdict (a scaffold
    count, "incomplete", a bare "complete genome") means GetOrganelle did not
    close the circle, so it records circularised="false" and the row picks up the
    same not_circularised flag a MitoHiFi row would. Previously a scaffold count
    matched no branch and left circularised blank, which let GetOrganelle rows
    reach "complete" on weaker evidence than any other assembler needed.

    A missing or truncated log (no result-status line at all) leaves circularised
    blank -- unknown topology, which is exactly what MitoHiFi rows with no
    contig-stats / circularity-check sidecar do.
    """
    circularised = MISSING
    mean_coverage = None
    failed = False
    for path in files:
        if not path.name.endswith(".get_org.log.txt"):
            continue
        text = read_text(path).lower()
        coverage_match = re.search(r"average [a-z_ -]*base-coverage\s*=\s*([0-9.]+)", text)
        if coverage_match:
            mean_coverage = parse_number(coverage_match.group(1))
        # Last verdict wins: GetOrganelle restates the result status as it retries
        # disentangling strategies, and only the final line reflects what it wrote.
        status_lines = [line for line in text.splitlines() if "result status" in line]
        if status_lines:
            circularised = "true" if "circular" in status_lines[-1] else "false"
        # GetOrganelle emits benign INFO-level "Disentangling failed:" messages
        # while it tries successive disentangling strategies before succeeding;
        # those are not run failures. Only genuine ERROR-level log lines (e.g.
        # "ERROR: Assembling failed.", "ERROR: No animal_mt seed reads found!")
        # or a Python traceback indicate an actual failure.
        if " - error:" in text or "traceback (most recent call last)" in text:
            failed = True
    return circularised, mean_coverage, failed


def count_getorganelle_candidates(files: Iterable[Path]) -> int | None:
    paths = [
        path
        for path in files
        if path.suffix.lower() in FASTA_EXTENSIONS
        and ("path_sequence" in path.name.lower() or ".complete.graph" in path.name.lower())
    ]
    if paths:
        return len(paths)
    selected_graphs = [path for path in files if "selected_graph" in path.name.lower()]
    if selected_graphs:
        return len(selected_graphs)
    return None


def getorganelle_graph_ambiguous(files: Iterable[Path]) -> bool:
    """True only when GetOrganelle could not resolve a *unique* mitogenome path.

    The authoritative signal is how many resolved sequences GetOrganelle actually
    wrote: it emits one ``*.path_sequence.fasta`` per equally-supported path, so
    more than one distinct path sequence (or more than one selected graph) means a
    genuine ambiguity a human must arbitrate.

    A single resolved path whose ``selected_graph.gfa`` happens to contain several
    ``S`` segments is NOT ambiguous: the standard fish mitogenome graph carries
    3-6 segments (the control-region / tRNA repeats) joined into one circular
    path. Keying off the raw segment count -- as this function previously did --
    flagged essentially every clean circular assembly for manual review, which was
    the single largest source of false-positive review flags in the audit run.
    """
    def distinct(names_iter):
        return {p.name for p in names_iter}

    path_sequences = distinct(
        path for path in files
        if path.suffix.lower() in FASTA_EXTENSIONS and "path_sequence" in path.name.lower()
    )
    if len(path_sequences) > 1:
        return True
    selected_graphs = distinct(path for path in files if "selected_graph" in path.name.lower())
    if len(selected_graphs) > 1:
        return True
    return False


def parse_getorganelle_run(run: RunFiles, thresholds: Thresholds, all_files: Iterable[Path]) -> dict[str, str]:
    row = {column: MISSING for column in COLUMNS}
    row.update({"sample_id": run.sample_id, "assembly_prefix": run.prefix, "assembler": "GetOrganelle"})

    final_fasta = choose_final_fasta(run.files, "GetOrganelle", run.prefix)
    if final_fasta:
        length, count = parse_fasta(final_fasta)
        row["final_length_bp"] = format_number(length)
        row["num_final_contigs"] = format_number(count)
    else:
        warn(f"No GetOrganelle final FASTA detected for {run.prefix}")

    candidate_count = count_getorganelle_candidates(run.files)
    if candidate_count is not None:
        row["num_candidate_contigs"] = str(candidate_count)

    circularised, mean_coverage, getorg_failed = getorganelle_evidence_from_log(run.files)
    row["circularised"] = circularised
    # A single scaffold the reference test confirms circular is recorded as
    # circularised even though GetOrganelle's log said "N scaffold(s)".
    circ_override = getorg_circular_override(run)
    if circ_override:
        row["circularised"] = circ_override
    if mean_coverage is not None:
        row["mean_coverage"] = format_number(mean_coverage)
    # The uniform remap depth supersedes GetOrganelle's own log figure wherever it
    # exists. The log value stays as the fallback for runs that predate the remap:
    # it is graph base-coverage over the reduced read set GetOrganelle selected, so
    # it is not comparable to the other assemblers the way mean_depth is.
    remap_cov, remap_cv = parse_coverage(run.files)
    if remap_cov is not None:
        row["mean_coverage"] = format_number(remap_cov)
    if remap_cv is not None:
        row["coverage_cv"] = format_number(remap_cv)

    row.update({key: value for key, value in parse_annotation_stats(all_files, run.prefix).items() if value})
    row["numt_flag"] = "true" if has_numt_signal(run.files) else "false"
    row.update(anomaly_for_run(run))
    row.update(collapse_child_evidence(run, all_files))
    apply_collapse_provenance(row, run)

    # The GetOrganelle reseed reference (relabelled <prefix>.reference.gb, or the
    # published mtdna/reference_seed/NC_*.gb) is not tied to a run by classify_file,
    # so look it up scoped to this run to avoid picking another sample's seed when
    # everything is staged into one flat summary directory.
    species, accession = reference_for_run(run, all_files)
    if species:
        row["reference_species"] = species
    if accession:
        row["reference_accession"] = accession
    row["reference_relevance"] = reference_relevance_for_run(run)

    apply_qc(row, thresholds)
    if getorganelle_graph_ambiguous(run.files):
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "ambiguous_getorganelle_graph")
        row[BLOCKING_KEY] = add_reason(row.get(BLOCKING_KEY, ""), "ambiguous_getorganelle_graph")
    finalise_status(row, failed=getorg_failed)
    return row


def oatk_circular_from_gfa(files: Iterable[Path]) -> str:
    """Oatk marks a circular contig with a self-link in its GFA (<prefix>.gfa): an
    ``L`` line whose from-segment and to-segment are the same unitig. Returns
    'true' / 'false' / '' (no GFA / unknown)."""
    for path in files:
        if path.name.endswith(".gfa"):
            saw_segment = False
            for line in read_text(path, max_chars=5_000_000).splitlines():
                if line.startswith("S\t"):
                    saw_segment = True
                elif line.startswith("L\t"):
                    parts = line.split("\t")
                    if len(parts) >= 4 and parts[1] == parts[3]:
                        return "true"
            return "false" if saw_segment else ""
    return ""


def parse_oatk_run(run: RunFiles, thresholds: Thresholds, all_files: Iterable[Path]) -> dict[str, str]:
    row = {column: MISSING for column in COLUMNS}
    row.update({"sample_id": run.sample_id, "assembly_prefix": run.prefix, "assembler": "Oatk"})

    final_fasta = choose_final_fasta(run.files, "Oatk", run.prefix)
    if final_fasta:
        length, count = parse_fasta(final_fasta)
        row["final_length_bp"] = format_number(length)
        row["num_final_contigs"] = format_number(count)
    else:
        warn(f"No Oatk final FASTA detected for {run.prefix}")

    # Circularity from the Oatk graph self-link; gene counts from the shared
    # annotation stats (the Oatk contig is annotated by the same EMMA/MITOS2 path).
    row["circularised"] = oatk_circular_from_gfa(run.files)
    # Oatk emitted no coverage of any kind before the uniform remap, so this is the
    # first time an oatk row can carry one. It is the same measurement the other two
    # assemblers now get, so the QC thresholds apply to it on equal terms.
    mean_cov, cov_cv = parse_coverage(run.files)
    if mean_cov is not None:
        row["mean_coverage"] = format_number(mean_cov)
    if cov_cv is not None:
        row["coverage_cv"] = format_number(cov_cv)
    row.update({key: value for key, value in parse_annotation_stats(all_files, run.prefix).items() if value})
    row["numt_flag"] = "false"
    # Oatk is reference-free: no reference species/accession/divergence, so those
    # columns and the no_congeneric_reference gate stay empty by construction.
    row.update(anomaly_for_run(run))
    row.update(collapse_child_evidence(run, all_files))
    apply_collapse_provenance(row, run)

    apply_qc(row, thresholds)
    finalise_status(row)
    return row


def write_rows(rows: list[dict[str, str]], output: Path) -> None:
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, MISSING) for column in COLUMNS})


def build_summary(inputs: list[Path], thresholds: Thresholds) -> list[dict[str, str]]:
    discover_sample_dirs(inputs)
    runs = discover_assembler_runs(inputs)
    if not runs:
        warn("No MitoHiFi or GetOrganelle assembly outputs were detected")
    all_files = collect_input_files(inputs)
    rows = []
    for run in runs:
        if run.assembler == "MitoHiFi":
            rows.append(parse_mitohifi_run(run, thresholds, all_files))
        elif run.assembler == "GetOrganelle":
            rows.append(parse_getorganelle_run(run, thresholds, all_files))
        elif run.assembler == "Oatk":
            rows.append(parse_oatk_run(run, thresholds, all_files))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", type=Path, default=[Path(".")], help="Files or directories to scan")
    parser.add_argument("--output", type=Path, default=Path("mitogenome_assembly_summary_mqc.tsv"))
    parser.add_argument("--min-mean-coverage", type=float, default=None)
    parser.add_argument("--max-coverage-cv", type=float, default=None)
    parser.add_argument("--min-length", type=int, default=None)
    parser.add_argument("--max-length", type=int, default=None)
    parser.add_argument("--expected-gene-count", type=int, default=None)
    parser.add_argument("--trna-tolerance", dest="trna_tolerance", type=int,
                        default=TRNA_TOLERANCE,
                        help="tRNA-only shortfall tolerated on an otherwise complete "
                             "assembly before missing_genes becomes a blocking reason. "
                             "Must match annotation_stats.py --trna-tolerance (both are "
                             "wired to params.annotation_trna_tolerance).")
    parser.add_argument("--expected-pcg-count", type=int, default=None,
                        help="Flag missing_protein_coding_genes when the CDS count is "
                             "below this (e.g. 13 for a vertebrate mitogenome).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = Thresholds(
        min_mean_coverage=args.min_mean_coverage,
        max_coverage_cv=args.max_coverage_cv,
        min_length=args.min_length,
        max_length=args.max_length,
        expected_gene_count=args.expected_gene_count,
        expected_pcg_count=args.expected_pcg_count,
        trna_tolerance=args.trna_tolerance,
    )
    rows = build_summary(args.input, thresholds)
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
