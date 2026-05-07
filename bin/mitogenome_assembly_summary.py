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
    "sample_id",
    "assembler",
    "status",
    "final_length_bp",
    "circularised",
    "num_candidate_contigs",
    "num_final_contigs",
    "num_genes",
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
        ".coverage.tsv",
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
    return strip_known_suffix(name)


def classify_file(path: Path) -> str | None:
    name = path.name.lower()
    parts = {part.lower() for part in path.parts}
    ignored_parts = {"emma", "lca", "genbank", "cds", "proteins", "genes"}
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
    if parts.intersection(ignored_parts) or name.startswith(ignored_prefixes):
        return None
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


def discover_assembler_runs(inputs: Iterable[Path]) -> list[RunFiles]:
    runs: dict[tuple[str, str], RunFiles] = {}
    files = collect_input_files(inputs)

    for path in files:
        assembler = classify_file(path)
        if not assembler:
            continue
        prefix = infer_prefix(path, assembler)
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
        else:
            if path.name == f"{prefix}.fasta":
                value += 50
            if "mitohifi" in name:
                value += 20
        return (-value, len(path.name), path.name)

    return sorted(candidates, key=score)[0]


def parse_annotation_stats(files: Iterable[Path], prefix: str, sample_id: str) -> dict[str, str]:
    rows = []
    for path in files:
        if path.name.endswith(".annotation_stats.csv"):
            rows.extend(read_table(path))

    for row in rows:
        row_sample = first_value(row, ["sample", "sample_id", "og_id"])
        row_code = first_value(row, ["code", "assembly", "assembler", "mt_assembly_prefix"])
        if row_sample and row_sample != sample_id and row_sample not in prefix:
            continue
        if row_code and row_code not in prefix:
            continue

        num_cds = first_numeric(row, ["num_cds"])
        num_trna = first_numeric(row, ["num_trna"])
        num_rrna = first_numeric(row, ["num_rrna"])
        num_genes = first_numeric(row, ["num_genes", "gene_count"])
        if num_genes is None and None not in (num_cds, num_trna, num_rrna):
            num_genes = (num_cds or 0) + (num_trna or 0) + (num_rrna or 0)

        return {
            "num_genes": format_number(num_genes),
            "missing_genes": first_value(row, ["missing_genes", "num_missing"]),
            "frameshift_flag": parse_bool(first_value(row, ["frameshift", "frameshift_flag", "frameshifts"])),
        }
    return {}


def parse_coverage(files: Iterable[Path]) -> tuple[float | None, float | None]:
    for path in files:
        name = path.name.lower()
        if not (name.endswith(".coverage.tsv") or "coverage" in name and path.suffix.lower() in {".tsv", ".txt", ".csv"}):
            continue
        rows = read_table(path)
        if rows:
            mean = first_numeric(rows[0], ["mean_coverage", "avg_coverage", "average_coverage", "coverage"])
            cv = first_numeric(rows[0], ["coverage_cv", "cv"])
            if mean is not None or cv is not None:
                return mean, cv

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
            continue
        if depths:
            mean = statistics.fmean(depths)
            cv = statistics.pstdev(depths) / mean if mean else None
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


def find_reference_genbank(files: Iterable[Path]) -> tuple[str, str]:
    genbanks = [path for path in files if path.suffix.lower() in {".gb", ".gbk", ".gbf"}]

    def score(path: Path) -> tuple[int, str]:
        lowered = [part.lower() for part in path.parts]
        name = path.name.lower()
        value = 0
        if "mitoreference" in lowered:
            value -= 100
        if name.startswith(("nc_", "nw_", "nz_")):
            value -= 50
        if "final" in name or "mitohifi" in name:
            value += 20
        return value, str(path)

    for path in sorted(genbanks, key=score):
        species, accession = parse_genbank_reference(path)
        if not is_missing(species) and not is_missing(accession):
            return species, accession
    return MISSING, MISSING


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
    genes = first_numeric(final_row, ["num_genes", "gene_count", "genes_count", "genes", "number_of_genes"])
    circular = parse_bool(first_value(final_row, ["circularised", "circularized", "circular", "is_circular", "was_circular"]))
    frameshift = parse_bool(first_value(final_row, ["frameshift_flag", "frameshift", "frameshifts", "frameshifts_found"]))
    mean_cov = first_numeric(final_row, ["mean_coverage", "avg_coverage", "average_coverage"])

    if length is not None:
        stats["final_length_bp"] = format_number(length)
    if genes is not None:
        stats["num_genes"] = format_number(genes)
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
    patterns = [
        "numt",
        "nuclear mitochondrial",
        "rejected",
        "low-confidence",
        "low confidence",
        "low_coverage",
        "divergent",
        "contigs_filtering",
    ]
    for path in files:
        text = str(path).lower()
        if path.suffix.lower() in {".tsv", ".txt", ".log"}:
            text += "\n" + read_text(path, max_chars=200_000).lower()
        if any(pattern in text for pattern in patterns):
            return True
    return False


def apply_qc(row: dict[str, str], thresholds: Thresholds) -> None:
    reasons = []

    length = parse_number(row.get("final_length_bp"))
    mean_cov = parse_number(row.get("mean_coverage"))
    cov_cv = parse_number(row.get("coverage_cv"))
    num_candidate = parse_number(row.get("num_candidate_contigs"))
    num_final = parse_number(row.get("num_final_contigs"))
    num_genes = parse_number(row.get("num_genes"))

    if not row.get("final_length_bp"):
        reasons.append("missing_final_fasta")
    if row.get("circularised") == "false":
        reasons.append("not_circularised")
    if num_candidate is not None and num_candidate > 1:
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
    if thresholds.expected_gene_count is not None and num_genes is not None and num_genes < thresholds.expected_gene_count:
        reasons.append("missing_genes")

    deduped = []
    for reason in reasons:
        if reason not in deduped:
            deduped.append(reason)
    row["manual_review_reason"] = ";".join(deduped)


def parse_mitohifi_run(run: RunFiles, thresholds: Thresholds, all_files: Iterable[Path]) -> dict[str, str]:
    row = {column: MISSING for column in COLUMNS}
    row.update({"sample_id": run.sample_id, "assembler": "MitoHiFi"})
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

    species, accession = find_reference_genbank([*run.files, *all_files])
    if species:
        row["reference_species"] = species
    if accession:
        row["reference_accession"] = accession

    row.update({key: value for key, value in parse_annotation_stats(all_files, run.prefix, run.sample_id).items() if value})
    row["numt_flag"] = "true" if has_numt_signal(run.files) else "false"

    apply_qc(row, thresholds)
    if not row["final_length_bp"]:
        row["status"] = "failed"
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "failed_run")
    elif row["manual_review_reason"]:
        row["status"] = "manual_review"
    else:
        row["status"] = "complete"
    return row


def add_reason(existing: str, reason: str) -> str:
    reasons = [item for item in existing.split(";") if item]
    if reason not in reasons:
        reasons.append(reason)
    return ";".join(reasons)


def getorganelle_status_from_log(files: Iterable[Path]) -> tuple[str, str, float | None]:
    status = "unknown"
    circularised = MISSING
    mean_coverage = None
    for path in files:
        if not path.name.endswith(".get_org.log.txt"):
            continue
        text = read_text(path).lower()
        coverage_match = re.search(r"average [a-z_ -]*base-coverage\s*=\s*([0-9.]+)", text)
        if coverage_match:
            mean_coverage = parse_number(coverage_match.group(1))
        if "result status" in text:
            status_lines = [line for line in text.splitlines() if "result status" in line]
            status_text = " ".join(status_lines)
            if "circular" in status_text:
                status = "circular"
                circularised = "true"
            elif "complete" in status_text:
                status = "complete"
            elif "incomplete" in status_text:
                status = "incomplete"
                circularised = "false"
        if "traceback" in text or "error:" in text or "failed" in text:
            status = "failed"
    return status, circularised, mean_coverage


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
    selected_graphs = [path for path in files if "selected_graph" in path.name.lower()]
    if len(selected_graphs) > 1:
        return True
    for path in selected_graphs:
        text = read_text(path, max_chars=300_000)
        segment_count = sum(1 for line in text.splitlines() if line.startswith("S\t"))
        if segment_count > 1:
            return True
    for path in files:
        if path.name.endswith(".get_org.log.txt"):
            text = read_text(path, max_chars=500_000).lower()
            if any(term in text for term in ("ambiguous", "multiple path", "path2", "more than one")):
                return True
    return False


def parse_getorganelle_run(run: RunFiles, thresholds: Thresholds, all_files: Iterable[Path]) -> dict[str, str]:
    row = {column: MISSING for column in COLUMNS}
    row.update({"sample_id": run.sample_id, "assembler": "GetOrganelle"})

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

    status, circularised, mean_coverage = getorganelle_status_from_log(run.files)
    row["status"] = status
    row["circularised"] = circularised
    if mean_coverage is not None:
        row["mean_coverage"] = format_number(mean_coverage)

    row.update({key: value for key, value in parse_annotation_stats(all_files, run.prefix, run.sample_id).items() if value})
    row["numt_flag"] = "true" if has_numt_signal(run.files) else "false"

    apply_qc(row, thresholds)
    if getorganelle_graph_ambiguous(run.files):
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "ambiguous_getorganelle_graph")
    if not row["final_length_bp"]:
        row["status"] = "failed"
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "failed_run")
    elif row["status"] == "failed":
        row["manual_review_reason"] = add_reason(row["manual_review_reason"], "failed_run")
    elif row["status"] in {"incomplete", "unknown"} and row["manual_review_reason"]:
        row["status"] = "manual_review"
    elif row["manual_review_reason"]:
        row["status"] = "manual_review"
    elif row["status"] == "unknown":
        row["status"] = "complete"
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = Thresholds(
        min_mean_coverage=args.min_mean_coverage,
        max_coverage_cv=args.max_coverage_cv,
        min_length=args.min_length,
        max_length=args.max_length,
        expected_gene_count=args.expected_gene_count,
    )
    rows = build_summary(args.input, thresholds)
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
