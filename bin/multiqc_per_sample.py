#!/usr/bin/env python3

"""Prepare and run one MultiQC report per detected sample."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shlex
import shutil
import subprocess
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - exercised on minimal login Python installs
    yaml = None


SUMMARY_NAME = "mitogenome_assembly_summary_mqc.tsv"
UPLOAD_SUMMARY_NAME = "upload_results_summary_mqc.tsv"
# Global MultiQC tables that carry one row per assembly/sample and must be
# filtered down to the current target before they enter a per-sample report.
FILTERED_SUMMARY_NAMES = {SUMMARY_NAME, UPLOAD_SUMMARY_NAME}
# Column names that hold the sample/OG id across the different summary tables.
SAMPLE_ID_FIELDS = ("sample_id", "og_id")
SAMPLE_RE = re.compile(r"(?<![A-Za-z0-9_-])(OG[0-9][A-Za-z0-9_-]*)(?![A-Za-z0-9_-])")
ASSEMBLY_RE = re.compile(
    r"(?<![A-Za-z0-9_-])"
    r"(?P<prefix>(?P<sample>OG[0-9][A-Za-z0-9_-]*)\.[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\."
    r"(?:getorg[0-9A-Za-z_-]+|(?:v[0-9A-Za-z_-]+)?mitohifi[0-9A-Za-z_-]*))"
    r"(?:\.emma[0-9A-Za-z_-]+)?"
)
TEXT_SUFFIXES = {
    ".csv",
    ".fa",
    ".fasta",
    ".fastg",
    ".fna",
    ".gfa",
    ".html",
    ".json",
    ".log",
    ".txt",
    ".tsv",
    ".yaml",
    ".yml",
}
COMMON_NAMES = {
    "workflow_summary_mqc.yaml",
    "methods_description_mqc.yaml",
    "software_versions_mqc_versions.yml",
    "software_versions_mqc_versions.yaml",
}
ASSEMBLY_SUFFIXES = (
    ".contigs_stats.with_coverage.tsv",
    ".contigs_stats.tsv",
    ".coverage.tsv",
    ".annotation_stats.csv",
    ".qc_summary.tsv",
    ".hifiasm.log",
    ".log",
    ".gb",
    ".gff",
    ".tbl",
    ".fasta",
    ".fa",
    ".fna",
)

TOOL_TO_GROUP = {
    "getorganelleconfig": "GETORGANELLE_CONFIG",
    "getorganellefromreads": "GETORGANELLE_FROMREADS",
    "catfastq": "CAT_FASTQ",
    "mitohififindreference": "MITOHIFI_FINDMITOREFERENCE",
    "mitohifi": "MITOHIFI_MITOHIFI",
    "mitohifiaveragecoverage": "MITOHIFI_AVERAGE_COVERAGE",
    "emma": "EMMA",
    "blastblastn": "BLAST_BLASTN",
    "lca": "LCA",
    "pushmtdnaassemblyresults": "PUSH_MTDNA_ASSM_RESULTS",
    "speciesvalidation": "SPECIES_VALIDATION",
    "pushmtdnaannotationresults": "PUSH_MTDNA_ANNOTATION_RESULTS",
    "pushlcablastresults": "PUSH_LCA_BLAST_RESULTS",
    "evaluateqcconditions": "EVALUATE_QC_CONDITIONS",
    "formatfiles": "FORMAT_FILES",
    "buildsourcemodifiers": "BUILD_SOURCE_MODIFIERS",
    "extractgenesgff": "EXTRACT_GENES_GFF",
    "translategenes": "TRANSLATE_GENES",
    "table2asn": "GEN_FILES_TABLE2ASN",
}

TOOL_ORDER = {
    "GETORGANELLE_CONFIG": 10,
    "GETORGANELLE_FROMREADS": 20,
    "CAT_FASTQ": 30,
    "MITOHIFI_FINDMITOREFERENCE": 40,
    "MITOHIFI_MITOHIFI": 50,
    "MITOHIFI_AVERAGE_COVERAGE": 60,
    "EMMA": 70,
    "BLAST_BLASTN": 80,
    "LCA": 90,
    "PUSH_MTDNA_ASSM_RESULTS": 100,
    "SPECIES_VALIDATION": 110,
    "PUSH_MTDNA_ANNOTATION_RESULTS": 120,
    "PUSH_LCA_BLAST_RESULTS": 130,
    "EVALUATE_QC_CONDITIONS": 140,
    "FORMAT_FILES": 150,
    "BUILD_SOURCE_MODIFIERS": 160,
    "EXTRACT_GENES_GFF": 170,
    "TRANSLATE_GENES": 180,
    "GEN_FILES_TABLE2ASN": 190,
}

ROW_PATTERN = re.compile(
    r"<tr><td>(?P<tool>.*?)</td><td>(?P<params>.*?)</td><td>(?P<notes>.*?)</td></tr>",
    re.IGNORECASE | re.DOTALL,
)


@dataclass(frozen=True)
class AssemblyTarget:
    sample: str
    prefix: str


def warn(message: str) -> None:
    print(f"[multiqc_per_sample] WARNING: {message}", file=sys.stderr)


def to_builtin(obj):
    if isinstance(obj, dict):
        return {key: to_builtin(value) for key, value in obj.items()}
    return obj


def simplify_group(group_name: str) -> str:
    group_name = str(group_name).strip().strip("'").strip('"')
    return group_name.rsplit(":", 1)[-1]


def strip_yaml_scalar(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def load_simple_versions_yaml(text: str) -> OrderedDict:
    """Parse the simple group/tool/version YAML emitted by nf-core version collation."""
    parsed = OrderedDict()
    current_group = None
    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        if raw_line[:1].isspace():
            if current_group is None or ":" not in raw_line:
                continue
            tool_name, version = raw_line.strip().split(":", 1)
            parsed[current_group][strip_yaml_scalar(tool_name)] = strip_yaml_scalar(version)
            continue
        if raw_line.rstrip().endswith(":"):
            current_group = simplify_group(raw_line.rstrip()[:-1])
            parsed[current_group] = OrderedDict()
    return parsed


def dump_versions_yaml(data: OrderedDict) -> str:
    if yaml is not None:
        return yaml.safe_dump(to_builtin(data), sort_keys=False)

    lines = []
    for group_name, tools in data.items():
        lines.append(f"{json.dumps(str(group_name))}:")
        for tool_name, version in tools.items():
            lines.append(f"  {json.dumps(str(tool_name))}: {json.dumps(str(version))}")
    return "\n".join(lines) + "\n"


def load_version_file(path: Path) -> OrderedDict:
    try:
        text = path.read_text()
    except Exception:
        return OrderedDict()

    if yaml is None:
        return load_simple_versions_yaml(text)

    try:
        data = yaml.safe_load(text) or {}
    except Exception:
        return OrderedDict()

    parsed = OrderedDict()
    if not isinstance(data, dict):
        return parsed

    for group_name, tools in data.items():
        simple_group = simplify_group(group_name)
        if not isinstance(tools, dict):
            continue
        parsed[simple_group] = OrderedDict(
            (str(tool_name), str(version).strip())
            for tool_name, version in tools.items()
        )

    return parsed


def load_multiqc_args(path: Path | None) -> list[str]:
    if not path:
        return []

    try:
        chunks = json.loads(path.read_text())
    except Exception as exc:
        raise SystemExit(f"Could not read MultiQC argument chunks from {path}: {exc}") from exc

    args: list[str] = []
    for chunk in chunks:
        if chunk:
            args.extend(shlex.split(str(chunk)))
    return args


def text_for_matching(path: Path) -> str:
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return ""
    try:
        if path.stat().st_size > 5_000_000:
            return ""
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def assembly_prefixes_from_text(text: str) -> set[AssemblyTarget]:
    targets = set()
    for match in ASSEMBLY_RE.finditer(text):
        prefix = match.group("prefix")
        targets.add(AssemblyTarget(match.group("sample"), prefix))
    return targets


def summary_targets(path: Path) -> tuple[set[AssemblyTarget], set[str]]:
    targets = set()
    samples = set()
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample_id = (row.get("sample_id") or "").strip()
                assembly_prefix = (
                    row.get("assembly_prefix")
                    or row.get("mt_assembly_prefix")
                    or row.get("prefix")
                    or ""
                ).strip()
                if sample_id:
                    samples.add(sample_id)
                if sample_id and assembly_prefix:
                    targets.add(AssemblyTarget(sample_id, assembly_prefix))
    except Exception as exc:
        warn(f"could not parse {path}: {exc}")
    return targets, samples


def discover_assembly_targets(files: list[Path]) -> list[AssemblyTarget]:
    targets = set()
    samples = set()
    for path in files:
        if path.name == SUMMARY_NAME:
            summary_file_targets, summary_samples = summary_targets(path)
            targets.update(summary_file_targets)
            samples.update(summary_samples)

        targets.update(assembly_prefixes_from_text(str(path)))
        targets.update(assembly_prefixes_from_text(text_for_matching(path)))
        for match in SAMPLE_RE.finditer(str(path)):
            samples.add(match.group(1))

    if targets:
        return sorted(targets, key=lambda target: (target.sample, target.prefix))

    return [AssemblyTarget(sample, sample) for sample in sorted(samples)]


def text_contains_sample(path: Path, sample: str) -> bool:
    return sample in text_for_matching(path)


def text_contains_assembly(path: Path, target: AssemblyTarget) -> bool:
    return target.prefix in text_for_matching(path)


def has_any_sample_id(path: Path) -> bool:
    if SAMPLE_RE.search(str(path)):
        return True
    return SAMPLE_RE.search(text_for_matching(path)) is not None


def has_any_assembly_prefix(path: Path) -> bool:
    if ASSEMBLY_RE.search(str(path)):
        return True
    return ASSEMBLY_RE.search(text_for_matching(path)) is not None


def is_common_file(path: Path) -> bool:
    return path.name in COMMON_NAMES or path.name.endswith(("_mqc_versions.yml", "_mqc_versions.yaml"))


def belongs_to_sample(path: Path, sample: str) -> bool:
    if path.name in FILTERED_SUMMARY_NAMES:
        return True
    if is_common_file(path):
        return True
    if path.name.endswith(".tool_params_mqcrow.html"):
        return sample in str(path) or text_contains_sample(path, sample) or not has_any_sample_id(path)
    return sample in str(path) or text_contains_sample(path, sample)


def belongs_to_assembly(path: Path, target: AssemblyTarget) -> bool:
    if path.name in FILTERED_SUMMARY_NAMES:
        return True
    if is_common_file(path):
        return True
    if target.prefix in str(path):
        return True
    if text_contains_assembly(path, target):
        return True
    if path.name.endswith(".tool_params_mqcrow.html"):
        return (
            target.sample in str(path)
            or text_contains_sample(path, target.sample)
            or not has_any_sample_id(path)
        ) and not has_any_assembly_prefix(path)
    return (target.sample in str(path) or text_contains_sample(path, target.sample)) and not has_any_assembly_prefix(path)


def safe_name(sample: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", sample).strip("._")
    return cleaned or "sample"


def strip_assembly_suffix(name: str) -> str:
    for suffix in ASSEMBLY_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def infer_assembly_prefix(files: list[Path], sample: str) -> str:
    candidates: set[str] = set()
    for path in files:
        name = path.name
        if not name.startswith(f"{sample}."):
            continue
        candidate = strip_assembly_suffix(name)
        if candidate and candidate != sample:
            candidates.add(candidate)

    if not candidates:
        return sample

    def score(candidate: str) -> tuple[int, int, str]:
        lower = candidate.lower()
        value = 0
        if "mitohifi" in lower:
            value += 50
        if "getorg" in lower:
            value += 40
        if ".emma" in lower:
            value -= 20
        return (-value, len(candidate), candidate)

    return sorted(candidates, key=score)[0]


def unique_destination(directory: Path, source_name: str) -> Path:
    candidate = directory / source_name
    if not candidate.exists():
        return candidate
    stem = candidate.stem
    suffix = candidate.suffix
    counter = 2
    while True:
        candidate = directory / f"{stem}.{counter}{suffix}"
        if not candidate.exists():
            return candidate
        counter += 1


def write_filtered_summary(source: Path, destination: Path, target: AssemblyTarget) -> bool:
    try:
        with source.open(newline="") as in_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            sample_field = next((field for field in SAMPLE_ID_FIELDS if field in fieldnames), None)
            has_assembly_prefix = "assembly_prefix" in fieldnames
            rows = [
                row
                for row in reader
                if (sample_field is None or (row.get(sample_field) or "").strip() == target.sample)
                and (
                    not has_assembly_prefix
                    or (row.get("assembly_prefix") or "").strip() == target.prefix
                )
            ]
    except Exception as exc:
        warn(f"could not filter {source} for {target.prefix}: {exc}")
        return False

    with destination.open("w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return True


def collect_versions(sample_dir: Path) -> OrderedDict:
    for filename in ("software_versions_mqc_versions.yml", "software_versions_mqc_versions.yaml"):
        for path in sample_dir.rglob(filename):
            path.unlink()

    version_files = []
    for pattern in (
        "software_versions.yml",
        "software_versions.yaml",
        "versions.yml",
        "versions.yaml",
        "versions_*.yml",
        "versions_*.yaml",
        "*_mqc_versions.yml",
        "*_mqc_versions.yaml",
    ):
        version_files.extend(
            sorted(
                path for path in sample_dir.rglob(pattern)
                if path.name not in {"software_versions_mqc_versions.yml", "software_versions_mqc_versions.yaml"}
            )
        )

    version_files = sorted(
        set(version_files),
        key=lambda path: (0 if path.name.startswith("software_versions") or "_software_mqc_versions" in path.name else 1, str(path)),
    )

    versions_by_group = OrderedDict()
    for version_file in version_files:
        for group_name, tools in load_version_file(version_file).items():
            merged_tools = versions_by_group.setdefault(group_name, OrderedDict())
            for tool_name, version in tools.items():
                merged_tools.setdefault(tool_name, version)

    if versions_by_group:
        (sample_dir / "software_versions_mqc_versions.yml").write_text(
            dump_versions_yaml(versions_by_group),
            encoding="utf-8",
        )

    return versions_by_group


def prepare_tool_parameters(sample_dir: Path, versions_by_group: OrderedDict) -> None:
    tool_rows = []
    seen_rows = set()
    for row_file in sorted(sample_dir.rglob("*.tool_params_mqcrow.html")):
        match = ROW_PATTERN.search(row_file.read_text(errors="ignore").strip())
        if not match:
            continue

        tool_html = match.group("tool").strip()
        params_html = match.group("params").strip()
        notes_html = match.group("notes").strip()
        tool_key = re.sub(r"[^a-z0-9]+", "", html.unescape(tool_html).lower())

        version_html = "NA"
        group_name = TOOL_TO_GROUP.get(tool_key)
        if group_name and group_name in versions_by_group:
            group_versions = versions_by_group[group_name]
            if len(group_versions) == 1:
                version_html = f"<samp>{html.escape(next(iter(group_versions.values())))}</samp>"
            else:
                formatted_versions = [
                    f"{html.escape(tool_name)} {html.escape(version)}"
                    for tool_name, version in group_versions.items()
                ]
                version_html = f"<samp>{'; '.join(formatted_versions)}</samp>"

        row_html = f"<tr><td>{tool_html}</td><td>{version_html}</td><td>{params_html}</td><td>{notes_html}</td></tr>"
        row_key = (tool_html, params_html, notes_html)
        if row_key in seen_rows:
            continue
        seen_rows.add(row_key)
        tool_rows.append((TOOL_ORDER.get(group_name, 999), html.unescape(tool_html).lower(), row_html))

    if not tool_rows:
        return

    tool_params_lines = [
        "id: 'nf-core-oceangenomesmitogenomes-tool-parameters'",
        "description: 'Exact tool parameters captured from the module command scripts.'",
        "section_name: 'Tool Parameters Used'",
        "plot_type: 'html'",
        "data: |",
        '    <table class="table table-condensed">',
        "    <thead><tr><th>Tool</th><th>Version</th><th>Effective Parameters</th><th>Notes</th></tr></thead>",
        "    <tbody>",
    ]
    tool_params_lines.extend(
        f"    {row}"
        for _order, _tool_name, row in sorted(tool_rows, key=lambda item: (item[0], item[1]))
    )
    tool_params_lines.extend(["    </tbody>", "    </table>"])
    (sample_dir / "tool_parameters_mqc.yaml").write_text("\n".join(tool_params_lines) + "\n", encoding="utf-8")


def prepare_multiqc_custom_content(sample_dir: Path) -> None:
    prepare_tool_parameters(sample_dir, collect_versions(sample_dir))


def prepare_sample_inputs(
    input_root: Path,
    sample_input_root: Path,
    output_root: Path,
    html_output_root: Path | None = None,
) -> list[tuple[str, str, Path, Path, Path]]:
    files = sorted(path for path in input_root.rglob("*") if path.is_file())
    targets = discover_assembly_targets(files)
    output_root.mkdir(parents=True, exist_ok=True)
    sample_input_root.mkdir(parents=True, exist_ok=True)

    if not targets:
        warn("no sample or assembly IDs found in MultiQC inputs")
        (output_root / "NO_SAMPLES_FOUND.txt").write_text(
            "No sample or assembly IDs were found in the MultiQC inputs.\n",
            encoding="utf-8",
        )
        return []

    prepared = []
    for target in targets:
        sample_label = safe_name(target.sample)
        assembly_label = safe_name(target.prefix)
        sample_input_dir = sample_input_root / sample_label / assembly_label
        sample_files = [source for source in files if belongs_to_assembly(source, target)]
        sample_report_dir = output_root / sample_label / assembly_label
        sample_html_dir = (html_output_root / sample_label / assembly_label) if html_output_root else sample_report_dir
        sample_input_dir.mkdir(parents=True, exist_ok=True)
        sample_report_dir.mkdir(parents=True, exist_ok=True)
        sample_html_dir.mkdir(parents=True, exist_ok=True)

        matched = 0
        for source in sample_files:
            destination = unique_destination(sample_input_dir, source.name)
            if source.name in FILTERED_SUMMARY_NAMES:
                if not write_filtered_summary(source, destination, target):
                    continue
            else:
                shutil.copy2(source, destination)
            matched += 1

        if matched == 0:
            warn(f"no MultiQC inputs matched assembly {target.prefix}")
            continue

        prepare_multiqc_custom_content(sample_input_dir)
        prepared.append((sample_label, assembly_label, sample_input_dir, sample_report_dir, sample_html_dir))

    return prepared


def run_multiqc(prepared_samples: list[tuple[str, str, Path, Path, Path]], multiqc_args: list[str]) -> None:
    for _sample_label, assembly_label, sample_input_dir, sample_report_dir, sample_html_dir in prepared_samples:
        report_name = f"{assembly_label}_multiqc_report.html"
        command = [
            "multiqc",
            "--force",
            *multiqc_args,
            "--filename",
            report_name,
            "--outdir",
            str(sample_report_dir),
            str(sample_input_dir),
        ]
        subprocess.run(command, check=True)
        shutil.copy2(sample_report_dir / report_name, sample_html_dir / report_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("multiqc_per_sample_inputs"))
    parser.add_argument("--sample-input-root", type=Path, default=Path("per_sample_inputs"))
    parser.add_argument("--output-root", type=Path, default=Path("per_sample"))
    parser.add_argument("--html-output-root", type=Path)
    parser.add_argument("--multiqc-arg-chunks-file", type=Path)
    parser.add_argument("--prepare-only", action="store_true", help="Prepare per-sample inputs without running MultiQC.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    prepared = prepare_sample_inputs(args.input_root, args.sample_input_root, args.output_root, args.html_output_root)
    if not args.prepare_only:
        run_multiqc(prepared, load_multiqc_args(args.multiqc_arg_chunks_file))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
