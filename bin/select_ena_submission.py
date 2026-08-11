#!/usr/bin/env python3
"""Compare finished ENA candidate packages and record selections per technology.

Every technology with a viable mitogenome is published, each to its own ENA
child study, so selection runs within a technology rather than across them: a
specimen can appear once as HiFi, once as Hi-C and once as Illumina, but never
twice within one of those.  A SELECTED status therefore means selected within
that technology, not the single record for the specimen.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import glob
import hashlib
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import psycopg2
except ImportError:
    psycopg2 = None

ROOT = Path(__file__).resolve().parent
PACKAGE_SPEC = importlib.util.spec_from_file_location(
    "ena_package_for_selection", ROOT / "ena_package.py"
)
ENA_PACKAGE = importlib.util.module_from_spec(PACKAGE_SPEC)
PACKAGE_SPEC.loader.exec_module(ENA_PACKAGE)


# Grouping vocabulary, not a ranking: no technology outranks another now that
# all three are published. The order is only used to keep report output stable.
TECHNOLOGIES = ("hifi", "hic", "ilmn")
LOCUS_TAG_PATTERN = re.compile(r"^[A-Z][A-Z0-9]{2,11}_([0-9]{6})([0-9]{3})$")
REPORT_FIELDS = [
    "og_id",
    "tech",
    "full_seqid",
    "package_path",
    "metadata_path",
    "package_status",
    "local_validation_status",
    "normalised_circular_sha256",
    "selection_status",
    "selected",
    "selection_reason",
]
SELECTED_FIELDS = [
    "og_id",
    "tech",
    "full_seqid",
    "package_path",
    "metadata_path",
    "study",
    "biosample_accession",
    "selection_reason",
]


def expand_paths(patterns: list[str]) -> list[Path]:
    found: set[Path] = set()
    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        if not matches and Path(pattern).is_file():
            matches = [pattern]
        found.update(Path(match).resolve() for match in matches)
    return sorted(found)


def load_packages(paths: list[Path]) -> list[dict[str, object]]:
    packages: list[dict[str, object]] = []
    seen: set[str] = set()
    for path in paths:
        data = json.loads(path.read_text())
        required = {
            "full_seqid",
            "og_id",
            "normalised_circular_sha256",
            "package_status",
            "local_validation_status",
        }
        missing = required - data.keys()
        if missing:
            raise ValueError(f"{path} lacks fields: {', '.join(sorted(missing))}")
        seqid = str(data["full_seqid"])
        if seqid in seen:
            raise ValueError(f"Duplicate package metadata for {seqid}")
        seen.add(seqid)
        data["_metadata_path"] = str(path)
        # Under Nextflow the metadata JSON is staged flat into a task work
        # directory, so path.parent is ephemeral and useless to anything that
        # reads package_path back out of the database later.  build recorded
        # where publishDir would put the package; trust that when present.
        published = str(data.get("published_package_path") or "").strip()
        if published:
            data["_package_path"] = published
        else:
            data["_package_path"] = str(path.parent)
            sys.stderr.write(
                f"WARNING: {path} has no published_package_path; recording the staged "
                f"location {path.parent}, which will not survive work-directory cleanup.\n"
            )
        for service in ("test", "production"):
            status_paths = sorted(
                Path(data["_package_path"]).parent.glob(
                    f"validation/webin_{service}/{seqid}.webin_{service}_status.tsv"
                )
            )
            if status_paths:
                with status_paths[-1].open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                if len(rows) == 1 and rows[0].get("status"):
                    data[f"webin_{service}_status"] = rows[0]["status"].strip()
        packages.append(data)
    return packages


def eligible(package: dict[str, object]) -> bool:
    return (
        package.get("package_status") == "READY"
        and package.get("local_validation_status") == "PASS"
    )


def tech_of(full_seqid: str) -> str:
    """Technology from position 1 of the SeqID: OG<n>.<tech>.<date>.<code>...

    Strict, because the technology chooses the ENA study: an unrecognised value
    would otherwise be grouped and submitted as if it were a fourth technology.
    """
    parts = str(full_seqid).split(".")
    tech = parts[1] if len(parts) > 1 else ""
    if tech not in TECHNOLOGIES:
        raise ValueError(
            f"{full_seqid} has technology {tech!r}, expected one of {', '.join(TECHNOLOGIES)}"
        )
    return tech


def group_key(package: dict[str, object]) -> tuple[str, str]:
    return str(package["og_id"]), tech_of(str(package["full_seqid"]))


def group_sort_key(key: tuple[str, str]) -> tuple[int, int]:
    return int(key[0][2:]), TECHNOLOGIES.index(key[1])


def candidate_sort_key(package: dict[str, object]) -> tuple[object, ...]:
    """Newest date first, then SeqID. Only ever applied within one technology."""
    parts = str(package["full_seqid"]).split(".")
    date = parts[2] if len(parts) > 2 else ""
    date_rank = -int(date) if date.isdigit() else 0
    return date_rank, str(package["full_seqid"])


def select_group(packages: list[dict[str, object]]) -> dict[str, object]:
    passing = [package for package in packages if eligible(package)]
    if not passing:
        return {
            "status": "NO_PASSING_CANDIDATE",
            "selected": None,
            "reason": "no_candidate_has_ready_package_and_local_pass",
        }
    if len(passing) == 1:
        return {
            "status": "SELECTED",
            "selected": passing[0]["full_seqid"],
            "reason": "unique_passing_candidate",
        }
    hashes = {str(package["normalised_circular_sha256"]) for package in passing}
    if len(hashes) == 1:
        selected = min(passing, key=candidate_sort_key)
        return {
            "status": "SELECTED",
            "selected": selected["full_seqid"],
            "reason": "equivalent_candidates_date_seqid_tiebreak",
        }
    return {
        "status": "MANUAL_REVIEW_REQUIRED",
        "selected": None,
        "reason": "multiple_distinct_passing_sequences",
    }


def read_decisions(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    """Manual overrides, one per specimen per technology.

    The technology comes from selected_seqid, so a reviewer can override the
    HiFi and Hi-C choices for one specimen independently by writing two rows.
    """
    if path is None:
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"og_id", "selected_seqid", "reviewer", "reason"}
    decisions: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if required - row.keys():
            raise ValueError(f"Decision file lacks columns: {', '.join(sorted(required - row.keys()))}")
        if not all(row[field].strip() for field in required):
            raise ValueError(f"Incomplete manual decision for {row['og_id']}")
        seqid = row["selected_seqid"].strip()
        if not seqid.startswith(f"{row['og_id'].strip()}."):
            raise ValueError(
                f"Manual selection {seqid} does not belong to specimen {row['og_id']}"
            )
        key = (row["og_id"].strip(), tech_of(seqid))
        if key in decisions:
            raise ValueError(f"Duplicate manual decision for {key[0]} {key[1]}")
        decisions[key] = row
    return decisions


def build_report(
    packages: list[dict[str, object]], decisions: dict[tuple[str, str], dict[str, str]]
) -> tuple[list[dict[str, str]], dict[tuple[str, str], dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for package in packages:
        grouped[group_key(package)].append(package)
    results: dict[tuple[str, str], dict[str, object]] = {}
    report: list[dict[str, str]] = []
    for key in sorted(grouped, key=group_sort_key):
        og_id, tech = key
        group = grouped[key]
        result = select_group(group)
        decision = decisions.get(key)
        if decision:
            available = {str(package["full_seqid"]) for package in group if eligible(package)}
            if decision["selected_seqid"] not in available:
                raise ValueError(
                    f"Manual selection {decision['selected_seqid']} is not a passing "
                    f"{tech} package for {og_id}"
                )
            result = {
                "status": "SELECTED",
                "selected": decision["selected_seqid"],
                "reason": f"manual:{decision['reason']}",
                "selected_by": decision["reviewer"],
            }
        results[key] = result
        for package in sorted(group, key=candidate_sort_key):
            report.append(
                {
                    "og_id": og_id,
                    "tech": tech,
                    "full_seqid": str(package["full_seqid"]),
                    "package_path": str(package.get("_package_path") or ""),
                    "metadata_path": str(package.get("_metadata_path") or ""),
                    "package_status": str(package["package_status"]),
                    "local_validation_status": str(package["local_validation_status"]),
                    "normalised_circular_sha256": str(
                        package["normalised_circular_sha256"]
                    ),
                    "selection_status": str(result["status"]),
                    "selected": str(package["full_seqid"] == result["selected"]).lower(),
                    "selection_reason": str(result["reason"]),
                }
            )
    return report, results


def write_report(path: Path, rows: list[dict[str, str]]) -> str:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_selected_packages(
    path: Path,
    packages: list[dict[str, object]],
    results: dict[tuple[str, str], dict[str, object]],
) -> None:
    by_seqid = {str(package["full_seqid"]): package for package in packages}
    rows: list[dict[str, str]] = []
    for key in sorted(results, key=group_sort_key):
        og_id, tech = key
        result = results[key]
        selected = result.get("selected")
        if selected is None:
            continue
        package = by_seqid[str(selected)]
        rows.append(
            {
                "og_id": og_id,
                "tech": tech,
                "full_seqid": str(selected),
                "package_path": str(package["_package_path"]),
                "metadata_path": str(package["_metadata_path"]),
                "study": str(package.get("study") or ""),
                "biosample_accession": str(
                    package.get("biosample_accession") or ""
                ),
                "selection_reason": str(result["reason"]),
            }
        )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SELECTED_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_db_config(path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser()
    parser.read(path)
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def refresh_packages_from_database(paths: list[Path], db_config: Path) -> list[Path]:
    """Refresh BioSample/depth manifests in place without rerunning science.

    "In place" has to mean the published package, not the staged copy.  Under
    Nextflow the inputs are symlinked into a task work directory that is
    discarded at the end of the run, so refreshing those would leave the
    package the submitter actually reads still carrying the stale manifest.

    Returns the metadata paths to load, which are the refreshed ones.
    """
    if psycopg2 is None:
        raise RuntimeError("psycopg2 is required to refresh candidate packages")
    refreshed: list[Path] = []
    connection = psycopg2.connect(**load_db_config(db_config))
    try:
        with connection.cursor() as cursor:
            for path in paths:
                metadata = json.loads(path.read_text())
                og_id = str(metadata["og_id"])
                assembly_prefix = str(metadata["assembly_prefix"])
                parts = assembly_prefix.split(".")
                if len(parts) != 4:
                    raise ValueError(
                        f"Invalid four-part assembly prefix: {assembly_prefix}"
                    )
                # Must match the source of truth PREPARE_ENA_METADATA reads, or
                # a refresh would disagree with the package it is refreshing.
                cursor.execute(
                    """
                    SELECT ncbi_biosample_id
                    FROM sample
                    WHERE og_id = %s
                    """,
                    (og_id,),
                )
                sample = cursor.fetchone()
                biosample = (sample[0] or "").strip() if sample else ""
                cursor.execute(
                    """
                    SELECT mean_depth
                    FROM mitogenome_data
                    WHERE og_id = %s AND tech = %s
                      AND seq_date = %s AND code = %s
                    """,
                    tuple(parts),
                )
                depth = cursor.fetchone()
                published = str(metadata.get("published_package_path") or "").strip()
                target = Path(published) if published else path.parent
                if not target.is_dir():
                    sys.stderr.write(
                        f"WARNING: published package {target} is missing; refreshing the "
                        f"staged copy instead, which will not reach the submitter.\n"
                    )
                    target = path.parent
                ENA_PACKAGE.refresh_package(
                    target,
                    {
                        "biosample_accession": biosample or None,
                        "mean_depth": depth[0] if depth else None,
                    },
                )
                refreshed.append(
                    path if target == path.parent
                    else sorted(target.glob("*.package_metadata.json"))[0]
                )
    finally:
        connection.close()
    return refreshed


def register_and_apply(
    packages: list[dict[str, object]],
    results: dict[tuple[str, str], dict[str, object]],
    report_digest: str,
    db_config: Path,
    selected_by: str,
) -> None:
    if psycopg2 is None:
        raise RuntimeError("psycopg2 is required for --mode apply")
    connection = psycopg2.connect(**load_db_config(db_config))
    try:
        with connection.cursor() as cursor:
            for package in packages:
                og_id = str(package["og_id"])
                numeric = int(og_id[2:])
                cursor.execute(
                    """
                    INSERT INTO ena_specimen_accessions (og_id, og_numeric, ena_biosample_accession)
                    VALUES (%s, %s, %s)
                    ON CONFLICT (og_id) DO UPDATE SET
                        ena_biosample_accession = COALESCE(
                            EXCLUDED.ena_biosample_accession,
                            ena_specimen_accessions.ena_biosample_accession
                        ),
                        updated_at = CURRENT_TIMESTAMP
                    """,
                    (og_id, numeric, package.get("biosample_accession")),
                )
                parts = str(package["full_seqid"]).split(".")
                cursor.execute(
                    """
                    INSERT INTO ena_candidate_packages (
                        full_seqid, og_id, assembly_prefix, annotation_version,
                        ena_study_accession, package_path, package_digest,
                        sequence_sha256, normalised_circular_sha256,
                        package_status, local_validation_status, mean_depth,
                        assembly_program, platform, biosample_accession,
                        webin_test_status, webin_production_status
                    ) VALUES (
                        %s, %s, %s, %s, %s, %s, %s, %s, %s,
                        %s, %s, %s, %s, %s, %s, %s, %s
                    )
                    ON CONFLICT (full_seqid) DO UPDATE SET
                        package_path = EXCLUDED.package_path,
                        package_digest = EXCLUDED.package_digest,
                        package_status = EXCLUDED.package_status,
                        local_validation_status = EXCLUDED.local_validation_status,
                        -- Selection only ever observes Webin statuses second
                        -- hand, from files beside the package.  Those are
                        -- absent whenever selection runs somewhere the
                        -- validation outputs are not, so an incoming NOT_RUN
                        -- means "no news", never "reset".  Only
                        -- RECORD_ENA_PACKAGE_VALIDATION clears a status.
                        webin_test_status = CASE
                            WHEN EXCLUDED.webin_test_status = 'NOT_RUN'
                            THEN ena_candidate_packages.webin_test_status
                            ELSE EXCLUDED.webin_test_status END,
                        webin_production_status = CASE
                            WHEN EXCLUDED.webin_production_status = 'NOT_RUN'
                            THEN ena_candidate_packages.webin_production_status
                            ELSE EXCLUDED.webin_production_status END,
                        updated_at = CURRENT_TIMESTAMP
                    """,
                    (
                        package["full_seqid"],
                        og_id,
                        package.get("assembly_prefix") or ".".join(parts[:4]),
                        package.get("annotation_version") or parts[-1],
                        package.get("study") or "",
                        package["_package_path"],
                        package.get("package_digest"),
                        package["sequence_sha256"],
                        package["normalised_circular_sha256"],
                        package["package_status"],
                        package["local_validation_status"],
                        package.get("mean_depth"),
                        package.get("program"),
                        package.get("platform"),
                        package.get("biosample_accession"),
                        package.get("webin_test_status") or "NOT_RUN",
                        package.get("webin_production_status") or "NOT_RUN",
                    ),
                )
                mapping_path = (
                    Path(str(package["_package_path"]))
                    / f"{package['full_seqid']}.locus_tag_mapping.tsv"
                )
                if not mapping_path.is_file():
                    raise ValueError(f"Missing locus-tag mapping: {mapping_path}")
                with mapping_path.open(newline="") as handle:
                    loci = list(csv.DictReader(handle, delimiter="\t"))
                for locus in loci:
                    tag = locus.get("locus_tag", "")
                    match = LOCUS_TAG_PATTERN.fullmatch(tag)
                    if not match:
                        raise ValueError(
                            f"Invalid locus tag in {mapping_path}: {tag!r}"
                        )
                    # The prefix is technology specific and validated by the
                    # database check, but the OG number embedded in the tag must
                    # match the specimen the package claims to be.
                    if int(match.group(1)) != numeric:
                        raise ValueError(
                            f"Locus tag {tag} in {mapping_path} does not belong to {og_id}"
                        )
                    cursor.execute(
                        """
                        INSERT INTO ena_candidate_loci (
                            full_seqid, og_id, gene_serial, locus_tag,
                            feature_key, start_coordinate, end_coordinate, strand
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT (full_seqid, feature_key) DO UPDATE SET
                            gene_serial = EXCLUDED.gene_serial,
                            locus_tag = EXCLUDED.locus_tag,
                            start_coordinate = EXCLUDED.start_coordinate,
                            end_coordinate = EXCLUDED.end_coordinate,
                            strand = EXCLUDED.strand
                        """,
                        (
                            package["full_seqid"],
                            og_id,
                            int(tag[-3:]),
                            tag,
                            locus["feature_key"],
                            int(locus["start"]) if locus.get("start") else None,
                            int(locus["end"]) if locus.get("end") else None,
                            locus.get("strand") or None,
                        ),
                    )
            by_group = defaultdict(list)
            for package in packages:
                by_group[group_key(package)].append(package)
            for key, result in results.items():
                og_id, tech = key
                group = by_group[key]
                # Every candidate in a technology group must name the same child
                # study; a disagreement means one of them was packaged under a
                # different study accession and would be submitted to the wrong
                # place. The primary key below would silently keep only one.
                studies = {str(package.get("study") or "") for package in group}
                if len(studies) > 1:
                    raise ValueError(
                        f"{og_id} {tech} packages disagree on ENA study: "
                        f"{', '.join(sorted(studies))}"
                    )
                study = studies.pop()
                if not study:
                    raise ValueError(f"No ENA study recorded for {og_id} {tech}")
                cursor.execute(
                    """
                    INSERT INTO ena_submission_selections (
                        ena_study_accession, og_id, biosample_accession,
                        selected_full_seqid, selection_status, selection_reason,
                        selection_report_digest, selected_by
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (ena_study_accession, og_id) DO UPDATE SET
                        biosample_accession = EXCLUDED.biosample_accession,
                        selected_full_seqid = EXCLUDED.selected_full_seqid,
                        selection_status = EXCLUDED.selection_status,
                        selection_reason = EXCLUDED.selection_reason,
                        selection_report_digest = EXCLUDED.selection_report_digest,
                        selected_by = EXCLUDED.selected_by,
                        selected_at = CURRENT_TIMESTAMP,
                        updated_at = CURRENT_TIMESTAMP
                    WHERE ena_submission_selections.archive_status = 'NOT_SUBMITTED'
                    """,
                    (
                        study,
                        og_id,
                        group[0].get("biosample_accession"),
                        result["selected"],
                        result["status"],
                        result["reason"],
                        report_digest,
                        result.get("selected_by") or selected_by,
                    ),
                )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-metadata", action="append", required=True)
    parser.add_argument("--mode", choices=("report", "apply"), default="report")
    parser.add_argument("--decision-file")
    parser.add_argument("--output", required=True)
    parser.add_argument("--selected-output")
    parser.add_argument("--db-config")
    parser.add_argument("--selected-by", default="pipeline")
    args = parser.parse_args()
    try:
        paths = expand_paths(args.package_metadata)
        if not paths:
            raise ValueError("No package metadata files matched")
        if args.mode == "apply":
            if not args.db_config:
                raise ValueError("--db-config is required for --mode apply")
            paths = refresh_packages_from_database(paths, Path(args.db_config))
        packages = load_packages(paths)
        decisions = read_decisions(Path(args.decision_file) if args.decision_file else None)
        rows, results = build_report(packages, decisions)
        report_digest = write_report(Path(args.output), rows)
        if args.selected_output:
            write_selected_packages(Path(args.selected_output), packages, results)
        if args.mode == "apply":
            register_and_apply(
                packages,
                results,
                report_digest,
                Path(args.db_config),
                args.selected_by,
            )
        return 0
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
