#!/usr/bin/env python3
"""Fetch candidate-specific ENA metadata from the OceanOmics database."""

from __future__ import annotations

import argparse
import configparser
import json
import re
import sys
from pathlib import Path

try:
    import psycopg2
except ImportError:
    psycopg2 = None


def load_db_config(path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser()
    parser.read(path)
    if not parser.has_section("postgres"):
        raise ValueError(f"Missing [postgres] section in {path}")
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def dotted_version(digits: str) -> str:
    return ".".join(digits) if digits else ""


def assembly_program(code: str) -> str:
    lowered = code.lower()
    digits = "".join(re.findall(r"[0-9]+", lowered))
    if "getorg" in lowered:
        return f"GetOrganelle {dotted_version(digits)}".strip()
    if "mitohifi" in lowered:
        return f"MitoHiFi {dotted_version(digits)}".strip()
    if "oatk" in lowered:
        return f"Oatk {dotted_version(digits)}".strip()
    raise ValueError(f"Unsupported assembly code for ENA PROGRAM: {code}")


def platform_for_tech(tech: str) -> str:
    try:
        return {"hifi": "PACBIO_SMRT", "ilmn": "ILLUMINA", "hic": "ILLUMINA"}[
            tech.lower()
        ]
    except KeyError as error:
        raise ValueError(f"Unsupported sequencing technology for ENA: {tech}") from error


def fetch_metadata(connection, args: argparse.Namespace) -> dict[str, object]:
    with connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT mean_depth
            FROM mitogenome_data
            WHERE og_id = %s AND tech = %s AND seq_date = %s AND code = %s
            """,
            (args.og_id, args.tech, args.seq_date, args.code),
        )
        depth_row = cursor.fetchone()
        mean_depth = depth_row[0] if depth_row else None
    # The BioSample and the run accessions are not fetched: the submission
    # pipeline registers the sample and owns the raw-read submissions, so the
    # SAMPLE and RUN_REF manifest keys are its to supply and nothing recorded
    # here could be more than a guess at them.
    return {
        # 2: dropped biosample_accession, biosample_source and run_accessions.
        "schema_version": 2,
        "og_id": args.og_id,
        "assembly_prefix": args.assembly_prefix,
        "annotation_version": args.annotation_version,
        "full_seqid": args.full_seqid,
        "validation_study": args.study,
        "mean_depth": mean_depth,
        "program": assembly_program(args.code),
        "platform": platform_for_tech(args.tech),
        "scientific_name": args.scientific_name,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--og-id", required=True)
    parser.add_argument("--assembly-prefix", required=True)
    parser.add_argument("--annotation-version", required=True)
    parser.add_argument("--full-seqid", required=True)
    parser.add_argument("--tech", required=True)
    parser.add_argument("--seq-date", required=True)
    parser.add_argument("--code", required=True)
    parser.add_argument(
        "--study",
        required=True,
        help="Study used for sequence-context validation, recorded as validation_study.",
    )
    parser.add_argument("--scientific-name", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        if psycopg2 is None:
            raise RuntimeError("psycopg2 is required to fetch ENA metadata")
        with psycopg2.connect(**load_db_config(Path(args.config))) as connection:
            metadata = fetch_metadata(connection, args)
        Path(args.output).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
