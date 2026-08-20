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
        # sample.ncbi_biosample_id is the source of truth for the specimen's
        # BioSample. ena_specimen_accessions is a derived cache that goes stale
        # whenever sample is edited, so it is not read here.
        cursor.execute(
            """
            SELECT ncbi_biosample_id
            FROM sample
            WHERE og_id = %s
            """,
            (args.og_id,),
        )
        sample_row = cursor.fetchone()
        # Unregistered specimens hold '' rather than NULL; an empty string would
        # otherwise read as "present but malformed" and mislabel the package
        # BLOCKED_METADATA instead of WAITING_FOR_BIOSAMPLE.
        biosample = (sample_row[0] or "").strip() if sample_row else ""
        biosample = biosample or None
    return {
        "schema_version": 1,
        "og_id": args.og_id,
        "assembly_prefix": args.assembly_prefix,
        "annotation_version": args.annotation_version,
        "full_seqid": args.full_seqid,
        "study": args.study,
        "biosample_accession": biosample,
        "biosample_source": "sample.ncbi_biosample_id",
        "mean_depth": mean_depth,
        "program": assembly_program(args.code),
        "platform": platform_for_tech(args.tech),
        "scientific_name": args.scientific_name,
        # Run accessions belong to the downstream submission pipeline, which
        # owns the raw-read submissions to PRJEB123419/420/421. Nothing in this
        # schema records them, so an empty list is emitted and RUN_REF is left
        # out of the manifest for the submitter to add.
        "run_accessions": [],
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
    parser.add_argument("--study", required=True)
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
