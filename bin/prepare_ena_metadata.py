#!/usr/bin/env python3
"""Fetch candidate-specific ENA metadata from the OceanOmics database.

Coverage is the exception: it is read from this run's own depth TSV when one is
supplied, and only falls back to mitogenome_data.mean_depth when it is not. The
pipeline measures the depth, push_mtdna_assm_results.py writes that same number
to the row, and this script used to SELECT it straight back -- a round trip that
bought nothing and cost two things. It made submission prep wait on the committed
row (which is why prep could not run with --skip_upload_results), and on a re-run
that skipped the write it would silently return the PREVIOUS assembly's depth for
a molecule that had just been reassembled. Reading the file removes both.

The database fallback is not dead code: --skip_mitogenome_depth and runs with
precomputed assemblies produce no measurement, and there a previously stored
value is the best available answer. mean_depth_source records which one won, so a
package can be audited after the fact rather than guessed at.
"""

from __future__ import annotations

import argparse
import configparser
import json
import re
import sys
from pathlib import Path

# bin/ is on PATH at runtime and sys.path[0] is this file's directory, so the
# sibling import resolves the same way push_species_validation.py imports
# load_db_config.
from depth_tsv import read_mean_depth

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


def resolve_mean_depth(connection, args: argparse.Namespace):
    """This run's measured depth, else the stored one, else nothing.

    Returns (mean_depth, source) where source is 'pipeline', 'database' or
    'none'. The database is queried only when the file yielded nothing, so on
    the normal path this opens no cursor for depth at all.
    """
    measured = read_mean_depth(getattr(args, "depth_tsv", None))
    if measured is not None:
        return measured, "pipeline"

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
    stored = depth_row[0] if depth_row else None
    if stored is None:
        return None, "none"
    return stored, "database"


def fetch_metadata(connection, args: argparse.Namespace) -> dict[str, object]:
    mean_depth, mean_depth_source = resolve_mean_depth(connection, args)
    if mean_depth_source == "none":
        # Not fatal: manifest_fields() omits COVERAGE rather than emitting a
        # wrong one, and the package is still worth handing over. Said out loud
        # because in a --skip_upload_results run this is the expected shape of a
        # sample that has never been uploaded, and a silent omission there is
        # indistinguishable from a measurement that genuinely came back empty.
        print(
            f"WARNING: no mean_depth for {args.full_seqid} from either the depth "
            "TSV or mitogenome_data; the ENA manifest will carry no COVERAGE",
            file=sys.stderr,
        )
    # The BioSample and the run accessions are not fetched: the submission
    # pipeline registers the sample and owns the raw-read submissions, so the
    # SAMPLE and RUN_REF manifest keys are its to supply and nothing recorded
    # here could be more than a guess at them.
    return {
        # This is the ena_input_metadata.json schema, which versions separately
        # from the package_metadata.json schema in ena_package.py -- the two share
        # a key name and nothing else, and build_package copies values across by
        # name without ever reading this number.
        # 2: dropped biosample_accession, biosample_source and run_accessions.
        # 3: added mean_depth_source.
        "schema_version": 3,
        "og_id": args.og_id,
        "assembly_prefix": args.assembly_prefix,
        "annotation_version": args.annotation_version,
        "full_seqid": args.full_seqid,
        "validation_study": args.study,
        "mean_depth": mean_depth,
        "mean_depth_source": mean_depth_source,
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
    parser.add_argument(
        "--depth-tsv",
        help=(
            "This run's <prefix>.mito_depth.tsv. Its mean_depth wins over the "
            "stored mitogenome_data value. Omit it, or pass the header-only "
            "assets/placeholders/empty_mito_depth.tsv, to fall back to the database."
        ),
    )
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
