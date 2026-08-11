#!/usr/bin/env python3
"""Allocate OG-derived locus tags and inject them into an NCBI feature table."""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import importlib.util
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import psycopg2
except ImportError:  # unit tests and file-backed package preparation
    psycopg2 = None


ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("ena_package", ROOT / "ena_package.py")
ENA_PACKAGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ENA_PACKAGE)

TAGGABLE = {"gene", "CDS", "tRNA", "rRNA", "mRNA", "ncRNA", "misc_RNA", "intron", "exon"}
UNTAGGED = {"source", "control_region", "repeat_region", "misc_feature"}


@dataclass
class Feature:
    lines: list[str]
    key: str
    start: int
    end: int
    qualifiers: list[tuple[str, str]] = field(default_factory=list)
    canonical_gene: str = ""
    occurrence: int = 1
    locus_tag: str = ""

    @property
    def coordinate_key(self) -> tuple[int, int]:
        return min(self.start, self.end), max(self.start, self.end)

    @property
    def strand(self) -> str:
        return "+" if self.start <= self.end else "-"


def feature_start(line: str) -> tuple[int, int, str] | None:
    if line.startswith("\t"):
        return None
    columns = line.split("\t")
    if len(columns) < 3 or not columns[2]:
        return None
    try:
        start = int(columns[0].lstrip("<>"))
        end = int(columns[1].lstrip("<>"))
    except ValueError:
        return None
    return start, end, columns[2]


def qualifier(line: str) -> tuple[str, str] | None:
    if not line.startswith("\t"):
        return None
    values = line.strip("\t").split("\t", 1)
    if not values or not values[0]:
        return None
    return values[0], values[1] if len(values) > 1 else ""


def parse_tbl(path: Path) -> tuple[list[str], list[Feature]]:
    """Read a feature table into coordinate-ordered features.

    Serials are handed out by walking the returned list, and inject_tags and
    write_mapping emit in that order too, so the order is part of what gets
    published rather than an internal detail.  Sorting here rather than in the
    caller means an annotator cannot reach through and scramble the tags:
    Emma sorts its feature starts as strings, so its tables arrive as
    1, 10053, 1027, 10343, 1099 ...  bin/process_files.py already normalises
    that upstream; this is the guarantee that holds if it ever stops.

    The sort is stable and keyed on the lowest coordinate only, which keeps a
    gene adjacent to the CDS/mRNA/tRNA emitted beneath it.
    """
    preamble: list[str] = []
    features: list[Feature] = []
    current: Feature | None = None
    for raw in path.read_text().splitlines():
        start = feature_start(raw)
        if start:
            if current is not None:
                features.append(current)
            current = Feature(lines=[raw], key=start[2], start=start[0], end=start[1])
        elif current is None:
            preamble.append(raw)
        else:
            current.lines.append(raw)
            parsed = qualifier(raw)
            if parsed:
                current.qualifiers.append(parsed)
    if current is not None:
        features.append(current)
    features.sort(key=lambda feature: feature.coordinate_key[0])
    return preamble, features


def normalize_gene(value: str) -> str:
    value = value.strip().replace("MT-", "")
    value = re.sub(r"\s+", "", value)
    return value.upper()


def named_gene(feature: Feature) -> str:
    values = dict(feature.qualifiers)
    for key in ("gene", "gene_synonym", "product"):
        if values.get(key):
            return normalize_gene(values[key])
    return ""


def assign_locus_identities(features: list[Feature]) -> list[Feature]:
    roots = [feature for feature in features if feature.key == "gene"]
    for root in roots:
        root.canonical_gene = named_gene(root)
        if not root.canonical_gene:
            raise ValueError(
                f"Gene feature at {root.start}..{root.end} has no gene/product identity"
            )
    by_name: dict[str, list[Feature]] = {}
    for root in sorted(roots, key=lambda item: (item.coordinate_key, item.strand)):
        by_name.setdefault(root.canonical_gene, []).append(root)
    for same_name in by_name.values():
        for occurrence, root in enumerate(same_name, 1):
            root.occurrence = occurrence

    for feature in features:
        if feature.key in UNTAGGED or feature.key not in TAGGABLE or feature.key == "gene":
            continue
        own_name = named_gene(feature)
        exact = [root for root in roots if root.coordinate_key == feature.coordinate_key]
        named = [root for root in exact if not own_name or root.canonical_gene == own_name]
        candidates = named or exact
        if len(candidates) != 1 and own_name:
            candidates = [root for root in roots if root.canonical_gene == own_name]
        if len(candidates) == 1:
            feature.canonical_gene = candidates[0].canonical_gene
            feature.occurrence = candidates[0].occurrence
        elif len(candidates) > 1:
            raise ValueError(
                f"Ambiguous locus for {feature.key} at {feature.start}..{feature.end}"
            )
        else:
            raise ValueError(
                f"No gene feature matches {feature.key} at {feature.start}..{feature.end}"
            )
    return roots


def read_registry(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        (row["canonical_gene"], int(row["gene_occurrence"])): row
        for row in rows
    }


def write_registry(path: Path, registry: dict[tuple[str, int], dict[str, str]]) -> None:
    fields = [
        "og_id",
        "gene_serial",
        "canonical_gene",
        "gene_occurrence",
        "feature_type",
        "strand",
        "coordinate_snapshot",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(registry.values(), key=lambda item: int(item["gene_serial"])):
            writer.writerow({field: row.get(field, "") for field in fields})


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


def allocate_file_registry(
    roots: list[Feature], og_id: str, registry_path: Path
) -> dict[tuple[str, int], dict[str, str]]:
    registry = read_registry(registry_path)
    next_serial = max((int(row["gene_serial"]) for row in registry.values()), default=0) + 1
    for root in roots:
        key = root.canonical_gene, root.occurrence
        if key not in registry:
            if next_serial > 999:
                raise ValueError(f"{og_id} exceeds 999 allocated mitochondrial loci")
            registry[key] = {
                "og_id": og_id,
                "gene_serial": str(next_serial),
                "canonical_gene": root.canonical_gene,
                "gene_occurrence": str(root.occurrence),
                "feature_type": root.key,
                "strand": root.strand,
                "coordinate_snapshot": f"{root.start}..{root.end}",
            }
            next_serial += 1
    write_registry(registry_path, registry)
    return registry


def allocate_db_registry(
    roots: list[Feature], og_id: str, db_config: Path
) -> dict[tuple[str, int], dict[str, str]]:
    if psycopg2 is None:
        raise RuntimeError("psycopg2 is required for database-backed locus allocation")
    numeric = ENA_PACKAGE.og_numeric(og_id)
    connection = psycopg2.connect(**load_db_config(db_config))
    try:
        with connection.cursor() as cursor:
            cursor.execute("SELECT pg_advisory_xact_lock(%s)", (numeric,))
            cursor.execute(
                """
                SELECT og_id
                FROM ena_specimen_accessions
                WHERE og_numeric = %s AND og_id <> %s
                """,
                (numeric, og_id),
            )
            collision = cursor.fetchone()
            if collision:
                raise ValueError(
                    f"OG numeric collision: {og_id} and {collision[0]} both map to {numeric}"
                )
            cursor.execute(
                """
                INSERT INTO ena_specimen_accessions (og_id, og_numeric)
                VALUES (%s, %s)
                ON CONFLICT (og_id) DO UPDATE SET updated_at = CURRENT_TIMESTAMP
                """,
                (og_id, numeric),
            )
            cursor.execute(
                """
                SELECT gene_serial, canonical_gene, gene_occurrence,
                       feature_type, strand, coordinate_snapshot
                FROM ena_locus_registry
                WHERE og_id = %s
                ORDER BY gene_serial
                """,
                (og_id,),
            )
            registry = {
                (row[1], row[2]): {
                    "og_id": og_id,
                    "gene_serial": str(row[0]),
                    "canonical_gene": row[1],
                    "gene_occurrence": str(row[2]),
                    "feature_type": row[3],
                    "strand": row[4] or "",
                    "coordinate_snapshot": row[5] or "",
                }
                for row in cursor.fetchall()
            }
            next_serial = max(
                (int(row["gene_serial"]) for row in registry.values()), default=0
            ) + 1
            for root in roots:
                key = root.canonical_gene, root.occurrence
                if key in registry:
                    continue
                if next_serial > 999:
                    raise ValueError(f"{og_id} exceeds 999 allocated mitochondrial loci")
                cursor.execute(
                    """
                    INSERT INTO ena_locus_registry (
                        og_id, gene_serial, canonical_gene,
                        gene_occurrence, feature_type, strand, coordinate_snapshot
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s)
                    """,
                    (
                        og_id,
                        next_serial,
                        root.canonical_gene,
                        root.occurrence,
                        root.key,
                        root.strand,
                        f"{root.start}..{root.end}",
                    ),
                )
                registry[key] = {
                    "og_id": og_id,
                    "gene_serial": str(next_serial),
                    "canonical_gene": root.canonical_gene,
                    "gene_occurrence": str(root.occurrence),
                    "feature_type": root.key,
                    "strand": root.strand,
                    "coordinate_snapshot": f"{root.start}..{root.end}",
                }
                next_serial += 1
        connection.commit()
        return registry
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def render_tags(
    registry: dict[tuple[str, int], dict[str, str]], og_id: str, prefix: str
) -> dict[tuple[str, int], dict[str, str]]:
    """Resolve each allocated serial to a tag under this candidate's prefix.

    The registry is specimen level and prefix free: OG910 gene 1 is serial 1 for
    every technology.  The prefix belongs to the technology's ENA study, so the
    same gene is OGMTHIFI_000910001 in the HiFi record and OGMTHIC_000910001 in
    the Hi-C record.  Both are published, and ENA forbids sharing a tag between
    two records, so the rendering has to happen here rather than at allocation.
    """
    for row in registry.values():
        row["locus_tag"] = ENA_PACKAGE.locus_tag(og_id, int(row["gene_serial"]), prefix)
    return registry


def inject_tags(
    preamble: list[str],
    features: list[Feature],
    registry: dict[tuple[str, int], dict[str, str]],
    output_path: Path,
) -> None:
    output = list(preamble)
    for feature in features:
        if feature.canonical_gene:
            row = registry[(feature.canonical_gene, feature.occurrence)]
            feature.locus_tag = row["locus_tag"]
        lines = [
            line for line in feature.lines
            if not (line.startswith("\t") and qualifier(line) and qualifier(line)[0] == "locus_tag")
        ]
        if feature.locus_tag:
            lines.append(f"\t\t\tlocus_tag\t{feature.locus_tag}")
        output.extend(lines)
    output_path.write_text("\n".join(output) + "\n")


def write_mapping(path: Path, full_seqid: str, features: list[Feature]) -> None:
    fields = [
        "full_seqid",
        "feature_key",
        "feature_type",
        "canonical_gene",
        "gene_occurrence",
        "start",
        "end",
        "strand",
        "locus_tag",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, feature in enumerate(features, 1):
            if not feature.locus_tag:
                continue
            writer.writerow(
                {
                    "full_seqid": full_seqid,
                    "feature_key": f"{feature.key}:{index}",
                    "feature_type": feature.key,
                    "canonical_gene": feature.canonical_gene,
                    "gene_occurrence": feature.occurrence,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand,
                    "locus_tag": feature.locus_tag,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--og-id", required=True)
    parser.add_argument("--full-seqid", required=True)
    parser.add_argument("--input-tbl", required=True)
    parser.add_argument("--output-tbl", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument(
        "--locus-prefix",
        required=True,
        help="INSDC prefix registered to this candidate's technology study",
    )
    allocation = parser.add_mutually_exclusive_group(required=True)
    allocation.add_argument("--registry-tsv")
    allocation.add_argument("--db-config")
    args = parser.parse_args()
    try:
        if ENA_PACKAGE.full_seqid_og_id(args.full_seqid) != args.og_id:
            raise ValueError("Full SeqID does not belong to --og-id")
        preamble, features = parse_tbl(Path(args.input_tbl))
        roots = assign_locus_identities(features)
        if not roots:
            raise ValueError("Feature table contains no gene features")
        if args.db_config:
            registry = allocate_db_registry(roots, args.og_id, Path(args.db_config))
        else:
            registry = allocate_file_registry(roots, args.og_id, Path(args.registry_tsv))
        render_tags(registry, args.og_id, args.locus_prefix)
        for feature in features:
            if feature.canonical_gene:
                feature.locus_tag = registry[
                    (feature.canonical_gene, feature.occurrence)
                ]["locus_tag"]
        inject_tags(preamble, features, registry, Path(args.output_tbl))
        write_mapping(Path(args.mapping), args.full_seqid, features)
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
