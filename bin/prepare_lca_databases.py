#!/usr/bin/env python3
"""Build and verify the persistent database cache used by LCA tasks."""

import argparse
import fcntl
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from calculateLCA import Config, DatabaseManager, FISHBASE_DATABASES


CACHE_SCHEMA_VERSION = 1
PARQUET_REQUIREMENTS = {
    "fishbase_species.parquet": {"SpecCode", "Genus", "Species", "FamCode"},
    "fishbase_families.parquet": {"FamCode", "Family", "Order", "Class"},
    "fishbase_synonyms.parquet": {"SpecCode", "SynGenus", "SynSpecies"},
}
TAXDUMP_FILES = ("taxdump/nodes.dmp", "taxdump/names.dmp")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_cache(cache_dir: Path):
    verified = []
    for filename, required_columns in PARQUET_REQUIREMENTS.items():
        path = cache_dir / filename
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty LCA cache file: {path}")
        frame = pd.read_parquet(path)
        missing_columns = required_columns.difference(frame.columns)
        if frame.empty or missing_columns:
            raise RuntimeError(
                f"Invalid LCA cache file {path}: empty={frame.empty}, "
                f"missing_columns={sorted(missing_columns)}"
            )
        verified.append(path)

    for relative_path in TAXDUMP_FILES:
        path = cache_dir / relative_path
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty LCA cache file: {path}")
        verified.append(path)
    return verified


def write_json_atomic(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def cache_signature(files) -> str:
    """Return a stable fingerprint of the cache contents and schema.

    Runtime audit fields such as ``verified_at`` are deliberately excluded so
    validating an unchanged cache does not invalidate downstream Nextflow tasks.
    """
    payload = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "files": {
            filename: {
                "bytes": metadata["bytes"],
                "sha256": metadata["sha256"],
            }
            for filename, metadata in sorted(files.items())
        },
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def prepare(cache_dir: Path, manifest_path: Path):
    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    lock_path = cache_dir / ".prepare.lock"

    with lock_path.open("a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        manager = DatabaseManager(cache_dir)
        manager.load_fishbase_data()
        manager.load_ncbi_taxdump()
        verified_files = verify_cache(cache_dir)

        files = {
            str(path.relative_to(cache_dir)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in verified_files
        }
        manifest = {
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "cache_signature": cache_signature(files),
            "status": "verified",
            "cache_dir": str(cache_dir),
            "verified_at": datetime.now(timezone.utc).isoformat(),
            "sources": {
                **FISHBASE_DATABASES,
                "taxdump.tar.gz": Config.NCBI_TAXDUMP_URL,
            },
            "files": files,
        }
        write_json_atomic(cache_dir / "manifest.json", manifest)
        write_json_atomic(manifest_path, manifest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.cache_dir, args.manifest)


if __name__ == "__main__":
    main()
