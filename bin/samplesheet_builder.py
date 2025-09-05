#!/usr/bin/env python3
"""
Build a samplesheet CSV for an nf-core pipeline from a directory of FASTQ files
covering three data types:

1) Illumina PE named like:  OG764.ilmn.240716.R1.fastq.gz / R2
   -> sample=OG764, type=ilmn, date=240716, assembly_prefix=OG764.ilmn.240716,
      fastq_1=/abs/path/R1, fastq_2=/abs/path/R2 (exactly one pair expected)

2) Hi-C PE named like:      OG765W-3_HICL_S5_L001_R1_001.fastq.gz (plus mates/lanes)
   -> sample=OG765 (first token like OG\d+), type=hic
      Search in s3://<bucket>/<prefix>/NOVA_250101_AA/OG765/ for matching files
      (e.g., OG765W-3_HICL_S5_L001_R{1,2}_001.fastq.gz) to get the *date* 250101
      from the run folder name (e.g., "NOVA_250101_AA").
      There may be multiple R1/R2 pairs; emit one CSV row per pair with
      assembly_prefix=OG765.hic.250101 and s3 URIs as fastq_1/fastq_2.

3) PacBio HiFi SE like:     OG785_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz
   -> sample=OG785, type=hifi, date=241004, assembly_prefix=OG785.hifi.241004,
      fastq_1=/abs/path/file (fastq_2 empty)

Output CSV columns: sample,type,date,assembly_prefix,fastq_1,fastq_2

Notes
-----
- Requires AWS credentials (env/instance profile) to search S3 for Hi-C data.
- Uses boto3 if available; otherwise falls back to the AWS CLI (`aws s3 ls`).
- Designed to be idempotent and safe to run within an nf-core module.

Example
-------
python build_samplesheet.py \
  --input-dir /data/reads \
  --out samplesheet.csv \
  --s3-bucket oceanomics \
  --s3-prefix OceanGenomes/illumina-hic
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import subprocess
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# -----------------------------
# Helpers: S3 listing utilities
# -----------------------------

def _boto3_available() -> bool:
    try:
        import boto3  # type: ignore
        return True
    except Exception:
        return False


def list_s3_keys(bucket: str, prefix: str) -> Iterable[str]:
    """Yield object keys under s3://bucket/prefix using boto3 if available,
    else using `aws s3 ls --recursive`.
    """
    if _boto3_available():
        import boto3  # type: ignore
        s3 = boto3.client("s3")
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            for obj in page.get("Contents", []) or []:
                yield obj["Key"]
        return

    # Fallback: AWS CLI
    cmd = [
        "aws",
        "s3",
        "ls",
        f"s3://{bucket}/{prefix}",
        "--recursive",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        raise RuntimeError(
            "Neither boto3 nor aws CLI is available. Install boto3 or the AWS CLI."
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"aws s3 ls failed: {e.stderr.strip()}")

    # Output lines look like: '2025-01-01 12:34:56      12345 OceanGenomes/illumina-hic/.../file.fastq.gz'
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 4:
            key = " ".join(parts[3:])
            yield key


# ---------------------------------
# Parsing and grouping local datasets
# ---------------------------------
ILMN_RE = re.compile(
    r"^(?P<sample>[^.]+)\.(?P<type>ilmn)\.(?P<date>\d{6})\.(?P<read>R[12])\.fastq\.gz$"
)

# Example: OG765W-3_HICL_S5_L001_R1_001.fastq.gz (or ..._HIC_...)
HIC_RE = re.compile(
    r"^(?P<idprefix>(?P<sample>OG\d+)[^_]*)_(?P<hic>HIC\w*)_S\d+_L\d{3}_R[12]_001\.fastq\.gz$",
    re.IGNORECASE,
)

# Example: OG785_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz
HIFI_RE = re.compile(
    r"^(?P<sample>OG\d+)_.*_(?P<date>\d{6})_.*\.hifi_reads.*\.fastq\.gz$",
    re.IGNORECASE,
)


def scan_local(input_dir: Path) -> Tuple[Dict[str, Dict[str, str]], List[Tuple[str, str]], List[Path]]:
    """Scan local directory recursively for FASTQ files.

    Returns a tuple of:
      - ilmn_groups: mapping from assembly_prefix -> {"R1": path, "R2": path, "sample":, "date":}
      - hic_bases: list of tuples (sample, base_id) derived from local HIC filenames
      - hifi_files: list of Paths to HiFi files
    """
    ilmn_groups: Dict[str, Dict[str, str]] = {}
    hic_seen: set = set()
    hic_bases: List[Tuple[str, str]] = []
    hifi_files: List[Path] = []

    for p in input_dir.rglob("*.fastq.gz"):
        name = p.name
        m = ILMN_RE.match(name)
        if m:
            d = m.groupdict()
            sample = d["sample"]
            date = d["date"]
            read = d["read"]
            assembly_prefix = f"{sample}.ilmn.{date}"
            grp = ilmn_groups.setdefault(
                assembly_prefix,
                {"sample": sample, "date": date, "R1": "", "R2": ""},
            )
            grp[read] = str(p.resolve())
            continue

        m = HIC_RE.match(name)
        if m:
            gd = m.groupdict()
            sample = gd["sample"]
            # Base ID to search with on S3: up to _S (exclude lane/read parts)
            # Use the piece before _S in the full filename
            base_before_S = name.split("_S", 1)[0]  # e.g., OG765W-3_HICL
            # Also consider a variant trimmed to just '..._HIC' (strip trailing letters after HIC)
            trimmed = re.sub(r"(HIC)[A-Za-z]*$", r"\1", base_before_S)
            for base_id in {base_before_S, trimmed}:
                key = (sample, base_id)
                if key not in hic_seen:
                    hic_seen.add(key)
                    hic_bases.append((sample, base_id))
            continue

        m = HIFI_RE.match(name)
        if m:
            hifi_files.append(p.resolve())
            continue

        # Non-matching fastq.gz files are ignored (by design)

    # Validate ilmn pairs
    for ap, grp in list(ilmn_groups.items()):
        if not grp["R1"] or not grp["R2"]:
            logging.warning(
                "Incomplete Illumina pair for %s (R1: %s, R2: %s) — entry will be skipped.",
                ap,
                bool(grp["R1"]),
                bool(grp["R2"]),
            )
            ilmn_groups.pop(ap, None)

    return ilmn_groups, hic_bases, hifi_files


# ---------------------------------
# Hi-C: find pairs and date from S3
# ---------------------------------
RUN_DIR_DATE_RE = re.compile(r"/[A-Za-z]+_(?P<date>\d{6})_[A-Za-z]+/")


def find_hic_pairs_from_s3(
    bucket: str,
    root_prefix: str,
    sample: str,
    base_id_variants: List[str],
) -> Tuple[str, List[Tuple[str, str]]]:
    """Search S3 for Hi-C FASTQs for a given sample and base-id variants.

    Returns (date, pairs) where pairs is a list of (s3://...R1..., s3://...R2...).
    Raises RuntimeError if nothing found.
    """
    prefix = root_prefix.rstrip("/") + "/"
    keys = list(list_s3_keys(bucket, prefix))

    # Filter to this sample and base-id variants
    variant_set = set(base_id_variants)
    cand = [
        k for k in keys
        if f"/{sample}/" in k
        and any(v in os.path.basename(k) for v in variant_set)
        and k.endswith(".fastq.gz")
        and re.search(r"_R[12]_001\.fastq\.gz$", k)
    ]
    if not cand:
        raise RuntimeError(
            f"No Hi-C FASTQs found in s3://{bucket}/{prefix} for sample {sample} "
            f"with bases {sorted(variant_set)}."
        )

    # Decide date by majority vote based on run directory segment
    dates = []
    for k in cand:
        m = RUN_DIR_DATE_RE.search("/" + k)  # ensure leading slash for regex
        if m:
            dates.append(m.group("date"))
    if not dates:
        raise RuntimeError(
            "Could not infer date from S3 run directories (expected like NOVA_250101_AA)."
        )
    date = Counter(dates).most_common(1)[0][0]

    # Keep only keys from the chosen date
    cand = [k for k in cand if f"_{date}_" in k]

    # Group into R1/R2 pairs by collapsing the _R1_/_R2_ token
    groups: Dict[str, Dict[str, str]] = defaultdict(dict)
    for k in cand:
        fname = os.path.basename(k)
        pair_key = fname.replace("_R1_", "__").replace("_R2_", "__")
        if "_R1_" in fname:
            groups[pair_key]["R1"] = f"s3://{bucket}/{k}"
        elif "_R2_" in fname:
            groups[pair_key]["R2"] = f"s3://{bucket}/{k}"

    pairs: List[Tuple[str, str]] = []
    for pk, d in groups.items():
        if "R1" in d and "R2" in d:
            pairs.append((d["R1"], d["R2"]))
        else:
            logging.warning("Skipping incomplete Hi-C pair for key %s", pk)

    if not pairs:
        raise RuntimeError("No complete R1/R2 Hi-C pairs found after grouping.")

    return date, sorted(pairs)


# -----------------------------
# Main orchestration
# -----------------------------

def build_samplesheet(
    input_dir: Path,
    out_csv: Path,
    s3_bucket: str,
    s3_prefix: str,
) -> None:
    ilmn_groups, hic_bases, hifi_files = scan_local(input_dir)

    rows: List[Dict[str, str]] = []

    # Illumina
    for ap, grp in sorted(ilmn_groups.items()):
        rows.append(
            {
                "sample": grp["sample"],
                "type": "ilmn",
                "date": grp["date"],
                "assembly_prefix": ap,
                "fastq_1": grp["R1"],
                "fastq_2": grp["R2"],
            }
        )

    # Hi-C (one row per R1/R2 pair found in S3)
    hic_handled: set = set()
    for sample, base_id in hic_bases:
        if (sample, base_id) in hic_handled:
            continue
        # Collect variant patterns to be robust (with/without trailing letters after HIC)
        trimmed = re.sub(r"(HIC)[A-Za-z]*$", r"\1", base_id)
        variants = sorted({base_id, trimmed})
        try:
            date, pairs = find_hic_pairs_from_s3(s3_bucket, s3_prefix, sample, variants)
        except Exception as e:
            logging.error("Hi-C lookup failed for %s %s: %s", sample, base_id, e)
            continue
        assembly_prefix = f"{sample}.hic.{date}"
        for r1, r2 in pairs:
            rows.append(
                {
                    "sample": sample,
                    "type": "hic",
                    "date": date,
                    "assembly_prefix": assembly_prefix,
                    "fastq_1": r1,
                    "fastq_2": r2,
                }
            )
        hic_handled.add((sample, base_id))

    # HiFi
    for p in sorted(hifi_files):
        m = HIFI_RE.match(p.name)
        if not m:
            continue
        gd = m.groupdict()
        sample = gd["sample"]
        date = gd["date"]
        assembly_prefix = f"{sample}.hifi.{date}"
        rows.append(
            {
                "sample": sample,
                "type": "hifi",
                "date": date,
                "assembly_prefix": assembly_prefix,
                "fastq_1": str(p),
                "fastq_2": "",
            }
        )

    # Write CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["sample", "type", "date", "assembly_prefix", "fastq_1", "fastq_2"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    logging.info("Wrote %d rows to %s", len(rows), out_csv)


def main():
    p = argparse.ArgumentParser(
        description="Build nf-core samplesheet CSV for ilmn/hic/hifi datasets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-dir", required=True, type=Path, help="Directory to scan")
    p.add_argument("--out", required=True, type=Path, help="Output CSV path")
    p.add_argument(
        "--s3-bucket",
        default="oceanomics",
        help="S3 bucket name containing Hi-C datasets",
    )
    p.add_argument(
        "--s3-prefix",
        default="OceanGenomes/illumina-hic",
        help="S3 key prefix under the bucket where Hi-C runs are stored",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase log verbosity (-v, -vv)",
    )

    args = p.parse_args()

    logging.basicConfig(
        level=logging.WARNING - (10 * min(args.verbose, 2)),
        format="[%(levelname)s] %(message)s",
    )

    build_samplesheet(args.input_dir, args.out, args.s3_bucket, args.s3_prefix)


if __name__ == "__main__":
    main()
