#!/usr/bin/env python3
"""Correct stale [mgcode=N] tags in already-published extraction output.

Gene and CDS extraction used to write a literal `[mgcode=2]` into every FASTA
header it emitted (bin/extract_genes_gff.py, bin/extract_cds_from_tbl.py), while
the genome record, the .tbl and table2asn all carried the sample's real
translation table. Any non-code-2 assembly published before that fix therefore has
extracted genes, CDS and proteins whose header contradicts its own genome record.

Only the LABEL is wrong. The sequences were always translated with the correct
code -- both protein-producing paths already took meta.genetic_code, and nothing
read the header tag back to decide a translation. Re-running the pipeline was
verified to produce output byte-identical to the published files once this tag is
normalised, so rewriting the token in place lands on exactly the same bytes a
re-run would, without re-invoking table2asn or ENA's Webin validation service.

The authoritative code for an assembly is read from that assembly's OWN genome
record (genbank/processed/*.fa, written by process_files.py from meta.genetic_code)
and feature table (annotation/*.tbl `transl_table`). It is never assumed: an
assembly whose sources are missing or disagree is skipped and reported, so a
future code-5 or code-9 sample cannot be silently relabelled 4.

Usage:
    # inventory only (default; writes nothing but the manifest)
    fix_published_mgcode.py --root /scratch/pawsey1348/tpeirce/batch-20 \
        --manifest /tmp/mgcode_manifest.tsv

    # apply, after taking a backup of exactly the manifest's files
    fix_published_mgcode.py --root ... --manifest ... --backup /scratch/.../bk.tar.gz --apply

    # confirm afterwards: re-inventory should find nothing left to do
    fix_published_mgcode.py --root ... --manifest /tmp/after.tsv
"""

import argparse
import hashlib
import re
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

# Sibling import: Nextflow bind-mounts the whole bin/ dir onto PATH, so this
# resolves via sys.path[0] (same pattern as bin/process_files.py).
from orf_utils import SUPPORTED_CODES

MGCODE = re.compile(rb"\[mgcode=(\d+)\]")
TRANSL_TABLE = re.compile(rb"transl_table[=\t ]+(\d+)")

# Where the authoritative per-assembly code is recorded. These are pipeline
# outputs that always carried meta.genetic_code correctly.
AUTHORITATIVE_GLOBS = (
    "genbank/processed/*.fa",
    "genbank/processed/*.fasta",
    "annotation/*.tbl",
    "genbank/*.tbl",
)

# The extraction outputs that carried the hardcoded tag. `ena/package/*.genes.fa`
# is the copy BUILD_ENA_CANDIDATE_PACKAGE stages beside the genome record, which
# is where a code-4 genome shipped next to code-2 gene records.
EXTRACTION_GLOBS = (
    "genbank/genes/*.fa",
    "genbank/proteins/*.fa",
    "genbank/cds/*.fa",
    "ena/package/*.genes.fa",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def authoritative_code(assembly: Path):
    """(code, note). code is None when it cannot be established beyond doubt."""
    codes = set()
    seen_any = False
    for pattern in AUTHORITATIVE_GLOBS:
        for path in sorted(assembly.glob(pattern)):
            try:
                data = path.read_bytes()
            except OSError:
                continue
            seen_any = True
            codes |= {int(m) for m in MGCODE.findall(data)}
            codes |= {int(m) for m in TRANSL_TABLE.findall(data)}
    if not codes:
        return None, "no genome record or feature table states a code" if seen_any \
            else "no genome record or feature table found"
    if len(codes) > 1:
        return None, f"sources disagree: {sorted(codes)}"
    code = codes.pop()
    if code not in SUPPORTED_CODES:
        return None, f"code {code} is not a supported mitochondrial table"
    return code, ""


def find_assemblies(root: Path):
    """An assembly dir is one holding an annotation/ subdirectory."""
    return sorted({p.parent for p in root.rglob("annotation") if p.is_dir()})


def scan(roots):
    """Return (todo, skipped). todo: (file, old_tag, new_code) needing a rewrite."""
    todo, skipped = [], []
    for root in roots:
        for assembly in find_assemblies(root):
            code, note = authoritative_code(assembly)
            if code is None:
                # Only worth reporting if this assembly actually has tagged output.
                if any(assembly.glob(p) for p in EXTRACTION_GLOBS):
                    skipped.append((assembly, note))
                continue
            for pattern in EXTRACTION_GLOBS:
                for path in sorted(assembly.glob(pattern)):
                    try:
                        data = path.read_bytes()
                    except OSError:
                        continue
                    tags = {int(m) for m in MGCODE.findall(data)}
                    if not tags or tags == {code}:
                        continue  # untagged, or already correct -> idempotent
                    todo.append((path, sorted(tags), code))
    return todo, skipped


def write_manifest(path: Path, todo, skipped):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"# generated {datetime.now(timezone.utc).isoformat()}\n")
        fh.write("file\told_tags\tnew_code\tsha256_before\n")
        for f, old, new in todo:
            fh.write(f"{f}\t{','.join(str(o) for o in old)}\t{new}\t{sha256(f)}\n")
        for assembly, note in skipped:
            fh.write(f"# SKIPPED\t{assembly}\t{note}\n")


def backup(todo, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(dest, "w:gz") as tar:
        for f, _old, _new in todo:
            tar.add(f, arcname=str(f))
    return dest


def rewrite(path: Path, code: int) -> bool:
    """Replace only the bracketed tag. Everything else is untouched by construction."""
    data = path.read_bytes()
    new = MGCODE.sub(b"[mgcode=%d]" % code, data)
    if new == data:
        return False
    path.write_bytes(new)
    return True


def refresh_checksums(package_dirs):
    """Regenerate ena/package/checksums.sha256 where a packaged file changed.

    The file verifies clean before this script runs, so patching a packaged
    genes.fa without refreshing it would leave a package failing its own
    integrity check. Rewrites only entries already listed, preserving order.
    """
    refreshed = []
    for pkg in sorted(package_dirs):
        manifest = pkg / "checksums.sha256"
        if not manifest.exists():
            continue
        lines = []
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            _old_hash, _, name = line.partition("  ")
            target = pkg / name
            lines.append(f"{sha256(target)}  {name}" if target.exists() else line)
        manifest.write_text("\n".join(lines) + "\n")
        refreshed.append(manifest)
    return refreshed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", action="append", required=True, type=Path,
                    help="published tree to scan (repeatable), e.g. .../batch-20")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="TSV record of every file that needs (or got) a new tag")
    ap.add_argument("--backup", type=Path,
                    help="tar.gz of exactly the manifest's files; required with --apply")
    ap.add_argument("--apply", action="store_true",
                    help="actually rewrite. Without it this only inventories.")
    ap.add_argument("--expect-files", type=int, default=None,
                    help="refuse to --apply unless exactly this many files need a change "
                         "(guards against the inventory drifting under you)")
    args = ap.parse_args()

    for root in args.root:
        if not root.is_dir():
            sys.exit(f"root not found: {root}")

    todo, skipped = scan(args.root)
    write_manifest(args.manifest, todo, skipped)

    assemblies = sorted({f.parent.parent.parent for f, _o, _n in todo})
    by_code = {}
    for _f, _o, new in todo:
        by_code[new] = by_code.get(new, 0) + 1

    print(f"assemblies needing a tag change : {len(assemblies)}")
    print(f"files needing a tag change      : {len(todo)}")
    print(f"target codes                    : {by_code or '-'}")
    print(f"manifest                        : {args.manifest}")
    if skipped:
        print(f"\nSKIPPED (code not established, nothing changed): {len(skipped)}")
        for assembly, note in skipped[:10]:
            print(f"  {assembly}: {note}")

    if todo:
        f, old, new = todo[0]
        print(f"\nexample: {f}")
        print(f"  {old} -> {new}")

    if not args.apply:
        print("\nDRY RUN -- nothing written. Re-run with --apply --backup <path> to change files.")
        return 0

    if args.expect_files is not None and len(todo) != args.expect_files:
        sys.exit(f"\nREFUSING: expected {args.expect_files} files, found {len(todo)}. "
                 "The inventory drifted -- re-check before applying.")
    if not todo:
        print("\nNothing to do.")
        return 0
    if not args.backup:
        sys.exit("--backup is required with --apply")

    print(f"\nbacking up {len(todo)} files -> {args.backup}")
    backup(todo, args.backup)
    print(f"backup written ({args.backup.stat().st_size / 1e6:.1f} MB)")

    changed = 0
    package_dirs = set()
    for f, _old, new in todo:
        if rewrite(f, new):
            changed += 1
            if f.parent.name == "package":
                package_dirs.add(f.parent)
    print(f"rewrote {changed} files")

    refreshed = refresh_checksums(package_dirs)
    print(f"refreshed {len(refreshed)} ena/package/checksums.sha256")
    return 0


if __name__ == "__main__":
    sys.exit(main())
