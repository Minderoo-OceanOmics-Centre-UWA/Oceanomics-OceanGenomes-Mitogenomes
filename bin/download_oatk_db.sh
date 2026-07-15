#!/usr/bin/env bash
# Download an OatkDB mitochondrial profile-HMM database (a .fam plus its pressed
# nhmmer index files) for use as --oatk_mito_db in the Oatk reference-free fallback.
#
# Usage:
#   bash bin/download_oatk_db.sh <db_name> <dest_dir> [oatkdb_version]
#
# Examples:
#   bash bin/download_oatk_db.sh actinopterygii_mito /scratch/$USER/oatk_db
#   bash bin/download_oatk_db.sh aves_mito           /scratch/$USER/oatk_db v20230921
#
# Common vertebrate mito databases: actinopterygii_mito (ray-finned fish),
# amphibia_mito, aves_mito, crocodylia_mito. Browse all names at
# https://github.com/c-zhou/OatkDB/tree/main/<version>
#
# Then point the pipeline at the .fam:
#   --oatk_mito_db <dest_dir>/<db_name>.fam
set -euo pipefail

DB_NAME="${1:?db name required, e.g. actinopterygii_mito}"
DEST_DIR="${2:?destination directory required}"
VERSION="${3:-v20230921}"
BASE="https://raw.githubusercontent.com/c-zhou/OatkDB/main/${VERSION}"

mkdir -p "$DEST_DIR"
echo "[download_oatk_db] fetching ${DB_NAME} (${VERSION}) -> ${DEST_DIR}"
for ext in fam fam.h3f fam.h3i fam.h3m fam.h3p; do
    f="${DB_NAME}.${ext}"
    echo "  - ${f}"
    curl -fsSL -o "${DEST_DIR}/${f}" "${BASE}/${f}"
done

echo "[download_oatk_db] done. Use:"
echo "    --oatk_mito_db ${DEST_DIR}/${DB_NAME}.fam"
