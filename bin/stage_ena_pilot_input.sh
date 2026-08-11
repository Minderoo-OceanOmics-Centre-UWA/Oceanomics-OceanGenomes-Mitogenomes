#!/usr/bin/env bash
set -euo pipefail

readonly OUTPUT_DIR="${1:-/scratch/pawsey1348/${USER}/ena-pilot-input/OG111}"
readonly MAX_ATTEMPTS="${ENA_STAGE_MAX_ATTEMPTS:-144}"
readonly RETRY_SECONDS="${ENA_STAGE_RETRY_SECONDS:-600}"
readonly DATASET_L001="ds.35abcef6e88e495bbb1028a865279181"
readonly DATASET_L002="ds.6be108c0f599433b978193c5fff77175"

mkdir -p "${OUTPUT_DIR}/${DATASET_L001}" "${OUTPUT_DIR}/${DATASET_L002}"

fastq_count() {
    find "${OUTPUT_DIR}" -type f -name '*.fastq.gz' | wc -l
}

if (( $(fastq_count) >= 4 )); then
    echo "OG111 pilot input is already staged: ${OUTPUT_DIR}"
    exit 0
fi

bs unarchive dataset -i "${DATASET_L001}" --retry || true
bs unarchive dataset -i "${DATASET_L002}" --retry || true

for (( attempt=1; attempt<=MAX_ATTEMPTS; attempt++ )); do
    echo "BaseSpace staging attempt ${attempt}/${MAX_ATTEMPTS}"
    bs download dataset -i "${DATASET_L001}" \
        -o "${OUTPUT_DIR}/${DATASET_L001}" --summary --retry || true
    bs download dataset -i "${DATASET_L002}" \
        -o "${OUTPUT_DIR}/${DATASET_L002}" --summary --retry || true

    if (( $(fastq_count) >= 4 )); then
        touch "${OUTPUT_DIR}/.staged"
        echo "OG111 pilot input staged: ${OUTPUT_DIR}"
        find "${OUTPUT_DIR}" -type f -name '*.fastq.gz' -print
        exit 0
    fi

    if (( attempt < MAX_ATTEMPTS )); then
        sleep "${RETRY_SECONDS}"
    fi
done

echo "BaseSpace did not expose all four OG111 FASTQs within the retry window." >&2
exit 1
