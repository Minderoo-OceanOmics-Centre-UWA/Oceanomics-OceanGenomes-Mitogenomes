#!/usr/bin/env bash
set -euo pipefail

module load nextflow/25.04.6
module load singularity/4.1.0-nompi

if [[ $# -gt 2 ]]; then
    echo "Usage: $0 ['<OG111 FASTQ glob>'] [output-directory]" >&2
    exit 2
fi

readonly DEFAULT_INPUT_ROOT="/scratch/pawsey1348/${USER}/ena-pilot-input/OG111"
readonly INPUT_GLOB="${1:-${DEFAULT_INPUT_ROOT}/**/*.fastq.gz}"
readonly REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly PILOT_OUT="${2:-/scratch/pawsey1348/${USER}/ena-genome-pilot-OG111}"

if [[ "${INPUT_GLOB}" != *OG111* ]]; then
    echo "Pilot input must target OG111, the current SAMEA-backed specimen: ${INPUT_GLOB}" >&2
    exit 2
fi

if [[ "${INPUT_GLOB}" == "${DEFAULT_INPUT_ROOT}/**/*.fastq.gz" ]] &&
   (( $(find "${DEFAULT_INPUT_ROOT}" -type f -name '*.fastq.gz' 2>/dev/null | wc -l) < 4 )); then
    echo "OG111 pilot FASTQs are not fully staged yet." >&2
    echo "Check the staging job or run: ${REPO_DIR}/bin/stage_ena_pilot_input.sh" >&2
    exit 2
fi

mkdir -p "${PILOT_OUT}"
cd "${PILOT_OUT}"

nextflow -log "${PILOT_OUT}/nextflow.log" run "${REPO_DIR}/main.nf" \
    -work-dir "${PILOT_OUT}/work" \
    -c "${REPO_DIR}/pawsey_profile.config" \
    -profile singularity \
    -resume \
    -with-report \
    --input_dir "${INPUT_GLOB}" \
    --outdir "${PILOT_OUT}" \
    --blast_db_dir /scratch/pawsey1348/tpeirce/blast_dbs \
    --taxonkit_db_dir /scratch/pawsey1348/tpeirce \
    --curated_blast_db /software/projects/pawsey0964/curated_db/OceanGenomes.CuratedNT.NBDLTranche1and2and3.CuratedBOLD.NoDuplicate.fasta \
    --nt_blast_db /scratch/references/blastdb_update/blast-2026-07-01/db/mito \
    --mitos_refdb /software/projects/pawsey0964/mitos_refdb \
    --mitos_refseq_ver refseq89m \
    --organelle_type animal_mt \
    --kvalue 21 \
    --sql_config /home/tpeirce/postgresql_details/oceanomics.cfg \
    --enable_oatk_fallback true \
    --oatk_mito_db /software/projects/pawsey0964/oatk_db/actinopterygii_mito.fam \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    --template_sbt "${REPO_DIR}/bin/template.sbt" \
    --force_db_overwrite true \
    --ena_validate_webin_test true \
    --ena_validate_webin_production false \
    --ena_webin_validate false \
    --ena_validation_attempt "og111-pilot-$(date +%Y%m%d)"

echo "Pilot complete: ${PILOT_OUT}"
echo "Candidate packages: ${PILOT_OUT}/mitogenomes/OG111/*/ena/package"
