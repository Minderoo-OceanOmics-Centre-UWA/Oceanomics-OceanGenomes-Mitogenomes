module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Per-run overrides. All optional; the defaults below reproduce a plain manual
# run of this script from the repo, so `bash nextflow_run.sh` is unchanged.
#   RUN_DIR        path to this repo            (default: this script's directory)
#   OG_OUT_DIR     run output directory
#   OG_WORK_DIR    nextflow work dir, relative to OG_OUT_DIR
#   OG_LOG_NAME    nextflow log filename, inside OG_OUT_DIR
#   OG_INPUT_DIR   fastp read glob passed to --input_dir
#   OG_BS_CONFIG   basespace config
#   OG_SQL_CONFIG  postgres config
# The Draft-Genomes pipeline's TRIGGER_MITOGENOME module generates
# run_<RUN>_mitogenomes_<timestamp>.sh by prepending those exports to a verbatim
# copy of this file. This is therefore the only place launcher changes need to be
# made: triggered runs pick them up automatically.

# Get the absolute path to this repo
RUN_DIR="${RUN_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
BASE="/scratch/pawsey1348/$USER"
OUT_DIR="${OG_OUT_DIR:-$BASE/mitogenomes-missing-audit-3}"

mkdir -p "$OUT_DIR"

OUT_DIR="$(cd "$OUT_DIR" && pwd -P)"

# Record what actually launched this run: a byte-exact copy of this script (for a
# triggered run that is the generated run_*.sh, overrides included), the repo
# commit, and the OG_* values in effect. This is the run's launcher provenance.
mkdir -p "$OUT_DIR/pipeline_info/launcher"
command cp "${BASH_SOURCE[0]}" "$OUT_DIR/pipeline_info/launcher/launcher.snapshot.sh"
git -C "$RUN_DIR" rev-parse HEAD > "$OUT_DIR/pipeline_info/launcher/mito_repo_commit.txt" 2>/dev/null \
    || echo "not a git checkout" > "$OUT_DIR/pipeline_info/launcher/mito_repo_commit.txt"
env | grep '^OG_' | sort > "$OUT_DIR/pipeline_info/launcher/og_overrides.env" || true

# Stage the backup scripts alongside the results and bake this run's outdir into
# them, so the backup is a no-argument `sbatch $OUT_DIR/backup_scripts/backup.sh`
# once the pipeline finishes.
mkdir -p "$OUT_DIR/backup_scripts"
command cp -r "$RUN_DIR/backup_scripts/." "$OUT_DIR/backup_scripts/"
sed -i "s|^RUNDIR_DEFAULT=.*|RUNDIR_DEFAULT=\"$OUT_DIR\"|" "$OUT_DIR/backup_scripts/backup.sh"

if ! grep -q "^RUNDIR_DEFAULT=\"$OUT_DIR\"$" "$OUT_DIR/backup_scripts/backup.sh"; then
    echo "WARNING: could not bake RUNDIR into $OUT_DIR/backup_scripts/backup.sh;" \
         "run it with -r $OUT_DIR"
fi

# Change to output directory to run Nextflow there
cd $OUT_DIR

nextflow -log "$OUT_DIR/${OG_LOG_NAME:-nextflow.log}" \
    run $RUN_DIR/main.nf \
    -work-dir "./${OG_WORK_DIR:-work}" \
    -c $RUN_DIR/pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --input_dir "${OG_INPUT_DIR:-/scratch/pawsey1348/$USER/path/to/fastp/*.fastq.gz}" \
    --outdir "$OUT_DIR" \
    --blast_db_dir "$(realpath ../blast_dbs)" \
    --taxonkit_db_dir "$(realpath ../)" \
    --curated_blast_db /software/projects/pawsey0964/curated_db/OceanGenomes.CuratedNT.NBDLTranche1and2and3.CuratedBOLD.NoDuplicate.fasta \
    --nt_blast_db /scratch/references/blastdb_update/blast-2026-07-01/db/mito \
    --mitos_refdb /software/projects/pawsey0964/mitos_refdb \
    --mitos_refseq_ver refseq89m  \
    --organelle_type "animal_mt" \
    --kvalue "21" \
    --bs_config "${OG_BS_CONFIG:-$HOME/.basespace/default.cfg}" \
    --sql_config "${OG_SQL_CONFIG:-$HOME/postgresql_details/oceanomics.cfg}" \
    --enable_oatk_fallback true \
    --oatk_mito_db /software/projects/pawsey0964/oatk_db/actinopterygii_mito.fam \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/$USER/tmp \
    --refresh-modules \
    --skip_mitogenome_assembly_getorg false \
    --skip_mitogenome_assembly_hifi false \
    --skip_mitogenome_annotation false \
    --skip_upload_results false \
    --samplesheet_prefix "samplesheet" \
    --template_sbt "$RUN_DIR/bin/template.sbt" \
    --force_db_overwrite true \
    --translation_table "2" \
    --ena_webin_validate true \
    --ena_study "PRJEB110568" \
    --ena_validation_attempt "initial"
    
    # --getorganelle_fromreads_args "-R 20 -w 75 -k 21,45,65,85,105 --max-extending-len inf --max-n-words 1000000000 --continue" ## Include this line if you want to customise the getorganelle fromreads args.
    # nextflow -log test_nextflow.log  ### replace the top line with this if you want to define the log file, if youre running multiple runs of the nf-core
    # -work-dir /scratch/pawsey0964/$USER/directory \. ### include work dir if you want to run this nf-core on multiple occasions and keep the work files separate.
    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.

# --- Per-genome compute cost (best-effort, self-contained) -----------------
# Repo-local script: writes pipeline_info/cost_per_sample.csv (SU per OG sample)
# for this run. Skipped when the poller runs this (it records cost, incl. the
# central ledger, via post_mito.py). The `|| echo` keeps a cost failure from
# affecting the run's exit status. Inherits this script's NXF_HOME + modules.
if [ -z "${OCEANOMICS_SKIP_COST:-}" ]; then
    COST_SCRIPT="$RUN_DIR/compute-audit/nf_workflow_cost.sh"

    if [ -f "$COST_SCRIPT" ]; then
        mkdir -p "$OUT_DIR/pipeline_info"

        COST_LOG="$OUT_DIR/pipeline_info/compute_cost.log"

        if bash "$COST_SCRIPT" \
            "$OUT_DIR" \
            "$OUT_DIR/pipeline_info/compute_usage.csv" \
            >"$COST_LOG" 2>&1
        then
            echo "compute cost: accounting complete"
        else
            echo "compute cost: accounting failed (non-fatal); see $COST_LOG"
        fi
    fi
fi
