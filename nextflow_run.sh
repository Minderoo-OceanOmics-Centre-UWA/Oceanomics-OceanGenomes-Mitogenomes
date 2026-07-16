module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Get the absolute path to the current directory
RUN_DIR="$(pwd -P)"
BASE="/scratch/pawsey1348/$USER"
OUT_DIR="$BASE/mitogenomes-missing-audit-3"

mkdir -p "$OUT_DIR"

OUT_DIR="$(cd "$OUT_DIR" && pwd -P)"
# Change to output directory to run Nextflow there
cd $OUT_DIR

nextflow -log $OUT_DIR/nextflow.log \
    run $RUN_DIR/main.nf \
    -work-dir ./work \
    -c $RUN_DIR/pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --input_dir "/scratch/pawsey1348/$USER/path/to/fastp/*.fastq.gz" \
    --outdir "$OUT_DIR" \
    --blast_db_dir "$(realpath ../blast_dbs)" \
    --taxonkit_db_dir "$(realpath ../)" \
    --curated_blast_db /software/projects/pawsey0964/curated_db/OceanGenomes.CuratedNT.NBDLTranche1and2and3.CuratedBOLD.NoDuplicate.fasta \
    --nt_blast_db /scratch/references/blastdb_update/blast-2026-02-01/db/mito \
    --mitos_refdb /software/projects/pawsey0964/mitos_refdb \
    --mitos_refseq_ver refseq89m  \
    --organelle_type "animal_mt" \
    --kvalue "21" \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
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
    --template_sbt "bin/template.sbt" \
    --force_db_overwrite false \
    --translation_table "2" \
    --ena_webin_validate true \
    --ena_study "PRJEB110568" \
    --ena_validation_attempt "initial"
    
    # --getorganelle_fromreads_args "-R 20 -w 75 -k 21,45,65,85,105 --max-extending-len inf --max-n-words 1000000000 --continue" ## Include this line if you want to customise the getorganelle fromreads args.
    # nextflow -log test_nextflow.log  ### replace the top line with this if you want to define the log file, if youre running multiple runs of the nf-core
    # -work-dir /scratch/pawsey0964/$USER/directory \. ### include work dir if you want to run this nf-core on multiple occasions and keep the work files separate.
    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
