module load nextflow/25.04.6
module load singularity/4.1.0-nompi

nextflow \
    run main.nf \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --outdir "$(realpath ../outdir)" \
    --blast_db_dir "$(realpath ../blast_dbs)" \
    --taxonkit_db_dir "$(realpath ../)" \
    --curated_blast_db /software/projects/pawsey0964/curated_db/OceanGenomes.CuratedNT.NBDLTranche1and2and3.CuratedBOLD.NoDuplicate.fasta \
    --nt_blast_db /scratch/references/blastdb_update/blast-2026-02-01/db/core_nt \
    --organelle_type "animal_mt" \
    --kvalue "21" \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/$USER/tmp \
    --refresh-modules \
    --skip_mitogenome_assembly_getorg false \
    --skip_mitogenome_assembly_hifi false \
    --skip_mitogenome_annotation false \
    --skip_upload_results true \
    --samplesheet_prefix "samplesheet" \
    --translation_table "2"\
    --input assets/samplesheet.csv \
    
    
    # --getorganelle_fromreads_args "-R 20 -w 75 -k 21,45,65,85 --max-extending-len inf" ## Include this line if you want to customise the getorganelle fromreads args.
    # nextflow -log test_nextflow.log  ### replace the top line with this if you want to define the log file, if youre running multiple runs of the nf-core
    # -work-dir /scratch/pawsey0964/$USER/directory \. ### include work dir if you want to run this nf-core on multiple occasions and keep the work files separate.
    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    # --input_dir "/scratch/pawsey0964/tpeirce/_NFCORE/OG1288/*.fastq.gz" \ # include the input dir if you are working with OceanGenomes files directly rather than a samplesheet. The create samplesheet script uses the OceanGenomes database to create the samplesheet from the fastq files.
