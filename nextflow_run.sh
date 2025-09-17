module load nextflow/24.10.0
module load singularity/4.1.0-nompi

/scratch/pawsey0964/tpeirce/nextflow-24.10.5-dist \
    run main.nf \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --input_dir "/scratch/pawsey0964/tpeirce/_NFCORE/test_dir/*.fastq.gz" \
    --outdir "$(realpath ../_outdir)" \
    --blast_db_dir "$(realpath ../blast_dbs)" \
    --taxonkit_db_dir "$(realpath ../)" \
    --curated_blast_db /scratch/pawsey0964/pbayer/OceanGenomes.CuratedNT.NBDLTranche1and2.CuratedBOLD.fasta \
    --organelle_type "animal_mt" \
    --kvalue "21" \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    --refresh-modules \
    --skip_mitogenome_assembly_getorg false \
    --skip_mitogenome_assembly_hifi false \
    --skip_mitogenome_annotation false \
    --skip_upload_results false \
    --samplesheet_prefix "samplesheet" \
    --template_sbt "/home/tpeirce/template.sbt" \
    --translation_table "2"
    
    
    # /scratch/pawsey0964/tpeirce/nextflow-24.10.5-dist -log test_nextflow.log  ### replace the top line with this if you want to define the log file, if youre running multiple runs of the nf-core
    # -work-dir /scratch/pawsey0964/tpeirce/_NFCORE/work_mito_codexupdated \. ### include work dir if you want to run this nf-core on multiple occasions and keep the work files separate.
    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
