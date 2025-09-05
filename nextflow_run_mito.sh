module load nextflow/24.10.0
module load singularity/4.1.0-nompi

/scratch/pawsey0964/tpeirce/nextflow-24.10.5-dist -log test_nextflow.log \
    run main.nf \
    -work-dir /scratch/pawsey0964/tpeirce/_NFCORE/work_mito \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --input_dir "/scratch/pawsey0964/tpeirce/_NFCORE/test_dir/" \
    --outdir /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mito_test \
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
    --samplesheet_prefix "samplesheet"
    
    
   
    #--input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
