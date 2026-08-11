module load nextflow/25.04.6
module load singularity/4.1.0-nompi

nextflow run qc_only_from_annotations.nf \
  -profile singularity \
  --annotation_files "/scratch/pawsey1348/tpeirce/OG910/*/annotation/*.{fa,fasta,gff,tbl,gb}" \
  --sql_config ~/postgresql_details/oceanomics.cfg \
  --template_sbt "bin/template.sbt" \
  --outdir "/scratch/pawsey1348/tpeirce/OG910" \
  --ena_validate_webin_test true \
  --ena_validation_attempt "initial"
