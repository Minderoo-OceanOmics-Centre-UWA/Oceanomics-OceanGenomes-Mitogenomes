module load nextflow/25.04.6
module load singularity/4.1.0-nompi

nextflow run qc_only_from_annotations.nf \
  -profile singularity \
  --annotation_files "/scratch/pawsey0964/tpeirce/__mount/OG2016/OG2016.ilmn.251215.getorg1770/emma/*.{fa,fasta,gff,tbl,gb}" \
  --sql_config ~/postgresql_details/oceanomics.cfg \
  --template_sbt "bin/template.sbt" \
  --outdir "$(realpath ../qc_outdir)" \
