nextflow run ena.nf \
      -c pawsey_profile.config \
      -profile singularity \
      -resume \
      -with-report \
      --ena_mode convert_validate \
      --ena_input "/path/to/ena_gbf_inputs.csv" \
      --ena_study "PRJEB110568" \
      --ena_validation_attempt "initial" \
      --outdir "$(realpath ../ena_results)"