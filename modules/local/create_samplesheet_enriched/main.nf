process CREATE_SAMPLESHEET_ENRICHED {
    tag "Creating enriched samplesheet from ${input_files}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    path(input_files)
    val output_name
    path sql_config
    path taxdump

    output:
    path "${output_name}"             , emit: samplesheet
    path "taxonomy_resolution.tsv"    , emit: taxonomy_resolution
    path "versions.yml"               , emit: versions

    script:
    // The species table only holds curated taxa. Without the taxdump fallback a
    // miss leaves class='unknown', which downstream reads as "vertebrate".
    def taxdump_arg = taxdump ? "--taxdump-dir ${taxdump}" : ''
    """
    create_samplesheet.py \
        --output "${output_name}" \
        --sql-config "${sql_config}" \
        ${taxdump_arg} \
        --input-files ${input_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9"
        psycopg2: "2.9.5"
    END_VERSIONS
    """

    stub:
    // Header-only rather than empty: the sheet is parsed by samplesheetToList
    // downstream, which needs a header row to parse at all.
    """
    printf 'sample,sequencing_type,single_end,original_id,completion_date,date,assembly_prefix,nominal_species_id,reference_species_id,class,family,order,invertebrates,genetic_code,fastq_1,fastq_2\\n' > ${output_name}
    printf 'sample\\tnominal_species_id\\tclass\\tfamily\\torder\\treference_species_id\\tsource\\n' > taxonomy_resolution.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        psycopg2: "stub"
    END_VERSIONS
    """
}
