process SELECT_ENA_SUBMISSION {
    tag "selection:${selection_mode}"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    path package_metadata
    val selection_mode
    path decision_file
    path db_config
    val selected_by

    output:
    path "ena_selection_report.tsv", emit: report
    path "ena_selected_packages.tsv", emit: selected
    path "versions.yml", emit: versions

    script:
    def metadata_args = package_metadata.collect { "--package-metadata '${it}'" }.join(' ')
    def decision_arg = decision_file ? "--decision-file '${decision_file}'" : ''
    def db_arg = db_config ? "--db-config '${db_config}'" : ''
    """
    select_ena_submission.py \
        ${metadata_args} \
        --mode '${selection_mode}' \
        ${decision_arg} \
        ${db_arg} \
        --selected-by '${selected_by}' \
        --output ena_selection_report.tsv \
        --selected-output ena_selected_packages.tsv
    printf '"%s":\n    python: "%s"\n    select_ena_submission: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf 'og_id\ttech\tfull_seqid\tpackage_path\tmetadata_path\tpackage_status\tlocal_validation_status\tnormalised_circular_sha256\tselection_status\tselected\tselection_reason\n' > ena_selection_report.tsv
    printf 'og_id\ttech\tfull_seqid\tpackage_path\tmetadata_path\tstudy\tbiosample_accession\tselection_reason\n' > ena_selected_packages.tsv
    printf '"%s":\n    python: "stub"\n    select_ena_submission: "stub"\n' "${task.process}" > versions.yml
    """
}
