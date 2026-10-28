// Store each assembly's QC verdict from the run's mitogenome_assembly_summary_mqc.tsv.
//
// One task per run, after MITOGENOME_ASSEMBLY_SUMMARY. push_assembly_qc.py updates the
// assembly rows push_mtdna_assm_results.py created and never creates one (genomes_db_review
// migration 068).
process PUSH_ASSEMBLY_QC {
    tag "assembly_qc"
    label 'process_upload'
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    path summary_tsv
    path config

    output:
    path "assembly_qc.upload.txt", emit: upload
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    push_assembly_qc.py \\
        ${config} \\
        ${summary_tsv} \\
        > assembly_qc.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_assembly_qc: "1.0.0"
    END_VERSIONS
    """

    stub:
    """
    : > assembly_qc.upload.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        push_assembly_qc: "stub"
    END_VERSIONS
    """
}
