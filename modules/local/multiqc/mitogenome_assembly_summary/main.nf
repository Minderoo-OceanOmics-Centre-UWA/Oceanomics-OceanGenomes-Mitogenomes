process MITOGENOME_ASSEMBLY_SUMMARY {
    tag "mitogenome_assembly_summary"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path assembly_files, stageAs: "mitogenome_summary_inputs/*"

    output:
    path "mitogenome_assembly_summary_mqc.tsv", emit: table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def minCoverage = params.mitogenome_summary_min_mean_coverage != null ? "--min-mean-coverage ${params.mitogenome_summary_min_mean_coverage}" : ''
    def maxCoverageCv = params.mitogenome_summary_max_coverage_cv != null ? "--max-coverage-cv ${params.mitogenome_summary_max_coverage_cv}" : ''
    def minLength = params.mitogenome_summary_min_length != null ? "--min-length ${params.mitogenome_summary_min_length}" : ''
    def maxLength = params.mitogenome_summary_max_length != null ? "--max-length ${params.mitogenome_summary_max_length}" : ''
    def expectedGenes = params.mitogenome_summary_expected_gene_count != null ? "--expected-gene-count ${params.mitogenome_summary_expected_gene_count}" : ''
    """
    mitogenome_assembly_summary.py \\
        --input mitogenome_summary_inputs \\
        --output mitogenome_assembly_summary_mqc.tsv \\
        ${minCoverage} \\
        ${maxCoverageCv} \\
        ${minLength} \\
        ${maxLength} \\
        ${expectedGenes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitogenome_assembly_summary: "1.0.0"
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf "sample_id\\tassembler\\tstatus\\tfinal_length_bp\\tcircularised\\tnum_candidate_contigs\\tnum_final_contigs\\tnum_genes\\tmissing_genes\\tframeshift_flag\\tmean_coverage\\tcoverage_cv\\treference_species\\treference_accession\\tnumt_flag\\tmanual_review_reason\\n" > mitogenome_assembly_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitogenome_assembly_summary: "stub"
    END_VERSIONS
    """
}
