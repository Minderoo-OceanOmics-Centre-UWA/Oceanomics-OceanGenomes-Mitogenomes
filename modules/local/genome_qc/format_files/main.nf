process FORMAT_FILES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), val(species_name), val(proceed_qc), path(input_files)

    output:
    tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), emit: processed_files
    tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.gff"), emit: gff
    tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), emit: gb
    val meta        , emit: meta
    path "versions.yml", emit: versions

    // when:
    // proceed_qc == true

    /*
     * Notes:
     * - The helper script lives in ${moduleDir}/assets/process_files.py
     * - We pass metadata via CLI args (safer than embedding into heredocs).
     * - All staged input files are in CWD; the script scans the directory.
     */
    script:
    def args = task.ext.args ?: ''
    def species = species_name ?: 'Unknown species'

    """
    mkdir -p processed

    # All input files are staged in the current working directory
    # Pass the current directory as input-dir since files are staged here
    process_files.py \\
        $args \\
        --og-id "${meta.id}" \\
        --species "${species}" \\
        --input-dir . \\
        --outdir processed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        process_files: "1.0.0"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p processed
    touch processed/dummy.fa
    touch processed/dummy.gff
    touch processed/dummy.tbl
    touch processed/dummy.cmt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        process_files: "1.0.0"
    END_VERSIONS
    """
}