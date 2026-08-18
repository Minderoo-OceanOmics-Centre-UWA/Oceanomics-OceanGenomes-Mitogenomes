process BUILD_ENA_CANDIDATE_PACKAGE {
    tag "$meta.full_seqid"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(sample_fa), path(sample_tbl), path(sample_gff), path(embl_file), path(input_metadata), path(flatfile_status)

    output:
    tuple val(meta), path("package"), emit: package_dir
    tuple val(meta), path("package/*.package_metadata.json"), emit: metadata
    tuple val(meta), path("22_ena_candidate_package.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    script:
    // Optional: absent when --ena_webin_validate is off, in which case the
    // package records the verdict as NOT_REQUESTED rather than inventing one.
    def flatfile_arg = flatfile_status ? "--flatfile-status '${flatfile_status}'" : ''
    """
    ena_package.py build \
        --metadata-input '${input_metadata}' \
        --fasta '${sample_fa}' \
        --embl '${embl_file}' \
        --tbl '${sample_tbl}' \
        --gff '${sample_gff}' \
        ${flatfile_arg} \
        --outdir package
    printf '%s\n' '<tr><td>ENA candidate package</td><td><samp>ena_package.py build</samp></td><td>Builds a self-contained genome-context candidate package for ${meta.full_seqid}.</td></tr>' > 22_ena_candidate_package.tool_params_mqcrow.html
    printf '"%s":\n    python: "%s"\n    ena_package: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    mkdir -p package
    cp '${embl_file}' 'package/${meta.full_seqid}.embl.gz'
    cp '${sample_tbl}' 'package/${meta.full_seqid}.tbl'
    cp '${sample_fa}' 'package/${meta.full_seqid}.fa'
    cp '${sample_gff}' 'package/${meta.full_seqid}.gff'
    printf '{"full_seqid":"${meta.full_seqid}","og_id":"${meta.id}","biosample_accession":"SAMEA1","sequence_sha256":"stub","normalised_circular_sha256":"stub","flatfile_validation":{"status":"PASS","reason":"validated","error_count":0,"warning_count":0,"webin_cli_version":"stub"}}\n' > 'package/${meta.full_seqid}.package_metadata.json'
    printf '%s\n' '<tr><td>ENA candidate package</td><td><samp>stub</samp></td><td>Stub candidate package for ${meta.full_seqid}.</td></tr>' > 22_ena_candidate_package.tool_params_mqcrow.html
    printf '"%s":\n    python: "stub"\n    ena_package: "stub"\n' "${task.process}" > versions.yml
    """
}
