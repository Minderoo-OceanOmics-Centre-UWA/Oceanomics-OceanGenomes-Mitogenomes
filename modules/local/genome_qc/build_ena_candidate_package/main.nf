// Where publishDir will place this process's `package` output, resolved to an
// absolute path so it stays meaningful to a reader that is neither this run nor
// on this machine. Kept next to the process so the two cannot drift apart
// without someone editing this file.
def publishedPackageDir(outdir, meta) {
    if (!outdir) {
        return ''
    }
    def root = file(outdir).toAbsolutePath().normalize()
    return "${root}/mitogenomes/${meta.id}/${meta.mt_assembly_prefix}/ena/package"
}

process BUILD_ENA_CANDIDATE_PACKAGE {
    tag "$meta.full_seqid"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(sample_fa), path(sample_gff), path(embl_file), path(locus_map), path(tagged_tbl), path(input_metadata)

    output:
    tuple val(meta), path("package"), emit: package_dir
    tuple val(meta), path("package/*.package_metadata.json"), emit: metadata
    tuple val(meta), path("package/*.local_validation.tsv"), emit: validation
    tuple val(meta), path("22_ena_candidate_package.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    script:
    // Must track the publishDir target for this process in conf/modules.config.
    // The build runs in a work directory that gets cleaned, so the durable
    // location has to be recorded rather than discovered later.  outdir has no
    // default, so an unset one records nothing rather than failing the build.
    def publish_root = publishedPackageDir(params.outdir, meta)
    def publish_arg = publish_root ? "--publish-root '${publish_root}'" : ''
    """
    ena_package.py build \
        --metadata-input '${input_metadata}' \
        --fasta '${sample_fa}' \
        --embl '${embl_file}' \
        --locus-map '${locus_map}' \
        --tagged-tbl '${tagged_tbl}' \
        --gff '${sample_gff}' \
        ${publish_arg} \
        --outdir package
    printf '%s\n' '<tr><td>ENA candidate package</td><td><samp>ena_package.py build</samp></td><td>Builds a self-contained genome-context candidate package for ${meta.full_seqid}.</td></tr>' > 22_ena_candidate_package.tool_params_mqcrow.html
    printf '"%s":\n    python: "%s"\n    ena_package: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    def publish_root = publishedPackageDir(params.outdir, meta)
    """
    mkdir -p package
    cp '${embl_file}' 'package/${meta.full_seqid}.embl.gz'
    cp '${locus_map}' 'package/${meta.full_seqid}.locus_tag_mapping.tsv'
    cp '${tagged_tbl}' 'package/${meta.full_seqid}.tbl'
    cp '${sample_fa}' 'package/${meta.full_seqid}.fa'
    cp '${sample_gff}' 'package/${meta.full_seqid}.gff'
    printf '{"full_seqid":"${meta.full_seqid}","og_id":"${meta.id}","package_status":"READY","local_validation_status":"PASS","sequence_sha256":"stub","normalised_circular_sha256":"stub","published_package_path":"${publish_root}"}\n' > 'package/${meta.full_seqid}.package_metadata.json'
    printf 'full_seqid\tstatus\terrors\n${meta.full_seqid}\tPASS\t\n' > 'package/${meta.full_seqid}.local_validation.tsv'
    printf '%s\n' '<tr><td>ENA candidate package</td><td><samp>stub</samp></td><td>Stub candidate package for ${meta.full_seqid}.</td></tr>' > 22_ena_candidate_package.tool_params_mqcrow.html
    printf '"%s":\n    python: "stub"\n    ena_package: "stub"\n' "${task.process}" > versions.yml
    """
}
