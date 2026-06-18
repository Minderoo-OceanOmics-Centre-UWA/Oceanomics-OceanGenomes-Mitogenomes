// Relabel the reference GenBank from MITOHIFI_FINDMITOREFERENCE to a per-sample
// filename (<assembly_prefix>.reference.gb) so it can be staged into
// MITOGENOME_ASSEMBLY_SUMMARY's flat collect() without colliding with other
// samples' identically-named reference records (every HiFi sample emits an
// identically-named `MitoReference` dir, and two samples can share the same NCBI
// accession). The assembly summary reads reference species/accession from this file.
process RELABEL_REFERENCE_GB {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(reference_gb)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference.gb"), emit: gb

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp -L ${reference_gb} ${meta.mt_assembly_prefix}.reference.gb
    """

    stub:
    """
    touch ${meta.mt_assembly_prefix}.reference.gb
    """
}
