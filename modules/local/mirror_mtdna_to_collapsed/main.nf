// Mirror the assembly-stage mtdna files (plus the collapse artefacts) of a
// genuinely collapsed sample into its <prefix>_collapsed/mtdna publish dir, so the
// curated variant directory carries the full provenance of the mitogenome
// (original assembly -> circularity check found the concatemer -> collapse decision
// -> corrected assembly) right next to the annotation of the collapsed sequence.
//
// Files are staged as real process inputs (not copied out of the publish dir), so
// Nextflow guarantees every upstream mtdna-writing task finished first and -resume
// stays correct. The caller sets meta.mt_assembly_prefix to the _collapsed prefix
// so the publishDir (see conf/modules.config) resolves to the collapsed variant dir.
process MIRROR_MTDNA_TO_COLLAPSED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(mtdna_files)

    output:
    tuple val(meta), path("mtdna_mirror/*"), emit: files

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p mtdna_mirror
    for f in ${mtdna_files}; do
        cp -L "\$f" mtdna_mirror/
    done
    """

    stub:
    """
    mkdir -p mtdna_mirror
    for f in ${mtdna_files}; do
        cp -L "\$f" mtdna_mirror/
    done
    """
}
