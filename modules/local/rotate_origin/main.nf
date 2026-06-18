process ROTATE_ORIGIN {
    tag "$meta.id"
    label 'process_single'

    // cox1 search (tblastn) + FASTA rotation both run with Biopython; reuse the
    // pinned MITOS BioContainer (it ships blast+python) so the invert path needs
    // no extra image.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path cox1_ref

    output:
    // Re-origined assembly handed to MITOS2. cox1-anchored rotation lifts the
    // coral nad5 group I intron (and cox3) off the linearisation point so MITOS
    // annotates them in order without splitting across position 1.
    tuple val(meta), path("rotated/${fasta.baseName}.fa"), emit: fasta
    tuple val(meta), path("06_rotate.tool_params_mqcrow.html"), emit: tool_params
    path "versions_rotate.yml", emit: versions

    script:
        """
        mkdir -p rotated

        rotate_to_cox1.py \\
            --genome ${fasta} \\
            --cox1-ref ${cox1_ref} \\
            --out rotated/${fasta.baseName}.fa

        cat <<-END_TOOL_PARAMS > 06_rotate.tool_params_mqcrow.html
        <tr><td>Rotate origin</td><td><samp>rotate_to_cox1.py --genome ${fasta} --cox1-ref ${cox1_ref}</samp></td><td>Re-origins the circular assembly for ${meta.id} to the cox1 start (tblastn vs a coral cox1 panel) so the nad5 group I intron and cox3 no longer straddle position 1 before MITOS2. Passes through unrotated if no confident cox1 hit.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions_rotate.yml
        "${task.process}":
            tblastn: \$(tblastn -version 2>/dev/null | head -n1 | sed 's/^tblastn: //')
            python: \$(python --version | sed 's/^Python //')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p rotated
        cp ${fasta} rotated/${fasta.baseName}.fa

        cat <<-END_TOOL_PARAMS > 06_rotate.tool_params_mqcrow.html
        <tr><td>Rotate origin</td><td><samp>rotate_to_cox1.py --genome ${fasta} --cox1-ref ${cox1_ref}</samp></td><td>Re-origins the circular assembly for ${meta.id} to the cox1 start before MITOS2.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions_rotate.yml
        "${task.process}":
            tblastn: "2.1.10"
            python: "3.10.0"
        END_VERSIONS
        """
}
