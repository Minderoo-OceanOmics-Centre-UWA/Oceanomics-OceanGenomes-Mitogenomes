// Decide whether an EMMA annotation is a candidate for the ND4L / ATP8 rescue.
// EMMA's rationalise_matches! step drops a short CDS when its computed circular
// overlap with a longer neighbour exceeds half the shorter feature's length,
// which routinely loses ND4L (vs ND4) and ATP8 (vs ATP6) from otherwise-complete
// vertebrate mitogenomes. Emits a one-line decision file (FIX\t<targets> | PASS\t-)
// that the annotation subworkflow branches on, so only the recoverable cases go
// to EMMA_GENE_RESCUE and everything else passes through untouched.
process EMMA_GENE_RESCUE_GATE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // The annotation bundle EMMA produced (EMMA.out.results); only the *.gff is read.
    tuple val(meta), path(annotation, stageAs: 'emma_in/*')

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.emma_rescue_qc.txt"), emit: decision
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Taxonomy for the curated gene-order variant lookup. Without it a clade whose
    // real gene order is non-canonical is judged out-of-order by the gate and
    // silently declined for rescue, so an assembly missing a rescuable gene stays
    // held for a reason this gate could have fixed. Same derivation as the emma
    // upload module and modules/local/reference_divergence.
    def klass       = (meta.class ?: '').toString().trim()
    def family      = (meta.family ?: '').toString().trim()
    def taxon_order = (meta.order ?: '').toString().trim()
    def genus       = (meta.nominal_species_id ?: '').toString().trim().split(/\s+/)[0] ?: ''
    def taxon_args  = [
        klass       ? "--class '${klass}'"       : '',
        family      ? "--family '${family}'"     : '',
        taxon_order ? "--order '${taxon_order}'" : '',
        genus       ? "--genus '${genus}'"       : '',
    ].findAll { it }.join(' ')
    """
    gff=\$(find emma_in -maxdepth 1 -name '*.gff' | head -n1)
    if [ -z "\$gff" ]; then
        printf 'PASS\\t-\\n' > ${meta.mt_assembly_prefix}.emma_rescue_qc.txt
    else
        emma_rescue_gate.py --gff "\$gff" ${taxon_args} \\
            --out ${meta.mt_assembly_prefix}.emma_rescue_qc.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'PASS\\t-\\n' > ${meta.mt_assembly_prefix}.emma_rescue_qc.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
