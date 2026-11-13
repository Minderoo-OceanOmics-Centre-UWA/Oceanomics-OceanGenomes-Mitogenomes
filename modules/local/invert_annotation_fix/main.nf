// Repair a broken MITOS2 invertebrate annotation using the assembly alone, then
// re-run the EMMA adapter (mitos_to_emma.py) so the gff / cds/ / proteins/
// standardisation, splice-join, translation and trnM re-origin are all reused.
//
// This is the LINEAGE-AGNOSTIC half of what used to be coral-only. The coral
// fixer's ORF-boundary snap needs no reference, no taxonomy and no coral, but it
// was reachable only down the coral branch -- so in NOVA_260909_LA all 17
// assemblies carrying a genuinely broken CDS (13 Porifera, plus Malacostraca and
// Asteroidea) had no repair path at all, while five Demospongiae shared one
// systematic MITOS boundary miscall between them.
//
// Runs BEFORE CORAL_ANNOTATION_FIX, which now does only what is genuinely coral:
// BLAST transfer of rrnL / nad5 / cox1 / mtMutS from an Anthozoa reference.
//
// Emits the SAME channels as MITOS2 so the annotation subworkflow can mix
// repaired and untouched samples without rewiring.
//
// Fail-safe: invert_fix_bed.py consults no reference, so it cannot import a gene
// the assembly does not carry -- the worst it can do is decline to move a
// boundary, which it records with a reason. It always exits 0, so a sample it
// cannot repair reproduces the original MITOS2 output and is held downstream by
// the annotation QC step rather than breaking the run.
process INVERT_ANNOTATION_FIX {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // genome = the rotated fasta MITOS2 annotated (ROTATE_ORIGIN.out.fasta),
    // bed    = MITOS2's raw result.bed.
    // origin_gene = the MITOS feature key the published genome is re-origined to,
    // resolved per taxon in mitogenome_annotation_lca. Must be the same anchor
    // MITOS2 was given for this sample, or a repaired assembly would be published
    // on a different origin from an unrepaired batch-mate with nothing to flag it.
    // mitos_to_emma.py makes --origin-gene REQUIRED so that divergence fails loudly.
    tuple val(meta), path(genome), path(bed), val(origin_gene)

    output:
    tuple val(meta), path("annotation/cds/*CO1*.fa"),  emit: co1_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR1*.fa"), emit: s12_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR2*.fa"), emit: s16_sequences, optional: true
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("annotation/mitos_fix/result.fixed.bed"), emit: bed
    tuple val(meta), path("annotation/*.gff"), path("annotation/proteins"), path("annotation/cds"), emit: gff_proteins
    tuple val(meta), path("annotation/mitos_fix/*.invert_fix.status.tsv"), emit: status
    tuple val(meta), path("annotation/mitos_conversion_audit.tsv"), emit: conversion_audit, optional: true
    tuple val(meta), path("07_invert_fix.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        def prefix   = task.ext.prefix ?: meta.mt_assembly_prefix
        def gcode    = task.ext.code ?: meta.genetic_code
        def species  = meta.species ?: ''
        def base_args = (task.ext.args ?: '').toString().trim()
        def mitos_tag = '2110'
        def topology_arg = (meta.circular == false) ? '--linear' : ''
        def mitos_prefix = "${prefix}.mitos${mitos_tag}"
        // Class only WIDENS the accepted initiation codons (orf_utils.
        // LINEAGE_ALT_STARTS): decapod CO1 initiates on CGA, which no translation
        // table lists. Omitting it can make the fixer attempt a repair it need
        // not, never the reverse, so a blank class is safe.
        def lineage_arg = meta.class ? "--lineage '${meta.class}'" : ''
        def origin_arg = "--origin-gene ${origin_gene}"
        def effective_args = ["invert_fix_bed.py --bed ${bed} --genome ${genome} ${topology_arg} --code ${gcode} ${lineage_arg} ${base_args}".replaceAll(/ +/, ' ').trim(),
                              "mitos_to_emma.py --bed result.fixed.bed --genome ${genome} --code ${gcode} ${topology_arg} ${origin_arg}".replaceAll(/ +/, ' ').trim()].join('; ')
        """
        mkdir -p annotation

        # 1) Snap each PCG's boundaries to the nearest clean ORF under this
        #    sample's own genetic code. Assembly only -- no reference is read.
        invert_fix_bed.py \\
            --bed ${bed} \\
            --genome ${genome} \\
            --code ${gcode} \\
            ${lineage_arg} \\
            --out-bed result.fixed.bed \\
            --status ${prefix}.invert_fix.status.tsv \\
            ${topology_arg} \\
            ${base_args}

        # 2) Re-run the EMMA adapter on the patched BED (joins, translation,
        #    cds/proteins extraction and the trnM re-origin are all reused).
        mitos_to_emma.py \\
            --bed result.fixed.bed \\
            --genome ${genome} \\
            --prefix "${mitos_prefix}" \\
            --outdir annotation \\
            --code ${gcode} \\
            --species "${species}" \\
            ${topology_arg} \\
            ${origin_arg}

        # Fix artefacts go in their own subdir so MITOS2's raw output
        # (annotation/mitos_raw/, published separately for this same sample) is
        # left untouched -- same convention as CORAL_ANNOTATION_FIX.
        mkdir -p annotation/mitos_fix
        cp result.fixed.bed annotation/mitos_fix/result.fixed.bed
        cp ${prefix}.invert_fix.status.tsv annotation/mitos_fix/${mitos_prefix}.invert_fix.status.tsv

        fix_status=\$(head -1 ${prefix}.invert_fix.status.tsv | cut -f1)
        cat <<-END_TOOL_PARAMS > 07_invert_fix.tool_params_mqcrow.html
        <tr><td>Invertebrate Annotation Fix</td><td><samp>${effective_args}</samp></td><td>Reference-free ORF boundary repair of the MITOS2 annotation for ${meta.id} (status \${fix_status}), then re-running the EMMA adapter.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitos: 2.1.10
            python: \$(python --version | sed 's/^Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
        def mitos_prefix = "${prefix}.mitos2110"
        """
        mkdir -p annotation/cds annotation/proteins annotation/mitos_fix
        touch annotation/${mitos_prefix}.fa
        touch annotation/${mitos_prefix}.gff
        touch annotation/mitos_conversion_audit.tsv
        touch annotation/cds/CO1.${mitos_prefix}.fa
        touch annotation/cds/RNR1.${mitos_prefix}.fa
        touch annotation/cds/RNR2.${mitos_prefix}.fa
        touch annotation/proteins/MT-CO1.${mitos_prefix}.fa
        touch annotation/mitos_fix/result.fixed.bed
        printf 'FIXED\\tstub\\n' > ${prefix}.invert_fix.status.tsv
        cp ${prefix}.invert_fix.status.tsv annotation/mitos_fix/${mitos_prefix}.invert_fix.status.tsv
        : > 07_invert_fix.tool_params_mqcrow.html
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitos: "2.1.10"
            python: "3.10.0"
        END_VERSIONS
        """
}
