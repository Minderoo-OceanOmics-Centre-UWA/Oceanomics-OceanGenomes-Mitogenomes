// Repair a deficient anthozoan MITOS2 annotation by transferring the 16S rRNA
// and the intron-split nad5 from a close coral reference (GenBank) at the BED
// level, then re-running the existing EMMA adapter (mitos_to_emma.py) so the
// gff / cds/ / proteins/ standardisation, splice-join, translation and trnM
// re-origin are all reused unchanged.
//
// Emits the SAME channels as MITOS2 so the annotation subworkflow can mix the
// repaired (FIX) and untouched (PASS) anthozoans without rewiring. Runs in the
// MITOS2 BioContainer, which already provides blastn + biopython.
//
// Fail-safe: coral_fix_bed.py guards every edit with BLAST coverage/identity and
// a nad5 ORF check, and writes the (possibly unchanged) BED plus a status line;
// it always exits 0. So a poor or distant reference reproduces the original
// MITOS2 output rather than breaking the run.
process CORAL_ANNOTATION_FIX {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // genome = the cox1-rotated fasta MITOS2 annotated (ROTATE_ORIGIN.out.fasta),
    // bed = MITOS2's raw result.bed, reference_gb = the resolved coral reference.
    tuple val(meta), path(genome), path(bed), path(reference_gb)

    output:
    tuple val(meta), path("annotation/cds/*CO1*.fa"),  emit: co1_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR1*.fa"), emit: s12_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR2*.fa"), emit: s16_sequences, optional: true
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("annotation/mitos_fix/result.fixed.bed"), emit: bed
    tuple val(meta), path("annotation/*.gff"), path("annotation/proteins"), emit: gff_proteins
    // Fix artefacts published under annotation/mitos_fix/ (declared explicitly so
    // publishDir copies them out of the work dir, not just left as work-dir copies).
    tuple val(meta), path("annotation/mitos_fix/*.coral_fix.status.txt"), emit: status
    tuple val(meta), path("annotation/mitos_fix/*.reference.gb"), emit: reference
    tuple val(meta), path("08_coral_fix.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        def prefix   = task.ext.prefix ?: meta.mt_assembly_prefix
        def gcode    = task.ext.code ?: (meta.genetic_code ?: 4)
        def species  = meta.species ?: ''
        def base_args = (task.ext.args ?: '').toString().trim()
        def mitos_tag = '2110'
        def topology_arg = (meta.circular == false) ? '--linear' : ''
        def mitos_prefix = "${prefix}.mitos${mitos_tag}"
        def effective_args = ["coral_fix_bed.py --bed ${bed} --genome ${genome} --ref-gb ${reference_gb} --code ${gcode} ${base_args}".replaceAll(/ +/, ' ').trim(),
                              "mitos_to_emma.py --bed result.fixed.bed --genome ${genome} --code ${gcode} ${topology_arg}".replaceAll(/ +/, ' ').trim()].join('; ')
        """
        mkdir -p annotation

        # 1) Patch the MITOS BED: BLAST-transfer 16S + nad5 exons from the reference.
        coral_fix_bed.py \\
            --bed ${bed} \\
            --genome ${genome} \\
            --ref-gb ${reference_gb} \\
            --code ${gcode} \\
            --out-bed result.fixed.bed \\
            --status ${prefix}.coral_fix.status.txt \\
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
            ${topology_arg}

        # Provenance: collect the fix artefacts in their OWN subdir (mitos_fix/) so
        # MITOS2's raw output (annotation/mitos_raw/, published separately for this
        # same sample) is left untouched. Previously this module wrote the patched
        # BED to annotation/mitos_raw/result.bed, which clobbered MITOS2's whole raw
        # dir on publish (copy mode, same dest path) — losing result.faa/fas/gff/etc.
        # mitos_fix/ holds the patched BED, the fix status and the reference used so
        # the repair stays auditable alongside the untouched raw output.
        mkdir -p annotation/mitos_fix
        cp result.fixed.bed annotation/mitos_fix/result.fixed.bed
        cp ${prefix}.coral_fix.status.txt annotation/mitos_fix/${mitos_prefix}.coral_fix.status.txt
        cp -L ${reference_gb} annotation/mitos_fix/${mitos_prefix}.reference.gb

        ref_used=\$(basename ${reference_gb})
        fix_status=\$(cut -f1 ${prefix}.coral_fix.status.txt)
        cat <<-END_TOOL_PARAMS > 08_coral_fix.tool_params_mqcrow.html
        <tr><td>Coral Annotation Fix</td><td><samp>${effective_args}</samp></td><td>Repairs the MITOS2 anthozoan annotation for ${meta.id} (status \${fix_status}) by BLAST-transferring 16S + nad5 from reference \${ref_used}, then re-running the EMMA adapter.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitos: 2.1.10
            blast: \$(blastn -version 2>/dev/null | sed -n 's/^blastn: //p')
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
        touch annotation/cds/CO1.${mitos_prefix}.fa
        touch annotation/cds/RNR1.${mitos_prefix}.fa
        touch annotation/cds/RNR2.${mitos_prefix}.fa
        touch annotation/proteins/MT-CO1.${mitos_prefix}.fa
        touch annotation/mitos_fix/result.fixed.bed
        printf 'FIXED\\tstub\\n' > ${prefix}.coral_fix.status.txt
        cp ${prefix}.coral_fix.status.txt annotation/mitos_fix/${mitos_prefix}.coral_fix.status.txt
        touch annotation/mitos_fix/${mitos_prefix}.reference.gb
        : > 08_coral_fix.tool_params_mqcrow.html
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitos: "2.1.10"
            python: "3.10.0"
        END_VERSIONS
        """
}
