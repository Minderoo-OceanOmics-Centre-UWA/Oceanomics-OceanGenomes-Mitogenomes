// Stub harness for the assembly-upload row construction in
// workflows/oceangenomesmitogenomes.nf.
//
// Every finished assembly must produce its SQL upload row as soon as THAT assembly is done.
// This used to be a bare groupTuple(by: mt_assembly_prefix) reconciling the canonical
// post-curation row against the GetOrganelle first-pass variant, which emits nothing until
// its source channel closes -- so one queued assembly task left PUSH_MTDNA_ASSM_RESULTS with
// zero tasks, which emptied qc_ready (it gates on the upload receipt) and stopped every
// finished sample from reaching QC. Run mitogenomes-missing-audit-5 lost 141 samples that way.
//
// Also covers the rows that are never measured: a failed (empty) assembly, an under-length
// one, and a GetOrganelle provenance variant must all still be uploaded, with the placeholder
// depth, and must not be held to the end of the run either.

nextflow.enable.dsl = 2

// Import the real implementation rather than copying it, so this test cannot silently pass
// against a stale duplicate of the logic it is meant to protect.
include { buildAssemblyUploadRows } from '../../workflows/oceangenomesmitogenomes.nf'

// Stands in for an assembler, with a per-sample delay so the test can prove a fast sample is
// not held behind a slow assembly. bases controls whether the result is an empty (failed)
// assembly, an under-length one, or a full-length one.
process ASSEMBLE {
    input:
    tuple val(meta), val(delay_s), val(bases)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.fasta"), path("${meta.mt_assembly_prefix}.log")

    script:
    """
    sleep ${delay_s}
    if [ ${bases} -gt 0 ]; then
        printf '>%s\\n' "${meta.mt_assembly_prefix}"          > ${meta.mt_assembly_prefix}.fasta
        head -c ${bases} /dev/zero | tr '\\0' 'A'            >> ${meta.mt_assembly_prefix}.fasta
        printf '\\n'                                         >> ${meta.mt_assembly_prefix}.fasta
    else
        : > ${meta.mt_assembly_prefix}.fasta
    fi
    printf 'assembly log\\n' > ${meta.mt_assembly_prefix}.log
    """
}

// Stands in for MITOGENOME_COVERAGE: only ever runs on the molecule that reached annotation.
process MEASURE_DEPTH {
    input:
    tuple val(prefix), path(fasta), path(log)

    output:
    tuple val(prefix), path("${prefix}.mito_depth.tsv")

    script:
    """
    printf 'mean_depth\\n42\\n' > ${prefix}.mito_depth.tsv
    """
}

// Stands in for PUSH_MTDNA_ASSM_RESULTS: records when each row reached the SQL upload, and
// which depth file it arrived with (real measurement vs the placeholder).
process CONSUME_UPLOAD {
    input:
    tuple val(meta), path(fasta), path(log), path(depth)

    output:
    tuple val(meta), path("uploaded.${meta.mt_assembly_prefix}.txt")

    script:
    """
    printf '%s\\t%s\\t%s\\n' \\
        "${meta.mt_assembly_prefix}" "${depth.name}" "\$(date +%s%N)" \\
        > uploaded.${meta.mt_assembly_prefix}.txt
    """
}

workflow ASSEMBLY_UPLOAD_STREAMING {

    take:
    canonical_work // [ meta, delay_s, bases ] -- canonical assemblies (incl. failed / under-length)
    variant_work   // [ meta, delay_s, bases ] -- GetOrganelle provenance variants
    min_length     // under this many bases an assembly never reaches annotation, so never a depth

    main:
    ASSEMBLE(canonical_work)
    def canonical_rows = ASSEMBLE.out.map { meta, fasta, log -> [ meta.mt_assembly_prefix, meta, fasta, log ] }

    // Depth is measured only on assemblies that reach annotation, exactly as in the parent
    // workflow (MITOGENOME_COVERAGE sits on the post-length-filter channel).
    MEASURE_DEPTH(
        canonical_rows
            .filter { _prefix, _meta, fasta, _log ->
                fasta.size() > 0 && fasta.text.readLines().findAll { !it.startsWith('>') }.sum { it.trim().size() } >= min_length
            }
            .map { prefix, _meta, fasta, log -> [ prefix, fasta, log ] }
    )

    ASSEMBLE_VARIANTS(variant_work)
    def variant_rows = ASSEMBLE_VARIANTS.out.map { meta, fasta, log -> [ meta.mt_assembly_prefix, meta, fasta, log ] }

    def no_depth_file = file("${projectDir}/assets/empty_mito_depth.tsv", checkIfExists: true)

    rows = buildAssemblyUploadRows(
        canonical_rows,
        variant_rows,
        MEASURE_DEPTH.out,
        no_depth_file,
        min_length
    )

    CONSUME_UPLOAD(rows)

    emit:
    upload_rows = rows
    uploaded    = CONSUME_UPLOAD.out
}

// Separate instance so a slow canonical assembly and a slow variant can be timed independently.
process ASSEMBLE_VARIANTS {
    input:
    tuple val(meta), val(delay_s), val(bases)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.fasta"), path("${meta.mt_assembly_prefix}.log")

    script:
    """
    sleep ${delay_s}
    printf '>%s\\n' "${meta.mt_assembly_prefix}"  > ${meta.mt_assembly_prefix}.fasta
    head -c ${bases} /dev/zero | tr '\\0' 'A'    >> ${meta.mt_assembly_prefix}.fasta
    printf '\\n'                                 >> ${meta.mt_assembly_prefix}.fasta
    printf 'variant log\\n' > ${meta.mt_assembly_prefix}.log
    """
}
