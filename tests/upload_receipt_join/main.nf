// Stub harness for the assembly-upload-receipt gate in
// subworkflows/local/upload_results_mito.
//
// A QC-ready sample must be released once its OWN assembly upload row is committed. The two
// sides reach this point down different lineages -- the QC side passed through
// MITOGENOME_ANNOTATION, which rewrites meta.mt_assembly_prefix to the annotated FASTA's
// basename, while the upload side carries the assembly stage's meta -- and on a -resume each
// is restored from its own cache entry. While this was a whole-meta join it matched nothing
// and dropped every sample silently: run mitogenomes-missing-audit-5 put 123 samples through
// the gate, had an upload row for all 123, and still ran MITOGENOME_QC zero times.

nextflow.enable.dsl = 2

// Import the real implementation rather than copying it, so this test cannot silently pass
// against a stale duplicate of the logic it is meant to protect.
include { attachAssemblyUploadReceipt } from '../../subworkflows/local/upload_results_mito/main.nf'

// Stands in for PUSH_MTDNA_ASSM_RESULTS, with a per-sample delay so the test can prove a fast
// sample is not held behind a slow one.
process COMMIT_UPLOAD {
    input:
    tuple val(meta), val(delay_s)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.mtdna.upload.txt")

    script:
    """
    sleep ${delay_s}
    printf 'uploaded\\n' > ${meta.mt_assembly_prefix}.mtdna.upload.txt
    """
}

// Stands in for MITOGENOME_QC: records which samples were released, and when.
process CONSUME_QC {
    input:
    tuple val(meta), val(species_name), val(proceed_qc), val(circular)

    output:
    tuple val(meta), path("qc.${meta.mt_assembly_prefix}.txt")

    script:
    """
    printf '%s\\t%s\\t%s\\n' \\
        "${meta.mt_assembly_prefix}" "${species_name}" "\$(date +%s%N)" \\
        > qc.${meta.mt_assembly_prefix}.txt
    """
}

workflow UPLOAD_RECEIPT_JOIN {

    take:
    qc_rows       // [ meta, species_name, proceed_qc, circular ] -- the QC-gate lineage's meta
    upload_work   // [ meta, delay_s ]                           -- the assembly lineage's meta

    main:
    COMMIT_UPLOAD(upload_work)

    released = attachAssemblyUploadReceipt(qc_rows, COMMIT_UPLOAD.out)

    CONSUME_QC(released)

    emit:
    qc_ready = released
    consumed = CONSUME_QC.out
}
