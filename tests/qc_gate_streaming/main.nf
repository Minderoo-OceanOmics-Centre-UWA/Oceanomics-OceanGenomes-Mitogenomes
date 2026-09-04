// Stub harness for the circularity-evidence attachment at the QC gate in
// subworkflows/local/upload_results_mito.
//
// The gate must release each sample as soon as its OWN species validation, annotation stats
// and circularity evidence are in, rather than waiting for every assembly in the run. It must
// also cope with an assembly that never reached a circularity check (failed / no-contig /
// precomputed), and with the collapse path rewriting meta between the two sides of the join.
// All three properties are exercised here without the SQL processes the real subworkflow runs.

nextflow.enable.dsl = 2

// Import the real implementation rather than copying it, so this test cannot silently pass
// against a stale duplicate of the logic it is meant to protect.
include { attachCircularityEvidence } from '../../subworkflows/local/upload_results_mito/main.nf'

// Stands in for the assembly stage's circularity check, with a per-sample delay so the test
// can prove a fast sample is not held behind a slow assembly.
process ASSEMBLY_CHECK {
    input:
    tuple val(meta), val(delay_s)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.circularity_check.tsv")

    script:
    """
    sleep ${delay_s}
    printf 'anomaly_type\\tfinal_verdict_circular\\nnone\\ttrue\\n' \\
        > ${meta.mt_assembly_prefix}.circularity_check.tsv
    """
}

// Stands in for SPECIES_VALIDATION joined with PUSH_MTDNA_ANNOTATION_RESULTS: the per-sample
// side of the gate, which is ready long before a slow assembly elsewhere in the run.
process SAMPLE_READY {
    input:
    tuple val(meta), val(delay_s)

    output:
    tuple val(meta), path("blast.${meta.mt_assembly_prefix}.tsv"), path("stats.${meta.mt_assembly_prefix}.csv")

    script:
    """
    sleep ${delay_s}
    printf 'Found_in_blast_YN\\nYes\\n' > blast.${meta.mt_assembly_prefix}.tsv
    printf 'passed\\nyes\\n'           > stats.${meta.mt_assembly_prefix}.csv
    """
}

// Stands in for EVALUATE_QC_CONDITIONS: records when the sample reached the gate, and which
// evidence file it arrived with (real check vs the empty-check stand-in).
process CONSUME_GATE {
    input:
    tuple val(meta), path(blast), path(stats), path(evidence)

    output:
    tuple val(meta), path("gated.${meta.mt_assembly_prefix}.txt")

    script:
    """
    printf '%s\\t%s\\t%s\\n' \\
        "${meta.mt_assembly_prefix}" "${evidence.name}" "\$(date +%s%N)" \\
        > gated.${meta.mt_assembly_prefix}.txt
    """
}

workflow QC_GATE_STREAMING {

    take:
    sample_work   // [ meta, delay_s ] -- species validation + annotation stats for this sample
    assembly_work // [ meta, delay_s ] -- assemblies that DO reach a circularity check
    no_check      // [ meta ]          -- assemblies that never reach one

    main:
    ASSEMBLY_CHECK(assembly_work)
    SAMPLE_READY(sample_work)

    def no_circularity_evidence = file("${projectDir}/assets/placeholders/empty_circularity_check.tsv", checkIfExists: true)

    // Mirrors the totality contract the assembly subworkflows now guarantee: exactly one
    // evidence row per emitted assembly, the empty-check stand-in where no check ran, keyed
    // by mt_assembly_prefix.
    evidence = ASSEMBLY_CHECK.out
        .mix(no_check.map { meta -> [ meta, no_circularity_evidence ] })
        .map { meta, ev -> [ meta.mt_assembly_prefix, ev ] }

    gated = attachCircularityEvidence(SAMPLE_READY.out, evidence)

    CONSUME_GATE(gated)

    emit:
    conditions = gated
    released   = CONSUME_GATE.out
}
