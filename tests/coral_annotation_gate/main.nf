// Stub harness for the annotator routing in subworkflows/local/mitogenome_annotation_lca.
//
// Exists because of a whole-session abort that no test saw. MITOS2.out.gff_proteins gained a
// fourth element (annotation/cds, for the QC gate's PCG ORF check) while the filter feeding
// ANNOTATION_QC_GATE still destructured three, so Groovy could not spread the tuple and threw
// a MissingMethodException from inside the filter operator. An operator failure is not a task
// failure -- errorStrategy does not apply -- so the run aborted and killed every queued task
// the moment the first MITOS2 sample emitted. Both invert panels died that way on 2026-09-07
// with every task green.
//
// The subworkflow's own module tests could not see it: they build the gate's 4-element input
// by hand and never exercise the wiring. So this runs the REAL subworkflow with every module
// stubbed, which is enough to bind every channel and catch any tuple-arity drift between a
// module's emit and the closure that consumes it.

nextflow.enable.dsl = 2

include { MITOGENOME_ANNOTATION } from '../../subworkflows/local/mitogenome_annotation_lca/main.nf'

workflow CORAL_ANNOTATION_GATE {

    take:
    mito_assembly
    reference_gb

    main:

    MITOGENOME_ANNOTATION (
        mito_assembly,
        params.curated_blast_db,
        params.nt_blast_db,
        params.mitos_refdb,
        reference_gb
    )

    emit:
    annotation_results = MITOGENOME_ANNOTATION.out.annotation_results
    region_counts      = MITOGENOME_ANNOTATION.out.region_counts
}
