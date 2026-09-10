// Harness for the DB-free QC chain.
//
// Every step that answers "is this mitogenome any good?" used to sit inside
// `if (!params.skip_upload_results && params.sql_config)` in
// workflows/oceangenomesmitogenomes.nf, because two of them -- SPECIES_VALIDATION and
// the annotation statistics half of PUSH_MTDNA_ANNOTATION_RESULTS -- were each both a
// QC step and a database writer. A run without a database therefore produced no gene
// counts, no missing-gene set, no completeness verdict, no held-samples report, and a
// mitogenome_assembly_summary_mqc.tsv whose num_genes / num_cds / missing_genes /
// frameshift_flag columns were empty for every sample. The five-sample invert
// regression panel ran green and reported nothing.
//
// This pins the property that fixes it: MITOGENOME_QC reaches a verdict and emits the
// annotation statistics with NO sql_config anywhere. If a database dependency creeps
// back into any of these four processes, this test is what says so -- the symptom
// otherwise is silence, which is exactly how the original defect survived.

nextflow.enable.dsl = 2

// MITOGENOME_QC resolves its zero-region placeholders from ${projectDir}/assets, and
// nextflow sets projectDir to the directory of the script it launches -- this one.
// assets/ is a symlink back to the repo's own, so the subworkflow loads exactly the
// placeholders it would in a real run rather than test-only copies of them.

// The real subworkflow, not a copy, so this cannot pass against stale logic.
include { MITOGENOME_QC } from '../../subworkflows/local/mitogenome_qc/main.nf'

workflow {
    def meta = [
        id: 'OG001',
        mt_assembly_prefix: 'OG001.ilmn.240101.getorg1770',
        nominal_species_id: 'Genus species',
        class: 'Actinopteri',
        genetic_code: 2,
        circular: true,
    ]

    // One annotated region, so region_counts is 1 and the grouping closes immediately.
    ch_annotation = Channel.of([ meta, [ file("${projectDir}/../test_data/ena_stub.gff") ] ])
    ch_blast      = Channel.of([ meta, file("${projectDir}/../test_data/local_qc_no_db/blast.12s.tsv") ])
    ch_lca        = Channel.of([ meta, file("${projectDir}/../test_data/local_qc_no_db/lca.12s.tsv") ])
    ch_evidence   = Channel.of([ meta.mt_assembly_prefix,
                                 file("${projectDir}/assets/placeholders/empty_circularity_check.tsv") ])
    ch_regions    = Channel.of([ meta, 1 ])

    MITOGENOME_QC(ch_annotation, ch_blast, ch_lca, ch_evidence, ch_regions)

    // Prove the channel that feeds MITOGENOME_ASSEMBLY_SUMMARY is populated. Its being
    // empty is precisely what left the annotation columns blank.
    MITOGENOME_QC.out.assembly_summary_files.view { "ASSEMBLY_SUMMARY_FILE: ${it.name}" }
    MITOGENOME_QC.out.qc_ready.view { m, s, p, c -> "QC_READY: ${m.id} proceed=${p}" }
}
