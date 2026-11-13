nextflow.enable.dsl = 2

// BLAST_BLASTN emits `summary` as a PLAIN path and `tool_params` as a
// [meta, path] tuple, on adjacent lines of the same module:
//
//     path "filtered_summary.${gene_type}.txt"  , emit: summary
//     tuple val(meta), path("...mqcrow.html")   , emit: tool_params
//
// The MultiQC collection copied the tuple idiom onto the plain one:
//
//     ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect { it[1] })
//
// Indexing a Path does not fail -- Path implements the list-ish access Groovy
// wants, and returns its NAME ELEMENT at that index. For a summary under
// /scratch/pawsey1348/... element 1 is the literal string "pawsey1348", so
// Nextflow tried to hash and stage a directory component as a MultiQC input,
// 48 times in one run:
//
//     WARN: [HashBuilder] Unable to get file attributes file: /<outdir>/pawsey1348
//
// Every real filtered_summary.*.txt was silently dropped in the process. This
// workflow runs both idioms over the same channel so the fix is pinned by what
// it produces, not by the shape of the source line.
workflow MULTIQC_INPUT_SHAPES {
    take:
    summaries          // plain Path channel, as BLAST_BLASTN.out.summary is

    main:
    fixed = Channel.empty().mix(summaries)
    // The defect, kept executable so the test proves what it prevents. Emitted as a
    // STRING because emitting it as a path reproduces the bug too faithfully to
    // assert on: nf-test stages workflow outputs, so it dies resolving the bogus
    // component with the same NoSuchFileException Nextflow raised in the real run.
    buggy = Channel.empty().mix(summaries.map { it[1].toString() })

    emit:
    fixed
    buggy
}
