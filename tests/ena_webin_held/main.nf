nextflow.enable.dsl = 2

// Import the real implementations, so this cannot pass against a copy of the parsing it
// exists to protect. These reasons are read by COLUMN INDEX out of a status TSV, which is
// exactly the kind of thing that breaks silently when a column is inserted.
include { webinHeldReason    } from '../../subworkflows/local/ena_submission_prep/main.nf'
include { flatfileFailReason } from '../../subworkflows/local/ena_submission_prep/main.nf'

workflow HELD_REASONS {
    take:
    webin_statuses
    flatfile_statuses

    main:
    webin = webin_statuses.map { label, f -> [ label ] + webinHeldReason(f) }
    flatfile = flatfile_statuses.map { label, f -> [ label, flatfileFailReason(f) ] }

    emit:
    webin
    flatfile
}
