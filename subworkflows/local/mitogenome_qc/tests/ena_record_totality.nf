// Test-only harness that mirrors the ENA validation-record grouping in
// MITOGENOME_QC (subworkflows/local/mitogenome_qc/main.nf) to guard its
// fixed-totality streaming invariant:
//
//   Every sample -- fully passing (E), table2asn-passed-but-conversion-failed (F),
//   or quarantined (Q) -- must contribute EXACTLY `ena_record_slots` files, so a
//   constant groupKey releases each record the moment its own files arrive, with NO
//   bare groupTuple and NO close-time remainder.
//
// The real subworkflow cannot be driven through the F/Q branches under `-stub`
// (every process stub emits status=PASS), so this harness reproduces the exact
// contributor mix + placeholder totality + groupKey construction with synthetic
// per-category channels. If the subworkflow's contributor set changes, update this
// harness and its slot count together.
//
// The stand-in "real output" file is arbitrary: the invariant is the per-sample
// ROW count feeding groupKey, not file content, so one shared file is reused for
// every real contributor row.

workflow ENA_RECORD_TOTALITY {

    main:

    def stand_in = file("${projectDir}/assets/placeholders/NO_REFERENCE.gb", checkIfExists: true)

    def notrun_flatfile_status  = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_status.not_run",  checkIfExists: true)
    def notrun_flatfile_checks  = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_checks.not_run",  checkIfExists: true)
    def notrun_flatfile_embl    = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_embl.not_run",    checkIfExists: true)
    def notrun_package_metadata = file("${projectDir}/assets/placeholders/ena_not_run/package_metadata.not_run", checkIfExists: true)
    def notrun_webin_status     = file("${projectDir}/assets/placeholders/ena_not_run/webin_status.not_run",     checkIfExists: true)
    def notrun_webin_manifest   = file("${projectDir}/assets/placeholders/ena_not_run/webin_manifest.not_run",   checkIfExists: true)

    // E = pass + conversion PASS (embl/package/webin run); F = pass + conversion FAIL;
    // Q = quarantined (never reaches ENA_FLATFILE).
    def mE = [ id: 'OGE', mt_assembly_prefix: 'E', full_seqid: 'E.emma102' ]
    def mF = [ id: 'OGF', mt_assembly_prefix: 'F', full_seqid: 'F.emma102' ]
    def mQ = [ id: 'OGQ', mt_assembly_prefix: 'Q', full_seqid: 'Q.emma102' ]

    parse_status = Channel.of( [mE, stand_in], [mF, stand_in], [mQ, stand_in] )  // always, all samples
    ff_status    = Channel.of( [mE, stand_in], [mF, stand_in] )                  // table2asn PASS set
    ff_checks    = Channel.of( [mE, stand_in], [mF, stand_in] )                  // table2asn PASS set
    ff_embl      = Channel.of( [mE, stand_in] )                                  // conversion PASS set
    pkg_metadata = Channel.of( [mE, stand_in] )                                  // conversion PASS set
    webin_status   = Channel.of( [mE, stand_in] )
    webin_manifest = Channel.of( [mE, stand_in] )

    ch_ena_quarantined = Channel.of( mQ )         // no ENA_FLATFILE / package / Webin
    ch_ena_no_embl     = Channel.of( mQ, mF )     // quarantined + conversion-failed

    ch_ena_validation_files = parse_status
        .mix( ff_status )
        .mix( ch_ena_quarantined.map { m -> [ m, notrun_flatfile_status ] } )
        .mix( ff_checks )
        .mix( ch_ena_quarantined.map { m -> [ m, notrun_flatfile_checks ] } )
        .mix( ff_embl )
        .mix( ch_ena_no_embl.map { m -> [ m, notrun_flatfile_embl ] } )
        .mix( pkg_metadata )
        .mix( ch_ena_no_embl.map { m -> [ m, notrun_package_metadata ] } )
    def ena_record_slots = (params.ena_webin_validate as boolean) ? 7 : 5
    if (params.ena_webin_validate) {
        ch_ena_validation_files = ch_ena_validation_files
            .mix( webin_status )
            .mix( ch_ena_no_embl.map { m -> [ m, notrun_webin_status ] } )
            .mix( webin_manifest )
            .mix( ch_ena_no_embl.map { m -> [ m, notrun_webin_manifest ] } )
    }

    group_sizes = ch_ena_validation_files
        .map { m, f -> tuple( groupKey(m.full_seqid ?: m.mt_assembly_prefix, ena_record_slots), m, f ) }
        .groupTuple(by: 0)
        .map { _key, metas, files -> [ metas[0].full_seqid, files.flatten().size() ] }

    emit:
    group_sizes = group_sizes
}
