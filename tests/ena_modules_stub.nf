#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ENA_FLATFILE } from '../modules/local/genome_qc/ena_flatfile/main'
include { WEBIN_VALIDATE } from '../modules/local/genome_qc/webin_validate/main'
include { GEN_FILES_TABLE2ASN } from '../modules/local/genome_qc/gen_files_table2asn/main'

workflow {
    meta = [id: 'ena_stub', mt_assembly_prefix: 'ena_stub']
    gbf = file("${projectDir}/test_data/ena_stub.gbf", checkIfExists: true)
    fa = file("${projectDir}/test_data/ena_stub.fa", checkIfExists: true)
    tbl = file("${projectDir}/test_data/ena_stub.tbl", checkIfExists: true)
    cmt = file("${projectDir}/test_data/ena_stub.cmt", checkIfExists: true)
    src = file("${projectDir}/test_data/ena_stub.src", checkIfExists: true)
    sbt = file("${projectDir}/test_data/ena_stub.sbt", checkIfExists: true)
    table2asn_input = channel.of(tuple(meta, fa, tbl, cmt, src, true))
    GEN_FILES_TABLE2ASN(table2asn_input, sbt)
    ENA_FLATFILE(channel.of(tuple(meta, gbf)))
    WEBIN_VALIDATE(ENA_FLATFILE.out.embl_file, 'PRJEB000000', 'stub-attempt')
}
