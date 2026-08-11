#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALLOCATE_ENA_LOCUS_TAGS } from '../modules/local/genome_qc/allocate_ena_locus_tags/main'
include { PREPARE_ENA_METADATA } from '../modules/local/genome_qc/prepare_ena_metadata/main'
include { BUILD_ENA_CANDIDATE_PACKAGE } from '../modules/local/genome_qc/build_ena_candidate_package/main'
include { WEBIN_VALIDATE_GENOME } from '../modules/local/genome_qc/webin_validate_genome/main'
include { RECORD_ENA_PACKAGE_VALIDATION } from '../modules/local/genome_qc/record_ena_package_validation/main'

workflow {
    meta = [
        id: 'OG910',
        mt_assembly_prefix: 'OG910.hifi.250101.v3mitohifi',
        full_seqid: 'OG910.hifi.250101.v3mitohifi.emma102',
        annotation_version: 'emma102',
        scientific_name: 'Test species',
        sequencing_type: 'hifi',
        date: '250101',
        // Resolved from the technology by EnaTargets in the real subworkflow.
        ena_study: 'PRJEB123419',
        ena_locus_prefix: 'OGMTHIFI'
    ]
    fa = file("${projectDir}/test_data/ena_stub.fa", checkIfExists: true)
    tbl = file("${projectDir}/test_data/ena_stub.tbl", checkIfExists: true)
    gff = file("${projectDir}/test_data/ena_stub.gff", checkIfExists: true)
    embl = file("${projectDir}/test_data/ena_stub.embl.gz", checkIfExists: true)
    config = file("${projectDir}/test_data/ena_stub_db.cfg", checkIfExists: true)

    ALLOCATE_ENA_LOCUS_TAGS(channel.of(tuple(meta, fa, tbl)), config)
    PREPARE_ENA_METADATA(channel.of(meta), config)

    package_input = channel.of(tuple(meta, fa, gff, embl))
        .join(ALLOCATE_ENA_LOCUS_TAGS.out.mapping, by: 0)
        .join(ALLOCATE_ENA_LOCUS_TAGS.out.tagged_tbl, by: 0)
        .join(PREPARE_ENA_METADATA.out.metadata, by: 0)
    BUILD_ENA_CANDIDATE_PACKAGE(package_input)
    WEBIN_VALIDATE_GENOME(BUILD_ENA_CANDIDATE_PACKAGE.out.package_dir, 'test', 'stub')
    RECORD_ENA_PACKAGE_VALIDATION(WEBIN_VALIDATE_GENOME.out.status, config)
}
