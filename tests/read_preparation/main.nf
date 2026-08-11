nextflow.enable.dsl = 2

include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'

//
// Stand-in for the first per-sample consumer of prepared reads (the reference join
// feeding MITOHIFI_MITOHIFI). It records the epoch-second window in which it ran, so a
// test can prove a singleton sample was assembled while a multi-file sample was still
// concatenating rather than waiting behind it.
//
process DOWNSTREAM_MARKER {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.marker.tsv"), emit: marker

    script:
    def readList = reads instanceof List ? reads : [reads]
    def names = readList.collect { it.name }.join(',')
    """
    started=\$(date +%s)
    sleep 1
    printf '%s\\t%s\\t%s\\t%s\\n' '${meta.id}' '${names}' "\${started}" "\$(date +%s)" \\
        > ${meta.id}.marker.tsv
    """
}

//
// The read-preparation routing under test: one canonical prefixed channel, two
// mutually exclusive branches, CAT_FASTQ on the multi-file branch only, and `mix` to
// recombine without a channel-close barrier.
//
// The three statements below (map, branch, mix) are a verbatim copy of the
// corresponding block in subworkflows/local/mitogenome_assembly/mitohifi/main.nf --
// only the hard-coded MitoHiFi version tag differs. MITOGENOME_ASSEMBLY_MITOHIFI
// itself cannot be exercised here: it needs live NCBI lookups and the MitoHiFi
// container. Keep the two in step when either is edited.
//
workflow READ_PREPARATION {
    take:
    fastp_reads

    main:
    ch_reads_prefixed = fastp_reads
        .map { meta, reads ->
            def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v323mitohifi"
            [ meta + [ mt_assembly_prefix: mt_assembly_prefix ], reads ]
        }

    ch_reads_routed = ch_reads_prefixed.branch { meta, reads ->
        def readList = reads instanceof List ? reads : [ reads ]
        needs_concat: meta.single_end ? readList.size() > 1 : readList.size() > 2
        passthrough: true
    }

    CAT_FASTQ (
        ch_reads_routed.needs_concat
    )

    final_reads = ch_reads_routed.passthrough.mix(CAT_FASTQ.out.reads)

    DOWNSTREAM_MARKER (
        final_reads
    )

    emit:
    reads  = final_reads
    marker = DOWNSTREAM_MARKER.out.marker
}

//
// Stand-in for MITOHIFI_FINDMITOREFERENCE: emits one reference file per sample and hands
// the meta straight back, exactly as the real module does.
//
process REFERENCE_STUB {
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.reference.txt"), emit: reference

    script:
    """
    printf 'reference for %s\\n' '${meta.id}' > ${meta.id}.reference.txt
    """
}

//
// The reference join under test: keyed on mt_assembly_prefix, not on the whole meta map.
// A copy of the corresponding block in
// subworkflows/local/mitogenome_assembly/mitohifi/main.nf -- keep the two in step.
//
// `params.drop_meta_key` models what a RESUMED run does. On resume the process side hands
// back the meta stored in the cache database, which can carry a different set of keys from
// the live one: nf-schema materialises a samplesheet column the sheet does not carry as an
// empty list, and Nextflow's task hash cannot see an empty collection, so the task still
// resumes and returns the older meta. A join keyed on the whole meta map then matches
// nothing at all, with no error and no warning. A prefix-keyed join does not care.
//
workflow REFERENCE_JOIN {
    take:
    fastp_reads

    main:
    ch_reads_prefixed = fastp_reads
        .map { meta, reads ->
            def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v323mitohifi"
            [ meta + [ mt_assembly_prefix: mt_assembly_prefix ], reads ]
        }

    REFERENCE_STUB (
        ch_reads_prefixed.map { meta, _reads -> meta }
    )

    ch_reference_outcomes = REFERENCE_STUB.out.reference
        .map { meta, reference ->
            def returned = params.drop_meta_key
                ? meta.findAll { key, _value -> key != params.drop_meta_key }
                : meta
            [ returned.mt_assembly_prefix.toString(), reference ]
        }

    ch_reference_joined = ch_reads_prefixed
        .map { meta, reads -> [ meta.mt_assembly_prefix.toString(), meta, reads ] }
        .join(ch_reference_outcomes, by: 0)
        .map { _prefix, meta, reads, reference -> [ meta, reads, reference ] }

    DOWNSTREAM_MARKER (
        ch_reference_joined.map { meta, reads, _reference -> [ meta, reads ] }
    )

    emit:
    joined = ch_reference_joined
    marker = DOWNSTREAM_MARKER.out.marker
}

workflow {
    READ_PREPARATION(
        Channel.fromList(
            (params.samples ?: []).collect { sample ->
                [
                    [
                        id: sample.id,
                        sequencing_type: 'hifi',
                        date: '250101',
                        single_end: true,
                        assembly_prefix: "${sample.id}.hifi.250101",
                    ],
                    sample.reads.collect { file(it, checkIfExists: true) },
                ]
            }
        )
    )
    READ_PREPARATION.out.marker.view { meta, marker -> "MARKER\t${meta.id}\t${marker.text.trim()}" }
}
