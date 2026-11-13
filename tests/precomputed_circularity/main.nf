nextflow.enable.dsl = 2

// Import the real implementations rather than copying them, so this cannot pass against a
// stale duplicate of the logic it protects.
include { assemblyRunPrefixOf          } from '../../workflows/oceangenomesmitogenomes.nf'
include { attachPrecomputedCircularity } from '../../workflows/oceangenomesmitogenomes.nf'
include { selectPrecomputedAssemblies  } from '../../workflows/oceangenomesmitogenomes.nf'
include { precomputedGlobCount        } from '../../workflows/oceangenomesmitogenomes.nf'

// One workflow, three emissions. Deliberately not three nextflow_workflow blocks: only the
// last block in a .nf.test file is executed, so the earlier ones would silently never run.
workflow PRECOMPUTED_SELECTION {
    take:
    prefixes
    fasta_ch
    evidence_ch
    placeholder
    glob_patterns

    main:
    // The run a variant belongs to. On the live path this is written once and is identical
    // across variants; the precomputed path used to set it per variant, so nothing grouped
    // them and every variant was annotated as its own assembly.
    resolved = prefixes.map { p -> [ p, assemblyRunPrefixOf(p) ] }

    // meta.circular must come back from the check TSV the original run left on disk.
    stamped = attachPrecomputedCircularity(fasta_ch, evidence_ch, placeholder)
        .map { meta, _fasta, _ev, real ->
            [ meta.mt_assembly_prefix, meta.containsKey('circular') ? meta.circular.toString() : 'absent', real ]
        }

    // ...and the variants of one run must reduce to a single assembly.
    winners = selectPrecomputedAssemblies('test', fasta_ch, evidence_ch, placeholder, 1.5d)
        .map { meta, _fasta, _ev -> [ meta.mt_assembly_run_prefix, meta.mt_assembly_prefix ] }

    // The counts the launch guard decides on. A precomputed run whose FASTA glob matches but
    // whose evidence glob matches NOTHING is the exact shape of the incident this guards: every
    // assembly silently took the NA placeholder, so circular genomes went to table2asn as
    // linear and the QC gate stopped holding on not_circular.
    counts = glob_patterns.map { label, pattern -> [ label, precomputedGlobCount(label, pattern) ] }

    emit:
    resolved
    stamped
    winners
    counts
}
