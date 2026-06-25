// Test-only harness for MITOGENOME_ASSEMBLY_GETORG.
//
// The subworkflow's `multiqc_files` / `summary_files` channels deliberately mix
// bare paths with .collect()-ed lists, which nf-test's channel sorter cannot
// order (ClassCastException String vs ArrayList) when it loads ALL named
// outputs. This harness exposes only the clean `assembly_fasta` channel so the
// stub tests can run; the HiC->fastp->GetOrganelle wiring is still verified via
// the execution trace (which captures every task regardless of emitted channels).

include { MITOGENOME_ASSEMBLY_GETORG } from '../main.nf'

workflow GETORG_TRIM_HARNESS {

    take:
    reads          // tuple(meta, reads)
    organelle_type // val

    main:
    MITOGENOME_ASSEMBLY_GETORG ( reads, organelle_type )

    emit:
    assembly_fasta = MITOGENOME_ASSEMBLY_GETORG.out.assembly_fasta
}
