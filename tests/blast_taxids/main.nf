// Harness for the two properties that let the -taxids breakage run for a whole
// cohort unnoticed (NOVA_260909_LA_mitogenomes-3):
//
//  1. CARDINALITY -- the taxdb channel must reach EVERY blastn task. It is built by
//     mixing two DOWNLOAD_BLAST_DB outputs, and a bare .mix() yields a queue channel
//     of two emissions that a per-sample process consumes once, starving every task
//     but the first. VERIFIED by deliberately dropping .collect().first() from the
//     wiring below: the workflow HANGS rather than failing, so in CI this guard
//     surfaces as a test timeout, not a tidy assertion failure. Either way it stops
//     the regression reaching a cohort run, which is the point.
//  2. PREFLIGHT -- when ext.args asks for -taxids and taxonomy4blast.sqlite3 is not
//     staged, the task must FAIL. blastn's own behaviour is to warn, write an empty
//     table and exit 0, which downstream reads as "no hits".
//
// Note this file lives under tests/ rather than modules/nf-core/blast/blastn/tests/:
// nf-test.config ignores 'modules/nf-core/**/tests/*', so a test placed there is
// silently skipped.

nextflow.enable.dsl = 2

include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn'

// Builds a tiny real BLAST db so the module's actual search path runs. Uses the same
// container as BLAST_BLASTN so makeblastdb and blastn are the same build.
process MAKE_TINY_DB {
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    path subject

    output:
    path "tinydb.*", emit: db_files

    script:
    """
    makeblastdb -in ${subject} -dbtype nucl -out tinydb
    """
}

// Exercises the module's REAL script path end to end: the set +e / exit-code capture,
// the fatal-stderr-signature check, the empty-output warning, and the awk filter chain.
// The preflight test below never reaches blastn, so without this the restructured body
// would only ever be executed for the first time on a live cohort run.
workflow BLAST_REAL_SEARCH {

    take:
    subject   // fasta to build the tiny db from
    query     // [ meta, fasta, gene_type, annotation_name ]
    taxdb     // taxdb files (no -taxids here, so the sqlite3 is not required)

    main:
    MAKE_TINY_DB(subject)

    // The tiny db lives in MAKE_TINY_DB's work dir, which singularity does NOT mount
    // into the blastn task (only the task's own work dir and its staged inputs are
    // visible) -- reaching for it by absolute path gives "Database memory map file
    // error". So stage the db files through the module's own path(taxdb) input, which
    // puts them in the task CWD, and pass the db as a bare relative prefix. The
    // module's `export BLASTDB="$PWD"` is what makes that resolve, so this also
    // exercises that line.
    db_bundle = taxdb
        .mix(MAKE_TINY_DB.out.db_files.flatten())
        .collect()
        .first()

    queries = query.map { meta, fa, gt, an -> [ meta, fa, gt, an, 'tinydb' ] }

    BLAST_BLASTN(queries, db_bundle)

    emit:
    filtered = BLAST_BLASTN.out.filtered
    summary  = BLAST_BLASTN.out.summary
}

workflow BLAST_TAXDB_FANOUT {

    take:
    queries  // [ meta, fasta, gene_type, annotation_name, blast_db ]
    taxdb_a  // the taxdb.* files (stands in for DOWNLOAD_BLAST_DB.out.db_files)
    taxdb_b  // taxonomy4blast.sqlite3 (stands in for .out.taxdb_file)

    main:
    // Built exactly as subworkflows/local/mitogenome_annotation_lca does, so this
    // test fails if that wiring regresses to a bare .mix().
    ch_blast_db = taxdb_a
        .mix(taxdb_b)
        .collect()
        .first()

    BLAST_BLASTN(queries, ch_blast_db)

    emit:
    filtered = BLAST_BLASTN.out.filtered
    summary  = BLAST_BLASTN.out.summary
}
