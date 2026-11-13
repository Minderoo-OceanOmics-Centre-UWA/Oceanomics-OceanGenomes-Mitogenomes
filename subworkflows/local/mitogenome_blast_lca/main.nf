/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_BLAST_DB      } from '../../../modules/local/download_blast_db'
include { DOWNLOAD_TAXONKIT_DB   } from '../../../modules/local/download_taxonkit_db'
include { BLAST_BLASTN           } from '../../../modules/nf-core/blast/blastn'
include { LCA                    } from '../../../modules/local/LCA'
include { PREPARE_LCA_DATABASES  } from '../../../modules/local/lca/prepare_databases'

// Extract the annotation name from a per-gene FASTA filename: strip the trailing
// .fa and return everything after the first dot (e.g. CO1.<prefix>.fa -> <prefix>).
// Defined at file scope so it resolves inside the .map closures below (a closure
// assigned to a local `def` is not in scope for nested operator closures).
def getAnnotationName(filename) {
    def name = filename.toString().replaceAll(/\.fa$/, '')
    def parts = name.split('\\.', 2)
    return parts.size() > 1 ? parts[1] : name
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN BLAST AND LCA ON ANNOTATED CO1 / 12S / 16S REGIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Split out of MITOGENOME_ANNOTATION so BLAST/LCA can also run on annotations this run
    did not produce: with --skip_mitogenome_annotation and only
    --precomputed_mitogenome_annotation_results, the top-level workflow extracts the
    regions from the supplied GFF + FASTA (REBUILD_ANNOTATION_BUNDLE) and calls this directly.
    MITOGENOME_ANNOTATION calls it on its own annotator output, so both paths share one
    BLAST/LCA implementation and the same database selection and gates.
*/

workflow MITOGENOME_BLAST_LCA {

    take:
    co1_sequences    // tuple val(meta), path(CO1.<annotation_name>.fa)
    s12_sequences    // tuple val(meta), path(RNR1.<annotation_name>.fa)
    s16_sequences    // tuple val(meta), path(RNR2.<annotation_name>.fa)
    curated_blast_db // params.curated_blast_db
    nt_blast_db      // params.nt_blast_db

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Download taxonomy database
    //

    // storeDir is the shared central taxdb cache, so a null path would silently
    // disable caching and re-download 62 MB per run. Mirrors the --taxonkit_db_dir
    // assert later in this subworkflow.
    if (!params.blast_db_dir) {
        error "--blast_db_dir is required for the shared BLAST taxdb cache"
    }

    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))

    // taxdb.btd/bti give blastn its sscinames/scomnames columns; taxonomy4blast.sqlite3
    // is what `-taxids` needs (see BLAST_BLASTN's invertebrate ext.args branch in
    // conf/modules.config). db_files globs "taxdb*", which does NOT match the sqlite3,
    // and staging taxdb.* into the task CWD makes BLAST resolve taxonomy from the CWD
    // instead of from the -db directory -- so leaving the sqlite3 out silently disabled
    // -taxids for all 184 core_nt searches of NOVA_260909_LA (empty table, exit 0).
    //
    // .collect() is load-bearing, not decoration: db_files is a value channel today
    // only by the all-value-inputs singleton rule, and a bare .mix() yields a QUEUE
    // channel of two emissions, which a per-sample process consumes ONCE -- starving
    // every blastn task but the first and hanging the run. .collect() also flattens
    // the nested emission, so path(taxdb) receives a flat 3-file list and the module
    // signature is unchanged. .first() pins it as a value channel explicitly.
    ch_blast_db = DOWNLOAD_BLAST_DB.out.db_files
        .mix(DOWNLOAD_BLAST_DB.out.taxdb_file)
        .collect()
        .first()

    //
    // MODULE: Download taxonkit database
    //

    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))

    //
    // Use mix() to process CO1, 12s and 16s sequences through blast
    //

    regions = co1_sequences
        .map { meta, file ->
            def annotation_name = getAnnotationName(file.name)
            [meta, file, 'CO1', annotation_name]
        }
        .mix(
            s12_sequences.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '12s', annotation_name]
            },
            s16_sequences.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '16s', annotation_name]
            }
        )

    //
    // Select the BLAST database per sample (curated for vertebrates, nt for invertebrates)
    //

    // --nt_blast_db is a BLAST database PREFIX, not a file: a multi-volume database
    // such as core_nt has core_nt.nal + core_nt.00.nin..., a single-volume one such as
    // mito has mito.nin, and NEITHER has a file at the bare prefix. So this cannot be a
    // schema "exists": true check the way curated_blast_db (a real fasta) is -- it has
    // to look for the index that actually exists.
    //
    // Worth the check: NOVA_260909_LA ran its entire invertebrate cohort against
    // .../db/mito (RefSeq mitochondrion, 18,058 seqs) while every comment and doc in
    // this repo said core_nt. A wrong-but-readable database produces empty BLAST output
    // and no error at all, so nothing downstream noticed.
    combined_sequences_with_db = regions
        .map { meta, file_, gene_type, annotation_name ->
            if (meta.invertebrates && !nt_blast_db) {
                error "Sample ${meta.id} is marked invertebrates=true, but --nt_blast_db was not provided"
            }
            if (meta.invertebrates) {
                def prefix = nt_blast_db.toString()
                if (!file("${prefix}.nal").exists() && !file("${prefix}.nin").exists()) {
                    error "--nt_blast_db '${prefix}' is not a BLAST nucleotide database: expected " +
                          "<prefix>.nal (multi-volume, e.g. core_nt) or <prefix>.nin (single-volume) " +
                          "beside it. Pass the database PREFIX, e.g. /path/to/db/core_nt"
                }
            }
            def blast_db = (meta.invertebrates ? nt_blast_db : curated_blast_db)
            [meta, file_, gene_type, annotation_name, blast_db]
        }

    //
    // MODULE: Using CO1,12s and 16s run BLAST and filter the results to provide matches for the calculation of the LCA
    //

    BLAST_BLASTN (
        combined_sequences_with_db, // tuple val(meta), path(fasta), val(gene_type), val(annotation_name), val(blast_db)
        ch_blast_db // path(db)
    )


    //
    // MODULE: Calculate the Lowest Common Ancestor (LCA) from the filtered BLAST results
    //
    ch_worms = channel.fromPath("${projectDir}/assets/taxonomy/worms_species.txt.gz", checkIfExists: true)
    if (!params.taxonkit_db_dir) {
        error "--taxonkit_db_dir is required for the persistent LCA taxonomy cache"
    }
    ch_lca_cache_dir = Channel.value("${params.taxonkit_db_dir}/lca_cache")
    ch_lca_prepare_script = Channel.value(file("${projectDir}/bin/prepare_lca_databases.py", checkIfExists: true))
    ch_lca_script = Channel.value(file("${projectDir}/bin/calculateLCA.py", checkIfExists: true))
    PREPARE_LCA_DATABASES(ch_lca_cache_dir, ch_lca_prepare_script)
    // PREPARE_LCA_DATABASES deliberately runs on every launch, but its audit
    // manifest contains a fresh verified_at timestamp and lives in a new work
    // directory. Convert it to the deterministic content signature before using
    // it as the LCA dependency, so unchanged databases preserve task hashes.
    ch_lca_cache_signature = PREPARE_LCA_DATABASES.out.manifest
        .map { manifest ->
            def payload = new groovy.json.JsonSlurper().parse(manifest.toFile())
            def signature = payload.cache_signature?.toString()
            if (!(signature ==~ /[0-9a-f]{64}/)) {
                throw new IllegalStateException("Invalid or missing cache_signature in ${manifest}")
            }
            signature
        }

    LCA (
        BLAST_BLASTN.out.filtered,
        ch_worms.first(),
        ch_lca_cache_dir,
        ch_lca_cache_signature,
        ch_lca_script
    )

    // summary is a plain Path channel, not [meta, path]. Path[1] returns the second
    // component (for example "pawsey1348"), creating bogus MultiQC inputs.
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary)
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(LCA.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())
    ch_versions = ch_versions.mix(PREPARE_LCA_DATABASES.out.versions)

    emit:
    blast_filtered_results  = BLAST_BLASTN.out.validation   // channel: [ meta, path(blast.<gene>.<name>.filtered.tsv) ]
    lca_results             = LCA.out.lca
    lca_raw_results         = LCA.out.lca_raw
    multiqc_files           = ch_multiqc_files              // channel: [ path(multiqc_files) ]
    versions                = ch_versions                   // channel: [ path(versions.yml) ]
}
