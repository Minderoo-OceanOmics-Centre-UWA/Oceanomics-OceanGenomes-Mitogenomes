/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mitogenome modules
include { DOWNLOAD_BLAST_DB      } from '../../../modules/local/download_blast_db'
include { DOWNLOAD_TAXONKIT_DB   } from '../../../modules/local/download_taxonkit_db'
include { EMMA                   } from '../../../modules/local/EMMA'
include { EMMA_GENE_RESCUE_GATE  } from '../../../modules/local/emma_gene_rescue_gate'
include { EMMA_GENE_RESCUE       } from '../../../modules/local/emma_gene_rescue'
include { TRNA_RESCUE_GATE       } from '../../../modules/local/trna_rescue_gate'
include { TRNA_SCAN              } from '../../../modules/local/trna_scan'
include { TRNA_RESCUE            } from '../../../modules/local/trna_rescue'
include { ROTATE_ORIGIN          } from '../../../modules/local/rotate_origin'
include { MITOS2                 } from '../../../modules/local/mitos2'
include { ANNOTATION_QC_GATE     } from '../../../modules/local/annotation_qc_gate'
include { CORAL_ANNOTATION_FIX   } from '../../../modules/local/coral_annotation_fix'
include { REFERENCE_RELEVANCE    } from '../../../modules/local/reference_relevance'
include { SELECT_CORAL_REFERENCE } from '../../../modules/local/select_coral_reference'
include { BLAST_BLASTN           } from '../../../modules/nf-core/blast/blastn'
include { LCA                    } from '../../../modules/local/LCA'
include { PREPARE_LCA_DATABASES  } from '../../../modules/local/lca/prepare_databases'

// Helper functions
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'

// Extract the annotation name from a per-gene FASTA filename: strip the trailing
// .fa and return everything after the first dot (e.g. CO1.<prefix>.fa -> <prefix>).
// Defined at file scope so it resolves inside the .map closures below (a closure
// assigned to a local `def` is not in scope for nested operator closures).
def getAnnotationName(filename) {
    def name = filename.toString().replaceAll(/\.fa$/, '')
    def parts = name.split('\\.', 2)
    return parts.size() > 1 ? parts[1] : name
}

// How many of CO1 / 12S / 16S an annotation bundle yielded, i.e. how many items
// this sample contributes to ch_annot_co1/s12/s16 and therefore how many BLAST
// and LCA tasks it will spawn (0-3).
//
// This count is what lets upload_results_mito size each sample's result group
// with groupKey() and release the sample as soon as its own regions are done,
// instead of waiting for every LCA in the run to finish. It is derived from the
// annotation bundle rather than emitted by the annotators so that EMMA / MITOS2 /
// CORAL_ANNOTATION_FIX task hashes are untouched and -resume still works.
//
// Counts matched globs, not files: the annotators declare
// `path("annotation/cds/*CO1*.fa")`, which emits a single channel item even in
// the (unseen so far) case of a glob matching more than one file.
def countAnnotatedRegions(files) {
    try {
        def bundle = (files instanceof List) ? files : [files]
        def cds = bundle.find { it.name == 'cds' && it.isDirectory() }
        if (!cds) return 0
        def names = []
        cds.eachFile { names << it.name }
        return ['CO1', 'RNR1', 'RNR2'].count { region ->
            names.any { it.endsWith('.fa') && it.contains(region) }
        }
    } catch (ignored) {
        return 0
    }
}

// Did SELECT_CORAL_REFERENCE pick a reference for this sample? Its status line starts
// with SELECTED / SELECTED_LOW_CONFIDENCE when a DB record aligned (and a reference.gb
// was written), or NONE when nothing did. Read the always-emitted status so the fixer
// can branch on a per-sample VALUE rather than a remainder join against the optional
// reference output. Defined at file scope so it resolves inside the .map/.filter closures.
def coralReferenceSelected(statusFile) {
    try {
        def line = statusFile.text.readLines().find { it?.trim() }
        return line?.split('\t', -1)?.first()?.trim()?.toUpperCase()?.startsWith('SELECTED') ?: false
    } catch (ignored) {
        return false
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ANNOTATION AND LCA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ANNOTATION {

    take:
    mito_assembly //  tuple val(meta), path(fasta)
    curated_blast_db // params.curated_blast_db
    nt_blast_db // params.nt_blast_db
    mitos_refdb // params.mitos_refdb (MITOS2 RefSeq reference data dir, for invertebrates)
    reference_gb // tuple val(meta), path(reference.gb) - per-sample, partial (assembly stage)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Download taxonomy database
    //

    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))
    
    ch_blast_db = DOWNLOAD_BLAST_DB.out.db_files

    // 
    // MODULE: Download taxonkit database
    //

    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))
    
    //
    // Assert the assembly's identity rather than assigning it.
    //
    // This used to OVERWRITE meta.mt_assembly_prefix with the FASTA basename, which is where
    // the pipeline's two notions of "the prefix" collided: the assembly stage kept a
    // sample-level value while everything from here on used the curated basename, so every
    // join spanning this line silently missed. That dropped 23 of 168 assemblies at the QC
    // gate and filed a reseed's depth against the first-pass row.
    //
    // Identity is now stamped by whoever renames the molecule (the reseed/rgj resolution, the
    // collapse, SANITISE_FASTA), so by here the two must already agree. Fail loudly if they do
    // not: a silent re-assignment is exactly what made the original defect invisible.
    //
    fasta_with_mt_assembly_prefix = mito_assembly
    .map { meta, fasta ->
        assert meta.mt_assembly_prefix == fasta.baseName : \
            "Assembly identity out of sync with its FASTA: meta.mt_assembly_prefix=" +
            "'${meta.mt_assembly_prefix}' but the annotated FASTA is '${fasta.baseName}'. " +
            "Whichever stage produced this FASTA must stamp mt_assembly_prefix from its basename."
        [meta, fasta]
    }

    //
    // MODULE: Reference relevance check.
    // The reference is resolved from the sample's species label, so a wrong/coarse
    // label yields a wrong-family reference that silently degrades seeding and the
    // coral annotation fix. BLAST the resolved reference against the assembly and
    // record a PASS/MISMATCH review flag for every sample that has one. Label- and
    // taxonomy-DB-free; always exits 0.
    //
    // Key the reference reuse on the LINEAGE prefix, not on identity and not on the whole meta
    // map: one reference is resolved per assembly run, before curation, so reference_gb's meta
    // carries the assembly-stage name while mito_assembly's carries the curated one
    // (<prefix>reseed / _collapsed / _concat). Joining on either identity or the whole map
    // would silently drop every curated assembly from the relevance check. The output keeps
    // the assembly's identity so it publishes into the same dir as the annotation.
    ch_reference_gb_keyed = reference_gb.map { meta, ref -> [ meta.mt_assembly_run_prefix, ref ] }

    REFERENCE_RELEVANCE (
        mito_assembly
            .map { meta, fasta -> [ meta.mt_assembly_run_prefix, meta, fasta ] }
            .join(ch_reference_gb_keyed, by: 0)
            .map { _key, meta, fasta, ref -> [ meta, fasta, ref ] }
    )
    ch_reference_relevance = REFERENCE_RELEVANCE.out.flag

    //
    // MODULE: Annotate the mitogenome.
    // Route by meta.invertebrates: vertebrates -> EMMA, invertebrates -> MITOS2.
    // EMMA does not annotate inverts correctly, so coral/invert samples use
    // MITOS2 instead. MITOS2 emits the same output channels as EMMA so the rest
    // of this subworkflow is annotator-agnostic.
    //

    // Hard assertion before the EMMA/MITOS2 branch: every sample must carry an
    // integer meta.genetic_code in the supported set. prepare_samplesheet resolves
    // it once (from taxonomic class or the explicit genetic_code column); this is
    // the fail-fast that stops an unset/garbage code reaching annotation, where
    // MITOS2 would translate under one table and table2asn validate under another.
    // Mirrors the --mitos_refdb / --nt_blast_db asserts later in this subworkflow.
    def SUPPORTED_GENETIC_CODES = [2, 4, 5, 9, 13, 14, 21, 24, 33] as Set
    ch_annot_branched = fasta_with_mt_assembly_prefix
        .map { meta, fasta ->
            def gc = meta.genetic_code
            if (!(gc instanceof Integer) || !SUPPORTED_GENETIC_CODES.contains(gc)) {
                error "Sample ${meta.id}: genetic_code '${gc}' is not a supported " +
                      "mitochondrial translation table (${SUPPORTED_GENETIC_CODES.sort().join(', ')})"
            }
            [meta, fasta]
        }
        .branch { meta, _fasta ->
            invert: meta.invertebrates
            vert:   true
        }

    // MITOS2 needs its RefSeq reference data dir; fail clearly if an invert
    // sample is present but --mitos_refdb was not provided.
    ch_mitos_refdb = mitos_refdb
        ? Channel.fromPath(mitos_refdb, checkIfExists: true).first()
        : Channel.value([])

    EMMA (
        ch_annot_branched.vert // tuple val(meta), path(fasta)
    )

    //
    // MODULE: EMMA gene rescue (ND4L / ATP8).
    // EMMA's rationalise_matches! overlap filter drops these short PCGs against
    // their longer neighbour (ND4, ATP6) even though the gene is in the assembly.
    // The gate flags an EMMA bundle as FIX (only ND4L/ATP8 missing, flanks present
    // and the rest in order) or PASS; only FIX bundles are re-annotated, and only
    // the annotation `results` bundle is swapped. co1/12S/16S -> BLAST/LCA and the
    // sample's onward flow are untouched: a rescue that cannot recover the gene
    // re-emits EMMA's original bundle and the sample is held at the QC gate as now.
    //
    // Every join below uses `remainder: true`. The rescue processes run with
    // errorStrategy 'ignore' so a crashed rescue cannot kill the run -- but an
    // ignored task emits nothing, and a plain join would then drop that sample
    // out of the pipeline entirely: no annotation bundle, no hold row, no error.
    // With the remainder joins, a missing gate decision reads as PASS and a
    // missing rescue output falls back to the input bundle, so the worst a failed
    // rescue can do is leave the annotation exactly as EMMA produced it.
    //
    EMMA_GENE_RESCUE_GATE ( EMMA.out.results )

    ch_emma_gate = EMMA_GENE_RESCUE_GATE.out.decision
        .map { meta, dfile ->
            def parts = dfile.text.trim().split('\t')
            [ meta, parts[0].trim(), parts.size() > 1 ? parts[1].trim() : '-' ]  // [meta, FIX|PASS, targets]
        }

    // join(remainder: true) pads an UNMATCHED row with a single trailing null, not
    // one null per right-hand element, so the row length varies -- index into it
    // rather than destructuring, or the closure fails on the very rows this is
    // here to keep.
    ch_emma_branched = EMMA.out.results.join(ch_emma_gate, remainder: true)
        .filter { row -> row[1] != null }                    // drop decision-only rows
        .map { row ->                                        // no decision -> PASS
            [ row[0], row[1],
              row.size() > 2 && row[2] ? row[2] : 'PASS',
              row.size() > 3 && row[3] ? row[3] : '-' ]
        }
        .branch { _meta, _bundle, decision, _targets ->
            fix:  decision == 'FIX'
            pass: true
        }
    ch_emma_pass = ch_emma_branched.pass.map { meta, bundle, _d, _t      -> [meta, bundle] }
    ch_emma_fix  = ch_emma_branched.fix.map  { meta, bundle, _d, targets -> [meta, bundle, targets] }

    ch_rescue_ref = Channel.fromPath("${projectDir}/assets/rescue_pcg_refs.faa", checkIfExists: true).first()

    EMMA_GENE_RESCUE (
        ch_emma_fix.combine(ch_rescue_ref).map { meta, bundle, targets, ref -> [meta, bundle, targets, ref] }
    )

    // EMMA bundle after the ND4L/ATP8 rescue: PASS bundles unchanged, FIX patched,
    // and a FIX whose rescue task failed falls back to its original bundle.
    ch_emma_results = ch_emma_pass.mix(
        ch_emma_fix.map { meta, bundle, _t -> [meta, bundle] }
            .join(EMMA_GENE_RESCUE.out.results, remainder: true)
            .filter { row -> row[1] != null }
            .map { row -> [ row[0], row.size() > 2 && row[2] ? row[2] : row[1] ] }
    )

    //
    // MODULE: tRNA rescue.
    // EMMA's covariance model sometimes misses a tRNA that is in the assembly on
    // an otherwise complete, correctly ordered mitogenome. The gate flags a
    // bundle as FIX (only tRNAs missing, whole PCG+rRNA core present and ordered)
    // or PASS; only FIX bundles are re-annotated, and only the annotation
    // `results` bundle is swapped. co1/12S/16S -> BLAST/LCA are untouched. A
    // rescue that cannot confidently place the tRNA re-emits the original bundle
    // and the sample is held at the QC gate as now (unless the residual shortfall
    // is within annotation_trna_tolerance).
    //
    // The rescue is vertebrate-code only (tRNAscan-SE's -M vert model, and the
    // canonical REF gene order the gate checks). That is expressed as an explicit
    // branch here rather than `ext.when` on the processes: a process skipped by
    // ext.when emits nothing, which -- exactly like an ignored failure -- would
    // silently drop every non-code-2 sample at the join below. Branching keeps
    // them visible and passes them straight through.
    //
    ch_trna_eligible = ch_emma_results.branch { meta, _bundle ->
        code2: meta.genetic_code == 2
        other: true
    }

    TRNA_RESCUE_GATE ( ch_trna_eligible.code2 )

    ch_trna_gate = TRNA_RESCUE_GATE.out.decision
        .map { meta, dfile ->
            def parts = dfile.text.trim().split('\t')
            [ meta, parts[0].trim(), parts.size() > 1 ? parts[1].trim() : '-' ]  // [meta, FIX|PASS, targets]
        }

    ch_trna_branched = ch_trna_eligible.code2.join(ch_trna_gate, remainder: true)
        .filter { row -> row[1] != null }
        .map { row ->
            [ row[0], row[1],
              row.size() > 2 && row[2] ? row[2] : 'PASS',
              row.size() > 3 && row[3] ? row[3] : '-' ]
        }
        .branch { _meta, _bundle, decision, _targets ->
            fix:  decision == 'FIX'
            pass: true
        }
    ch_trna_pass = ch_trna_branched.pass.map { meta, bundle, _d, _t      -> [meta, bundle] }
    ch_trna_fix  = ch_trna_branched.fix.map  { meta, bundle, _d, targets -> [meta, bundle, targets] }

    // tRNAscan-SE runs in its own (Perl) container; TRNA_RESCUE parses + splices
    // its output in the stdlib psycopg2 container. A sample whose scan failed has
    // no tsv, so it never reaches TRNA_RESCUE and keeps its bundle below.
    TRNA_SCAN ( ch_trna_fix.map { meta, bundle, _targets -> [meta, bundle] } )

    TRNA_RESCUE (
        ch_trna_fix.join(TRNA_SCAN.out.tsv)
            .map { meta, bundle, targets, scan -> [meta, bundle, scan, targets] }
    )

    // Vertebrate annotation bundle used from here on: non-code-2 and PASS bundles
    // unchanged, FIX patched, and any FIX whose scan or rescue failed falls back
    // to the bundle it went in with.
    ch_emma_results_trna = ch_trna_eligible.other
        .mix(ch_trna_pass)
        .mix(
            ch_trna_fix.map { meta, bundle, _t -> [meta, bundle] }
                .join(TRNA_RESCUE.out.results, remainder: true)
                .filter { row -> row[1] != null }
                .map { row -> [ row[0], row.size() > 2 && row[2] ? row[2] : row[1] ] }
        )

    // Re-origin invert assemblies to the cox1 start before MITOS2 so intron-split
    // genes (e.g. the hexacoral nad5 group I intron) no longer straddle
    // GetOrganelle's linearisation point. This mirrors EMMA's `--rotate MT-TF`
    // for vertebrates; several invert groups lack tRNA-Phe (Cnidaria always;
    // Porifera in some lineages), so cox1 is the anchor instead. Each sample is
    // rotated against its own phylum-appropriate curated cox1 panel -- see
    // InvertTaxonGroups.cox1PanelGroup() in lib/. The original (un-rotated)
    // assembly is published separately by the assembly stage and left untouched.
    ch_cox1_panel_paths = [
        reduced_trna  : file("${projectDir}/assets/cox1_anthozoa.faa", checkIfExists: true),
        mollusca      : file("${projectDir}/assets/cox1_mollusca.faa", checkIfExists: true),
        arthropoda    : file("${projectDir}/assets/cox1_arthropoda.faa", checkIfExists: true),
        echinodermata : file("${projectDir}/assets/cox1_echinodermata.faa", checkIfExists: true),
    ]

    ROTATE_ORIGIN (
        ch_annot_branched.invert.map { meta, fasta ->
            if (!mitos_refdb) {
                error "Sample ${meta.id} is marked invertebrates=true, but --mitos_refdb was not provided"
            }
            def group = InvertTaxonGroups.cox1PanelGroup(meta.class)
            [meta, fasta, ch_cox1_panel_paths[group]]
        }
    )

    MITOS2 (
        ROTATE_ORIGIN.out.fasta,
        ch_mitos_refdb
    )

    //
    // Anthozoan annotation QC gate + reference-based fixer.
    // MITOS2 annotates coral PCGs + 12S correctly but routinely drops the
    // divergent 16S and one exon of the intron-split nad5 (a Hexacorallia-wide
    // group I intron trait). The gate flags each Cnidarian annotation as FIX
    // (deficient) or PASS, so only broken corals are re-annotated; correctly
    // annotated batch-mates pass through MITOS2 untouched.
    //
    // This gate is Cnidaria-only: it encodes a failure mode discovered from real
    // coral MITOS2 output, and guessing an equivalent for other invert phyla
    // (including Porifera, despite sharing genetic code/rotation handling with
    // Cnidaria -- see InvertTaxonGroups in lib/) would be exactly that, a guess.
    // Other invert phyla's MITOS2 output merges straight through below
    // (ch_mitos_*_passthrough); a phylum-specific fix gets designed later from
    // real evidence once this batch has run.
    //
    ANNOTATION_QC_GATE (
        MITOS2.out.gff_proteins.filter { meta, _gff, _proteins -> InvertTaxonGroups.isCoralFixEligible(meta.class) }
    )

    ch_gate = ANNOTATION_QC_GATE.out.decision
        .map { meta, dfile -> [meta, dfile.text.split('\t')[0].trim()] }   // [meta, FIX|PASS]

    // Non-Cnidarian invert phyla never enter the gate above, so joining any
    // MITOS2 output channel against ch_gate below only ever matches Cnidaria
    // samples (join() is an inner join on meta) -- these are their untouched
    // MITOS2 outputs, merged directly into the final annotation channels.
    ch_mitos_co1_passthrough     = MITOS2.out.co1_sequences.filter { meta, _f -> !InvertTaxonGroups.isCoralFixEligible(meta.class) }
    ch_mitos_s12_passthrough     = MITOS2.out.s12_sequences.filter { meta, _f -> !InvertTaxonGroups.isCoralFixEligible(meta.class) }
    ch_mitos_s16_passthrough     = MITOS2.out.s16_sequences.filter { meta, _f -> !InvertTaxonGroups.isCoralFixEligible(meta.class) }
    ch_mitos_results_passthrough = MITOS2.out.results.filter { meta, _f -> !InvertTaxonGroups.isCoralFixEligible(meta.class) }

    // PASS corals: keep MITOS2 output unchanged.
    ch_mitos_co1_pass     = MITOS2.out.co1_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_s12_pass     = MITOS2.out.s12_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_s16_pass     = MITOS2.out.s16_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_results_pass = MITOS2.out.results.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }

    // FIX corals: fixer base = the cox1-rotated genome MITOS2 annotated + the raw BED.
    ch_fix_base = ROTATE_ORIGIN.out.fasta
        .join(MITOS2.out.bed)
        .join(ch_gate)
        .filter { it[-1] == 'FIX' }
        .map { meta, genome, bed, _d -> [meta, genome, bed] }

    // Reference resolution -- label-free. Pick the annotation reference from the
    // curated Anthozoa mitogenome DB by SEQUENCE similarity to the assembly, so a
    // wrong/coarse species label can no longer hand a wrong-family reference to
    // the fixer. The bundled curated anthozoan reference is the only fallback (used
    // just for the rare sample no DB record aligns to).
    ch_coral_db = Channel.fromPath("${projectDir}/assets/coral_mito_refdb.gb", checkIfExists: true).first()

    SELECT_CORAL_REFERENCE ( ch_fix_base.map { meta, genome, _bed -> [meta, genome] }, ch_coral_db )
    ch_selected_ref = SELECT_CORAL_REFERENCE.out.reference   // [meta, reference.gb]

    ch_fix_selected = ch_fix_base.join(ch_selected_ref)
        .map { meta, genome, bed, ref -> [meta, genome, bed, ref] }

    // Fallback: a FIX sample for which the selector emitted no reference (status NONE:
    // no DB record aligned) gets the bundled curated anthozoan reference. The unselected
    // set is read from SELECT_CORAL_REFERENCE.out.status, which is emitted for EVERY FIX
    // sample, so this is a plain per-sample join (the sample falls back the moment its own
    // status arrives) rather than a remainder join that waits for the whole run to close.
    ch_coral_ref_unselected = SELECT_CORAL_REFERENCE.out.status
        .filter { _meta, status_file -> !coralReferenceSelected(status_file) }
        .map { meta, _status_file -> [meta, true] }
    ch_curated_ref = Channel.fromPath("${projectDir}/assets/anthozoa_reference.gb", checkIfExists: true).first()
    ch_fix_fallback = ch_fix_base.join(ch_coral_ref_unselected, by: 0)
        .map { meta, genome, bed, _flag -> [meta, genome, bed] }
        .combine(ch_curated_ref)
        .map { meta, genome, bed, ref -> [meta, genome, bed, ref] }

    ch_coral_fix_input = ch_fix_selected.mix(ch_fix_fallback)

    CORAL_ANNOTATION_FIX ( ch_coral_fix_input )

    // Merge annotators: verts (EMMA, ND4L/ATP8- then tRNA-rescued) + PASS corals
    // (MITOS2) + FIX corals (fixer) + untouched non-Cnidarian invert phyla
    // (MITOS2, no gate applied). co1/12S/16S come straight from EMMA for every
    // vertebrate sample -- neither rescue touches those, so BLAST/LCA is
    // unaffected by the rescue outcomes.
    ch_annot_co1     = EMMA.out.co1_sequences.mix(ch_mitos_co1_pass, CORAL_ANNOTATION_FIX.out.co1_sequences, ch_mitos_co1_passthrough)
    ch_annot_s12     = EMMA.out.s12_sequences.mix(ch_mitos_s12_pass, CORAL_ANNOTATION_FIX.out.s12_sequences, ch_mitos_s12_passthrough)
    ch_annot_s16     = EMMA.out.s16_sequences.mix(ch_mitos_s16_pass, CORAL_ANNOTATION_FIX.out.s16_sequences, ch_mitos_s16_passthrough)
    ch_annot_results = ch_emma_results_trna.mix(ch_mitos_results_pass, CORAL_ANNOTATION_FIX.out.results, ch_mitos_results_passthrough)

    // Per-sample region count, emitted once per annotated sample (including the
    // 0-region case). Downstream this both sizes the result groups and identifies
    // the samples that will never reach BLAST/LCA at all.
    ch_annot_region_counts = ch_annot_results
        .map { meta, files -> [ meta, countAnnotatedRegions(files) ] }
    ch_annot_params  = EMMA.out.tool_params.mix(EMMA_GENE_RESCUE.out.tool_params, TRNA_RESCUE.out.tool_params,
                                                ROTATE_ORIGIN.out.tool_params,
                                                MITOS2.out.tool_params, CORAL_ANNOTATION_FIX.out.tool_params)
    ch_annot_versions = EMMA.out.versions.mix(EMMA_GENE_RESCUE_GATE.out.versions, EMMA_GENE_RESCUE.out.versions,
                                              TRNA_RESCUE_GATE.out.versions, TRNA_SCAN.out.versions, TRNA_RESCUE.out.versions,
                                              ROTATE_ORIGIN.out.versions, MITOS2.out.versions,
                                              ANNOTATION_QC_GATE.out.versions, CORAL_ANNOTATION_FIX.out.versions,
                                              SELECT_CORAL_REFERENCE.out.versions)

    //
    // Use mix() to process CO1, 12s and 16s sequences through blast
    //

    combined_sequences = ch_annot_co1
        .map { meta, file ->
            def annotation_name = getAnnotationName(file.name)
            [meta, file, 'CO1', annotation_name]
        }
        .mix(
            ch_annot_s12.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '12s', annotation_name]
            },
            ch_annot_s16.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '16s', annotation_name]
            }
        )
    
    //
    // Select the BLAST database per sample (curated for vertebrates, nt for invertebrates)
    //

    combined_sequences_with_db = combined_sequences
        .map { meta, file, gene_type, annotation_name ->
            if (meta.invertebrates && !nt_blast_db) {
                error "Sample ${meta.id} is marked invertebrates=true, but --nt_blast_db was not provided"
            }
            def blast_db = (meta.invertebrates ? nt_blast_db : curated_blast_db)
            [meta, file, gene_type, annotation_name, blast_db]
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
    ch_worms = channel.fromPath("${projectDir}/assets/worms_species.txt.gz", checkIfExists: true)
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
        // valid_blast_results, // tuple val(meta), path(blast_filtered), val(gene_type), val(annotation_name)
        // DOWNLOAD_TAXONKIT_DB.out.db_files // path(db)
    )


    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    // Need to update this section to include everything
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(ch_annot_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(LCA.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(ch_annot_versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())
    ch_versions = ch_versions.mix(PREPARE_LCA_DATABASES.out.versions)
    ch_versions = ch_versions.mix(REFERENCE_RELEVANCE.out.versions.first())



    //
    // Emit outputs
    //

    emit:
    multiqc_files           = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    annotation_results      = ch_annot_results
    blast_filtered_results  = BLAST_BLASTN.out.validation
    lca_results             = LCA.out.lca
    lca_raw_results         = LCA.out.lca_raw
    region_counts           = ch_annot_region_counts   // channel: [ val(meta), val(0..3) ]
    reference_relevance     = ch_reference_relevance   // channel: [ meta, path(reference_relevance.txt) ] (PASS|MISMATCH|UNKNOWN)
    versions                = ch_versions              // channel: [ path(versions.yml) ]
}
