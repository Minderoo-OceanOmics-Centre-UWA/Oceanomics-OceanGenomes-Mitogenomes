// Build a custom GetOrganelle label (gene) database from a reference GenBank.
//
// GetOrganelle uses TWO databases for very different jobs (see the FAQ:
// https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference):
//   * the SEED database (-s) only recruits the initial overlapping reads, and
//   * the LABEL/gene database (--genes) classifies and disentangles contigs
//     during extension.
// For divergent animal mitogenomes the label database is the one that matters
// most. The reseed pass already supplies a closely-related mitogenome as the
// custom seed; this module turns the GenBank record that MITOHIFI_FINDMITOREFERENCE
// downloads alongside that seed into a matching custom label database, using
// GetOrganelle's own get_annotated_regions_from_gb.py utility (same container).
//
// The reseed is only ever run WITH a custom gene database (the user's rule:
// "always use --genes for the reseed"). To honour that without crashing samples
// whose nearest NCBI relative happens to be poorly annotated, the gene database
// is emitted ONLY when it clears a minimum CDS-count threshold. When it does
// not, the (optional) output is absent and the subworkflow falls the sample
// back to its first-pass result instead of reseeding without genes.
process GETORGANELLE_GENEDB {
    tag "$meta.id"
    label 'process_single'

    // MitoHiFi container (same image findMitoReference already uses) because it
    // ships Biopython. The GetOrganelle BioConda image does NOT, so its bundled
    // get_annotated_regions_from_gb.py prints "biopython not found" and silently
    // emits nothing. We extract the genes ourselves with bin/extract_getorganelle_genedb.py.
    conda "conda-forge::biopython=1.79"
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    tuple val(meta), path(gb)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}reseed.genedb.fasta"), emit: genes, optional: true
    tuple val(meta), path("02c_getorganelle_genedb.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def types     = task.ext.args ?: 'CDS,rRNA'
    def prefix    = "${meta.mt_assembly_prefix}reseed"
    def min_genes = params.getorganelle_genedb_min_genes ?: 10
    """
    set -uo pipefail

    # Extract the annotated gene regions (CDS + rRNA by default) from the
    # reference GenBank into a single FASTA whose headers follow the convention
    # of GetOrganelle's bundled animal_mt label database ('>gene gene - ACCESSION'),
    # which its --genes parser expects. Tolerate a non-zero exit (e.g. a malformed
    # record) and let the gene-count check below decide whether a usable database
    # was produced.
    extract_getorganelle_genedb.py \\
        $gb \\
        ${prefix}.genedb.candidate.fasta \\
        --types ${types} || echo "WARN: extract_getorganelle_genedb.py exited non-zero for ${meta.id}" >&2

    gene_count=0
    if [ -f ${prefix}.genedb.candidate.fasta ]; then
        gene_count=\$(grep -c '^>' ${prefix}.genedb.candidate.fasta || true)
    fi

    # Only publish the custom label database when it is well annotated enough to
    # be useful as a GetOrganelle gene database. A complete animal mitogenome has
    # 13 protein-coding genes; require at least ${min_genes}. Otherwise leave the
    # (optional) output absent so the reseed is skipped rather than run with a
    # sparse / empty gene database.
    if [ "\$gene_count" -ge ${min_genes} ]; then
        mv ${prefix}.genedb.candidate.fasta ${prefix}.genedb.fasta
        genedb_status="built (\${gene_count} genes)"
    else
        genedb_status="skipped (\${gene_count} genes < ${min_genes} threshold)"
        echo "WARN: ${meta.id} reference GenBank yielded only \${gene_count} genes (< ${min_genes}); custom gene database not built, reseed will fall back to first-pass." >&2
    fi

    cat <<-END_TOOL_PARAMS > 02c_getorganelle_genedb.tool_params_mqcrow.html
    <tr><td>GetOrganelle Gene Database</td><td><samp>extract_getorganelle_genedb.py ${gb} --types ${types}</samp></td><td>Builds a custom GetOrganelle label (gene) database from the reference for ${meta.id}: \${genedb_status}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def types  = task.ext.args ?: 'CDS,rRNA'
    def prefix = "${meta.mt_assembly_prefix}reseed"
    """
    touch ${prefix}.genedb.fasta

    cat <<-END_TOOL_PARAMS > 02c_getorganelle_genedb.tool_params_mqcrow.html
    <tr><td>GetOrganelle Gene Database</td><td><samp>extract_getorganelle_genedb.py ${gb} --types ${types}</samp></td><td>Builds a custom GetOrganelle label (gene) database from the reference for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
