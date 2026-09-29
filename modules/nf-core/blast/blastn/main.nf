process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    tuple val(meta), path(fasta), val(gene_type), val(annotation_name), val(blast_db)
    path(taxdb)

    output:
    tuple val(meta), path("blast.${gene_type}.${annotation_name}.filtered.tsv"), val(gene_type), val(annotation_name), emit: filtered
    tuple val(meta), path("blast.${gene_type}.${annotation_name}.filtered.tsv"), emit: validation
    path "filtered_summary.${gene_type}.txt"   , emit: summary
    tuple val(meta), path("08_blast_blastn.${gene_type}.${annotation_name}.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "blast.${gene_type}.${annotation_name}"
    def effective_args = ["-num_threads ${task.cpus}", "-db ${blast_db}", "-query ${fasta}", args, "-out ${prefix}.tsv"].findAll { it?.toString()?.trim() }.join(' ')
    def stderr_log = "blastn.${gene_type}.${annotation_name}.stderr"
    // Keyed on the resolved args string rather than on meta.invertebrates, so it stays
    // correct if conf/modules.config re-keys which branch passes -taxids.
    def taxid_preflight = args.contains('-taxids') ? """
    if [ ! -s taxonomy4blast.sqlite3 ]; then
        echo "ERROR: ext.args requests -taxids but taxonomy4blast.sqlite3 was not staged." >&2
        echo "       BLAST resolves taxonomy from the CWD because taxdb.* is staged here," >&2
        echo "       so the taxid filter would be silently disabled: blastn prints a" >&2
        echo "       warning, writes an empty result table, and still exits 0." >&2
        echo "       Check DOWNLOAD_BLAST_DB.out.taxdb_file is mixed into ch_blast_db" >&2
        echo "       (subworkflows/local/mitogenome_annotation_lca)." >&2
        exit 1
    fi
""" : ''

    """
    mkdir -p lca

    # Make the CWD taxonomy lookup explicit rather than incidental. Safe because every
    # -db value this module is passed is an absolute path.
    export BLASTDB="\$PWD"
${taxid_preflight}
    # blastn is run unguarded so its native exit code can be inspected: the process
    # shell is `bash -C -e -u -o pipefail` (see process.shell in nextflow.config), and
    # under -e a non-zero blastn would abort before the diagnostics below.
    set +e
    blastn \\
        -num_threads ${task.cpus} \\
        -db ${blast_db} \\
        -query ${fasta} \\
        ${args} \\
        -out ${prefix}.tsv 2> ${stderr_log}
    blastn_exit=\$?
    set -e
    cat ${stderr_log} >&2

    # Preserve the native exit code so conf/base.config's errorStrategy can classify a
    # resource kill (137/140) apart from a tool error.
    if [ "\$blastn_exit" -ne 0 ]; then
        echo "ERROR: blastn exited \$blastn_exit for ${gene_type} / ${annotation_name}" >&2
        exit "\$blastn_exit"
    fi

    # blastn can report a fatal condition on stderr and STILL EXIT 0. That is how the
    # missing taxonomy4blast.sqlite3 escaped an entire run despite -e being armed: 184
    # of 211 searches returned an empty table with exit 0 and were read downstream as
    # "no hits". Treat a known-fatal signature as the failure it is.
    if grep -Eq 'requires additional data files|BLAST Database error|BLAST engine error|Taxonomy.*not found' ${stderr_log}; then
        echo "ERROR: blastn reported a fatal condition but exited 0; treating as a failure." >&2
        echo "       See ${stderr_log} in the task work dir." >&2
        exit 1
    fi

    # An empty table is a LEGITIMATE result for a query with no hits in the restricted
    # search space, so this is a warning and never an assertion.
    if [ ! -s ${prefix}.tsv ]; then
        echo "WARN: blastn returned no hits for ${gene_type} / ${annotation_name}" >&2
    fi

    cat <<-END_TOOL_PARAMS > 08_blast_blastn.${gene_type}.${annotation_name}.tool_params_mqcrow.html
    <tr><td>BLAST BLASTN</td><td><samp>${effective_args}</samp></td><td>Searches ${gene_type} sequence ${annotation_name} against the selected BLAST database for ${meta.id}; output is filtered for length and identity.</td></tr>
    END_TOOL_PARAMS

    # Filter the results
    awk -F '\t' '{if ((\$15 - \$14 > 200) && (\$7 > 98)) print}' ${prefix}.tsv > ${prefix}.filtered.tsv
    
    # Add in a column for the days date and the gene type
    sed -i "s/\$/\\t\$(date +%y%m%d)/" ${prefix}.filtered.tsv
    wait
    sed -i "s/\$/\\t${gene_type}/" ${prefix}.filtered.tsv

    # Create a summary for multiqc
    awk -F '\t' '{combo[\$1"\t"\$4]++} END {print "Sample\tSpecies\tHits"; for (c in combo) print c, combo[c]}' ${prefix}.filtered.tsv > filtered_summary.${gene_type}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    // prefix must match the script block's default, and the stub must create every
    // declared output: it previously touched only "${meta.id}.txt", so -stub runs of
    // this module failed on missing outputs and no test could use them.
    def prefix = task.ext.prefix ?: "blast.${gene_type}.${annotation_name}"
    def args = task.ext.args ?: ''
    def effective_args = ["-num_threads ${task.cpus}", "-db ${blast_db}", "-query ${fasta}", args, "-out ${prefix}.tsv"].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch ${prefix}.tsv
    touch ${prefix}.filtered.tsv
    printf 'Sample\\tSpecies\\tHits\\n' > filtered_summary.${gene_type}.txt

    cat <<-END_TOOL_PARAMS > 08_blast_blastn.${gene_type}.${annotation_name}.tool_params_mqcrow.html
    <tr><td>BLAST BLASTN</td><td><samp>${effective_args}</samp></td><td>Searches ${gene_type} sequence ${annotation_name} against the selected BLAST database for ${meta.id}; output is filtered for length and identity.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
