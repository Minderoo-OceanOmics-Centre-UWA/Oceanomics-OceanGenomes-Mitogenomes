process GEN_FILES_TABLE2ASN {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::table2asn"

    input:
    tuple val(meta), path(sample_fa), path(sample_tbl), path(sample_cmt), path(sample_src)  
    path sample_sbt // generic sbt file for the oceanomics lab that you can generate from GenBank.

    output:
    tuple val(meta), path("*.sqn")    , emit: sqn_file
    tuple val(meta), path("*.val")    , emit: val_file
    tuple val(meta), path("*.gbf")    , emit: gbf_file

    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def LOGFILE = "${meta.id}.table2asn.log"
    def sample_out = "${meta.mt_assembly_prefix}.sqn"
    
    """
    # Run table2asn with correct syntax
    table2asn \\
        -indir . \\
        -euk \\
        -J \\
        -t "$sample_sbt" \\
        -i "$sample_fa" \\
        -f "$sample_tbl" \\
        -w "$sample_cmt" \\
        -src-file "$sample_src" \\
        -o "$sample_out" \\
        -M n \\
        -j "[mgcode=2] [location=mitochondrion] [topology=circular]" \\
        -V vb \\
        -Z \\
        -W \\
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}