process DOWNLOAD_TAXONKIT_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.taxonkit_db_dir // Cache the database
    
    input:
    val db_name
    
    output:
    path("taxonkit_dbs"), emit: db_files
   
    script:
    db_dir = params.taxonkit_db_dir 
    """
    mkdir -p taxonkit_dbs

    # Check if files already exist
    if [ -f "${db_dir}/taxonkit_dbs/citations.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/delnodes.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/division.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/gencode.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/images.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/merged.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/names.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/nodes.dmp" ] \\
    && [ -f "${db_dir}/taxonkit_dbs/gc.prt" ]; then
        echo "Files already exist — linking to work dir"
        for f in citations.dmp delnodes.dmp division.dmp gencode.dmp images.dmp merged.dmp names.dmp nodes.dmp gc.prt; do
            ln -s "${db_dir}/taxonkit_dbs/\$f" taxonkit_dbs/
        done
    else
        echo "Downloading fresh taxonomy database..."
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf taxdump.tar.gz -C taxonkit_dbs
        wait
        rm taxdump.tar.gz
    fi
    
    """
    
    stub:
    """
    mkdir -p ${params.taxonkit_db_dir}/taxonkit_db1s/

    taxonkit_dir = ${params.taxonkit_db_dir}/taxonkit_db1s/
    touch \$taxonkit_dir/citations.dmp \\
        \$taxonkit_dir/delnodes.dmp \\
        \$taxonkit_dir/division.dmp \\
        \$taxonkit_dir/gencode.dmp \\
        \$taxonkit_dir/images.dmp \\
        \$taxonkit_dir/merged.dmp \\
        \$taxonkit_dir/names.dmp \\
        \$taxonkit_dir/nodes.dmp gc.prt
    """
}