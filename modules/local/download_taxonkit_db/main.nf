process DOWNLOAD_TAXONKIT_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.taxonkit_db_dir ?: "${launchDir}/"  // Cache the database
    
    input:
    val db_name
    
    output:
    path("taxonkit_dbs"), emit: db_files
   
    script:
    """
    mkdir -p taxonkit_dbs

    # Check if files already exist
    if [ -f "${launchDir}/taxonkit_dbs/citations.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/delnodes.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/division.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/gencode.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/images.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/merged.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/names.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/nodes.dmp" ] \\
    && [ -f "${launchDir}/taxonkit_dbs/gc.prt" ]; then
        echo "Files already exist — linking to work dir"
        for f in citations.dmp delnodes.dmp division.dmp gencode.dmp images.dmp merged.dmp names.dmp nodes.dmp gc.prt; do
            ln -s "${launchDir}/taxonkit_dbs/\$f" taxonkit_dbs/
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
    mkdir -p taxonkit_dbs/
    touch taxonkit_dbs/citations.dmp \\
        taxonkit_dbs/delnodes.dmp \\
        taxonkit_dbs/division.dmp \\
        taxonkit_dbs/gencode.dmp \\
        taxonkit_dbs/images.dmp \\
        taxonkit_dbs/merged.dmp \\
        taxonkit_dbs/names.dmp \\
        taxonkit_dbs/nodes.dmp gc.prt
    """
}