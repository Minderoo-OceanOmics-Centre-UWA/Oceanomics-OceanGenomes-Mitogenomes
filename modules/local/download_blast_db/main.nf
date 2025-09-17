process DOWNLOAD_BLAST_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.blast_db_dir // Cache the database
    
    input:
    val db_name
    
    output:
    path "${db_name}*", emit: db_files
    
    script:
    db_path = params.blast_db_dir

    """
    echo "Downloading taxonomy database..."
    
    # Check if files already exist
    if [ -f "${db_path}/blast_dbs/taxdb.btd" ] && [ -f "${db_path}/blast_dbs/taxdb.bti" ]; then
        echo "Taxonomy files already exist, skipping download"
    else
        echo "Downloading fresh taxonomy database..."
        wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz

        wait
        
        tar -xzf taxdb.tar.gz
        rm taxdb.tar.gz
    fi
    
    # Verify files
    ls -la taxdb.*
    """
    
    stub:
    """
    touch taxdb.btd taxdb.bti
    """
}