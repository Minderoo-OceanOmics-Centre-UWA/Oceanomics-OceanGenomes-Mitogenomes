process DOWNLOAD_BLAST_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.blast_db_dir // Cache the database

    input:
    val db_name

    output:
    path "${db_name}*"           , emit: db_files
    path "taxonomy4blast.sqlite3", emit: taxdb_file, optional: true

    script:
    db_path = params.blast_db_dir

    """
    echo "Downloading taxonomy database..."

    # storeDir already short-circuits this task when the central cache is complete.
    # This branch covers a partially-populated cache: reuse what is there rather
    # than re-pulling 62 MB. It must populate the work dir, not merely skip, because
    # the declared taxdb* outputs are resolved from the work dir.
    if [ -f "${db_path}/taxdb.btd" ] && [ -f "${db_path}/taxdb.bti" ]; then
        echo "Reusing cached taxonomy database from ${db_path}"
        cp -f "${db_path}"/taxdb.btd "${db_path}"/taxdb.bti .
        if [ -f "${db_path}/taxonomy4blast.sqlite3" ]; then
            cp -f "${db_path}/taxonomy4blast.sqlite3" .
        fi
    else
        echo "Downloading fresh taxonomy database..."

        # https, not ftp: the FTP control/data channel split is what dropped
        # mid-transfer and left a full-length but corrupt file after resume.
        wget --continue --tries=5 --retry-connrefused --waitretry=10 --timeout=60 \\
            https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
        wget --tries=5 --retry-connrefused --timeout=60 \\
            https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.md5

        # Never hand tar an unverified archive: wget reporting a complete save is
        # not proof of integrity when the transfer was resumed.
        if ! md5sum --check taxdb.tar.gz.md5; then
            echo "ERROR: taxdb.tar.gz failed md5 verification; discarding partial download" >&2
            rm -f taxdb.tar.gz taxdb.tar.gz.md5
            exit 1
        fi

        tar -xzf taxdb.tar.gz
        rm -f taxdb.tar.gz taxdb.tar.gz.md5
    fi

    # Verify files
    ls -la taxdb.*
    """

    stub:
    """
    touch taxdb.btd taxdb.bti taxonomy4blast.sqlite3
    """
}
