process DOWNLOAD_BLAST_DB {
    tag "${db_name}"
    label 'process_low'
    storeDir params.blast_db_dir // Cache the database

    input:
    val db_name

    output:
    path "${db_name}*"           , emit: db_files
    // Required, NOT optional. blastn's -taxids filter resolves taxonomy from this
    // file, and because taxdb.* is staged into the task CWD, BLAST looks in the CWD
    // rather than beside the -db path -- so the sqlite3 must be staged with them.
    // Leaving it optional is what silently disabled -taxids for 184 of 211 searches
    // in NOVA_260909_LA_mitogenomes-3 (empty result table, exit 0). Making it
    // optional again would also starve BLAST_BLASTN's taxdb input on an empty
    // channel and hang the run, which is worse than the bug it hides.
    path "taxonomy4blast.sqlite3", emit: taxdb_file

    script:
    db_path = params.blast_db_dir

    """
    echo "Downloading taxonomy database..."

    # storeDir already short-circuits this task when the central cache is complete.
    # This branch covers a partially-populated cache: reuse what is there rather
    # than re-pulling 62 MB. It must populate the work dir, not merely skip, because
    # the declared taxdb* outputs are resolved from the work dir.
    if [ -f "${db_path}/taxdb.btd" ] && [ -f "${db_path}/taxdb.bti" ] && [ -f "${db_path}/taxonomy4blast.sqlite3" ]; then
        echo "Reusing cached taxonomy database from ${db_path}"
        cp -f "${db_path}"/taxdb.btd "${db_path}"/taxdb.bti "${db_path}"/taxonomy4blast.sqlite3 .
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

    # Verify all three files, not just the two that taxdb.* globs. This process is
    # upstream of every blastn task, so a cache missing the sqlite3 costs one failed
    # task here instead of a whole cohort of silently empty -taxids searches.
    missing=""
    for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
        [ -s "\$f" ] || missing="\$missing \$f"
    done
    if [ -n "\$missing" ]; then
        echo "ERROR: taxonomy database incomplete after staging; missing or empty:\$missing" >&2
        echo "       Cache dir (--blast_db_dir): ${db_path}" >&2
        echo "       All three files ship in taxdb.tar.gz; remove the cache dir to force a clean re-download." >&2
        exit 1
    fi
    ls -la taxdb.btd taxdb.bti taxonomy4blast.sqlite3
    """

    stub:
    """
    touch taxdb.btd taxdb.bti taxonomy4blast.sqlite3
    """
}
