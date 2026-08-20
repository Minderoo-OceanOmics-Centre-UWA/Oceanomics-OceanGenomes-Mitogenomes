-- Collapses ena_validation_attempts to one row per
-- (assembly_prefix, ena_study, validation_attempt) so pipeline reruns overwrite
-- the previous attempt instead of appending a new history row every time.
-- sql/004 plus push_ena_validation_results.py later supersede the original
-- freeze rule: rows remain overwritable until the corresponding selected
-- archive record is SUBMITTED or ACCESSION_ASSIGNED.
--
-- Run once, manually, the same way as 001_create_ena_validation_attempts.sql.
-- Connection details are in /home/tpeirce/postgresql_details/oceanomics.cfg (the
-- same file the pipeline passes as --sql_config). Note the database is
-- oceanomics_genomes; the "--dbname oceanomics" this file used to name does not exist:
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/002_ena_validation_attempts_single_row_per_attempt.sql

BEGIN;

-- Keep exactly one row per key: prefer a submission_ready row over a failed
-- one, then the most recently recorded attempt. Deletes every row that has
-- some other row in its group with a strictly "greater" (submission_ready,
-- recorded_at, id) tuple, which always leaves exactly one survivor per group.
DO $dedupe$
BEGIN
    -- 014 replaces assembly_prefix with full_seqid.  bin/apply_ena_migrations.py
    -- replays the whole chain on every run, so this legacy statement has to
    -- stand down on a database that is already past 014.
    IF EXISTS (
        SELECT 1 FROM pg_attribute
        WHERE attrelid = to_regclass('ena_validation_attempts')
          AND attname = 'assembly_prefix'
          AND NOT attisdropped
    ) THEN
        EXECUTE $dedupe_sql$
            DELETE FROM ena_validation_attempts a
            USING ena_validation_attempts b
            WHERE a.assembly_prefix = b.assembly_prefix
              AND a.ena_study = b.ena_study
              AND a.validation_attempt = b.validation_attempt
              AND (b.submission_ready, b.recorded_at, b.id) > (a.submission_ready, a.recorded_at, a.id)
        $dedupe_sql$;
    END IF;
END
$dedupe$;

DROP INDEX IF EXISTS ena_validation_exact_result_idx;

ALTER TABLE ena_validation_attempts
    ADD COLUMN IF NOT EXISTS attempt_count INTEGER NOT NULL DEFAULT 1;

DO $key_idx$
BEGIN
    -- 014 replaces assembly_prefix with full_seqid.  bin/apply_ena_migrations.py
    -- replays the whole chain on every run, so this legacy statement has to
    -- stand down on a database that is already past 014.
    IF EXISTS (
        SELECT 1 FROM pg_attribute
        WHERE attrelid = to_regclass('ena_validation_attempts')
          AND attname = 'assembly_prefix'
          AND NOT attisdropped
    ) THEN
        EXECUTE $key_sql$
            CREATE UNIQUE INDEX IF NOT EXISTS ena_validation_attempts_key_idx
                ON ena_validation_attempts (assembly_prefix, ena_study, validation_attempt)
        $key_sql$;
    END IF;
END
$key_idx$;

COMMIT;
