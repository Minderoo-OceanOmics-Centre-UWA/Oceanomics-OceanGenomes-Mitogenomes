-- Adds og_num to ena_validation_attempts, as the second column (right after
-- id), following the same convention as sample, draft_genomes, lca,
-- lca_raw_results, sequencing and mitogenome_data (migration 018): a stored
-- generated column parsing the numeric suffix off og_id (OG703 -> 703).
--
-- PostgreSQL cannot add a column at an arbitrary position or reorder columns
-- in place, so this is a table rebuild, same shape as 014 and 018: build the
-- table with columns in the wanted order, copy rows across, swap it in.
--
-- Migration 019 exists because an earlier rebuild-style migration silently
-- reset a column for every row and nothing caught it until the corruption
-- was found independently, well after the fact. This rebuild does not repeat
-- that: before the old table is dropped, the copy is verified row-for-row
-- against the source with a bidirectional EXCEPT diff over every carried
-- column, not just a row count. Any mismatch raises inside the transaction,
-- so the old table is never dropped and nothing is renamed.
--
-- Guarded like 014 and 018: bin/apply_ena_migrations.py replays the whole
-- chain on every run, so this is a no-op once og_num already sits at
-- column 2, or where the table does not exist.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/020_ena_validation_attempts_og_num.sql

BEGIN;

SET LOCAL lock_timeout = '5s';

DO $do$
DECLARE
    dependent_views TEXT[][] := ARRAY[]::TEXT[][];
    dependent   RECORD;
    idx         INT;
    orig        BIGINT;
    copied      BIGINT;
    mismatched  BIGINT;
    null_og_num BIGINT;
BEGIN
    IF to_regclass('ena_validation_attempts') IS NULL THEN
        RETURN;
    END IF;

    -- Already rebuilt: nothing to do.
    IF EXISTS (
        SELECT 1 FROM pg_attribute
        WHERE attrelid = to_regclass('ena_validation_attempts')
          AND attname = 'og_num' AND attnum = 2
          AND NOT attisdropped
    ) THEN
        RETURN;
    END IF;

    -- ena_validation_latest is SELECT DISTINCT ON (full_seqid) *, so it
    -- depends on every column and has to go first. Restored below, unchanged.
    DROP VIEW IF EXISTS ena_validation_latest;

    FOR dependent IN
        SELECT v.oid::regclass::text AS view_name,
               pg_get_viewdef(v.oid, true) AS view_def
        FROM pg_depend d
        JOIN pg_rewrite r ON r.oid = d.objid
        JOIN pg_class v ON v.oid = r.ev_class
        WHERE d.refobjid = to_regclass('ena_validation_attempts')
          AND v.relkind = 'v'
          AND v.oid <> to_regclass('ena_validation_attempts')
        GROUP BY v.oid
    LOOP
        dependent_views := dependent_views || ARRAY[ARRAY[dependent.view_name, dependent.view_def]];
        EXECUTE format('DROP VIEW %s', dependent.view_name);
    END LOOP;

    CREATE TABLE ena_validation_attempts_new (
        id BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
        og_num INTEGER GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::INTEGER) STORED,
        full_seqid TEXT NOT NULL,
        og_id TEXT NOT NULL,
        tech TEXT,
        seq_date TEXT,
        code TEXT,
        annotation TEXT,
        ena_study TEXT NOT NULL DEFAULT '',
        validation_mode TEXT NOT NULL,
        validation_attempt TEXT NOT NULL,
        table2asn_status TEXT NOT NULL,
        reject_count INTEGER,
        error_count INTEGER,
        warning_count INTEGER,
        info_count INTEGER,
        fatal_discrepancy_count INTEGER,
        nostop_count INTEGER,
        blocking_codes TEXT,
        warning_codes TEXT,
        conversion_status TEXT NOT NULL,
        conversion_reason TEXT,
        conversion_exit INTEGER,
        preflight_status TEXT NOT NULL,
        preflight_reason TEXT,
        preflight_exit INTEGER,
        webin_status TEXT NOT NULL,
        webin_reason TEXT,
        webin_exit INTEGER,
        submission_ready BOOLEAN NOT NULL DEFAULT FALSE,
        recorded_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
        attempt_count INTEGER NOT NULL DEFAULT 1,
        CONSTRAINT ena_validation_submission_ready_check
            CHECK (NOT submission_ready OR webin_status = 'PASS')
    );

    -- og_num is omitted: a generated column rejects an explicit value.
    INSERT INTO ena_validation_attempts_new (
        id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
        validation_mode, validation_attempt, table2asn_status,
        reject_count, error_count, warning_count, info_count,
        fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
        conversion_status, conversion_reason, conversion_exit,
        preflight_status, preflight_reason, preflight_exit,
        webin_status, webin_reason, webin_exit, submission_ready,
        recorded_at, attempt_count
    )
    OVERRIDING SYSTEM VALUE
    SELECT
        id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
        validation_mode, validation_attempt, table2asn_status,
        reject_count, error_count, warning_count, info_count,
        fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
        conversion_status, conversion_reason, conversion_exit,
        preflight_status, preflight_reason, preflight_exit,
        webin_status, webin_reason, webin_exit, submission_ready,
        recorded_at, attempt_count
    FROM ena_validation_attempts;

    -- Row count alone would not have caught 019's silent corruption, so
    -- verify every carried column, not just the count.
    SELECT count(*) INTO orig   FROM ena_validation_attempts;
    SELECT count(*) INTO copied FROM ena_validation_attempts_new;
    IF copied <> orig THEN
        RAISE EXCEPTION 'ena_validation_attempts rebuild copied % of % rows', copied, orig;
    END IF;

    SELECT count(*) INTO mismatched FROM (
        SELECT id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
               validation_mode, validation_attempt, table2asn_status,
               reject_count, error_count, warning_count, info_count,
               fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
               conversion_status, conversion_reason, conversion_exit,
               preflight_status, preflight_reason, preflight_exit,
               webin_status, webin_reason, webin_exit, submission_ready,
               recorded_at, attempt_count
        FROM ena_validation_attempts
        EXCEPT
        SELECT id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
               validation_mode, validation_attempt, table2asn_status,
               reject_count, error_count, warning_count, info_count,
               fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
               conversion_status, conversion_reason, conversion_exit,
               preflight_status, preflight_reason, preflight_exit,
               webin_status, webin_reason, webin_exit, submission_ready,
               recorded_at, attempt_count
        FROM ena_validation_attempts_new
    ) old_not_in_new;
    IF mismatched <> 0 THEN
        RAISE EXCEPTION 'ena_validation_attempts rebuild: % source rows have no matching copy', mismatched;
    END IF;

    SELECT count(*) INTO mismatched FROM (
        SELECT id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
               validation_mode, validation_attempt, table2asn_status,
               reject_count, error_count, warning_count, info_count,
               fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
               conversion_status, conversion_reason, conversion_exit,
               preflight_status, preflight_reason, preflight_exit,
               webin_status, webin_reason, webin_exit, submission_ready,
               recorded_at, attempt_count
        FROM ena_validation_attempts_new
        EXCEPT
        SELECT id, full_seqid, og_id, tech, seq_date, code, annotation, ena_study,
               validation_mode, validation_attempt, table2asn_status,
               reject_count, error_count, warning_count, info_count,
               fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
               conversion_status, conversion_reason, conversion_exit,
               preflight_status, preflight_reason, preflight_exit,
               webin_status, webin_reason, webin_exit, submission_ready,
               recorded_at, attempt_count
        FROM ena_validation_attempts
    ) new_not_in_old;
    IF mismatched <> 0 THEN
        RAISE EXCEPTION 'ena_validation_attempts rebuild: % copied rows do not match their source', mismatched;
    END IF;

    -- og_id is NOT NULL, so the generated og_num should never be NULL. A
    -- NULL here means some og_id does not match the OG<digits> pattern the
    -- generation expression assumes, and the rebuild should not proceed
    -- silently on top of that.
    SELECT count(*) INTO null_og_num FROM ena_validation_attempts_new WHERE og_num IS NULL;
    IF null_og_num <> 0 THEN
        RAISE EXCEPTION 'ena_validation_attempts rebuild: % rows have a NULL og_num (og_id not OG<digits>)', null_og_num;
    END IF;

    DROP TABLE ena_validation_attempts;
    ALTER TABLE ena_validation_attempts_new RENAME TO ena_validation_attempts;

    IF to_regclass('ena_validation_attempts_new_pkey') IS NOT NULL THEN
        ALTER INDEX ena_validation_attempts_new_pkey
            RENAME TO ena_validation_attempts_pkey;
    END IF;

    -- Carried-over ids were inserted around the identity sequence, so move it
    -- past them before the next insert.
    PERFORM setval(
        pg_get_serial_sequence('ena_validation_attempts', 'id'),
        COALESCE((SELECT max(id) FROM ena_validation_attempts), 1),
        (SELECT count(*) > 0 FROM ena_validation_attempts)
    );

    -- The ON CONFLICT target used by push_ena_validation_results.py.
    CREATE UNIQUE INDEX ena_validation_attempts_key_idx
        ON ena_validation_attempts (full_seqid, ena_study, validation_attempt);

    CREATE INDEX ena_validation_attempts_identity_idx
        ON ena_validation_attempts (og_id, tech, seq_date, code, annotation);

    CREATE INDEX ena_validation_attempts_recorded_at_idx
        ON ena_validation_attempts (recorded_at DESC);

    CREATE VIEW ena_validation_latest AS
    SELECT DISTINCT ON (full_seqid)
        *
    FROM ena_validation_attempts
    ORDER BY full_seqid, recorded_at DESC, id DESC;

    -- Dropped in dependency order, so restore in reverse.
    FOR idx IN REVERSE COALESCE(array_length(dependent_views, 1), 0) .. 1 LOOP
        EXECUTE format('CREATE VIEW %s AS %s', dependent_views[idx][1], dependent_views[idx][2]);
    END LOOP;

    COMMENT ON COLUMN ena_validation_attempts.og_num IS
        'Numeric part of og_id, maintained by the database. Generated, not writable: '
        'do not include it in any INSERT or UPDATE column list.';
END
$do$;

COMMIT;

-- Sanity check after running:
--
--   SELECT attnum, attname, attgenerated FROM pg_attribute
--   WHERE attrelid = 'ena_validation_attempts'::regclass AND attname = 'og_num';
--   -- expect: 2 | og_num | s
