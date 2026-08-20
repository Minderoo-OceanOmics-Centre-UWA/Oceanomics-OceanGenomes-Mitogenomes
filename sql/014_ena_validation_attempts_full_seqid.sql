-- Key ena_validation_attempts on the full seq id instead of the assembly prefix.
--
-- What this pipeline validates is a flatfile built from one annotation of one
-- assembly: OG82.ilmn.240313.getorg1770.emma102, not OG82.ilmn.240313.getorg1770.
-- Every other ENA artifact (the EMBL, the manifest, the package, the metadata
-- JSON) already carries that full seq id; only the validation record was keyed
-- on the 4-field prefix, so re-annotating an assembly overwrote the record for
-- the annotation validated before it, and no row could say which annotation
-- earned submission_ready.
--
-- assembly_prefix is replaced by full_seqid, and the annotation version becomes
-- its own column sitting immediately after code.  Column order is part of the
-- request, so this is a table rebuild rather than an ADD COLUMN.
--
-- Existing rows are carried over with full_seqid = assembly_prefix and a NULL
-- annotation: the annotation they were validated with is not recoverable from
-- the table, and they are superseded the next time those assemblies are
-- validated (which inserts properly keyed rows, leaving the legacy placeholder
-- to be deleted).
--
-- Guarded like 009, 011 and 012: bin/apply_ena_migrations.py replays the whole
-- chain on every run and the integration tests build partial schemas, so this
-- is a no-op once full_seqid exists, or where the table does not.
--
-- Run after 013_drop_ena_candidate_runs.sql.

BEGIN;

DO $do$
DECLARE
    -- Views outside this repository also read the table (the live database has
    -- mitogenome_submission_view).  A rebuild has to drop them and put them back
    -- exactly as they were, so capture their definitions rather than assuming
    -- ena_validation_latest is the only dependent.
    dependent_views TEXT[][] := ARRAY[]::TEXT[][];
    dependent RECORD;
    idx INT;
BEGIN
    IF to_regclass('ena_validation_attempts') IS NULL THEN
        RETURN;
    END IF;

    -- Already rebuilt: nothing to do.
    IF EXISTS (
        SELECT 1 FROM pg_attribute
        WHERE attrelid = to_regclass('ena_validation_attempts')
          AND attname = 'full_seqid'
          AND NOT attisdropped
    ) THEN
        RETURN;
    END IF;

    -- ena_validation_latest is SELECT DISTINCT ON (assembly_prefix) *, so it
    -- depends on every column of the old table and has to go first.  014
    -- recreates it below, keyed on full_seqid.
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
        id, assembly_prefix, og_id, tech, seq_date, code, NULL, ena_study,
        validation_mode, validation_attempt, table2asn_status,
        reject_count, error_count, warning_count, info_count,
        fatal_discrepancy_count, nostop_count, blocking_codes, warning_codes,
        conversion_status, conversion_reason, conversion_exit,
        preflight_status, preflight_reason, preflight_exit,
        webin_status, webin_reason, webin_exit, submission_ready,
        recorded_at, attempt_count
    FROM ena_validation_attempts;

    DROP TABLE ena_validation_attempts;
    ALTER TABLE ena_validation_attempts_new RENAME TO ena_validation_attempts;

    -- The primary key came across under the build name.
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

    COMMENT ON COLUMN ena_validation_attempts.full_seqid IS
        'Assembly prefix plus annotation version (OG82.ilmn.240313.getorg1770.emma102): '
        'the id the validated flatfile, manifest and package all carry.';
    COMMENT ON COLUMN ena_validation_attempts.annotation IS
        'Annotation version, the fifth field of full_seqid. NULL on rows carried '
        'over from before the record was keyed on the full seq id.';
END
$do$;

COMMIT;
