CREATE TABLE IF NOT EXISTS ena_validation_attempts (
    id BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    assembly_prefix TEXT NOT NULL,
    og_id TEXT NOT NULL,
    tech TEXT,
    seq_date TEXT,
    code TEXT,
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
    flatfile_name TEXT,
    flatfile_sha256 CHAR(64),
    flatfile_size BIGINT,
    manifest_name TEXT,
    manifest_sha256 CHAR(64),
    manifest_size BIGINT,
    workflow_run_name TEXT,
    workflow_session_id TEXT,
    pipeline_revision TEXT,
    result_digest CHAR(64) NOT NULL,
    recorded_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_validation_submission_ready_check
        CHECK (NOT submission_ready OR webin_status = 'PASS')
);

-- Guarded on result_digest still existing.  bin/apply_ena_migrations.py replays
-- the whole chain on every run, and 012_drop_ena_selection_layer.sql drops that
-- column, so on a database already past 012 this index cannot be built.  Losing
-- nothing: 002_ena_validation_attempts_single_row_per_attempt.sql drops this
-- index one migration later and replaces it with ena_validation_attempts_key_idx.
DO $exact$
BEGIN
    -- to_regclass, not a bare relname: bare names match across every schema, so
    -- a migration applied into a throwaway schema would consult the public copy.
    IF EXISTS (
        SELECT 1 FROM pg_attribute
        WHERE attrelid = to_regclass('ena_validation_attempts')
          AND attname = 'result_digest'
          AND NOT attisdropped
    ) THEN
        CREATE UNIQUE INDEX IF NOT EXISTS ena_validation_exact_result_idx
            ON ena_validation_attempts (assembly_prefix, ena_study, validation_attempt, result_digest);
    END IF;
END
$exact$;

CREATE INDEX IF NOT EXISTS ena_validation_attempts_identity_idx
    ON ena_validation_attempts (og_id, tech, seq_date, code);

CREATE INDEX IF NOT EXISTS ena_validation_attempts_recorded_at_idx
    ON ena_validation_attempts (recorded_at DESC);

CREATE OR REPLACE VIEW ena_validation_latest AS
SELECT DISTINCT ON (assembly_prefix)
    *
FROM ena_validation_attempts
ORDER BY assembly_prefix, recorded_at DESC, id DESC;
