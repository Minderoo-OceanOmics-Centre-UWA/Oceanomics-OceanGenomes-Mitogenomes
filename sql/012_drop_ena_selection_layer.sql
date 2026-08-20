-- Retire the in-repo ENA selection layer and shrink ena_validation_attempts to
-- what this pipeline actually decides.
--
-- Selection and submission moved out of this repository.  ena_selection.nf,
-- select_ena_submission.py and record_ena_package_validation.py are gone, so
-- nothing writes ena_candidate_packages or ena_submission_selections any more,
-- and the ena_submission_queue view they fed has no producer.
--
-- submission_ready changes meaning with them.  It used to be earned only by a
-- production Webin validation on a package that was already SELECTED, which the
-- pipeline never ran, so the column was false on every row.  It now means what
-- this pipeline can actually attest: the flatfile was produced and cleared every
-- gate, so it is ready to hand to the submission pipeline.  That is the
-- constraint sql/001 originally carried, restored here.
--
-- Existing rows are deliberately NOT backfilled; they carry whatever verdict
-- they were written with and are corrected on their next pipeline run.
--
-- Run after 011_drop_local_package_validation.sql.
-- Guarded like 009 and 011: the integration tests build partial schemas, so
-- touch each relation only where it is already present.  Idempotent.

BEGIN;

DROP VIEW IF EXISTS ena_submission_queue;

-- ena_submission_selections references ena_candidate_packages, so it goes first.
DROP TABLE IF EXISTS ena_submission_selections;
DROP TABLE IF EXISTS ena_candidate_packages;

DO $do$
BEGIN
    IF to_regclass('ena_validation_attempts') IS NOT NULL THEN
        -- ena_validation_latest is SELECT DISTINCT ON (assembly_prefix) *, and
        -- Postgres expanded that * when the view was created.  It therefore
        -- depends on every column below and has to be dropped before them.
        -- Only when the table still has the pre-014 shape: after 014 the view is
        -- keyed on full_seqid, does not reference the columns touched below, and
        -- nothing here would recreate it if it were dropped.
        IF EXISTS (
            SELECT 1 FROM pg_attribute
            WHERE attrelid = to_regclass('ena_validation_attempts')
              AND attname = 'assembly_prefix'
              AND NOT attisdropped
        ) THEN
            DROP VIEW IF EXISTS ena_validation_latest;
        END IF;

        ALTER TABLE ena_validation_attempts
            DROP CONSTRAINT IF EXISTS ena_validation_submission_ready_check;

        ALTER TABLE ena_validation_attempts
            DROP CONSTRAINT IF EXISTS ena_validation_package_digest_check;

        -- Dropping result_digest also drops ena_validation_exact_result_idx.
        -- The ON CONFLICT target used by push_ena_validation_results.py is
        -- ena_validation_attempts_key_idx from sql/002, which is unaffected.
        ALTER TABLE ena_validation_attempts
            DROP COLUMN IF EXISTS package_status,
            DROP COLUMN IF EXISTS webin_production_status,
            DROP COLUMN IF EXISTS webin_error_count,
            DROP COLUMN IF EXISTS webin_warning_count,
            DROP COLUMN IF EXISTS webin_error_codes,
            DROP COLUMN IF EXISTS webin_cli_version,
            DROP COLUMN IF EXISTS webin_production_report_path,
            DROP COLUMN IF EXISTS overall_status,
            DROP COLUMN IF EXISTS package_digest,
            DROP COLUMN IF EXISTS flatfile_name,
            DROP COLUMN IF EXISTS flatfile_sha256,
            DROP COLUMN IF EXISTS flatfile_size,
            DROP COLUMN IF EXISTS manifest_name,
            DROP COLUMN IF EXISTS manifest_sha256,
            DROP COLUMN IF EXISTS manifest_size,
            DROP COLUMN IF EXISTS workflow_run_name,
            DROP COLUMN IF EXISTS workflow_session_id,
            DROP COLUMN IF EXISTS pipeline_revision,
            DROP COLUMN IF EXISTS result_digest;

        -- The collator has only ever written false, so this should hold on
        -- every stored row.  If a legacy row does not satisfy it, the ALTER
        -- fails loudly rather than silently rewriting history; fix that row, or
        -- add the constraint NOT VALID so it governs new writes only.
        ALTER TABLE ena_validation_attempts
            ADD CONSTRAINT ena_validation_submission_ready_check
                CHECK (NOT submission_ready OR webin_status = 'PASS');

        -- 014 replaces assembly_prefix with full_seqid and owns the view from
        -- then on; skip this once the rebuild has happened.
        IF EXISTS (
            SELECT 1 FROM pg_attribute
            WHERE attrelid = to_regclass('ena_validation_attempts')
              AND attname = 'assembly_prefix'
              AND NOT attisdropped
        ) THEN
            EXECUTE $view$
                CREATE VIEW ena_validation_latest AS
                SELECT DISTINCT ON (assembly_prefix)
                    *
                FROM ena_validation_attempts
                ORDER BY assembly_prefix, recorded_at DESC, id DESC
            $view$;
        END IF;

        COMMENT ON COLUMN ena_validation_attempts.submission_ready IS
            'True when the flatfile passed every gate this pipeline runs, ending at the '
            'Webin format check; it means ready to hand to the submission pipeline, not submitted.';
    END IF;
END
$do$;

COMMIT;
