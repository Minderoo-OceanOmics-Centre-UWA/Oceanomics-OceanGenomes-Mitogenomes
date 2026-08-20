-- Retire the local package validation gate.
--
-- ena_package.py used to re-check its own output (manifest fields, EMBL entry
-- name, chromosome list) and record the verdict as local_validation_status.
-- That check never told us anything package_status did not already say: a
-- package is either READY, or it is blocked with the reason named
-- (WAITING_FOR_BIOSAMPLE, BLOCKED_METADATA, ...).  The pipeline now ends at the
-- flatfile format check (WEBIN_VALIDATE, ena-webin-cli -context sequence), and
-- everything BioSample-dependent belongs to ena_selection.nf.
--
-- ena_candidate_packages.webin_test_status stays: RECORD_ENA_PACKAGE_VALIDATION
-- in the selection workflow still writes it.  The same-named column on
-- ena_validation_attempts goes, because the per-candidate collator no longer
-- runs a genome-context Webin check and can only ever write NOT_RUN there.
--
-- Run after 010_drop_ena_locus_tables.sql.
-- Guarded like 009: the integration tests build partial schemas, so touch
-- ena_validation_attempts and rebuild ena_submission_queue only where they are
-- already present.  Idempotent.

BEGIN;

DO $do$
DECLARE
    had_queue boolean := to_regclass('ena_submission_queue') IS NOT NULL;
BEGIN
    DROP VIEW IF EXISTS ena_submission_queue;

    ALTER TABLE ena_candidate_packages
        DROP CONSTRAINT IF EXISTS ena_candidate_local_status_check;

    ALTER TABLE ena_candidate_packages
        DROP COLUMN IF EXISTS local_validation_status;

    IF to_regclass('ena_validation_attempts') IS NOT NULL THEN
        -- ena_validation_latest is SELECT DISTINCT ON (assembly_prefix) *, and
        -- Postgres expanded that * when the view was created, so it depends on
        -- every column below and has to be dropped before them.
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

        -- submission_ready is now earned by production Webin validation alone;
        -- the constraint has to go before the column it references.
        ALTER TABLE ena_validation_attempts
            DROP CONSTRAINT IF EXISTS ena_validation_submission_ready_check;

        ALTER TABLE ena_validation_attempts
            DROP COLUMN IF EXISTS local_package_status,
            DROP COLUMN IF EXISTS webin_test_status,
            DROP COLUMN IF EXISTS webin_test_report_path;

        ALTER TABLE ena_validation_attempts
            ADD CONSTRAINT ena_validation_submission_ready_check
                CHECK (NOT submission_ready OR webin_production_status = 'PASS');

        -- Rebuilt against the narrowed table.  012 drops and recreates it
        -- again; this migration still has to leave the schema whole on its own.
        -- 014 replaces assembly_prefix with full_seqid and owns the view from
        -- then on, so skip this once the rebuild has happened.
        IF EXISTS (
            SELECT 1 FROM pg_attribute
            WHERE attrelid = to_regclass('ena_validation_attempts')
              AND attname = 'assembly_prefix'
              AND NOT attisdropped
        ) THEN
            EXECUTE $latest$
                CREATE VIEW ena_validation_latest AS
                SELECT DISTINCT ON (assembly_prefix)
                    *
                FROM ena_validation_attempts
                ORDER BY assembly_prefix, recorded_at DESC, id DESC
            $latest$;
        END IF;
    END IF;

    IF had_queue THEN
        -- As sql/008, minus local_validation_status in the projection and the
        -- filter.
        EXECUTE $view$
            CREATE VIEW ena_submission_queue AS
            SELECT
                s.og_id,
                split_part(p.assembly_prefix, '.', 2)   AS tech,
                p.full_seqid,                            -- also the Webin ASSEMBLYNAME
                p.assembly_prefix,
                p.annotation_version,
                s.ena_study_accession,
                p.biosample_accession,                   -- as written on the manifest
                p.package_path,
                p.package_digest,
                p.platform,
                p.assembly_program,
                p.mean_depth,
                p.package_status,
                p.webin_test_status,
                p.webin_production_status,
                s.selection_status,
                s.selection_reason,
                s.selected_at,
                s.archive_status,
                sa.embargo_status
            FROM ena_submission_selections s
            JOIN ena_candidate_packages p
                ON p.full_seqid = s.selected_full_seqid
            LEFT JOIN sample sa
                ON sa.og_id = s.og_id
            WHERE s.selection_status = 'SELECTED'
              AND s.archive_status = 'NOT_SUBMITTED'
              AND p.package_status = 'READY'
        $view$;

        EXECUTE $comment$
            COMMENT ON VIEW ena_submission_queue IS
                'Packages cleared for ENA submission but not yet submitted. Consumed by the '
                'ENA-mito-genomes submitter; see docs/ena_submission_handoff.md. Embargo is '
                'exposed, not applied.'
        $comment$;
    END IF;
END
$do$;

COMMIT;
