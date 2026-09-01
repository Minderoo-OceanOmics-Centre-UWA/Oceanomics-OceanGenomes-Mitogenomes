-- Recomputes submission_ready for rows migration 005 wrongly cleared when it
-- was replayed against a restored database.
--
-- Background. Migration 005 adds local_package_status / webin_production_status
-- with DEFAULT 'NOT_RUN', then clears submission_ready wherever those columns
-- are not 'PASS'. When 005 first ran against a live database, those columns
-- held real historical state, so the cleanup only touched rows that actually
-- lacked package validation. Replayed today from an empty migration ledger
-- against a restored ena_validation_attempts table, the columns came back
-- with nothing but the fresh 'NOT_RUN' default -- so the cleanup matched
-- every row that had submission_ready = TRUE and reset all of them. Migration
-- 011 later drops those columns for good, and nothing downstream ever
-- recomputed the flag, so the reset was silent and permanent.
--
-- Fix. Recompute submission_ready from each row's own current gate columns,
-- using the same rule bin/collate_ena_validation.py already applies when it
-- first sets the flag:
--
--   ready = webin_status = 'PASS' AND (
--       (validation_mode = 'validate' AND preflight_status = 'PASS')
--       OR (validation_mode <> 'validate' AND table2asn_status = 'PASS'
--           AND conversion_status = 'PASS')
--   )
--
-- This is a blanket recompute, not a date-scoped patch: a row that is
-- correctly false today (a failed gate, or the FAIL_INFRASTRUCTURE webin runs
-- since the restore) still fails this condition and stays false. Only rows
-- 005's replay wrongly reset come back to true. Scoping the UPDATE to
-- WHERE NOT submission_ready makes a re-run a no-op.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/019_ena_validation_attempts_recompute_submission_ready.sql

BEGIN;

UPDATE ena_validation_attempts
SET submission_ready = TRUE
WHERE NOT submission_ready
  AND webin_status = 'PASS'
  AND (
      (validation_mode = 'validate' AND preflight_status = 'PASS')
      OR (validation_mode <> 'validate'
          AND table2asn_status = 'PASS'
          AND conversion_status = 'PASS')
  );

COMMIT;

-- Sanity check after running:
--
--   SELECT submission_ready, count(*) FROM ena_validation_attempts
--   GROUP BY submission_ready;
