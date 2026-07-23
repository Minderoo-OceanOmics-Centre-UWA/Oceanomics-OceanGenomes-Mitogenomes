-- Collapses ena_validation_attempts to one row per
-- (assembly_prefix, ena_study, validation_attempt) so pipeline reruns overwrite
-- the previous attempt instead of appending a new history row every time.
-- push_ena_validation_results.py enforces the other half: once a row reaches
-- submission_ready = true it is frozen and later reruns under the same
-- validation_attempt label no longer overwrite it.
--
-- Run once, manually, the same way as 001_create_ena_validation_attempts.sql:
--   psql --dbname oceanomics --file sql/002_ena_validation_attempts_single_row_per_attempt.sql

BEGIN;

-- Keep exactly one row per key: prefer a submission_ready row over a failed
-- one, then the most recently recorded attempt. Deletes every row that has
-- some other row in its group with a strictly "greater" (submission_ready,
-- recorded_at, id) tuple, which always leaves exactly one survivor per group.
DELETE FROM ena_validation_attempts a
USING ena_validation_attempts b
WHERE a.assembly_prefix = b.assembly_prefix
  AND a.ena_study = b.ena_study
  AND a.validation_attempt = b.validation_attempt
  AND (b.submission_ready, b.recorded_at, b.id) > (a.submission_ready, a.recorded_at, a.id);

DROP INDEX IF EXISTS ena_validation_exact_result_idx;

ALTER TABLE ena_validation_attempts
    ADD COLUMN IF NOT EXISTS attempt_count INTEGER NOT NULL DEFAULT 1;

CREATE UNIQUE INDEX IF NOT EXISTS ena_validation_attempts_key_idx
    ON ena_validation_attempts (assembly_prefix, ena_study, validation_attempt);

COMMIT;
