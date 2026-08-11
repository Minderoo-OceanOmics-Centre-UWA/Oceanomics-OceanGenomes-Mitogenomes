-- Split Webin test/production results and distinguish package readiness from
-- submission readiness.  Production validation is not submission:
-- validation rows remain overwritable until the selected record is explicitly
-- marked SUBMITTED or ACCESSION_ASSIGNED.

BEGIN;

ALTER TABLE ena_validation_attempts
    ADD COLUMN IF NOT EXISTS package_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    ADD COLUMN IF NOT EXISTS local_package_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    ADD COLUMN IF NOT EXISTS webin_test_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    ADD COLUMN IF NOT EXISTS webin_production_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    ADD COLUMN IF NOT EXISTS webin_error_count INTEGER,
    ADD COLUMN IF NOT EXISTS webin_warning_count INTEGER,
    ADD COLUMN IF NOT EXISTS webin_error_codes TEXT,
    ADD COLUMN IF NOT EXISTS webin_cli_version TEXT,
    ADD COLUMN IF NOT EXISTS webin_test_report_path TEXT,
    ADD COLUMN IF NOT EXISTS webin_production_report_path TEXT,
    ADD COLUMN IF NOT EXISTS overall_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    ADD COLUMN IF NOT EXISTS package_digest CHAR(64);

ALTER TABLE ena_validation_attempts
    DROP CONSTRAINT IF EXISTS ena_validation_submission_ready_check;

-- Rows created by the retired sequence-context path have never passed the new
-- genome-package and production-Webin gates. Preserve their detailed legacy
-- statuses, but clear the derived readiness boolean before enforcing the new
-- definition.
UPDATE ena_validation_attempts
SET submission_ready = FALSE
WHERE submission_ready
  AND (
      local_package_status <> 'PASS'
      OR webin_production_status <> 'PASS'
  );

ALTER TABLE ena_validation_attempts
    ADD CONSTRAINT ena_validation_submission_ready_check
        CHECK (
            NOT submission_ready
            OR (
                local_package_status = 'PASS'
                AND webin_production_status = 'PASS'
            )
        );

ALTER TABLE ena_validation_attempts
    DROP CONSTRAINT IF EXISTS ena_validation_package_digest_check;

ALTER TABLE ena_validation_attempts
    ADD CONSTRAINT ena_validation_package_digest_check
        CHECK (package_digest IS NULL OR package_digest ~ '^[0-9a-f]{64}$');

COMMENT ON COLUMN ena_validation_attempts.table2asn_status IS
    'Supplementary NCBI table2asn validation gate used before ENA package conversion.';

COMMENT ON COLUMN ena_validation_attempts.submission_ready IS
    'True only for the package currently selected for this technology study, after local and production Webin validation; it does not mean submitted.';

COMMIT;
