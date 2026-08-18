-- One read-only surface for the ENA submitter.
--
-- Submission is a separate pipeline in a separate repository
-- (Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes), because specimens sit
-- under embargo for an arbitrary period after assembly and cannot be pushed to
-- ENA when the mitogenome happens to finish.  That pipeline needs to know which
-- package to submit and where it is; this view is the whole of what it needs, so
-- it never has to reproduce selection logic or guess at table joins.
--
-- Rows appear only once a package is genuinely submittable:
--   selection_status       = SELECTED   (not MANUAL_REVIEW_REQUIRED)
--   archive_status         = NOT_SUBMITTED
--   package_status         = READY
--   local_validation_status = PASS   (dropped again in 011; see that file)
--
-- A row therefore means: this exact package, at this exact path, is the one
-- candidate chosen for this specimen in this technology, and nobody has
-- submitted it yet.  Setting archive_status to SUBMITTED removes the row and
-- simultaneously engages the freeze rule in select_ena_submission.py, so
-- re-selection can no longer move a specimen out from under a live submission.
--
-- embargo_status is exposed but deliberately NOT filtered on.  Embargo is the
-- submitter's gate to apply, and keeping it there keeps one owner for that
-- decision.  It comes from sample, not mitogenome_data: embargo attaches to the
-- specimen, not to an individual assembly.
--
-- Run after 007_insdc_biosample_accessions.sql.
-- Idempotent: the view is dropped and recreated.

BEGIN;

DROP VIEW IF EXISTS ena_submission_queue;

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
    p.local_validation_status,
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
  AND p.local_validation_status = 'PASS';

COMMENT ON VIEW ena_submission_queue IS
    'Packages cleared for ENA submission but not yet submitted. Consumed by the '
    'ENA-mito-genomes submitter; see docs/ena_submission_handoff.md. Embargo is '
    'exposed, not applied.';

COMMIT;
