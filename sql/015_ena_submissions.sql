-- Record what happened to a validated sequence after it left this pipeline.
--
-- 014 made the database able to say which annotation of which assembly cleared
-- every gate.  It still cannot say whether that sequence was submitted, what
-- ENA gave back, or which BioSample went out on the manifest: that lives only
-- as receipt XML under receipts/<OG>/ in the downstream submitter, plus
-- mitogenome_data.genbank_accession, which is keyed without the annotation
-- version and named for an archive that does not mint ERZ accessions.
--
-- This is deliberately NOT columns on ena_validation_attempts.  That table is
-- upserted with DO UPDATE SET across every non-key column
-- (bin/push_ena_validation_results.py), so an accession stored there would be
-- erased by the next validation rerun.  A ledger has to survive revalidation,
-- so it gets its own table.
--
-- Ownership is unchanged: this pipeline builds and format-validates, and never
-- writes or reads ena_submissions.  Validation must not depend on submission
-- state.  The writer is the downstream submitter
-- (Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes); see
-- docs/ena_submission_handoff.md for the contract.
--
-- Run after 014_ena_validation_attempts_full_seqid.sql.  Idempotent.

BEGIN;

CREATE TABLE IF NOT EXISTS ena_submissions (
    full_seqid TEXT NOT NULL,
    -- Test submissions are answered by a separate service with throwaway
    -- accessions.  Keeping the mode in the key means a dry run can never
    -- overwrite the production record of the same sequence.
    webin_mode TEXT NOT NULL DEFAULT 'production',

    -- Identity is derived, never written, so it cannot drift from full_seqid.
    -- annotation is everything after the fourth field, so a dotted version like
    -- emma1.0.2 stays whole.  It is NULL for a bare 4-field id.
    og_id      TEXT GENERATED ALWAYS AS (split_part(full_seqid, '.', 1)) STORED,
    tech       TEXT GENERATED ALWAYS AS (split_part(full_seqid, '.', 2)) STORED,
    seq_date   TEXT GENERATED ALWAYS AS (split_part(full_seqid, '.', 3)) STORED,
    code       TEXT GENERATED ALWAYS AS (split_part(full_seqid, '.', 4)) STORED,
    annotation TEXT GENERATED ALWAYS AS
        (substring(full_seqid from '^(?:[^.]+[.]){4}(.+)$')) STORED,

    submission_status TEXT NOT NULL DEFAULT 'NOT_SUBMITTED',
    submitted_at TIMESTAMPTZ,
    submitted_by TEXT,
    error_message TEXT,
    receipt_path TEXT,
    receipt_sha256 CHAR(64),
    webin_cli_version TEXT,

    ena_study_accession    TEXT,   -- PRJEB..., the per-technology child study
    ena_sample_accession   TEXT,   -- ERS..., only when a Webin sample was registered
    ena_analysis_accession TEXT,   -- ERZ..., returned by webin-cli at submit time
    ena_assembly_accession TEXT,   -- GCA_..., assigned later
    ena_sequence_accession TEXT,   -- OU/LR..., assigned later

    -- What the manifest actually carried, not what the catalogue holds: the two
    -- differ whenever a specimen is registered at NCBI only.
    biosample_accession TEXT,
    biosample_source TEXT,

    locus_tag_prefix TEXT,
    run_accessions TEXT[],

    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,

    PRIMARY KEY (full_seqid, webin_mode),

    CONSTRAINT ena_submissions_full_seqid_check
        CHECK (full_seqid ~ '^OG[0-9]+[.][A-Za-z0-9._-]+$'),
    CONSTRAINT ena_submissions_webin_mode_check
        CHECK (webin_mode IN ('production', 'test')),
    CONSTRAINT ena_submissions_status_check
        CHECK (
            submission_status IN (
                'NOT_SUBMITTED',
                'SUBMITTED',
                'ACCESSION_ASSIGNED',
                'FAILED'
            )
        ),

    -- A claim of having submitted has to carry the evidence for it.
    CONSTRAINT ena_submissions_submitted_at_check
        CHECK (
            submission_status NOT IN ('SUBMITTED', 'ACCESSION_ASSIGNED')
            OR submitted_at IS NOT NULL
        ),
    CONSTRAINT ena_submissions_accession_presence_check
        CHECK (
            submission_status <> 'ACCESSION_ASSIGNED'
            OR ena_analysis_accession IS NOT NULL
        ),
    CONSTRAINT ena_submissions_failure_reason_check
        CHECK (submission_status <> 'FAILED' OR error_message IS NOT NULL),

    CONSTRAINT ena_submissions_study_accession_check
        CHECK (ena_study_accession IS NULL OR ena_study_accession ~ '^PRJEB[0-9]+$'),
    CONSTRAINT ena_submissions_sample_accession_check
        CHECK (ena_sample_accession IS NULL OR ena_sample_accession ~ '^ERS[0-9]+$'),
    CONSTRAINT ena_submissions_analysis_accession_check
        CHECK (ena_analysis_accession IS NULL OR ena_analysis_accession ~ '^ERZ[0-9]+$'),
    CONSTRAINT ena_submissions_assembly_accession_check
        CHECK (ena_assembly_accession IS NULL OR ena_assembly_accession ~ '^GCA_[0-9]{9}[.][0-9]+$'),
    CONSTRAINT ena_submissions_sequence_accession_check
        CHECK (ena_sequence_accession IS NULL OR ena_sequence_accession ~ '^[A-Z]{2}[0-9]{6,8}([.][0-9]+)?$'),

    -- Same namespace rule as ena_specimen_accessions (sql/007): the accession
    -- prefix records only which archive minted it, and the catalogue is SAMN.
    CONSTRAINT ena_submissions_biosample_check
        CHECK (biosample_accession IS NULL OR biosample_accession ~ '^SAM(EA|N|D)[0-9]+$'),

    CONSTRAINT ena_submissions_locus_tag_prefix_check
        CHECK (locus_tag_prefix IS NULL OR locus_tag_prefix ~ '^[A-Z][A-Z0-9]{2,11}$'),
    CONSTRAINT ena_submissions_receipt_digest_check
        CHECK (receipt_sha256 IS NULL OR receipt_sha256 ~ '^[0-9a-f]{64}$'),
    CONSTRAINT ena_submissions_run_accessions_check
        CHECK (
            run_accessions IS NULL
            OR array_to_string(run_accessions, ' ') ~ '^((ERR|SRR|DRR)[0-9]+( |$))*$'
        )
);

CREATE INDEX IF NOT EXISTS ena_submissions_identity_idx
    ON ena_submissions (og_id, tech, seq_date, code, annotation);

CREATE INDEX IF NOT EXISTS ena_submissions_status_idx
    ON ena_submissions (submission_status);

COMMENT ON TABLE ena_submissions IS
    'Submission ledger, written by the downstream ENA submission pipeline. This pipeline '
    'never writes or reads it: validation does not depend on submission state.';
COMMENT ON COLUMN ena_submissions.webin_mode IS
    'Which Webin service answered: production or test. Part of the key so a dry run '
    'cannot overwrite the production record.';
COMMENT ON COLUMN ena_submissions.biosample_accession IS
    'The BioSample the manifest actually carried. SAMN means the specimen is registered '
    'at NCBI only, which webin-cli cannot resolve (see sql/007).';
COMMENT ON COLUMN ena_submissions.biosample_source IS
    'Where that accession came from, e.g. sample.ncbi_biosample_id or webin_sample_receipt.';
COMMENT ON COLUMN ena_submissions.run_accessions IS
    'Raw-read runs the assembly came from, as ERR/SRR/DRR accessions.';

-- The one question both halves answer together: what is validated, and what has
-- happened to it since.  Defined on ena_validation_attempts with its own
-- DISTINCT ON rather than on ena_validation_latest: a view on that view would
-- make DROP VIEW ena_validation_latest fail inside 011, 012 and 014.
CREATE OR REPLACE VIEW ena_submission_status AS
WITH latest AS (
    SELECT DISTINCT ON (full_seqid) *
    FROM ena_validation_attempts
    ORDER BY full_seqid, recorded_at DESC, id DESC
)
SELECT
    v.full_seqid,
    v.og_id,
    v.tech,
    v.seq_date,
    v.code,
    v.annotation,
    v.ena_study,
    v.webin_status,
    v.submission_ready,
    v.recorded_at AS validated_at,
    COALESCE(s.submission_status, 'NOT_SUBMITTED') AS submission_status,
    s.ena_analysis_accession,
    s.ena_assembly_accession,
    s.ena_sequence_accession,
    s.biosample_accession,
    s.locus_tag_prefix,
    s.submitted_at,
    s.receipt_path
FROM latest v
LEFT JOIN ena_submissions s
       ON s.full_seqid = v.full_seqid
      AND s.webin_mode = 'production';

COMMIT;
