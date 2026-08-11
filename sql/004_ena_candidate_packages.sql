-- Candidate-level ENA package registry and specimen-derived locus tags.
--
-- Every scientifically viable assembly/annotation version can have a complete
-- local ENA package.  Submission selection is a later pointer to one of those
-- immutable packages; it is not a prerequisite for package construction.
--
-- Selection is per (specimen, technology): every technology with a viable
-- mitogenome is published to its own ENA child study, and no specimen appears
-- twice within one technology.
--
-- Gene serials are deterministic at specimen level and carry no prefix; the tag
-- is rendered per candidate under the prefix registered to its technology
-- (see 006_ena_tech_aware_locus_tags.sql):
--   OG910 + gene serial 1 -> OGMTHIFI_000910001 in the HiFi record
--
-- Run once manually after 003_mitogenome_data_uniform_depth.sql.
-- Idempotent: all objects use IF NOT EXISTS.

BEGIN;

CREATE TABLE IF NOT EXISTS ena_specimen_accessions (
    og_id TEXT PRIMARY KEY,
    og_numeric INTEGER NOT NULL UNIQUE,
    ena_biosample_accession TEXT UNIQUE,
    accession_source TEXT,
    verified_at TIMESTAMPTZ,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_specimen_og_id_check
        CHECK (og_id ~ '^OG[0-9]+$'),
    CONSTRAINT ena_specimen_og_numeric_check
        CHECK (og_numeric BETWEEN 0 AND 999999),
    CONSTRAINT ena_specimen_og_numeric_matches_check
        CHECK (og_numeric = substring(og_id FROM 3)::INTEGER),
    CONSTRAINT ena_specimen_biosample_check
        CHECK (
            ena_biosample_accession IS NULL
            OR ena_biosample_accession ~ '^SAMEA[0-9]+$'
        )
);

CREATE TABLE IF NOT EXISTS ena_related_assemblies (
    id BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    og_id TEXT NOT NULL REFERENCES ena_specimen_accessions (og_id),
    relationship_type TEXT NOT NULL,
    archive TEXT NOT NULL,
    accession TEXT NOT NULL,
    is_primary BOOLEAN NOT NULL DEFAULT FALSE,
    accession_source TEXT,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_related_relationship_check
        CHECK (relationship_type IN ('reference', 'draft')),
    CONSTRAINT ena_related_archive_check
        CHECK (archive IN ('ENA', 'NCBI', 'OTHER')),
    UNIQUE (og_id, relationship_type, archive, accession)
);

CREATE UNIQUE INDEX IF NOT EXISTS ena_related_one_primary_idx
    ON ena_related_assemblies (og_id, relationship_type)
    WHERE is_primary;

CREATE TABLE IF NOT EXISTS ena_candidate_runs (
    full_seqid TEXT NOT NULL,
    run_accession TEXT NOT NULL,
    platform TEXT NOT NULL,
    contribution_role TEXT NOT NULL DEFAULT 'primary',
    accession_source TEXT,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_candidate_runs_platform_check
        CHECK (platform IN ('PACBIO_SMRT', 'ILLUMINA')),
    CONSTRAINT ena_candidate_runs_role_check
        CHECK (contribution_role IN ('primary', 'polishing', 'supporting')),
    PRIMARY KEY (full_seqid, run_accession)
);

CREATE TABLE IF NOT EXISTS ena_candidate_packages (
    full_seqid TEXT PRIMARY KEY,
    og_id TEXT NOT NULL REFERENCES ena_specimen_accessions (og_id),
    assembly_prefix TEXT NOT NULL,
    annotation_version TEXT NOT NULL,
    ena_study_accession TEXT NOT NULL,
    package_path TEXT NOT NULL,
    package_digest CHAR(64),
    sequence_sha256 CHAR(64) NOT NULL,
    normalised_circular_sha256 CHAR(64) NOT NULL,
    package_status TEXT NOT NULL,
    local_validation_status TEXT NOT NULL,
    webin_test_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    webin_production_status TEXT NOT NULL DEFAULT 'NOT_RUN',
    mean_depth DOUBLE PRECISION,
    assembly_program TEXT,
    platform TEXT,
    biosample_accession TEXT,
    pipeline_revision TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_candidate_full_seqid_check
        CHECK (full_seqid ~ '^OG[0-9]+[.][A-Za-z0-9._-]+$'),
    CONSTRAINT ena_candidate_package_digest_check
        CHECK (package_digest IS NULL OR package_digest ~ '^[0-9a-f]{64}$'),
    CONSTRAINT ena_candidate_sequence_digest_check
        CHECK (sequence_sha256 ~ '^[0-9a-f]{64}$'),
    CONSTRAINT ena_candidate_circular_digest_check
        CHECK (normalised_circular_sha256 ~ '^[0-9a-f]{64}$'),
    CONSTRAINT ena_candidate_package_status_check
        CHECK (
            package_status IN (
                'READY',
                'WAITING_FOR_BIOSAMPLE',
                'BLOCKED_METADATA',
                'BLOCKED_LOCUS_REVIEW',
                'PACKAGE_BLOCKED'
            )
        ),
    CONSTRAINT ena_candidate_local_status_check
        CHECK (local_validation_status IN ('PASS', 'FAIL', 'NOT_RUN')),
    CONSTRAINT ena_candidate_mean_depth_check
        CHECK (mean_depth IS NULL OR mean_depth >= 0)
);

CREATE INDEX IF NOT EXISTS ena_candidate_packages_og_idx
    ON ena_candidate_packages (og_id, created_at DESC);

CREATE TABLE IF NOT EXISTS ena_locus_registry (
    og_id TEXT NOT NULL REFERENCES ena_specimen_accessions (og_id),
    gene_serial INTEGER NOT NULL,
    locus_tag TEXT NOT NULL UNIQUE,
    canonical_gene TEXT NOT NULL,
    gene_occurrence INTEGER NOT NULL DEFAULT 1,
    feature_type TEXT NOT NULL,
    strand CHAR(1),
    feature_sequence_sha256 CHAR(64),
    coordinate_snapshot TEXT,
    allocation_status TEXT NOT NULL DEFAULT 'ACTIVE',
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_locus_gene_serial_check
        CHECK (gene_serial BETWEEN 1 AND 999),
    CONSTRAINT ena_locus_tag_check
        CHECK (locus_tag ~ '^OGMT_[0-9]{9}$'),
    CONSTRAINT ena_locus_strand_check
        CHECK (strand IS NULL OR strand IN ('+', '-')),
    CONSTRAINT ena_locus_status_check
        CHECK (allocation_status IN ('ACTIVE', 'RETIRED', 'REVIEW_REQUIRED')),
    PRIMARY KEY (og_id, gene_serial),
    UNIQUE (og_id, canonical_gene, gene_occurrence)
);

CREATE TABLE IF NOT EXISTS ena_candidate_loci (
    full_seqid TEXT NOT NULL REFERENCES ena_candidate_packages (full_seqid),
    og_id TEXT NOT NULL,
    gene_serial INTEGER NOT NULL,
    locus_tag TEXT NOT NULL,
    feature_key TEXT NOT NULL,
    start_coordinate INTEGER,
    end_coordinate INTEGER,
    strand CHAR(1),
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (og_id, gene_serial)
        REFERENCES ena_locus_registry (og_id, gene_serial),
    FOREIGN KEY (locus_tag)
        REFERENCES ena_locus_registry (locus_tag),
    PRIMARY KEY (full_seqid, feature_key)
);

CREATE TABLE IF NOT EXISTS ena_submission_selections (
    ena_study_accession TEXT NOT NULL,
    og_id TEXT NOT NULL REFERENCES ena_specimen_accessions (og_id),
    biosample_accession TEXT,
    selected_full_seqid TEXT REFERENCES ena_candidate_packages (full_seqid),
    selection_status TEXT NOT NULL,
    selection_reason TEXT NOT NULL,
    selection_report_digest CHAR(64),
    selected_by TEXT NOT NULL,
    selected_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    archive_status TEXT NOT NULL DEFAULT 'NOT_SUBMITTED',
    ena_analysis_accession TEXT,
    ena_assembly_accession TEXT,
    submitted_at TIMESTAMPTZ,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ena_submission_status_check
        CHECK (
            selection_status IN (
                'SELECTED',
                'MANUAL_REVIEW_REQUIRED',
                'NO_PASSING_CANDIDATE'
            )
        ),
    CONSTRAINT ena_submission_presence_check
        CHECK (
            (selection_status = 'SELECTED' AND selected_full_seqid IS NOT NULL)
            OR (selection_status <> 'SELECTED')
        ),
    CONSTRAINT ena_submission_archive_status_check
        CHECK (
            archive_status IN (
                'NOT_SUBMITTED',
                'SUBMITTED',
                'ACCESSION_ASSIGNED'
            )
        ),
    PRIMARY KEY (ena_study_accession, og_id)
);

CREATE UNIQUE INDEX IF NOT EXISTS ena_submission_one_per_biosample_idx
    ON ena_submission_selections (ena_study_accession, biosample_accession)
    WHERE biosample_accession IS NOT NULL
      AND selection_status = 'SELECTED';

COMMIT;

-- Collision audit to run before loading specimen rows:
-- SELECT substring(og_id FROM 3)::INTEGER AS og_numeric, array_agg(og_id)
-- FROM sample
-- WHERE og_id ~ '^OG[0-9]+$'
-- GROUP BY 1 HAVING count(*) > 1;
