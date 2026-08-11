-- Record any INSDC BioSample accession, and name the ENA-registration gap.
--
-- BioSample is a single INSDC namespace mirrored across the three archives, and
-- the accession prefix records only which archive minted it:
--   SAMEA -> EBI/ENA
--   SAMN  -> NCBI
--   SAMD  -> DDBJ
--
-- OceanOmics registers its specimens at NCBI, so every accession the database
-- holds is SAMN: sample.ncbi_biosample_id, draft_genomes.biosample_accession and
-- ref_genomes_assembly_uploads.biosample are SAMN-only, with zero SAMEA between
-- them.  The original SAMEA-only CHECK therefore rejected the entire catalogue,
-- which is why ena_specimen_accessions held 2848 rows and one accession.
--
-- A mirrored SAMN is still not submittable.  webin-cli resolves the manifest
-- SAMPLE against ENA's own submission sample service, which knows only samples
-- registered through Webin -- not the EBI BioSamples mirror the ENA browser
-- serves.  Verified against ena-webin-cli 9.0.3, -context genome -validate -test:
--   SAMEA132129018 -> "Submission(s) validated successfully."
--   SAMN40589646   -> "Failed to initialise validator ... sample is null"
--   SAMEA40589646  -> no such accession (the archives mint from independent
--                     ranges, so the prefix is not a relabelling of one number)
--
-- So the registry must be able to store the SAMN we hold, while the candidate
-- packages must be able to say that a specimen is registered at NCBI and needs
-- an ENA-referenceable accession -- distinct both from having no BioSample at
-- all and from holding a malformed value.

BEGIN;

ALTER TABLE ena_specimen_accessions
    DROP CONSTRAINT IF EXISTS ena_specimen_biosample_check;

ALTER TABLE ena_specimen_accessions
    ADD CONSTRAINT ena_specimen_biosample_check
    CHECK (
        ena_biosample_accession IS NULL
        OR ena_biosample_accession ~ '^SAM(EA|N|D)[0-9]+$'
    );

ALTER TABLE ena_candidate_packages
    DROP CONSTRAINT IF EXISTS ena_candidate_package_status_check;

ALTER TABLE ena_candidate_packages
    ADD CONSTRAINT ena_candidate_package_status_check
    CHECK (
        package_status IN (
            'READY',
            'WAITING_FOR_BIOSAMPLE',
            'BLOCKED_NCBI_ONLY_BIOSAMPLE',
            'BLOCKED_METADATA',
            'BLOCKED_LOCUS_REVIEW',
            'PACKAGE_BLOCKED'
        )
    );

COMMIT;
