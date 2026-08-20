-- Retire ena_candidate_runs, the last remnant of the in-repo ENA selection layer.
--
-- 010, 011 and 012 moved locus-tag allocation, package validation and selection
-- to Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes.  ena_candidate_runs was
-- created alongside the rest of that layer in sql/004 and was simply overlooked
-- in the sweep.  Run accessions describe the raw-read submissions to
-- PRJEB123419/420/421, which the downstream submitter owns, exactly like the
-- locus tags 010 moved, so this pipeline has no source for them and no reason
-- to keep the table.
--
-- Nothing in this repository has ever written it, so there is nothing to
-- archive.  Run after 012_drop_ena_selection_layer.sql.  Idempotent.

BEGIN;

DROP TABLE IF EXISTS ena_candidate_runs;

COMMIT;
