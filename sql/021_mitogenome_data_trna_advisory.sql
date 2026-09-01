-- Record the tRNAs waived by the annotation QC gate.
--
-- annotation_stats.py now lets a vertebrate mitogenome pass when it carries the
-- whole conserved core (13 protein-coding genes + both rRNAs) in the correct
-- order but is short at most `annotation_trna_tolerance` tRNAs (default 2) -- an
-- EMMA tRNA-model limitation rather than an assembly defect. missing_genes still
-- lists every absent gene; trna_advisory names the tRNAs that were tolerated so a
-- pass driven by the allowance stays queryable and auditable.
--
-- '' / NULL and the literal 'no' both mean "nothing waived" (annotation_stats.py
-- writes 'no'); a non-empty ';'-joined list (e.g. 'TP' or 'TW;TA;TN') is the set
-- of tolerated tRNAs.

BEGIN;

ALTER TABLE mitogenome_data
    ADD COLUMN IF NOT EXISTS trna_advisory TEXT;

COMMIT;
