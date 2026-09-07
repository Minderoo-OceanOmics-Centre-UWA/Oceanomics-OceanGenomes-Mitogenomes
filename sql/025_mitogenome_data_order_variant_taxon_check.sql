-- Record whether a SECOND taxonomy source agreed with the rank a gene-order
-- variant rule matched on.
--
-- The ORDER_VARIANTS table (migration 023) relaxes the annotation QC gate, and it
-- is keyed on the samplesheet taxonomy in meta -- which is resolved by a fuzzy
-- trigram match and is known to mislabel samples. The obvious conservative design
-- is to require the samplesheet taxonomy and the BLAST-derived consensus lineage
-- to concur BEFORE a rule may fire.
--
-- That was rejected as a precondition, deliberately. The samples the table exists
-- to unblock include ones whose samplesheet taxon is unresolved or wrong, so
-- requiring agreement would have shipped the mechanism without releasing anything
-- it was built for. Recording the disagreement instead makes a rule applied to a
-- mislabelled sample queryable after the fact rather than invisible, which is the
-- part that actually matters.
--
-- Values:
--   'agree'               the rule matched at rank R and the LCA agrees at R
--   'disagree:<taxon>'    the LCA resolved R to something else. THE interesting
--                         one: a relaxed gate on a sample whose own BLAST hits do
--                         not support the taxon the rule keyed on. Query it
--                         periodically; it holds nothing.
--   'unresolved'          the LCA has no consensus at R, or the file was absent
--   'no'                  no rule fired, so there is nothing to check
--
-- This column must NEVER become a gate. Turning it into one re-introduces the
-- agreement-as-precondition design above, and would hold correctly-annotated
-- samples whose only problem is a bad database label.

BEGIN;

ALTER TABLE mitogenome_data
    ADD COLUMN IF NOT EXISTS order_variant_taxon_check TEXT;

COMMIT;
