-- Record why a non-canonical gene order was accepted, and how far off canonical a
-- rejected one was.
--
-- The annotation QC gate used to decide gene order with a single strict test
-- against one hard-coded vertebrate order, so a real, published, lineage-level
-- rearrangement failed exactly like a scrambled assembly. annotation_stats.py now
-- consults a curated per-taxon table (bin/mito_gene_order.py, ORDER_VARIANTS) and
-- records what it did.
--
-- order_variant
--   The id of the curated rule that accepted this assembly's order, e.g.
--   'scarine_imq'. Set only when order_correct = 'variant'. A relaxed gate stays
--   honest only if what was relaxed is recorded and queryable after the fact --
--   the same reasoning as trna_advisory in migration 021.
--
--   Note a variant rule ADDS an accepted order rather than replacing the canonical
--   one, so a canonical member of a rule-carrying taxon keeps order_correct='yes'
--   and order_variant='no'. A row with a non-'no' value here is an assembly that
--   actually showed the rearrangement.
--
-- order_deviation
--   When the order was rejected, the number of genes falling outside the longest
--   common subsequence of the observed order against the nearest accepted one.
--   order_correct='no' otherwise collapses "one tRNA out of place" and "half the
--   genome is inverted" into the same value, which makes a held pile impossible to
--   triage. It is also the cheapest signal for spotting the NEXT clade variant: a
--   group of samples sharing a low non-zero deviation is what a missing table row
--   looks like.
--
-- annotation_gaps
--   ';'-joined interior intergenic spans at or above params.annotation_gap_threshold
--   (default 50 bp), each rendered 'ND4:11768-TS1:11982(213)'. This is what
--   separates a real transposition (a tRNA moved) from a missed call (a tRNA-shaped
--   hole left at both the origin and the destination of the apparent move), which is
--   otherwise only visible by reading the GFF by hand. Unlike the two columns above
--   it is populated for EVERY completeness profile, including the invertebrate core
--   profile, where gene order is not evaluated at all and gene presence is the only
--   other signal.
--
-- All three are PURELY ADVISORY except order_variant's effect via order_correct:
-- order_deviation and annotation_gaps never influence `passed`.
--
-- NULL and the literal 'no' both mean "nothing to report" in all three columns
-- (annotation_stats.py writes 'no'), matching the convention already used by
-- missing_genes, extra_genes and trna_advisory.
--
-- order_deviation is TEXT rather than INTEGER so it can carry that same 'no'
-- convention rather than overloading NULL or 0 -- 0 would otherwise be ambiguous
-- with a correct order.
--
-- order_correct is already VARCHAR (sql/018_mitogenome_data_og_num_first.sql), so
-- its new 'variant' value needs no schema change, and nothing branches on it:
-- bin/evaluate_qc_conditions.py gates purely on `passed`, so the new value cannot
-- change routing.

BEGIN;

ALTER TABLE mitogenome_data
    ADD COLUMN IF NOT EXISTS order_variant TEXT,
    ADD COLUMN IF NOT EXISTS order_deviation TEXT,
    ADD COLUMN IF NOT EXISTS annotation_gaps TEXT;

COMMIT;
