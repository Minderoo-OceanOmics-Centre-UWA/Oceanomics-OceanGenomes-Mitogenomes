-- Record the taxonomic rank at which a sample's species ID was validated.
--
-- species_validation.py used to compare the nominal species ID against the BLAST
-- hits as a whole normalised string, which assumed every nominal ID is a
-- binomial. Many are not: informal field labels ('Diaphus sp 1'), tentative
-- determinations ('Squalus notocaudatus?'), family-level IDs ('Ophidiidae') and
-- bare genera ('Serrivomer') all occur, and invertebrate labels carry them more
-- often than fish ones. The comparison now matches at the rank the label actually
-- ASSERTS, using the lineage the LCA already computed, and records which rank
-- carried it.
--
-- validated_rank
--   'species'           the binomial was found in the BLAST species set
--   'genus'             the genus was found in the BLAST genus set AND the LCA
--                       row agreed at genus
--   'genus_downgraded'  as 'genus', but the label was a tentative species
--                       ('cf.', 'aff.', '?'). Recorded separately because the
--                       release DROPPED an uncertain claim rather than asserting
--                       one, which is safe in the direction that matters but is
--                       still a downgrade and should not be invisible.
--   'family'            the LCA row agreed at family
--   'unmatched'         no rank was supported; the sample stays held
--
--   A relaxed gate stays honest only if what was relaxed is recorded, which is the
--   same reasoning as trna_advisory (021) and order_variant (023). A genus-level
--   release and a species-level release must be distinguishable after the fact.
--
--   Note the matching rule is deliberately NOT "the loosest rank that succeeds".
--   A binomial whose genus is supported but whose species is not stays held, since
--   releasing it would convert "we could not confirm the species" into
--   "validated".
--
-- lca_genus
--   For a family-level release, the genus the LCA resolved. In the observed cases
--   the LCA is MORE specific than the label -- an 'Ophidiidae' label against an
--   LCA genus of Lamprogrammus -- so submitting 'Ophidiidae sp.' discards
--   something the pipeline already worked out. Recording it makes these
--   selectable as label-upgrade candidates for a curator, and means a later
--   policy could hold family-level matches for curation with no code change.
--   NULL for every other rank.

BEGIN;

ALTER TABLE lca_validation
    ADD COLUMN IF NOT EXISTS validated_rank TEXT,
    ADD COLUMN IF NOT EXISTS lca_genus TEXT;

COMMIT;
