-- Additive columns for the final annotation-integrity verdict, the lineage gene
-- expectations, the MITOS->EMMA conversion audit, graded gene order, and the two
-- non-vertebrate genes the pipeline now carries end to end (ATP9, mtMutS).
--
-- Purely additive: no existing column is altered or dropped, and every new column
-- is nullable, so an older annotation CSV that lacks them still uploads.
--
-- Freshness: annotation_stats.py always writes an explicit value ('ok'/'none')
-- into the status and issue columns on a successful run, never a blank, so NULL
-- here means "produced before this migration" rather than "clean". The uploader
-- refreshes all of them together whenever a complete row arrives, so a rerun that
-- fixes an annotation cannot leave a stale issue string behind.

ALTER TABLE mitogenome_data
    -- Graded gene order. order_correct (yes/no/NA) and order_deviation (migration
    -- 023, an integer count of genes outside the longest common subsequence) are
    -- both kept unchanged. order_status says WHICH KIND of deviation it was, so a
    -- single displaced tRNA is distinguishable from a real block rearrangement,
    -- and order_deviation_detail NAMES the displaced elements. Deliberately not
    -- reusing order_deviation: that column is an integer and this is a list.
    ADD COLUMN IF NOT EXISTS order_status                   TEXT,
    ADD COLUMN IF NOT EXISTS order_deviation_detail         TEXT,

    -- Genes expected for this lineage but outside the vertebrate set.
    ADD COLUMN IF NOT EXISTS expected_lineage_genes         TEXT,
    ADD COLUMN IF NOT EXISTS missing_expected_lineage_genes TEXT,
    ADD COLUMN IF NOT EXISTS lineage_gene_advisories        TEXT,

    -- What the MITOS -> EMMA conversion could not map, and what it kept but
    -- deliberately excluded from completeness scoring.
    ADD COLUMN IF NOT EXISTS mitos_unmapped_features        TEXT,
    ADD COLUMN IF NOT EXISTS mitos_known_auxiliary_features TEXT,
    ADD COLUMN IF NOT EXISTS duplicate_loci                 TEXT,

    -- Final integrity of the PUBLISHED annotation (post-repair), covering every
    -- exported CDS rather than a fixed 13-gene list.
    ADD COLUMN IF NOT EXISTS annotation_integrity_status    TEXT,
    ADD COLUMN IF NOT EXISTS annotation_integrity_issues    TEXT,

    -- Non-vertebrate genes: nucleotide span and translated length.
    ADD COLUMN IF NOT EXISTS atp9                           INTEGER,
    ADD COLUMN IF NOT EXISTS atp9_trans                     INTEGER,
    ADD COLUMN IF NOT EXISTS mtmuts                         INTEGER,
    ADD COLUMN IF NOT EXISTS mtmuts_trans                   INTEGER,

    -- Provenance, so "this row was not refreshed by the rerun" is detectable
    -- rather than invisible.
    ADD COLUMN IF NOT EXISTS annotation_stats_version       TEXT,
    ADD COLUMN IF NOT EXISTS annotation_updated_at          TIMESTAMPTZ;

COMMENT ON COLUMN mitogenome_data.annotation_integrity_status IS
    'ok | advisory | broken | not_evaluated. Computed on the published annotation after all repair branches merge.';
COMMENT ON COLUMN mitogenome_data.order_status IS
    'reference | trna_displacement | rearranged_block | rearranged_major | not_evaluated.';
COMMENT ON COLUMN mitogenome_data.mtmuts IS
    'Octocoral mitochondrial mismatch-repair gene (INSDC /gene=mtMutS). Expected in Malacalcyonacea/Scleralcyonacea, absent in Hexacorallia.';
