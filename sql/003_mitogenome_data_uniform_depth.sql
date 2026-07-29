-- Adds the uniform, cross-platform read-depth metric to mitogenome_data.
--
-- Background. avg_coverage held three different quantities depending on which
-- assembler produced the row, so comparing that column across platforms was
-- always wrong:
--
--   GetOrganelle  k-mer coverage off the assembly graph (~0.2x true depth), and
--                 measured on the reduced read set GetOrganelle selects by
--                 default. avg_base_coverage on those rows IS a base depth, but
--                 still only over the reduced set.
--   MitoHiFi      per-base depth of ONLY the reads recruited by mapping to a
--                 related-species reference, so a divergent reference silently
--                 depressed it. Both columns hold that same single value.
--   Oatk          NULL. Neither parser matched its log.
--
-- mean_depth replaces all three with one definition: mean per-base depth of the
-- sample's own reads remapped to the assembly that went to annotation, folded on
-- a doubled reference when the molecule is circular.
--
-- avg_coverage and avg_base_coverage are deliberately NOT altered or recomputed.
-- Historical rows keep exactly the values they were written with; depth_method
-- records which kind of number a given row carries, so the two generations can
-- never be silently mixed in a query.
--
-- Run once, manually. Connection details live in
-- /home/tpeirce/postgresql_details/oceanomics.cfg (the same file the pipeline
-- passes as --sql_config):
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/003_mitogenome_data_uniform_depth.sql
--
-- Idempotent: safe to re-run.

BEGIN;

ALTER TABLE mitogenome_data
    ADD COLUMN IF NOT EXISTS mean_depth             DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS median_depth           DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS depth_sd               DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS depth_cv               DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS breadth_1x             DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS breadth_10x            DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS mito_mapped_reads      BIGINT,
    ADD COLUMN IF NOT EXISTS total_reads            BIGINT,
    ADD COLUMN IF NOT EXISTS mito_read_fraction     DOUBLE PRECISION,
    ADD COLUMN IF NOT EXISTS depth_target_length_bp INTEGER,
    ADD COLUMN IF NOT EXISTS depth_target_fasta     TEXT,
    ADD COLUMN IF NOT EXISTS depth_method           TEXT,
    ADD COLUMN IF NOT EXISTS depth_measured_at      TIMESTAMPTZ;

COMMENT ON COLUMN mitogenome_data.mean_depth IS
    'Mean per-base read depth of the sample''s own reads remapped to this assembly. '
    'Comparable across GetOrganelle / MitoHiFi / Oatk and across Illumina / HiC / HiFi. '
    'Use this, not avg_coverage, for any cross-platform comparison.';

COMMENT ON COLUMN mitogenome_data.depth_method IS
    'remap_full_v1  = mean_depth measured by remapping the full post-QC read set to '
    'this assembly (circular molecules folded on a doubled reference). '
    'not_measured   = this assembly never reached annotation (failed, under-length, '
    'or a discarded assembly variant), or the depth step was skipped. '
    'legacy_*       = pre-dates the uniform measurement; only the assembler-specific '
    'avg_coverage / avg_base_coverage values exist for this row.';

COMMENT ON COLUMN mitogenome_data.avg_coverage IS
    'LEGACY, assembler-specific and NOT comparable across assemblers: k-mer coverage '
    'for GetOrganelle rows, reference-recruited read depth for MitoHiFi rows, NULL for '
    'Oatk. Retained for provenance. See mean_depth.';

-- Label-only backfill. Recomputes nothing and touches no row that already carries
-- a measured depth: both guards (depth_method IS NULL AND mean_depth IS NULL) must
-- hold, so re-running this after real measurements have landed is a no-op.
--
-- `code` is the 4th dot-separated field of the assembly prefix
-- (og_id.tech.seq_date.code), e.g. getorg1770 / v323mitohifi / v10oatk.
UPDATE mitogenome_data SET depth_method = 'legacy_getorg_kmer'
 WHERE depth_method IS NULL AND mean_depth IS NULL AND code LIKE 'getorg%';

UPDATE mitogenome_data SET depth_method = 'legacy_mitohifi_recruited'
 WHERE depth_method IS NULL AND mean_depth IS NULL AND code LIKE '%mitohifi%';

UPDATE mitogenome_data SET depth_method = 'legacy_none'
 WHERE depth_method IS NULL AND mean_depth IS NULL AND code LIKE '%oatk%';

-- Kept as its own bucket rather than folded into legacy_none: legacy_none means
-- "we know Oatk wrote no coverage", whereas this means "the assembler code did not
-- match any known pattern", which is a data-quality signal worth being able to see.
UPDATE mitogenome_data SET depth_method = 'legacy_unknown'
 WHERE depth_method IS NULL AND mean_depth IS NULL;

CREATE INDEX IF NOT EXISTS mitogenome_data_depth_method_idx
    ON mitogenome_data (depth_method);

COMMIT;

-- Sanity check after running:
--
--   SELECT depth_method, count(*), round(avg(mean_depth)::numeric, 1) AS avg_mean_depth
--   FROM mitogenome_data GROUP BY depth_method ORDER BY 2 DESC;
