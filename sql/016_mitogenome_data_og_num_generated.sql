-- Restores og_num on mitogenome_data as a database-generated column.
--
-- Background. og_num is the numeric part of og_id (OG703 -> 703). Every other
-- table in the schema has it as a stored generated column and every one of them
-- is 100% populated:
--
--   sample           GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::integer)
--   draft_genomes    GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::integer)
--   lca              GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::integer)
--   lca_raw_results  GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::integer)
--   sequencing       GENERATED ALWAYS AS ((substring(sequencing_id, '^OG([0-9]+)'))::integer)
--
-- mitogenome_data was the sole outlier: a plain nullable integer with no
-- default, no generation expression and no trigger, populated on 28 of 289 rows.
-- The snapshot table mitogenome_data_SS260818 still carries the generated column
-- (2181/2181 populated), so the live table almost certainly lost it when it was
-- recreated alongside the rename. Neither push_mtdna_assm_results.py nor
-- push_emma_annotation_results.py ever wrote the column, so every row inserted
-- since then arrived NULL.
--
-- PostgreSQL 14 cannot convert a plain column into a generated one in place, so
-- the column is dropped and re-added. That is safe here: og_num is not in the
-- primary key (og_id, tech, seq_date, code), not in the unique constraint, not
-- indexed, not part of any foreign key, and no view depends on it
-- (mitogenome_submission_view reads mitogenome_data_SS260818, not this table).
-- The 28 pre-existing values all equal SUBSTRING(og_id FROM 3), so the rewrite
-- reproduces them exactly. og_num moves to the last column position, which is
-- where it already sits in sample and draft_genomes.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/016_mitogenome_data_og_num_generated.sql
--
-- Idempotent: the guard makes a re-run a no-op, so the column is never
-- needlessly rewritten.

BEGIN;

DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'public'
          AND table_name   = 'mitogenome_data'
          AND column_name  = 'og_num'
          AND is_generated = 'ALWAYS'
    ) THEN
        ALTER TABLE public.mitogenome_data DROP COLUMN IF EXISTS og_num;
        ALTER TABLE public.mitogenome_data
            ADD COLUMN og_num INTEGER
            GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::INTEGER) STORED;
    END IF;
END $$;

COMMENT ON COLUMN mitogenome_data.og_num IS
    'Numeric part of og_id, maintained by the database. Generated, not writable: '
    'do not include it in any INSERT or UPDATE column list.';

COMMIT;

-- Sanity check after running:
--
--   SELECT count(*) AS rows,
--          count(og_num) AS populated,
--          count(*) FILTER (WHERE og_num <> SUBSTRING(og_id FROM 3)::int) AS mismatched
--   FROM mitogenome_data;
