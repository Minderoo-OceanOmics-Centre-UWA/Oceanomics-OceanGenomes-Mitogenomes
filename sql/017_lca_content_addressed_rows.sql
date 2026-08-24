-- Records the content-addressing rework of lca and lca_raw_results.
--
-- Background. This change was applied directly to the live database around
-- 2026-08-20 and never written down. The code half of it shipped in 55321df,
-- the same commit that added migrations 014 and 015, so the ENA side was
-- captured and this side was not. A database rebuilt from sql/ alone would fail
-- on the first LCA push with
--
--   constraint "lca_content_unique" for table "lca" does not exist
--
-- because bin/push_lca_blast_results.py names that constraint as its ON CONFLICT
-- target and bin/push_lca_raw_results.py names lca_raw_results_content_unique.
--
-- What it does. Both tables used to be keyed on lca_run_date, so every re-run
-- appended a near-duplicate row whether or not the result had changed. They are
-- now keyed on a content_hash maintained by a BEFORE INSERT OR UPDATE trigger,
-- computed over the whole row minus lca_run_date, og_num and content_hash
-- itself. A re-run that reproduces a hit exactly therefore conflicts and just
-- refreshes that row's lca_run_date; a re-run that produces a different result
-- for the same hit inserts a new row and the earlier one is kept. taxon_rank_db
-- lands alongside it: the rank reported by the reference database, which the
-- BLAST parser now carries through.
--
-- The FK constraints are also renamed. Both tables called theirs fk_mitogenome,
-- which is ambiguous once more than one table points at mitogenome_data.
--
-- Every step is guarded, so this is a no-op against the live database (where all
-- of it already exists) and correct against one rebuilt from sql/. Order
-- matters: the trigger has to exist before content_hash can be backfilled, and
-- the backfill has to finish before the column can take NOT NULL.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/017_lca_content_addressed_rows.sql
--
-- Idempotent: safe to re-run.

BEGIN;

-- 1. The hash function. og_num is excluded because it is generated from og_id
--    and so adds nothing; lca_run_date is excluded because "when we last saw
--    this" must not change the identity of the row.
CREATE OR REPLACE FUNCTION public.lca_set_content_hash()
RETURNS trigger
LANGUAGE plpgsql
AS $function$
BEGIN
    NEW.content_hash := md5(
        (to_jsonb(NEW) - 'lca_run_date' - 'og_num' - 'content_hash')::text
    );
    RETURN NEW;
END;
$function$;

-- 2/3. The new columns. content_hash starts nullable: a non-empty table cannot
--      take NOT NULL until step 5 has given every row a value.
ALTER TABLE public.lca
    ADD COLUMN IF NOT EXISTS taxon_rank_db TEXT,
    ADD COLUMN IF NOT EXISTS content_hash  TEXT;

ALTER TABLE public.lca_raw_results
    ADD COLUMN IF NOT EXISTS taxon_rank_db TEXT,
    ADD COLUMN IF NOT EXISTS content_hash  TEXT;

-- 4. The triggers.
DROP TRIGGER IF EXISTS lca_content_hash ON public.lca;
CREATE TRIGGER lca_content_hash
    BEFORE INSERT OR UPDATE ON public.lca
    FOR EACH ROW EXECUTE FUNCTION public.lca_set_content_hash();

DROP TRIGGER IF EXISTS lca_raw_results_content_hash ON public.lca_raw_results;
CREATE TRIGGER lca_raw_results_content_hash
    BEFORE INSERT OR UPDATE ON public.lca_raw_results
    FOR EACH ROW EXECUTE FUNCTION public.lca_set_content_hash();

-- 5. Backfill. The BEFORE UPDATE trigger computes the value; the assignment
--    itself is a no-op. The IS NULL guard means no row that already has a hash
--    is ever rewritten, so on the live database this matches nothing.
UPDATE public.lca            SET content_hash = content_hash WHERE content_hash IS NULL;
UPDATE public.lca_raw_results SET content_hash = content_hash WHERE content_hash IS NULL;

-- 6. Now the column can be required.
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema='public' AND table_name='lca'
          AND column_name='content_hash' AND is_nullable='YES'
    ) THEN
        ALTER TABLE public.lca ALTER COLUMN content_hash SET NOT NULL;
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema='public' AND table_name='lca_raw_results'
          AND column_name='content_hash' AND is_nullable='YES'
    ) THEN
        ALTER TABLE public.lca_raw_results ALTER COLUMN content_hash SET NOT NULL;
    END IF;
END $$;

-- 7. Swap the unique key from lca_run_date to content_hash. These are the
--    constraint names the push scripts name as ON CONFLICT targets, so they are
--    part of the pipeline's contract and cannot be renamed casually.
DO $$
BEGIN
    ALTER TABLE public.lca DROP CONSTRAINT IF EXISTS lca_results_unique;
    IF NOT EXISTS (
        SELECT 1 FROM pg_constraint
        WHERE conrelid='public.lca'::regclass AND conname='lca_content_unique'
    ) THEN
        ALTER TABLE public.lca ADD CONSTRAINT lca_content_unique
            UNIQUE (og_id, tech, seq_date, code, annotation, region, content_hash);
    END IF;

    ALTER TABLE public.lca_raw_results DROP CONSTRAINT IF EXISTS lca_raw_results_unique;
    IF NOT EXISTS (
        SELECT 1 FROM pg_constraint
        WHERE conrelid='public.lca_raw_results'::regclass
          AND conname='lca_raw_results_content_unique'
    ) THEN
        ALTER TABLE public.lca_raw_results ADD CONSTRAINT lca_raw_results_content_unique
            UNIQUE (og_id, tech, seq_date, code, annotation, sequence_region,
                    accession_id, content_hash);
    END IF;
END $$;

-- 8. Disambiguate the foreign key names.
DO $$
BEGIN
    IF EXISTS (SELECT 1 FROM pg_constraint
               WHERE conrelid='public.lca'::regclass AND conname='fk_mitogenome')
       AND NOT EXISTS (SELECT 1 FROM pg_constraint
               WHERE conrelid='public.lca'::regclass AND conname='fk_mitogenome_lca') THEN
        ALTER TABLE public.lca RENAME CONSTRAINT fk_mitogenome TO fk_mitogenome_lca;
    END IF;

    IF EXISTS (SELECT 1 FROM pg_constraint
               WHERE conrelid='public.lca_raw_results'::regclass AND conname='fk_mitogenome')
       AND NOT EXISTS (SELECT 1 FROM pg_constraint
               WHERE conrelid='public.lca_raw_results'::regclass
                 AND conname='fk_mitogenome_lca_raw_results') THEN
        ALTER TABLE public.lca_raw_results
            RENAME CONSTRAINT fk_mitogenome TO fk_mitogenome_lca_raw_results;
    END IF;
END $$;

-- 9. Say so on the columns themselves.
COMMENT ON COLUMN lca.content_hash IS
    'md5 of the row minus lca_run_date, og_num and content_hash, maintained by the '
    'lca_content_hash trigger. Do not include it in any INSERT or UPDATE column list.';

COMMENT ON COLUMN lca_raw_results.content_hash IS
    'md5 of the row minus lca_run_date, og_num and content_hash, maintained by the '
    'lca_raw_results_content_hash trigger. Do not include it in any INSERT or UPDATE '
    'column list.';

COMMENT ON COLUMN lca.taxon_rank_db IS
    'Taxonomic rank as reported by the reference database, as opposed to taxon_rank, '
    'which is the rank the LCA itself resolved to.';

COMMENT ON COLUMN lca_raw_results.taxon_rank_db IS
    'Taxonomic rank as reported by the reference database for this hit.';

COMMIT;

-- Sanity check after running:
--
--   SELECT 'lca' AS t, count(*) AS rows, count(content_hash) AS hashed,
--          count(DISTINCT content_hash) AS distinct_hashes FROM lca
--   UNION ALL
--   SELECT 'lca_raw_results', count(*), count(content_hash), count(DISTINCT content_hash)
--   FROM lca_raw_results;
