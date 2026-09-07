-- Widen the three LCA columns that were too narrow for real BLAST output.
--
-- Background. Three column types silently cost 87 assemblies their LCA rows
-- across batch-12 .. batch-20:
--
--   blast_filtered_lca.taxon_id (integer)
--     BLAST reports staxids, which is a ';'-joined list whenever an accession
--     is registered under more than one taxon. The invert/coral database
--     returns these routinely: 533 of 16,079 rows in batch-20 and 427 of
--     25,435 in batch-19 carry a ';'. Every one of them raised
--     'invalid input syntax for type integer' and, because the push script had
--     no per-row savepoint, took the whole sample's upload down with it. Held
--     as text the value survives intact, matching lca_raw_results.taxon_id,
--     which is already character varying.
--
--   lca.top_confidence_score / lca_raw_results.confidence_score (real)
--     real bottoms out around 1.18e-38. HiFi hits produce confidence values
--     like 5.27e-163, which raised 'is out of range for type real' and lost
--     the same way.
--
-- The USING clauses round-trip the floats through their shortest text form
-- rather than widening the float4 bit pattern: a stored 0.995 stays 0.995
-- instead of becoming 0.99500000476837158.
--
-- The rehash at the end matters. Rows in both tables are identified by
-- content_hash (migration 017), an md5 over to_jsonb of the row, and a column
-- type change rewrites the table without firing the BEFORE UPDATE trigger. Any
-- row whose float now renders differently would otherwise keep a stale hash
-- and be re-inserted as a duplicate the next time the pipeline reproduced it.
-- The no-op assignment makes the trigger recompute every hash under the new
-- types. If two rows collapse onto one hash the unique constraint aborts the
-- transaction, which is the outcome we want: it means the rows only ever
-- differed by float noise.
--
-- Nothing in the repo reads blast_filtered_lca.taxon_id -- it is written by
-- bin/push_lca_blast_results.py alone, filtered_lca_view does not reference it,
-- and blast_filtered_lca_new_pk does not include it.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/022_lca_widen_taxon_id_and_confidence.sql
--
-- Idempotent: safe to re-run.

BEGIN;

DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'public' AND table_name = 'blast_filtered_lca'
          AND column_name = 'taxon_id' AND data_type <> 'text'
    ) THEN
        ALTER TABLE public.blast_filtered_lca
            ALTER COLUMN taxon_id TYPE text USING taxon_id::text;
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'public' AND table_name = 'lca'
          AND column_name = 'top_confidence_score' AND data_type <> 'double precision'
    ) THEN
        ALTER TABLE public.lca
            ALTER COLUMN top_confidence_score TYPE double precision
            USING top_confidence_score::text::double precision;

        UPDATE public.lca SET content_hash = content_hash;
    END IF;

    IF EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_schema = 'public' AND table_name = 'lca_raw_results'
          AND column_name = 'confidence_score' AND data_type <> 'double precision'
    ) THEN
        ALTER TABLE public.lca_raw_results
            ALTER COLUMN confidence_score TYPE double precision
            USING confidence_score::text::double precision;

        UPDATE public.lca_raw_results SET content_hash = content_hash;
    END IF;
END $$;

COMMENT ON COLUMN blast_filtered_lca.taxon_id IS
    'BLAST staxids for the matched accession, verbatim. Usually a single NCBI '
    'taxon id, but a '';''-joined list when the accession is registered under '
    'more than one taxon -- cast it before using it as a number.';

COMMENT ON COLUMN lca.top_confidence_score IS
    'double precision, not real: LCA confidence values reach ~1e-163, far below '
    'the ~1.18e-38 floor of real.';

COMMENT ON COLUMN lca_raw_results.confidence_score IS
    'double precision, not real: LCA confidence values reach ~1e-163, far below '
    'the ~1.18e-38 floor of real.';

COMMIT;

-- Sanity check after running:
--
--   SELECT table_name, column_name, data_type
--   FROM information_schema.columns
--   WHERE (table_name, column_name) IN (
--             ('blast_filtered_lca', 'taxon_id'),
--             ('lca', 'top_confidence_score'),
--             ('lca_raw_results', 'confidence_score'));
