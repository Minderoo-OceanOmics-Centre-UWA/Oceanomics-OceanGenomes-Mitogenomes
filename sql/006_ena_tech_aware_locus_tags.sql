-- Make the specimen locus registry prefix-free so one specimen can be published
-- once per sequencing technology.
--
-- Each technology has its own ENA study and therefore its own registered INSDC
-- locus-tag prefix:
--   hifi -> PRJEB123419 / OGMTHIFI
--   hic  -> PRJEB123420 / OGMTHIC
--   ilmn -> PRJEB123421 / OGMTILMN
--
-- A prefix belongs to exactly one study and a locus tag cannot be shared between
-- two published records, so the same specimen gene must render differently in
-- each per-tech record.  The registry therefore stores the specimen-level serial
-- only (OG910 gene 1 is serial 1 in every technology) and the rendered tag
-- OGMTHIFI_000910001 / OGMTHIC_000910001 / OGMTILMN_000910001 is produced per
-- candidate at allocation time from the tech of its full_seqid.
--
-- ena_candidate_loci keeps the resolved tag: that table is per package, so it is
-- the only place the tech-specific tag is recorded.
--
-- Run after 005_ena_validation_attempts_genome_context.sql via
-- bin/apply_ena_migrations.py.  Idempotent.

BEGIN;

-- The rendered-tag foreign key cannot survive a prefix-free registry.  The
-- composite (og_id, gene_serial) foreign key stays and is the real integrity
-- guarantee.  Drop by lookup rather than by generated name so this works on
-- databases where 004 was applied under a different constraint name.
-- Every catalogue lookup below is scoped with to_regclass() rather than by bare
-- relname or conname. Bare names match across every schema, so a migration
-- applied into a throwaway schema (as the integration test does) would find and
-- try to alter the copy in public.
DO $$
DECLARE
    constraint_name TEXT;
BEGIN
    SELECT con.conname INTO constraint_name
    FROM pg_constraint con
    WHERE con.contype = 'f'
      AND con.conrelid = to_regclass('ena_candidate_loci')
      AND con.confrelid = to_regclass('ena_locus_registry')
      AND con.conkey = ARRAY[
          (SELECT attnum FROM pg_attribute
           WHERE attrelid = to_regclass('ena_candidate_loci')
             AND attname = 'locus_tag')
      ]::SMALLINT[];
    IF constraint_name IS NOT NULL THEN
        EXECUTE format(
            'ALTER TABLE ena_candidate_loci DROP CONSTRAINT %I', constraint_name
        );
    END IF;
END
$$;

-- Dropping the column also drops its UNIQUE constraint and the hardcoded
-- '^OGMT_[0-9]{9}$' check.  No CASCADE: if anything else still depends on the
-- rendered tag the migration must fail loudly rather than silently remove it.
ALTER TABLE ena_locus_registry DROP COLUMN IF EXISTS locus_tag;

COMMENT ON TABLE ena_locus_registry IS
    'Specimen-level gene serial allocation. Shared across every technology for a '
    'specimen; carries no prefix. Rendered tags live in ena_candidate_loci.';
COMMENT ON COLUMN ena_locus_registry.gene_serial IS
    'Stable 1-999 serial for this specimen gene. Rendered as '
    '<tech prefix>_<og numeric:06d><gene_serial:03d> per candidate.';

-- ena_candidate_loci keeps rendered tags, now under any registered prefix.
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_constraint
        WHERE conrelid = to_regclass('ena_candidate_loci')
          AND conname = 'ena_candidate_loci_locus_tag_check'
    ) THEN
        ALTER TABLE ena_candidate_loci
            DROP CONSTRAINT ena_candidate_loci_locus_tag_check;
    END IF;
    ALTER TABLE ena_candidate_loci
        ADD CONSTRAINT ena_candidate_loci_locus_tag_check
        CHECK (locus_tag ~ '^[A-Z][A-Z0-9]{2,11}_[0-9]{9}$');
END
$$;

COMMENT ON COLUMN ena_candidate_loci.locus_tag IS
    'Tag as submitted for this candidate, under the prefix registered to its '
    'technology study.';

-- Selection is per (specimen, technology). Because each technology has its own
-- study accession, the existing primary key (ena_study_accession, og_id) is
-- exactly the "never two records for one OG number within one technology" rule,
-- and the partial unique index on (ena_study_accession, biosample_accession)
-- enforces the same thing against the BioSample.
COMMENT ON TABLE ena_submission_selections IS
    'One selected candidate per specimen per technology study. SELECTED means '
    'the chosen version within that technology, not the only version of the '
    'specimen: every technology with a viable mitogenome is published.';

-- 004 used to create this table as ena_canonical_selections, with a
-- canonical_status column whose CANONICAL value implied one record per
-- specimen.  004 now creates the renamed table directly, so all that is left
-- here is removing the legacy one.  It must be empty: refuse rather than
-- destroy selections that were never carried across.
DO $$
DECLARE
    legacy_rows BIGINT;
BEGIN
    -- Unqualified on purpose: this resolves through search_path, so the
    -- migration behaves the same in public and in the throwaway schema the
    -- integration test applies it into.
    IF to_regclass('ena_canonical_selections') IS NOT NULL THEN
        EXECUTE 'SELECT count(*) FROM ena_canonical_selections' INTO legacy_rows;
        IF legacy_rows > 0 THEN
            RAISE EXCEPTION
                'ena_canonical_selections still holds % row(s); migrate them into '
                'ena_submission_selections before this migration can drop it',
                legacy_rows;
        END IF;
        DROP TABLE ena_canonical_selections;
    END IF;
END
$$;

COMMIT;

-- Post-migration checks:
-- SELECT count(*), count(DISTINCT og_id) FROM ena_locus_registry;
-- SELECT column_name FROM information_schema.columns
--  WHERE table_name = 'ena_locus_registry' ORDER BY ordinal_position;
