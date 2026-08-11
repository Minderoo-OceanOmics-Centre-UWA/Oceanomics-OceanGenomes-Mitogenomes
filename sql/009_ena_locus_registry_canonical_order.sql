-- Renumber existing gene serials into coordinate order.
--
-- Emma sorts its feature table by the start coordinate as a *string*, so a
-- molecule comes out as 1, 10053, 1027, 10343, 1099 ...  The locus-tag
-- allocator numbered loci by walking the file, so that ordering was baked into
-- every allocated serial: in OG5 the 16S gene at 71..1026 holds serial 29 while
-- the ND4L gene at 10053..10349 holds serial 2.  bin/process_files.py and
-- bin/allocate_ena_locus_tags.py now sort numerically, which fixes newly
-- allocated specimens only: ena_locus_registry is keyed on
-- (og_id, canonical_gene, gene_occurrence) and the allocators skip a gene that
-- already holds a serial, so everything allocated before this migration keeps
-- its scrambled numbering until it is rewritten here.
--
-- Ordering comes from coordinate_snapshot, which the allocator writes as
-- '<start>..<end>' at allocation time.  Minus-strand features are stored
-- start > end (ND6 as '14279..13755'), so the key is the lower of the two.
--
-- Run after 008_ena_submission_queue.sql via bin/apply_ena_migrations.py.
-- Idempotent: a second run computes new_serial = old_serial for every row.
--
-- Locus tags are frozen once published, so this refuses to run if any selection
-- has left NOT_SUBMITTED.  On-disk packages are not touched: the ENA packaging
-- path has to be re-run so each ena/package/ picks up the new tags.

BEGIN;

-- Nothing to renumber before 004 has created the registry.
DO $$
DECLARE
    published BIGINT;
BEGIN
    IF to_regclass('ena_locus_registry') IS NULL THEN
        RETURN;
    END IF;

    -- Unqualified on purpose: this resolves through search_path, so the
    -- migration behaves the same in public and in the throwaway schema the
    -- integration test applies it into.
    IF to_regclass('ena_submission_selections') IS NOT NULL THEN
        EXECUTE $q$
            SELECT count(*) FROM ena_submission_selections
            WHERE archive_status <> 'NOT_SUBMITTED'
        $q$ INTO published;
        IF published > 0 THEN
            RAISE EXCEPTION
                '% selection(s) have already been submitted to ENA; their locus '
                'tags are published and cannot be renumbered', published;
        END IF;
    END IF;

    -- The rendered tag in ena_candidate_loci carries the serial as its last
    -- three digits, so the two tables have to move together.
    --
    -- ON COMMIT DROP clears this on the normal path.  The drop is for the
    -- integration test, which strips the BEGIN/COMMIT so it can roll the whole
    -- thing back, and therefore applies the migration twice inside one
    -- transaction to prove it is idempotent.
    DROP TABLE IF EXISTS ena_locus_serial_map;
    CREATE TEMP TABLE ena_locus_serial_map ON COMMIT DROP AS
    SELECT
        og_id,
        gene_serial AS old_serial,
        CAST(row_number() OVER (
            PARTITION BY og_id
            ORDER BY
                CASE
                    WHEN coordinate_snapshot ~ '^[0-9]+\.\.[0-9]+$'
                    THEN LEAST(
                        split_part(coordinate_snapshot, '..', 1)::INTEGER,
                        split_part(coordinate_snapshot, '..', 2)::INTEGER
                    )
                END NULLS LAST,
                canonical_gene,
                gene_occurrence
        ) AS INTEGER) AS new_serial
    FROM ena_locus_registry;
END
$$;

-- The composite foreign key is not deferrable and the 1-999 check leaves no
-- parking range, so both come off for the duration of the renumber.  Drop the
-- foreign key by catalogue lookup rather than by generated name so this works
-- on databases where 004 was applied under a different constraint name.  Every
-- lookup is scoped with to_regclass(): a bare relname or conname matches across
-- every schema and would reach into public from the test's throwaway schema.
DO $$
DECLARE
    constraint_name TEXT;
BEGIN
    IF to_regclass('ena_candidate_loci') IS NULL THEN
        RETURN;
    END IF;
    SELECT con.conname INTO constraint_name
    FROM pg_constraint con
    WHERE con.contype = 'f'
      AND con.conrelid = to_regclass('ena_candidate_loci')
      AND con.confrelid = to_regclass('ena_locus_registry')
      AND con.conkey = ARRAY[
          (SELECT attnum FROM pg_attribute
           WHERE attrelid = to_regclass('ena_candidate_loci')
             AND attname = 'og_id'),
          (SELECT attnum FROM pg_attribute
           WHERE attrelid = to_regclass('ena_candidate_loci')
             AND attname = 'gene_serial')
      ]::SMALLINT[];
    IF constraint_name IS NOT NULL THEN
        EXECUTE format(
            'ALTER TABLE ena_candidate_loci DROP CONSTRAINT %I', constraint_name
        );
    END IF;
END
$$;

DO $$
DECLARE
    constraint_name TEXT;
BEGIN
    IF to_regclass('ena_locus_registry') IS NULL THEN
        RETURN;
    END IF;
    SELECT con.conname INTO constraint_name
    FROM pg_constraint con
    WHERE con.contype = 'c'
      AND con.conrelid = to_regclass('ena_locus_registry')
      AND con.conname = 'ena_locus_gene_serial_check';
    IF constraint_name IS NOT NULL THEN
        EXECUTE format(
            'ALTER TABLE ena_locus_registry DROP CONSTRAINT %I', constraint_name
        );
    END IF;
END
$$;

-- Two phases: PRIMARY KEY (og_id, gene_serial) means a direct update would
-- collide with a serial the same specimen still holds.  Negating first moves
-- every row out of the positive range in one statement.
DO $$
BEGIN
    IF to_regclass('ena_locus_registry') IS NULL THEN
        RETURN;
    END IF;

    UPDATE ena_locus_registry SET gene_serial = -gene_serial;

    UPDATE ena_locus_registry AS registry
    SET gene_serial = map.new_serial,
        updated_at = CURRENT_TIMESTAMP
    FROM ena_locus_serial_map AS map
    WHERE registry.og_id = map.og_id
      AND registry.gene_serial = -map.old_serial;

    -- Prefix and the six-digit OG numeric are unchanged; only the trailing
    -- three-digit serial moves, so the tag stays within
    -- ena_candidate_loci_locus_tag_check.
    IF to_regclass('ena_candidate_loci') IS NOT NULL THEN
        UPDATE ena_candidate_loci AS loci
        SET gene_serial = map.new_serial,
            locus_tag = regexp_replace(
                loci.locus_tag, '[0-9]{3}$', lpad(map.new_serial::TEXT, 3, '0')
            )
        FROM ena_locus_serial_map AS map
        WHERE loci.og_id = map.og_id
          AND loci.gene_serial = map.old_serial;
    END IF;
END
$$;

DO $$
BEGIN
    IF to_regclass('ena_locus_registry') IS NOT NULL THEN
        ALTER TABLE ena_locus_registry
            ADD CONSTRAINT ena_locus_gene_serial_check
            CHECK (gene_serial BETWEEN 1 AND 999);
    END IF;
    IF to_regclass('ena_candidate_loci') IS NOT NULL THEN
        ALTER TABLE ena_candidate_loci
            ADD CONSTRAINT ena_candidate_loci_og_id_gene_serial_fkey
            FOREIGN KEY (og_id, gene_serial)
            REFERENCES ena_locus_registry (og_id, gene_serial);
    END IF;
END
$$;

COMMENT ON COLUMN ena_locus_registry.gene_serial IS
    'Stable 1-999 serial for this specimen gene, ascending with the lower '
    'coordinate of the feature. Rendered as '
    '<tech prefix>_<og numeric:06d><gene_serial:03d> per candidate.';

COMMIT;

-- Post-migration checks:
-- Every specimen's serials ascend with position; expects zero rows.
-- SELECT og_id FROM (
--     SELECT og_id, gene_serial,
--            LEAST(split_part(coordinate_snapshot, '..', 1)::INTEGER,
--                  split_part(coordinate_snapshot, '..', 2)::INTEGER) AS low,
--            lag(gene_serial) OVER w AS prev_serial,
--            lag(LEAST(split_part(coordinate_snapshot, '..', 1)::INTEGER,
--                      split_part(coordinate_snapshot, '..', 2)::INTEGER))
--                OVER w AS prev_low
--     FROM ena_locus_registry
--     WHERE coordinate_snapshot ~ '^[0-9]+\.\.[0-9]+$'
--     WINDOW w AS (PARTITION BY og_id ORDER BY gene_serial)
-- ) ordered
-- WHERE prev_low IS NOT NULL AND low < prev_low;
--
-- Rendered tags still agree with the registry; expects zero rows.
-- SELECT loci.full_seqid, loci.locus_tag
-- FROM ena_candidate_loci loci
-- JOIN ena_locus_registry registry USING (og_id, gene_serial)
-- WHERE loci.locus_tag !~ '^[A-Z][A-Z0-9]{2,11}_[0-9]{9}$'
--    OR right(loci.locus_tag, 3) <> lpad(registry.gene_serial::TEXT, 3, '0');
