-- Retire this pipeline's locus-tag registry.
--
-- Locus tags are now allocated and injected by the downstream ENA submission
-- pipeline, so nothing here writes a /locus_tag any more: no allocator, no
-- mapping TSV, no tag in the .tbl, .gbf, .embl or .gff.  The two tables that
-- held the allocation are dead weight and, worse, a trap: leaving them in place
-- invites a future reader to treat their contents as the tags a record was
-- actually submitted under, which they no longer are.
--
-- Archived before dropping.  Records already submitted to ENA carry tags that
-- came out of ena_locus_registry, and the (og_id, canonical_gene,
-- gene_occurrence) -> gene_serial assignment is not reconstructible from the
-- flat files once it is gone.  The archives are frozen history, not an
-- interface: the downstream pipeline owns tags from here on and must not read
-- them as a source of truth.
--
-- Run after 009_ena_locus_registry_canonical_order.sql via
-- bin/apply_ena_migrations.py.  Idempotent: the archives are only created on the
-- first run, and dropping an already-dropped table is a no-op.

BEGIN;

-- Snapshot before the drop.  IF NOT EXISTS keeps a re-run from overwriting the
-- archive with an empty result once the source tables are gone.
DO $$
BEGIN
    IF to_regclass('ena_locus_registry') IS NOT NULL THEN
        CREATE TABLE IF NOT EXISTS ena_locus_registry_archive AS
            SELECT * FROM ena_locus_registry;
    END IF;

    IF to_regclass('ena_candidate_loci') IS NOT NULL THEN
        CREATE TABLE IF NOT EXISTS ena_candidate_loci_archive AS
            SELECT * FROM ena_candidate_loci;
    END IF;
END
$$;

-- ena_candidate_loci carries foreign keys onto ena_locus_registry
-- (see 004_ena_candidate_packages.sql), so it has to go first.
DROP TABLE IF EXISTS ena_candidate_loci;
DROP TABLE IF EXISTS ena_locus_registry;

DO $$
BEGIN
    IF to_regclass('ena_locus_registry_archive') IS NOT NULL THEN
        COMMENT ON TABLE ena_locus_registry_archive IS
            'Frozen copy of ena_locus_registry as of migration 010. Records the '
            'specimen gene serials this pipeline allocated before locus-tag '
            'assignment moved to the downstream submission pipeline. Historical '
            'reference only: do not write to it and do not treat it as the '
            'current tag assignment.';
    END IF;

    IF to_regclass('ena_candidate_loci_archive') IS NOT NULL THEN
        COMMENT ON TABLE ena_candidate_loci_archive IS
            'Frozen copy of ena_candidate_loci as of migration 010. Records the '
            'per-candidate locus tags this pipeline submitted before locus-tag '
            'assignment moved to the downstream submission pipeline. Historical '
            'reference only.';
    END IF;
END
$$;

COMMIT;
