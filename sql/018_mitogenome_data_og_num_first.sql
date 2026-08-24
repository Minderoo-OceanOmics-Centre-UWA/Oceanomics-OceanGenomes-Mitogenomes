-- Puts og_num back as the first column of mitogenome_data.
--
-- Background. og_num was deliberately column 1, so the table is readable in a
-- client without scrolling 90 columns. Migration 016 made it a generated column;
-- PostgreSQL cannot convert a column in place, so 016 had to drop and re-add it,
-- which moved it to position 91. The generated behaviour is right and stays; the
-- position is the regression this fixes.
--
-- PostgreSQL cannot reorder columns, so the table is rebuilt with the columns
-- declared in the wanted order and swapped in. That is the same operation done
-- by hand once before -- the live indexes all carry a _1 suffix because the
-- clean names were already taken by what is now mitogenome_data_SS260818 -- and
-- that earlier rebuild is what silently dropped og_num's generation expression.
-- Doing it as a migration means the generated column is declared explicitly
-- rather than lost, and the whole thing is written down.
--
-- Names are preserved exactly. mitogenome_data_pkey, mitogenome_data_unique and
-- mitogenome_data_depth_method_idx are NOT free to reclaim: they still belong to
-- the SS260818 snapshot, so the _1 suffixes are deliberate, not leftovers.
--
-- Three inbound foreign keys (lca, lca_raw_results, lca_validation) follow the
-- old table through the rename, so they are dropped and re-added against the new
-- one. The readonly role's SELECT grant is re-granted: a rebuilt table would
-- lose it silently.
--
-- Run this only when no pipeline is writing to mitogenome_data. It is one
-- transaction, but the foreign-key swap will block or fail against an in-flight
-- push.
--
-- Applied by bin/apply_ena_migrations.py, or manually:
--
--   psql -h 146.118.120.134 -p 5432 -U postgres -d oceanomics_genomes \
--        -f sql/018_mitogenome_data_og_num_first.sql
--
-- Idempotent: guarded on og_num not already being column 1, so a re-run is a
-- no-op rather than a needless rebuild.

BEGIN;

-- Fail fast rather than queue behind an in-flight push. The rebuild takes
-- ACCESS EXCLUSIVE on mitogenome_data and on the three tables whose foreign keys
-- are re-pointed; without this, a concurrent PUSH_MTDNA_ASSM_RESULTS would sit
-- behind us for as long as we hold them. With it, a busy database aborts the
-- migration cleanly instead, and it can simply be re-run later.
SET LOCAL lock_timeout = '5s';

DO $$
DECLARE
    copied  BIGINT;
    orig    BIGINT;
BEGIN
    IF to_regclass('public.mitogenome_data') IS NULL
       OR EXISTS (
           SELECT 1 FROM pg_attribute
           WHERE attrelid = 'public.mitogenome_data'::regclass
             AND attname = 'og_num' AND attnum = 1
       ) THEN
        RETURN;
    END IF;

    -- 1. Move the current table and its object names aside.
    ALTER TABLE public.mitogenome_data RENAME TO mitogenome_data_reorder_old;
    ALTER TABLE public.mitogenome_data_reorder_old
        RENAME CONSTRAINT mitogenome_data_pkey_1 TO mitogenome_data_pkey_reorder_old;
    ALTER TABLE public.mitogenome_data_reorder_old
        RENAME CONSTRAINT mitogenome_data_unique_1 TO mitogenome_data_unique_reorder_old;
    ALTER INDEX public.mitogenome_data_depth_method_idx_1
        RENAME TO mitogenome_data_depth_method_idx_reorder_old;

    -- 2. The inbound keys now point at the renamed table.
    ALTER TABLE public.lca             DROP CONSTRAINT fk_mitogenome_lca;
    ALTER TABLE public.lca_raw_results DROP CONSTRAINT fk_mitogenome_lca_raw_results;
    ALTER TABLE public.lca_validation  DROP CONSTRAINT fk_mitogenome_lca_validation;

    -- 3. The table as it should be: og_num first, generated.
    CREATE TABLE public.mitogenome_data (
        og_num                 INTEGER GENERATED ALWAYS AS ((SUBSTRING(og_id FROM 3))::INTEGER) STORED,
        og_id                  TEXT NOT NULL,
        tech                   TEXT NOT NULL,
        seq_date               TEXT NOT NULL,
        code                   TEXT NOT NULL,
        annotation             VARCHAR,
        stats                  TEXT,
        length                 INTEGER,
        length_emma            INTEGER,
        seqlength_12s          INTEGER,
        seqlength_16s          INTEGER,
        seqlength_co1          INTEGER,
        cds_no                 INTEGER,
        trna_no                INTEGER,
        rrna_no                INTEGER,
        status                 TEXT,
        genbank                TEXT,
        rrna12s                INTEGER,
        rrna16s                INTEGER,
        atp6                   INTEGER,
        atp8                   INTEGER,
        cox1                   INTEGER,
        cox2                   INTEGER,
        cox3                   INTEGER,
        cytb                   INTEGER,
        nad1                   INTEGER,
        nad2                   INTEGER,
        nad3                   INTEGER,
        nad4                   INTEGER,
        nad4l                  INTEGER,
        nad5                   INTEGER,
        nad6                   INTEGER,
        trna_phe               INTEGER,
        trna_val               INTEGER,
        trna_leuuag            INTEGER,
        trna_leuuaa            INTEGER,
        trna_ile               INTEGER,
        trna_met               INTEGER,
        trna_thr               INTEGER,
        trna_pro               INTEGER,
        trna_lys               INTEGER,
        trna_asp               INTEGER,
        trna_glu               INTEGER,
        trna_sergcu            INTEGER,
        trna_seruga            INTEGER,
        trna_tyr               INTEGER,
        trna_cys               INTEGER,
        trna_trp               INTEGER,
        trna_ala               INTEGER,
        trna_asn               INTEGER,
        trna_gly               INTEGER,
        trna_arg               INTEGER,
        trna_his               INTEGER,
        trna_gln               INTEGER,
        manual_curation_notes  TEXT,
        bankit                 VARCHAR,
        genbank_accession      VARCHAR,
        date_submitted_genbank DATE,
        atp6_trans             INTEGER,
        atp8_trans             INTEGER,
        cox1_trans             INTEGER,
        cox2_trans             INTEGER,
        cox3_trans             INTEGER,
        cytb_trans             INTEGER,
        nad1_trans             INTEGER,
        nad2_trans             INTEGER,
        nad3_trans             INTEGER,
        nad4_trans             INTEGER,
        nad4l_trans            INTEGER,
        nad5_trans             INTEGER,
        nad6_trans             INTEGER,
        extra_genes            VARCHAR,
        missing_genes          VARCHAR,
        order_correct          VARCHAR,
        passed                 VARCHAR,
        mean_depth             DOUBLE PRECISION,
        median_depth           DOUBLE PRECISION,
        depth_sd               DOUBLE PRECISION,
        depth_cv               DOUBLE PRECISION,
        breadth_1x             DOUBLE PRECISION,
        breadth_10x            DOUBLE PRECISION,
        mito_mapped_reads      BIGINT,
        total_reads            BIGINT,
        mito_read_fraction     DOUBLE PRECISION,
        depth_target_length_bp INTEGER,
        depth_target_fasta     TEXT,
        depth_method           TEXT,
        depth_measured_at      TIMESTAMPTZ,
        avg_coverage           REAL,
        avg_base_coverage      REAL,
        CONSTRAINT mitogenome_data_pkey_1   PRIMARY KEY (og_id, tech, seq_date, code),
        CONSTRAINT mitogenome_data_unique_1 UNIQUE (og_id, tech, seq_date, code, annotation)
    );

    -- 4. Copy. og_num is omitted: a generated column rejects an explicit value.
    INSERT INTO public.mitogenome_data (
            og_id, tech, seq_date, code, annotation, stats, length, length_emma, seqlength_12s,
            seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, status, genbank, rrna12s,
            rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, nad5,
            nad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, trna_thr,
            trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, trna_tyr,
            trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, trna_gln,
            manual_curation_notes, bankit, genbank_accession, date_submitted_genbank,
            atp6_trans, atp8_trans, cox1_trans, cox2_trans, cox3_trans, cytb_trans, nad1_trans,
            nad2_trans, nad3_trans, nad4_trans, nad4l_trans, nad5_trans, nad6_trans,
            extra_genes, missing_genes, order_correct, passed, mean_depth, median_depth,
            depth_sd, depth_cv, breadth_1x, breadth_10x, mito_mapped_reads, total_reads,
            mito_read_fraction, depth_target_length_bp, depth_target_fasta, depth_method,
            depth_measured_at, avg_coverage, avg_base_coverage
    )
    SELECT
            og_id, tech, seq_date, code, annotation, stats, length, length_emma, seqlength_12s,
            seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, status, genbank, rrna12s,
            rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, nad5,
            nad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, trna_thr,
            trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, trna_tyr,
            trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, trna_gln,
            manual_curation_notes, bankit, genbank_accession, date_submitted_genbank,
            atp6_trans, atp8_trans, cox1_trans, cox2_trans, cox3_trans, cytb_trans, nad1_trans,
            nad2_trans, nad3_trans, nad4_trans, nad4l_trans, nad5_trans, nad6_trans,
            extra_genes, missing_genes, order_correct, passed, mean_depth, median_depth,
            depth_sd, depth_cv, breadth_1x, breadth_10x, mito_mapped_reads, total_reads,
            mito_read_fraction, depth_target_length_bp, depth_target_fasta, depth_method,
            depth_measured_at, avg_coverage, avg_base_coverage
    FROM public.mitogenome_data_reorder_old;

    -- 5. Refuse to commit a short copy.
    SELECT count(*) INTO copied FROM public.mitogenome_data;
    SELECT count(*) INTO orig   FROM public.mitogenome_data_reorder_old;
    IF copied <> orig THEN
        RAISE EXCEPTION 'mitogenome_data rebuild copied % of % rows', copied, orig;
    END IF;

    -- 6. The non-key index.
    CREATE INDEX mitogenome_data_depth_method_idx_1
        ON public.mitogenome_data (depth_method);

    -- 7. Re-point the inbound keys at the new table, same names as before.
    ALTER TABLE public.lca ADD CONSTRAINT fk_mitogenome_lca
        FOREIGN KEY (og_id, tech, seq_date, code)
        REFERENCES public.mitogenome_data (og_id, tech, seq_date, code);
    ALTER TABLE public.lca_raw_results ADD CONSTRAINT fk_mitogenome_lca_raw_results
        FOREIGN KEY (og_id, tech, seq_date, code)
        REFERENCES public.mitogenome_data (og_id, tech, seq_date, code);
    ALTER TABLE public.lca_validation ADD CONSTRAINT fk_mitogenome_lca_validation
        FOREIGN KEY (og_id, tech, seq_date, code)
        REFERENCES public.mitogenome_data (og_id, tech, seq_date, code);

    -- 8. A rebuilt table starts with no grants.
    GRANT SELECT ON public.mitogenome_data TO readonly;

    -- 9. The old table has served its purpose.
    DROP TABLE public.mitogenome_data_reorder_old;
END $$;

-- Comments live outside the guard so they are restored whichever branch ran.
COMMENT ON COLUMN mitogenome_data.og_num IS
    'Numeric part of og_id, maintained by the database. Generated, not writable: '
    'do not include it in any INSERT or UPDATE column list.';

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

COMMIT;

-- Sanity check after running:
--
--   SELECT attnum, attname, attgenerated FROM pg_attribute
--   WHERE attrelid = 'mitogenome_data'::regclass AND attname = 'og_num';
--   -- expect: 1 | og_num | s
