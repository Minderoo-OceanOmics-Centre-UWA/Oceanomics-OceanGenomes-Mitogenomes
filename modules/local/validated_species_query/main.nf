process VALIDATED_SPECIES_QUERY {
    tag "Validated species query for ${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val(meta)
    path sql_config
    path taxdump        // NCBI taxdump dir, or [] when --taxonkit_db_dir is unset

    output:
    tuple val(meta), stdout, emit: species
    tuple val(meta), path("*.tax_class.txt"), emit: tax_class
    path "versions.yml", emit: versions

    script:
    def taxdump_dir = taxdump ? taxdump.toString() : ''
    """
    #!/usr/bin/env python3

    import psycopg2
    import os
    import sys
    import configparser
    from pathlib import Path

    config_file = "${sql_config}"
    og_id = "${meta.id}"
    tech = "${meta.sequencing_type}"
    seq_date = "${meta.date}"
    code = "${meta.code}"
    annotation = "${meta.annotation}"
    taxdump_dir = "${taxdump_dir}"

    def load_db_config(cfg):
        if not Path(cfg).exists():
            raise FileNotFoundError(f"Config file '{cfg}' does not exist.")

        parser = configparser.ConfigParser()
        parser.read(cfg)

        if not parser.has_section('postgres'):
            raise ValueError("Missing [postgres] section in config file.")

        required = ['dbname', 'user', 'password', 'host', 'port']
        for key in required:
            if not parser.has_option('postgres', key):
                raise ValueError(f"Missing '{key}' in [postgres] section.")

        return {
            'dbname': parser.get('postgres', 'dbname'),
            'user': parser.get('postgres', 'user'),
            'password': parser.get('postgres', 'password'),
            'host': parser.get('postgres', 'host'),
            'port': parser.getint('postgres', 'port')
        }

    species = 'unknown'
    tax_class = ''
    candidates = []
    conn = None
    cursor = None

    try:
        db = load_db_config(config_file)
        conn = psycopg2.connect(**db)
        cursor = conn.cursor()

        query = (
            "SELECT validated_species_name "
            "FROM lca_validation "
            "WHERE og_id = %s "
            "AND tech = %s "
            "AND seq_date = %s "
            "AND code = %s "
            "AND annotation = %s "
            "LIMIT 1"
        )
        cursor.execute(query, (og_id, tech, seq_date, code, annotation))
        result = cursor.fetchone()

        if result and result[0]:
            species = str(result[0]).strip()

        # Taxonomic class, for the mitochondrial genetic code. Resolved from the
        # VALIDATED species name first and only then from the sample's nominal one:
        # this entrypoint exists for samples whose species was corrected by hand, so
        # the nominal name in `sample` is exactly the value not to trust. Falling back
        # to it still beats no class at all -- a corrected name almost always agrees
        # with the nominal one at class rank, which is the only rank that matters here.
        #
        # Matching mirrors bin/create_samplesheet.py: exact species or synonym in
        # master_species, then genus -- a congener's row is safe here, since only the
        # class is taken and congeners share it. The fuzzy family/order tiers it also
        # has are deliberately left out; a fuzzy match is far too loose a basis for
        # choosing a translation table. master_species carries no invertebrates, so
        # those fall through to the taxdump below.
        class_query = (
            "SELECT sp.class FROM master_species sp "
            "WHERE sp.class IS NOT NULL AND btrim(sp.class) <> '' "
            "AND (lower(sp.species) = lower(%s) OR lower(sp.genus) = lower(%s) "
            "     OR EXISTS (SELECT 1 FROM unnest(string_to_array(sp.synonym, ',')) AS syn(name) "
            "                WHERE lower(btrim(syn.name)) = lower(%s))) "
            "ORDER BY (lower(sp.species) = lower(%s)) DESC, sp.species "
            "LIMIT 1"
        )

        if species and species.lower() != 'unknown':
            candidates.append(species)
        cursor.execute("SELECT trim(nominal_species_id) FROM sample WHERE og_id = %s", (og_id,))
        nominal = cursor.fetchone()
        if nominal and nominal[0]:
            candidates.append(str(nominal[0]).strip())

        for name in candidates:
            genus = name.split(' ')[0]
            cursor.execute(class_query, (name, genus, name, name))
            hit = cursor.fetchone()
            if hit and hit[0]:
                tax_class = str(hit[0]).strip()
                break
    except Exception as exc:
        print(f"Error querying lca_validation: {exc}", file=sys.stderr)
    finally:
        if cursor is not None:
            cursor.close()
        if conn is not None:
            conn.close()

    # NCBI taxdump fallback, as create_samplesheet.py has. Imported from bin/ on
    # PATH for the same reason species_name_utils is below. A phylum stands in for
    # a missing class ('Porifera'), which mito_genetic_codes.json also carries.
    if not tax_class and taxdump_dir and candidates:
        for _entry in os.environ.get("PATH", "").split(os.pathsep):
            if _entry and os.path.isfile(os.path.join(_entry, "taxdump_lineage.py")):
                sys.path.insert(0, _entry)
                break
        try:
            from taxdump_lineage import TaxdumpLineage
            resolver = TaxdumpLineage(taxdump_dir)
            if resolver.available:
                for name in candidates:
                    lineage = resolver.lineage_for_name(name)
                    found = lineage.get('class') or lineage.get('phylum')
                    if found:
                        tax_class = str(found).strip()
                        break
        except Exception as exc:
            print(f"[WARN] taxdump class lookup failed for {og_id}: {exc}", file=sys.stderr)

    # Normalise open nomenclature to the ENA-submittable 'Genus sp.' form. This
    # value becomes the /organism= in the flatfile on the qc-only path, and rows
    # stored before the normalisation existed still hold 'Genus sp' / 'Genus spp.',
    # which webin-cli rejects.
    #
    # IMPORT the canonical implementation rather than copying it. This heredoc runs
    # from the task work dir, so sys.path[0] is not bin/ and a plain sibling import
    # fails -- but Nextflow puts bin/ on PATH, so the directory can be found there.
    # The copy this replaces had drifted in scope: it handled only 'sp'/'spp' and
    # so left a BARE GENUS unnormalised, which ENA rejects with 'Organism is not
    # Submittable' after the sample has cleared every gate upstream. Two
    # implementations of one rule is how that happens.
    for _entry in os.environ.get("PATH", "").split(os.pathsep):
        if _entry and os.path.isfile(os.path.join(_entry, "species_name_utils.py")):
            sys.path.insert(0, _entry)
            break
    try:
        from species_name_utils import normalise_open_nomenclature
        species = normalise_open_nomenclature(species)
    except ImportError:
        # Never fail the task over a cosmetic field. Falling back to the raw value
        # reproduces the pre-normalisation behaviour rather than inventing a name.
        print("[WARN] species_name_utils not on PATH; organism left unnormalised",
              file=sys.stderr)

    print(species, end='')

    # Always written, empty when the class did not resolve: the caller joins on this
    # output, and an optional path would drop every unresolved sample from the run
    # rather than falling back to --translation_table as intended.
    with open("${meta.annotation_prefix}.tax_class.txt", 'w') as fh:
        fh.write(tax_class)

    with open('versions.yml', 'w') as fh:
        fh.write('"${task.process}":\\n')
        fh.write('    python: "3.9"\\n')
        fh.write('    psycopg2: "2.9.5"\\n')
    """

    stub:
    """
    echo -n "unknown"
    : > ${meta.annotation_prefix}.tax_class.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        psycopg2: "stub"
    END_VERSIONS
    """
}
