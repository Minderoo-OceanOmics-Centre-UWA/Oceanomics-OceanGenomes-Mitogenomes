process VALIDATED_SPECIES_QUERY {
    tag "Validated species query for ${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val(meta)
    path sql_config

    output:
    tuple val(meta), stdout, emit: species
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import psycopg2
    import sys
    import configparser
    from pathlib import Path

    config_file = "${sql_config}"
    og_id = "${meta.id}"
    tech = "${meta.sequencing_type}"
    seq_date = "${meta.date}"
    code = "${meta.code}"
    annotation = "${meta.annotation}"

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
    except Exception as exc:
        print(f"Error querying lca_validation: {exc}", file=sys.stderr)
    finally:
        if cursor is not None:
            cursor.close()
        if conn is not None:
            conn.close()

    print(species, end='')

    with open('versions.yml', 'w') as fh:
        fh.write('"${task.process}":\\n')
        fh.write('    python: "3.9"\\n')
        fh.write('    psycopg2: "2.9.5"\\n')
    """

    stub:
    """
    echo -n "unknown"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        psycopg2: "stub"
    END_VERSIONS
    """
}
