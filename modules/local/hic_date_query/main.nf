process HIC_DATE_QUERY {
    tag "Querying database for ${meta.id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"
    
    input:
    tuple val(sample_id), val(meta)
    path sql_config
    
    output:
    tuple val(meta), stdout, emit: date
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env python3
    
    import psycopg2
    import sys
    import os
    import configparser
    from pathlib import Path
    
    config_file = "${sql_config}"
    # Database connection parameters (adjust as needed)
    def load_db_config(config_file):
        if not Path(config_file).exists():
            raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
        
        config = configparser.ConfigParser()
        config.read(config_file)

        if not config.has_section('postgres'):
            raise ValueError("❌ Missing [postgres] section in config file.")

        required_keys = ['dbname', 'user', 'password', 'host', 'port']
        for key in required_keys:
            if not config.has_option('postgres', key):
                raise ValueError(f"❌ Missing '{key}' in [postgres] section of config file.")

        return {
            'dbname': config.get('postgres', 'dbname'),
            'user': config.get('postgres', 'user'),
            'password': config.get('postgres', 'password'),
            'host': config.get('postgres', 'host'),
            'port': config.getint('postgres', 'port')    
        }
    
    try:
        # Connect to database
        db_params = load_db_config(config_file)
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()
        
        # Query for sample date (adjust query as needed)
        query = "SELECT seq_date FROM sequencing WHERE hic_library_tube_id = %s"
        cursor.execute(query, ('${sample_id}',))
        
        result = cursor.fetchone()
        if result:
            # Format date as needed (assuming YYMMDD format)
            date_str = result[0].strftime('%y%m%d') if hasattr(result[0], 'strftime') else str(result[0])
            print(date_str, end='')
        else:
            print('unknown', end='')
            
    except Exception as e:
        print(f'Error querying database: {e}', file=sys.stderr)
        print('unknown', end='')
    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()
    
    # Create versions file
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.9"\\n')
        f.write('    psycopg2: "2.9.5"\\n')
    """
}