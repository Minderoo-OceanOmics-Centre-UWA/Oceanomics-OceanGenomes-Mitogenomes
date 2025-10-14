process HIFI_DATE_QUERY {
    tag "Querying database for ${meta.id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"
    
    input:
    tuple val(sample_id), val(completion_date), val(meta)
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
    from datetime import datetime, timedelta
    
    config_file = "${sql_config}"
    
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
        
        completion_date_str = "${completion_date}"
        
        if completion_date_str != "unknown" and completion_date_str.strip():
            try:
                # Parse the completion date (assuming YYMMDD format)
                completion_date = datetime.strptime(completion_date_str, '%y%m%d').date()
                
                # Calculate search window (e.g., 30 days before completion)
                search_start = completion_date - timedelta(days=30)
                
                # Convert dates back to YYMMDD text format for comparison with database
                search_start_str = search_start.strftime('%y%m%d')
                completion_date_str_formatted = completion_date.strftime('%y%m%d')
                
                # Query using text comparison since seq_date is stored as text
                query = '''
                SELECT seq_date 
                FROM sequencing 
                WHERE og_id = %s 
                  AND seq_date >= %s 
                  AND seq_date <= %s
                ORDER BY seq_date DESC 
                LIMIT 1
                '''
                cursor.execute(query, ('${sample_id}', search_start_str, completion_date_str_formatted))
                
            except ValueError as ve:
                print(f'Error parsing completion date "{completion_date_str}": {ve}', file=sys.stderr)
                # Fall back to original query without date filtering
                query = "SELECT seq_date FROM sequencing WHERE og_id = %s ORDER BY seq_date DESC LIMIT 1"
                cursor.execute(query, ('${sample_id}',))
                
        else:
            # If no completion date available, get the most recent sequencing run
            query = "SELECT seq_date FROM sequencing WHERE og_id = %s ORDER BY seq_date DESC LIMIT 1"
            cursor.execute(query, ('${sample_id}',))
        
        result = cursor.fetchone()
        if result:
            # seq_date is already in YYMMDD format (text), so just return it
            date_str = result[0]
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

    stub:
    """
    echo -n "000000"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        psycopg2: "stub"
    END_VERSIONS
    """
}