process MULTIQC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_1' :
        'biocontainers/multiqc:1.30--pyhdfd78af_1' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    python <<'PY'
from collections import OrderedDict
from pathlib import Path

import yaml


def to_builtin(obj):
    if isinstance(obj, dict):
        return {key: to_builtin(value) for key, value in obj.items()}
    return obj


def simplify_group(group_name):
    group_name = str(group_name).strip().strip("'").strip('"')
    return group_name.rsplit(':', 1)[-1]


def load_version_file(path):
    try:
        data = yaml.safe_load(path.read_text()) or {}
    except Exception:
        return OrderedDict()

    parsed = OrderedDict()
    if not isinstance(data, dict):
        return parsed

    for group_name, tools in data.items():
        simple_group = simplify_group(group_name)
        if not isinstance(tools, dict):
            continue

        parsed[simple_group] = OrderedDict(
            (str(tool_name), str(version).strip())
            for tool_name, version in tools.items()
        )

    return parsed


for filename in ('software_versions_mqc_versions.yml', 'software_versions_mqc_versions.yaml'):
    for path in Path('.').rglob(filename):
        path.unlink()

version_files = []
for pattern in (
    'software_versions.yml',
    'software_versions.yaml',
    'versions.yml',
    'versions.yaml',
    'versions_*.yml',
    'versions_*.yaml',
    '*_mqc_versions.yml',
    '*_mqc_versions.yaml',
):
    version_files.extend(
        sorted(
            path for path in Path('.').rglob(pattern)
            if path.name not in {'software_versions_mqc_versions.yml', 'software_versions_mqc_versions.yaml'}
        )
    )

version_files = sorted(
    set(version_files),
    key=lambda path: (0 if path.name.startswith('software_versions') or '_software_mqc_versions' in path.name else 1, str(path)),
)

versions_by_group = OrderedDict()
for version_file in version_files:
    for group_name, tools in load_version_file(version_file).items():
        merged_tools = versions_by_group.setdefault(group_name, OrderedDict())
        for tool_name, version in tools.items():
            merged_tools.setdefault(tool_name, version)

if versions_by_group:
    Path('software_versions_mqc_versions.yml').write_text(
        yaml.safe_dump(to_builtin(versions_by_group), sort_keys=False),
        encoding='utf-8',
    )

# Tool parameter rows (*.tool_params_mqcrow.html) are intentionally NOT aggregated
# into the main report; the "Tool Parameters Used" section is rendered only in the
# per-sample reports (see bin/multiqc_per_sample.py) to keep the main report small.
PY

    multiqc \\
        --force \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
