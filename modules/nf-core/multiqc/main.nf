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
import html
import re

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

tool_to_group = {
    'getorganelleconfig': 'GETORGANELLE_CONFIG',
    'getorganellefromreads': 'GETORGANELLE_FROMREADS',
    'catfastq': 'CAT_FASTQ',
    'mitohififindreference': 'MITOHIFI_FINDMITOREFERENCE',
    'mitohifi': 'MITOHIFI_MITOHIFI',
    'mitohifiaveragecoverage': 'MITOHIFI_AVERAGE_COVERAGE',
    'emma': 'EMMA',
    'blastblastn': 'BLAST_BLASTN',
    'lca': 'LCA',
    'pushmtdnaassemblyresults': 'PUSH_MTDNA_ASSM_RESULTS',
    'speciesvalidation': 'SPECIES_VALIDATION',
    'pushmtdnaannotationresults': 'PUSH_MTDNA_ANNOTATION_RESULTS',
    'pushlcablastresults': 'PUSH_LCA_BLAST_RESULTS',
    'evaluateqcconditions': 'EVALUATE_QC_CONDITIONS',
    'formatfiles': 'FORMAT_FILES',
    'buildsourcemodifiers': 'BUILD_SOURCE_MODIFIERS',
    'extractgenesgff': 'EXTRACT_GENES_GFF',
    'translategenes': 'TRANSLATE_GENES',
    'table2asn': 'GEN_FILES_TABLE2ASN',
}

tool_order = {
    'GETORGANELLE_CONFIG': 10,
    'GETORGANELLE_FROMREADS': 20,
    'CAT_FASTQ': 30,
    'MITOHIFI_FINDMITOREFERENCE': 40,
    'MITOHIFI_MITOHIFI': 50,
    'MITOHIFI_AVERAGE_COVERAGE': 60,
    'EMMA': 70,
    'BLAST_BLASTN': 80,
    'LCA': 90,
    'PUSH_MTDNA_ASSM_RESULTS': 100,
    'SPECIES_VALIDATION': 110,
    'PUSH_MTDNA_ANNOTATION_RESULTS': 120,
    'PUSH_LCA_BLAST_RESULTS': 130,
    'EVALUATE_QC_CONDITIONS': 140,
    'FORMAT_FILES': 150,
    'BUILD_SOURCE_MODIFIERS': 160,
    'EXTRACT_GENES_GFF': 170,
    'TRANSLATE_GENES': 180,
    'GEN_FILES_TABLE2ASN': 190,
}

row_pattern = re.compile(
    r'<tr><td>(?P<tool>.*?)</td><td>(?P<params>.*?)</td><td>(?P<notes>.*?)</td></tr>',
    re.IGNORECASE | re.DOTALL,
)

tool_rows = []
seen_rows = set()
for row_file in sorted(Path('.').rglob('*.tool_params_mqcrow.html')):
    match = row_pattern.search(row_file.read_text().strip())
    if not match:
        continue

    tool_html = match.group('tool').strip()
    params_html = match.group('params').strip()
    notes_html = match.group('notes').strip()
    tool_key = re.sub(r'[^a-z0-9]+', '', html.unescape(tool_html).lower())

    version_html = 'NA'
    group_name = tool_to_group.get(tool_key)
    if group_name and group_name in versions_by_group:
        group_versions = versions_by_group[group_name]
        if len(group_versions) == 1:
            version_html = f"<samp>{html.escape(next(iter(group_versions.values())))}</samp>"
        else:
            formatted_versions = [
                f"{html.escape(tool_name)} {html.escape(version)}"
                for tool_name, version in group_versions.items()
            ]
            version_html = f"<samp>{'; '.join(formatted_versions)}</samp>"

    row_html = f"<tr><td>{tool_html}</td><td>{version_html}</td><td>{params_html}</td><td>{notes_html}</td></tr>"
    row_key = (tool_html, params_html, notes_html)
    if row_key in seen_rows:
        continue
    seen_rows.add(row_key)
    tool_rows.append(
        (
            tool_order.get(group_name, 999),
            html.unescape(tool_html).lower(),
            row_html,
        )
    )

if tool_rows:
    tool_params_lines = [
        "id: 'nf-core-oceangenomesmitogenomes-tool-parameters'",
        "description: 'Exact tool parameters captured from the module command scripts.'",
        "section_name: 'Tool Parameters Used'",
        "plot_type: 'html'",
        "data: |",
        '    <table class="table table-condensed">',
        '    <thead><tr><th>Tool</th><th>Version</th><th>Effective Parameters</th><th>Notes</th></tr></thead>',
        "    <tbody>",
    ]
    tool_params_lines.extend(
        f"    {row}"
        for _order, _tool_name, row in sorted(tool_rows, key=lambda item: (item[0], item[1]))
    )
    tool_params_lines.extend(["    </tbody>", "    </table>"])
    newline = chr(10)
    Path('tool_parameters_mqc.yaml').write_text(newline.join(tool_params_lines) + newline, encoding='utf-8')
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
