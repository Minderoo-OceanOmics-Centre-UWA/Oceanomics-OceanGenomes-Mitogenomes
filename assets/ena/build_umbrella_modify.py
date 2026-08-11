#!/usr/bin/env python3
"""Build the umbrella MODIFY XML from the live umbrella record plus new child accessions.

    ./build_umbrella_modify.py <umbrella_current.xml> <out.xml> PRJEB1 PRJEB2 PRJEB3
    ./build_umbrella_modify.py --description-file desc.txt <current.xml> <out.xml> PRJEB1 ...

MODIFY replaces the project record wholesale, so anything dropped here is dropped at ENA.
The existing PARENT_PROJECT relation is preserved deliberately -- losing it would detach
the umbrella from its NCBI parent. ENA-generated blocks (IDENTIFIERS, PROJECT_LINKS,
PROJECT_ATTRIBUTES) and center_name are stripped because they are server-owned and are
rejected or ignored on resubmission.
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

SERVER_OWNED = ("IDENTIFIERS", "PROJECT_LINKS", "PROJECT_ATTRIBUTES")


def main(argv: list[str]) -> int:
    description = None
    if len(argv) > 2 and argv[1] == "--description-file":
        description = Path(argv[2]).read_text().strip()
        argv = [argv[0], *argv[3:]]

    if len(argv) < 4:
        return f"usage: {argv[0]} [--description-file f] <current.xml> <out.xml> <child> [...]"

    current, out, children = argv[1], argv[2], argv[3:]

    project = ET.parse(current).getroot().find("PROJECT")
    if project is None:
        return f"no PROJECT element in {current}"

    for tag in SERVER_OWNED:
        for node in project.findall(tag):
            project.remove(node)
    project.attrib.pop("center_name", None)

    if project.find("UMBRELLA_PROJECT") is None:
        return f"{current} is not an umbrella project; refusing to build a MODIFY"

    if description is not None:
        node = project.find("DESCRIPTION")
        if node is None:
            return f"no DESCRIPTION element in {current}"
        changed = node.text != description
        node.text = description
        print(f"  desc      : {'REPLACED' if changed else 'already matches'}")

    related = project.find("RELATED_PROJECTS")
    if related is None:
        related = ET.SubElement(project, "RELATED_PROJECTS")

    parents = [
        rp.find("PARENT_PROJECT").get("accession")
        for rp in related.findall("RELATED_PROJECT")
        if rp.find("PARENT_PROJECT") is not None
    ]
    existing = {
        rp.find("CHILD_PROJECT").get("accession")
        for rp in related.findall("RELATED_PROJECT")
        if rp.find("CHILD_PROJECT") is not None
    }

    for accession in children:
        if accession in existing:
            continue
        node = ET.SubElement(related, "RELATED_PROJECT")
        ET.SubElement(node, "CHILD_PROJECT", accession=accession)

    wrapper = ET.Element("PROJECT_SET")
    wrapper.append(project)
    ET.indent(wrapper, space="  ")
    ET.ElementTree(wrapper).write(out, encoding="UTF-8", xml_declaration=True)

    print(f"wrote {out}")
    print(f"  accession : {project.get('accession')}")
    print(f"  parents   : {parents or 'NONE -- investigate before submitting'}")
    print(f"  children  : {sorted(existing | set(children))}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
