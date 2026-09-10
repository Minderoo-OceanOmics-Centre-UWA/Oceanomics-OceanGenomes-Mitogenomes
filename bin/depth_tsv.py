#!/usr/bin/env python3
"""Read mean_depth out of a <prefix>.mito_depth.tsv without pandas.

A dependency-light sibling of parse_depth_tsv() in push_mtdna_assm_results.py.
That one is the canonical reader and pulls the whole DEPTH_COLUMNS set for the
mitogenome_data upsert; this one wants a single number and has to run inside
tylerpeirce/psycopg2:0.1, which carries psycopg2 but not pandas. Importing the
pandas version would drag numpy into a container that does not have it, so the
two are kept deliberately separate and only their SEMANTICS are shared:

    a missing file, an empty file, a header-only placeholder, an absent
    mean_depth column and a blank / NA / unparseable cell all return None.

None means "this assembly was never measured" and must never be coerced to 0.0
downstream: ena_package.manifest_fields() omits COVERAGE for None but would
happily emit "COVERAGE 0" for a zero, which is a false statement about the
assembly rather than a missing one.

Keep this in step with parse_depth_tsv() if the TSV layout changes.
"""

import csv


def read_mean_depth(path):
    """Return the mean_depth of the TSV's first data row as a float, or None.

    See the module docstring for what None covers. Any read or parse failure is
    also None: this is a best-effort enrichment of ENA metadata, and a corrupt
    depth file should leave COVERAGE unstated rather than fail the run.
    """
    if path is None:
        return None
    try:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                # First data row only, matching parse_depth_tsv's df.iloc[0].
                # mito_depth.py writes exactly one row per assembly.
                token = (row.get("mean_depth") or "").strip()
                if not token or token.upper() == "NA":
                    return None
                try:
                    return float(token)
                except ValueError:
                    return None
            # No data rows: the header-only placeholder.
            return None
    except OSError:
        return None
