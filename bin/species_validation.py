#!/usr/bin/env python3
import argparse
import csv
import re
from collections import namedtuple
import configparser
import sys

try:
    import psycopg2
except ImportError:  # Allows unit tests to inject a connection factory.
    psycopg2 = None

from species_name_utils import (
    UNRESOLVED_TAXON_TOKENS,
    genus_of,
    normalise_open_nomenclature,
    parse_nominal,
)

SPECIES_IN_LCA_COLUMN = "species_in_LCA"   # <--- NEW

def load_db_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }

def concatenate_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                outfile.writelines(infile)

def concatenate_lca_files(file_list, output_file):
    """
    Concatenate LCA files that each have a header line.
    Writes the header only once (from the first file).
    """
    first_file = True
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                for line_num, line in enumerate(infile):
                    # Always write the header from the first file
                    if first_file:
                        outfile.write(line)
                    else:
                        # Skip header line for subsequent files
                        if line_num == 0:
                            continue
                        outfile.write(line)
            first_file = False

def normalise_name(name):
    return name.strip().lower().replace('_', ' ')

def get_species_for_ogid(db_params, og_id):
    query = "SELECT nominal_species_id FROM sample WHERE og_id = %s"
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            cur.execute(query, (og_id,))
            result = cur.fetchone()
            return result[0] if result else None
    except Exception as e:
        print(f"[ERROR] Database access failed: {e}")
        return None
    finally:
        if conn:
            conn.close()

# The BLAST hits for one assembly, parsed rather than kept as a text blob.
BlastHits = namedtuple("BlastHits", ["species", "genera"])

# blast_combined.<prefix>.tsv is headerless; column 4 (index 3) is the hit's
# scientific name, taken from the NCBI taxid.
BLAST_SCINAME_INDEX = 3
_BLAST_NULL_NAMES = {"n/a", "na", "", "unknown", "none", "null"}

# ... but a hit against OceanOmics' OWN reference database has no NCBI taxid, so
# columns 3-6 are all 'N/A' and the species name is only in the subject title, as
# a structured '[organism=...]' tag. Those are real, identified hits -- often the
# only ones a sample has -- so they must be read, not discarded.
#
# This is still a PARSE, not a substring search: it extracts a named field, so a
# species name occurring incidentally in free text cannot validate a sample, which
# is the fail-open behaviour this whole function exists to remove.
BLAST_ORGANISM_RE = re.compile(r"\[organism=([^\]]+)\]", re.IGNORECASE)

# A BOLD subject title is '<taxid> <Genus species> <BOLD-ID>|<Genus species>|<MARKER>',
# whose second '|'-delimited field is the subject's own identification. It matters
# because it can differ from the NCBI name for the SAME taxid, which is how a
# synonym reaches us: taxid 334986 is 'Hydrolagus ogilbyi' to NCBI and
# 'Chimaera ogilbyi' in BOLD, and a sample labelled with the latter is the same
# animal. Dropping these would hold correctly-identified samples that had already
# reached ENA.
#
# Also a named field rather than free text, so this stays a parse. The binomial
# shape is required so a marker code or an accession cannot be mistaken for a name.
BLAST_BINOMIAL_RE = re.compile(r"^[A-Z][a-z-]+ [a-z][a-z-]+$")


def load_blast_species_set(filepath):
    """Parse the BLAST table into a species set and the genera derived from it.

    This used to return the whole file as one lowercased string, and the caller
    asked whether the nominal name appeared ANYWHERE in it. That failed in both
    directions. It failed OPEN, because a substring test matches inside a longer
    name and inside unrelated description fields: a bare-genus nominal ID of
    'Serrivomer' matched inside 'Serrivomer jesperseni', passed the QC gate, and
    was then rejected at Webin because ENA does not recognise a bare genus. And it
    failed CLOSED for every nominal ID that is not a binomial, since no amount of
    BLAST evidence makes 'Diaphus sp 1' a substring of 'Diaphus watasei'.

    Set membership on a parsed column fixes the first; matching at the asserted
    rank, in compare_lca_and_blast, fixes the second.

    An empty or missing file yields empty sets and no exception -- a zero-region
    sample's blast_combined is the empty placeholder, and samples with genuinely
    no BLAST hits exist.
    """
    species = set()
    try:
        with open(filepath, "r", newline="") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= BLAST_SCINAME_INDEX:
                    continue
                name = normalise_name(fields[BLAST_SCINAME_INDEX])
                have_ncbi_name = bool(name) and name not in _BLAST_NULL_NAMES
                if have_ncbi_name:
                    species.add(name)

                # The subject's own identification, which can be a SYNONYM of the
                # NCBI name for the same taxid. Read on every row, not only on
                # taxid-less ones, because a synonym appears precisely where an
                # NCBI name is also present.
                for field in fields[BLAST_SCINAME_INDEX + 1:]:
                    if "|" not in field:
                        continue
                    for part in field.split("|")[1:]:
                        part = part.strip()
                        if BLAST_BINOMIAL_RE.match(part):
                            species.add(normalise_name(part))
                            break
                    break

                if have_ncbi_name:
                    continue
                # No NCBI name: try the subject title's [organism=...] tag, which
                # is how an in-house reference hit carries its identification.
                for field in fields[BLAST_SCINAME_INDEX + 1:]:
                    tagged = BLAST_ORGANISM_RE.search(field)
                    if tagged:
                        name = normalise_name(tagged.group(1))
                        if name and name not in _BLAST_NULL_NAMES:
                            species.add(name)
                        break
    except FileNotFoundError:
        return BlastHits(species=frozenset(), genera=frozenset())
    return BlastHits(
        species=frozenset(species),
        # genus_of is the same helper the reference-divergence and
        # reference-relevance checks use, so "congeneric" means one thing.
        genera=frozenset(g for g in (genus_of(s) for s in species) if g),
    )

def parse_assembly_key_from_blast(blast_file):
    """
    Extract (og_id, tech, seq_date, code, annotation) from the first
    query_id field in a BLAST table. The query_id is expected to look like
    'og_id.tech.seq_date.code.annotation' (5 dot-separated parts).
    Returns (og_id, tech, seq_date, code, annotation) or None if the file
    is empty / malformed.
    """
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                first_field = line.split('\t', 1)[0]
                parts = first_field.split('.')
                if len(parts) >= 5:
                    return tuple(parts[:5])
                return None
    except Exception as e:
        print(f"[WARN] Could not parse assembly key from BLAST file '{blast_file}': {e}")
    return None


def upsert_lca_validation(
    db_params, key, validated_species_name, validator="nf-core", force=False,
    validated_rank=None, lca_genus=None
):
    """
    Insert/update the lca_validation row for this sample with the validated
    species name and validator tag. The composite key matches the assembly
    naming convention: (og_id, tech, seq_date, code, annotation).

    Guard: if an existing row carries a validator tag other than "nf-core",
    the row was probably set by a human reviewer or another pipeline and
    must not be silently overwritten by an automated re-run. Skip the
    upsert in that case unless ``force=True``.
    """
    og_id, tech, seq_date, code, annotation = key
    params = {
        "og_id": og_id,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "annotation": annotation,
        "validated_species_name": validated_species_name,
        "validator": validator,
        # The rank the evidence actually supported. A relaxed gate stays honest
        # only if what was relaxed is recorded, and a genus-level release must be
        # distinguishable from a species-level one after the fact.
        "validated_rank": validated_rank,
        # For a family-level release, the genus the LCA resolved. Submitting
        # 'Ophidiidae sp.' when the pipeline already resolved Lamprogrammus throws
        # away information, so record it: these surface as label-upgrade candidates
        # for a curator, and validated_rank makes them selectable if a later policy
        # wants to hold family matches for curation instead.
        "lca_genus": lca_genus,
    }
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            if not force:
                cur.execute(
                    """
                    SELECT validator, validated_species_name
                    FROM lca_validation
                    WHERE og_id = %(og_id)s
                      AND tech = %(tech)s
                      AND seq_date = %(seq_date)s
                      AND code = %(code)s
                      AND annotation = %(annotation)s
                    """,
                    params,
                )
                existing = cur.fetchone()
                if existing is not None:
                    existing_validator, existing_species = existing
                    if (
                        existing_validator is not None
                        and str(existing_validator).strip().lower() != "nf-core"
                    ):
                        print(
                            "⚠️ Existing values preserved for lca_validation "
                            f"{og_id}.{tech}.{seq_date}.{code}.{annotation}: "
                            f"validator='{existing_validator}' (not 'nf-core'); "
                            "pass --force to overwrite."
                        )
                        print(
                            "📌 Preserved stored values: "
                            f"{{'validator': '{existing_validator}', "
                            f"'validated_species_name': '{existing_species}'}}"
                        )
                        return True

            upsert_query = """
            INSERT INTO lca_validation (
                og_id, tech, seq_date, code, annotation,
                validated_species_name, validator, validated_rank, lca_genus
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s,
                %(validated_species_name)s, %(validator)s, %(validated_rank)s,
                %(lca_genus)s
            )
            ON CONFLICT (og_id, tech, seq_date, code, annotation)
            DO UPDATE SET
                validated_species_name = EXCLUDED.validated_species_name,
                validator = EXCLUDED.validator,
                validated_rank = EXCLUDED.validated_rank,
                lca_genus = EXCLUDED.lca_genus
            """
            cur.execute(upsert_query, params)
        conn.commit()
        if force:
            print(
                "✅ Success: lca_validation overwritten (--force) for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        else:
            print(
                "✅ Success: lca_validation upserted for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        return True
    except Exception as e:
        if conn is not None:
            conn.rollback()
        print(f"❌ Database error while writing lca_validation: {e}")
        return False
    finally:
        if conn is not None:
            conn.close()


# LCA lineage values that mean "the LCA declined to resolve this rank".
# calculateLCA.py writes 'dropped' where it declined, and 'Unknown' where the
# lineage lookup returned nothing. Both must count as absent, or a rank match
# would succeed against the literal string 'Unknown'.
_LCA_ABSENT = set(UNRESOLVED_TAXON_TOKENS) | {"dropped"}


def _lca_value(row, column):
    """A lineage value from an lca_combined row, or '' if the LCA declined."""
    value = (row.get(column) or "").strip()
    return "" if value.lower() in _LCA_ABSENT else value


def match_at_rank(nominal, blast_hits, row):
    """Decide whether the BLAST/LCA evidence supports the nominal ID's own claim.

    Returns (found_in_blast, validated_rank). The rule is: match at the rank the
    label ASSERTS, not at the loosest rank that happens to succeed.

      species            the binomial is in the BLAST species set    -> 'species'
      genus              the genus is in the BLAST genus set AND the
                         LCA row agrees at genus                     -> 'genus'
      species_uncertain  as for genus                     -> 'genus_downgraded'
      family             the LCA row agrees at family                -> 'family'
      None               never                                      -> 'unmatched'

    Not falling back to genus for a binomial is the point of the design, not a
    limitation of it. A nominal ID of 'Bathypterois parini' is a species-level
    claim; if the evidence supports only the genus, the species claim is
    UNVERIFIED, and releasing it anyway would quietly convert "we could not
    confirm the species" into "validated" -- the opposite of what this gate is
    for.

    A 'cf.'/'?' label is different: it asserts a TENTATIVE species, so resolving
    it at genus drops an uncertain claim rather than asserting one. That is safe
    in the direction that matters, and the downgrade is recorded rather than
    inferred.

    Requiring BOTH the BLAST genus set and the LCA row to agree for a genus-rank
    match is what keeps a swapped or contaminated sample held: when the hits and
    the LCA agree with each other and both disagree with the label, neither side
    can rescue it.
    """
    rank, value = nominal
    if rank is None or not value:
        return "No", "unmatched"

    if rank == "species":
        return (("Yes", "species") if normalise_name(value) in blast_hits.species
                else ("No", "unmatched"))

    if rank in ("genus", "species_uncertain"):
        genus = normalise_name(value)
        lca_genus = normalise_name(_lca_value(row, "genus"))
        if genus and genus in blast_hits.genera and genus == lca_genus:
            return "Yes", ("genus" if rank == "genus" else "genus_downgraded")
        return "No", "unmatched"

    if rank == "family":
        family = normalise_name(value)
        if family and family == normalise_name(_lca_value(row, "family")):
            return "Yes", "family"
        return "No", "unmatched"

    return "No", "unmatched"


def compare_lca_and_blast(config_path, og_id, lca_files, blast_files, output_file, assembly_prefix=None, force=False):
    """Write the per-region summary TSV and, when validated, upsert lca_validation.

    Returns True on success, False if a database write failed. The caller MUST
    propagate that into the exit status: a failed lca_validation write used to
    print and return normally, so the task exited 0 and Nextflow saw success --
    which is exactly how rows lost to the foreign-key write-ordering race stayed
    invisible until someone read a per-sample upload log. The sibling pushers
    (push_lca_raw_results.py, push_lca_blast_results.py) already fail loudly via
    bin/pg_row_guard.py; this brings the third writer into line.
    """
    db_params = load_db_config(config_path)
    # File-naming prefix only: an OG can have multiple assembly attempts, so
    # combined/summary filenames must be qualified with the unique per-assembly
    # prefix (meta.mt_assembly_prefix) to avoid colliding with other attempts
    # for the same OG. The DB lookups and output row content still use the
    # real og_id.
    prefix = assembly_prefix or og_id

    # Combine files
    concatenate_lca_files(lca_files, f"lca_combined.{prefix}.tsv")
    concatenate_files(blast_files, f"blast_combined.{prefix}.tsv")

    # Get nominal_species_id from DB
    db_species = get_species_for_ogid(db_params, og_id)
    has_nominal_species = db_species is not None

    if has_nominal_species:
        # This value becomes the /organism= in the ENA flatfile (via lca_results.tsv
        # -> evaluate_qc_conditions.py -> FORMAT_FILES --species -> process_files.py),
        # and is also what gets stored in lca_validation.validated_species_name. ENA
        # rejects 'Genus sp' and 'Genus spp.' as not submittable, so normalise to the
        # 'Genus sp.' form here rather than at the point of use.
        raw_db_species = db_species
        normalised_species = normalise_open_nomenclature(db_species)
        if normalised_species != db_species:
            print(f"[INFO] Normalised nominal species '{db_species}' -> '{normalised_species}'")
            db_species = normalised_species

        # Normalise once
        db_species_norm = normalise_name(db_species)
    else:
        raw_db_species = None
        print(
            f"[WARN] OG ID '{og_id}' nominal species not found in database — "
            "species match columns will be recorded as N/A."
        )

    # Load BLAST results as parsed sets, not a text blob.
    blast_hits = load_blast_species_set(f"blast_combined.{prefix}.tsv")

    # What rank does the label actually claim? Classified from the RAW value, so
    # the normalisation above cannot turn a bare genus into something that reads
    # as an explicit 'sp.' claim and change the answer.
    nominal = parse_nominal(raw_db_species) if has_nominal_species else (None, None)
    if has_nominal_species:
        print(f"[INFO] Nominal '{raw_db_species}' asserts rank "
              f"{nominal[0] or 'none'}: {nominal[1] or '-'}")

    sample_validated = False
    validated_ranks = set()
    family_release_genus = None

    # Process LCA and compare
    with open(f"lca_combined.{prefix}.tsv", newline='') as tsvfile, open(output_file, "w", newline='') as out:
        # Use DictReader so we can refer to columns by name
        reader = csv.DictReader(tsvfile, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        # validated_rank is APPENDED, never inserted. evaluate_qc_conditions.py
        # reads Found_in_blast_YN from this file, and inserting ahead of it would
        # silently shift what that reads.
        writer.writerow(["og_id", "LCA_result", "nom_species_id", "Match_YN",
                         "Found_in_blast_YN", "validated_rank"])

        for row in reader:
            if not row:
                continue

            # Get the comma-separated species list from the LCA file
            species_in_lca_raw = row.get(SPECIES_IN_LCA_COLUMN)
            if not species_in_lca_raw:
                # If the column isn't present or is empty, skip this row
                continue

            if not has_nominal_species:
                # Nothing to compare the LCA/BLAST hits against.
                writer.writerow([og_id, species_in_lca_raw, "N/A", "N/A", "N/A", "N/A"])
                continue

            # Split on commas and normalise each species name
            species_list_norm = [
                normalise_name(s)
                for s in species_in_lca_raw.split(',')
                if s.strip()
            ]

            # Check if DB nominal species is among the LCA species list
            match = "Yes" if db_species_norm in species_list_norm else "No"

            # Does the evidence support the label's own claim, at the rank it makes?
            in_blast, validated_rank = match_at_rank(nominal, blast_hits, row)
            if in_blast == "Yes":
                sample_validated = True
                validated_ranks.add(validated_rank)
                if validated_rank == "family" and not family_release_genus:
                    # The LCA is MORE specific than the label here. Keep it.
                    family_release_genus = _lca_value(row, "genus") or None

            # LCA_result column: store the raw species_in_LCA string
            writer.writerow([og_id, species_in_lca_raw, db_species, match, in_blast,
                             validated_rank])

    print(f"[INFO] Results written to: {output_file}")

    if not has_nominal_species:
        # No nominal species to compare against, but the sample was still
        # processed: record an lca_validation row with no species match so
        # it isn't silently absent from the table, and so PUSH_LCA_BLAST_RESULTS
        # (fed by this task's lca_combined/blast_combined outputs) still runs.
        print(
            f"[INFO] OG ID '{og_id}' has no nominal_species_id — recording "
            "lca_validation row with no species match."
        )
        key = parse_assembly_key_from_blast(f"blast_combined.{prefix}.tsv")
        if key is None:
            # Expected for a zero-region sample: its blast_combined is the empty
            # placeholder, so there is no sequence id to derive a key from. Not an
            # error, and deliberately not fatal.
            print(
                "⚠️ Could not derive (og_id, tech, seq_date, code, annotation) from "
                f"blast_combined.{prefix}.tsv — skipping lca_validation upsert."
            )
            return True
        return upsert_lca_validation(
            db_params, key, None, validator="nf-core", force=force
        )

    # Push the validation result to the lca_validation table. We only write a
    # row when the sample is validated (Found_in_blast_YN = Yes for at least
    # one region) so unvalidated samples leave the existing DB row untouched.
    if sample_validated:
        key = parse_assembly_key_from_blast(f"blast_combined.{prefix}.tsv")
        if key is None:
            print(
                "⚠️ Could not derive (og_id, tech, seq_date, code, annotation) from "
                f"blast_combined.{prefix}.tsv — skipping lca_validation upsert."
            )
            return True
        # One rank per assembly in practice, but be explicit if the regions
        # disagree rather than silently picking one.
        rank = ";".join(sorted(validated_ranks)) if validated_ranks else "unmatched"
        return upsert_lca_validation(
            db_params, key, db_species, validator="nf-core", force=force,
            validated_rank=rank,
            lca_genus=family_release_genus if rank == "family" else None,
        )

    print(
        f"ℹ️ Sample {og_id} not validated (no Found_in_blast_YN=Yes) — "
        "skipping lca_validation upsert."
    )
    return True


# ---------------------------
# Entry point
# ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compare per-region LCA / BLAST results against the nominal species "
            "stored in the OceanOmics DB, write a summary TSV, and (when "
            "validated) upsert the result into the lca_validation table."
        )
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite an existing lca_validation row even when its validator "
            "is something other than 'nf-core' (e.g. a human reviewer). Off "
            "by default."
        ),
    )
    parser.add_argument(
        "--assembly-prefix",
        default=None,
        help=(
            "Unique per-assembly filename prefix (meta.mt_assembly_prefix). "
            "Used only to name the combined/summary output files so multiple "
            "assembly attempts for the same OG don't collide; defaults to "
            "og_id when omitted."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("og_id")
    parser.add_argument(
        "lca_files",
        help="Comma-separated list of per-region LCA TSVs.",
    )
    parser.add_argument(
        "blast_files",
        help="Comma-separated list of filtered BLAST TSVs.",
    )
    args = parser.parse_args()

    lca_files = args.lca_files.split(',')
    blast_files = args.blast_files.split(',')

    prefix = args.assembly_prefix or args.og_id
    output_file = f"lca_results.{prefix}.tsv"
    ok = compare_lca_and_blast(
        args.config_file,
        args.og_id,
        lca_files,
        blast_files,
        output_file,
        assembly_prefix=args.assembly_prefix,
        force=args.force,
    )
    # Non-zero on a failed DB write so Nextflow sees the failure. The summary TSV is
    # still written either way, so a retry or a -resume has the same inputs.
    sys.exit(0 if ok else 1)
