#!/usr/bin/env python3
"""
Resolve a taxon name to its NCBI class / order / family from a taxdump directory.

Kept dependency-free (stdlib only) so any script in bin/ can import it, the same
way as species_name_utils: Nextflow bind-mounts bin/ into the task container and
puts it on PATH, so a sibling import resolves via sys.path[0].

The OceanOmics `species` table is the primary source of a sample's lineage, but
it only carries taxa someone has already curated. When it misses, the sample's
class falls back to 'unknown', which downstream is indistinguishable from
'vertebrate' -- wrong genetic code, wrong annotator, wrong BLAST database. This
resolver closes most of that gap from the NCBI taxdump the pipeline already
downloads and caches (modules/local/download_taxonkit_db).

Lookup is by scientific name, so merged.dmp (old taxid -> new taxid) is not
consulted: a retired taxid is only reachable by ID, never by name.
"""
import os
import sys


# Ranks worth indexing by name. A nominal_species_id is usually a binomial, but
# the species table's own fallbacks mean it can legitimately be a genus, family
# or order name, and each of those still pins down a class.
INDEXED_RANKS = frozenset({'species', 'genus', 'family', 'order', 'class'})

# Ranks read back out of a lineage walk.
WANTED_RANKS = ('class', 'order', 'family', 'genus', 'species')

# Open-nomenclature qualifiers. NCBI holds species-rank placeholder nodes such
# as 'Exocoetus sp.', one per submitter's unidentified organism. Matching one
# would resolve our sample to that particular record; drop to the genus instead.
OPEN_NOMENCLATURE = frozenset({'sp', 'sp.', 'spp', 'spp.', 'cf', 'cf.',
                               'aff', 'aff.', 'nr', 'nr.'})


class TaxdumpLineage:
    """Name -> lineage lookups over nodes.dmp / names.dmp.

    Parsing is lazy and happens once per instance: the two files are ~400 MB of
    text, so callers should construct one resolver and reuse it.
    """

    def __init__(self, taxdump_dir):
        self.taxdump_dir = taxdump_dir
        self._loaded = False
        self._parent = {}          # taxid -> parent taxid
        self._rank = {}            # taxid -> rank
        self._name = {}            # taxid -> scientific name (indexed ranks only)
        self._name_to_taxid = {}   # lowercased scientific name -> taxid
        self._ambiguous = set()    # names shared by more than one taxon

    # -- loading ---------------------------------------------------------

    @property
    def available(self):
        """True when both dmp files are present, without parsing them."""
        if not self.taxdump_dir:
            return False
        return (os.path.isfile(os.path.join(self.taxdump_dir, 'nodes.dmp')) and
                os.path.isfile(os.path.join(self.taxdump_dir, 'names.dmp')))

    def load(self):
        if self._loaded:
            return
        if not self.available:
            raise FileNotFoundError(
                f"nodes.dmp/names.dmp not found in taxdump dir: {self.taxdump_dir}")
        self._parse_nodes(os.path.join(self.taxdump_dir, 'nodes.dmp'))
        self._parse_names(os.path.join(self.taxdump_dir, 'names.dmp'))
        self._loaded = True

    def _parse_nodes(self, nodes_file):
        # Integer keys and interned rank strings: the full tree is ~2.5M nodes
        # and the naive str->str form costs well over a gigabyte.
        parent = self._parent
        rank = self._rank
        with open(nodes_file, 'r') as handle:
            for line in handle:
                parts = line.split('\t|\t', 3)
                if len(parts) < 3:
                    continue
                try:
                    taxid = int(parts[0])
                    parent[taxid] = int(parts[1])
                except ValueError:
                    continue
                rank[taxid] = sys.intern(parts[2].strip())

    def _parse_names(self, names_file):
        # nodes.dmp is parsed first, so ranks are known here and only the ranks
        # we actually look up need indexing.
        rank = self._rank
        with open(names_file, 'r') as handle:
            for line in handle:
                parts = line.split('\t|')
                if len(parts) < 4:
                    continue
                name_class = parts[3].strip().strip('|').strip()
                if name_class != 'scientific name':
                    continue
                try:
                    taxid = int(parts[0].strip())
                except ValueError:
                    continue
                if rank.get(taxid) not in INDEXED_RANKS:
                    continue
                name = parts[1].strip()
                self._name[taxid] = name
                key = name.lower()
                existing = self._name_to_taxid.get(key)
                if existing is not None and existing != taxid:
                    # Cross-kingdom homonyms are real (Morus the bird vs Morus
                    # the mulberry). No lineage is better than the wrong one.
                    self._ambiguous.add(key)
                else:
                    self._name_to_taxid[key] = taxid

    # -- lookups ---------------------------------------------------------

    def lineage_for_taxid(self, taxid):
        """Walk to the root collecting the wanted ranks. {} if the taxid is unknown."""
        if taxid is None or taxid not in self._parent:
            return {}
        wanted = set(WANTED_RANKS)
        found = {}
        current = taxid
        seen = set()
        while current and current != 1 and current not in seen:
            seen.add(current)
            rank = self._rank.get(current)
            if rank in wanted and rank not in found:
                name = self._name.get(current)
                if name:
                    found[rank] = name
            current = self._parent.get(current)
        return found

    def lineage_for_name(self, name):
        """
        Resolve a taxon name to {'class': ..., 'order': ..., 'family': ...}.

        Tries the name as given, then its genus (first token). Returns a dict
        with a 'matched_name' / 'matched_rank' pair describing what was hit, or
        {} when the name is absent from NCBI or is an unresolvable homonym.
        """
        self.load()
        for candidate in self._candidates(name):
            key = candidate.lower()
            if key in self._ambiguous:
                continue
            taxid = self._name_to_taxid.get(key)
            if taxid is None:
                continue
            lineage = self.lineage_for_taxid(taxid)
            if not lineage:
                continue
            lineage['matched_name'] = self._name.get(taxid, candidate)
            lineage['matched_rank'] = self._rank.get(taxid, '')
            return lineage
        return {}

    @staticmethod
    def _candidates(name):
        """Name forms to try, most specific first: binomial, then genus."""
        tokens = (name or '').strip().split()
        tokens = [t for t in tokens if t]
        if not tokens:
            return []
        out = []
        if (len(tokens) >= 2 and
                tokens[1].lower() not in OPEN_NOMENCLATURE and
                tokens[1].isalpha()):
            out.append(f"{tokens[0]} {tokens[1]}")
        if tokens[0] not in out:
            out.append(tokens[0])
        return out
