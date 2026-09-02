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
# or order name, and each of those still pins down a class. Phylum is indexed
# too, one rank coarser than useful: it does NOT pin down a class, but for a
# sample identified no further than 'Porifera' the phylum is still enough to pick
# the genetic code and the annotation route, and it is what the operator has.
INDEXED_RANKS = frozenset({'species', 'genus', 'family', 'order', 'class', 'phylum'})

# Ranks read back out of a lineage walk.
WANTED_RANKS = ('phylum', 'class', 'order', 'family', 'genus', 'species')

# Nodes used to tell an invertebrate from a vertebrate by ancestry rather than by
# a hand-maintained class list. Looked up by name during parsing because neither
# rank ('kingdom', 'clade') is in INDEXED_RANKS.
ANCESTRY_ANCHORS = ('Metazoa', 'Vertebrata')

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
        self._mito_code = {}       # taxid -> NCBI mitochondrial genetic code
        self._name = {}            # taxid -> scientific name (indexed ranks only)
        self._name_to_taxid = {}   # lowercased scientific name -> taxid
        self._name_to_taxids = {}  # lowercased scientific name -> all taxids
        self._ambiguous = set()    # names shared by more than one taxon
        self._anchors = {}         # 'Metazoa' / 'Vertebrata' -> taxid

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
        mito = self._mito_code
        with open(nodes_file, 'r') as handle:
            for line in handle:
                # Field 8 is the taxon's mitochondrial genetic code -- NCBI's own
                # assignment, and what ENA and table2asn validate a submission
                # against. Split far enough to reach it.
                parts = line.split('\t|\t', 9)
                if len(parts) < 3:
                    continue
                try:
                    taxid = int(parts[0])
                    parent[taxid] = int(parts[1])
                except ValueError:
                    continue
                rank[taxid] = sys.intern(parts[2].strip())
                if len(parts) > 8:
                    try:
                        code = int(parts[8].strip())
                    except ValueError:
                        continue
                    if code:            # 0 = unset (the root and a few stubs)
                        mito[taxid] = code

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
                name = parts[1].strip()
                # Anchors first: their ranks ('kingdom', 'clade') are outside
                # INDEXED_RANKS, so the filter below would drop them.
                if name in ANCESTRY_ANCHORS:
                    self._anchors.setdefault(name, []).append(taxid)
                if rank.get(taxid) not in INDEXED_RANKS:
                    continue
                self._name[taxid] = name
                key = name.lower()
                self._name_to_taxids.setdefault(key, []).append(taxid)
                existing = self._name_to_taxid.get(key)
                if existing is not None and existing != taxid:
                    # Cross-kingdom homonyms are real (Morus the bird vs Morus
                    # the mulberry). Recorded here, then narrowed to the animal
                    # candidate in _resolve_taxid() -- see the note there.
                    self._ambiguous.add(key)
                else:
                    self._name_to_taxid[key] = taxid

    # -- lookups ---------------------------------------------------------

    # -- ancestry --------------------------------------------------------

    def _anchor(self, name):
        """Resolve an ANCESTRY_ANCHORS name to its animal-kingdom taxid.

        'Vertebrata' is itself a homonym (the red algal genus Vertebrata), so the
        anchor is the candidate that sits under Metazoa. Metazoa is picked by
        rank instead, there being no Metazoa outside the animals.
        """
        cached = getattr(self, '_anchor_cache', None)
        if cached is None:
            cached = self._anchor_cache = {}
        if name in cached:
            return cached[name]
        candidates = self._anchors.get(name, [])
        chosen = None
        if name == 'Metazoa':
            for taxid in candidates:
                if self._rank.get(taxid) == 'kingdom':
                    chosen = taxid
                    break
        else:
            metazoa = self._anchor('Metazoa')
            for taxid in candidates:
                if metazoa is not None and self.has_ancestor(taxid, metazoa):
                    chosen = taxid
                    break
        if chosen is None and len(candidates) == 1:
            chosen = candidates[0]
        cached[name] = chosen
        return chosen

    def has_ancestor(self, taxid, ancestor):
        """True when `ancestor` is `taxid` or one of its parents."""
        if taxid is None or ancestor is None:
            return False
        current = taxid
        seen = set()
        while current and current != 1 and current not in seen:
            if current == ancestor:
                return True
            seen.add(current)
            current = self._parent.get(current)
        return False

    def is_animal(self, taxid):
        return self.has_ancestor(taxid, self._anchor('Metazoa'))

    def is_vertebrate(self, taxid):
        return self.has_ancestor(taxid, self._anchor('Vertebrata'))

    def _resolve_taxid(self, key):
        """Taxid for a lowercased name, narrowing homonyms to the animal one.

        A cross-kingdom homonym used to be dropped outright, on the grounds that
        no lineage beats the wrong one. That is too blunt here: this pipeline
        sequences animals and never plants, so 'Acanthella' is the sponge genus
        and 'Calantica' the barnacle, and refusing both cost two samples their
        whole lineage. When exactly one candidate is an animal it is the answer.
        Two animal candidates (the phylum Ctenophora and the crane-fly genus
        Ctenophora) are still genuinely ambiguous and still resolve to nothing.
        """
        taxid = self._name_to_taxid.get(key)
        if key not in self._ambiguous:
            return taxid
        animals = [t for t in self._name_to_taxids.get(key, []) if self.is_animal(t)]
        return animals[0] if len(animals) == 1 else None

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
            taxid = self._resolve_taxid(candidate.lower())
            if taxid is None:
                continue
            lineage = self.lineage_for_taxid(taxid)
            if not lineage:
                continue
            lineage['matched_name'] = self._name.get(taxid, candidate)
            lineage['matched_rank'] = self._rank.get(taxid, '')
            # Ancestry, not a class list, is what decides invertebrate status:
            # a class missing from a hand-maintained set reads as 'vertebrate'
            # and silently annotates the sample with the wrong code, annotator
            # and BLAST database. See is_invertebrate() in create_samplesheet.py.
            lineage['is_animal'] = self.is_animal(taxid)
            lineage['is_vertebrate'] = self.is_vertebrate(taxid)
            # NCBI's own mitochondrial code for this exact taxon. More precise
            # than any class-level map can be: Cephalodiscidae is code 33 while
            # its parent class Pterobranchia is 5, so the two cannot both be
            # expressed by a table keyed on class.
            lineage['mito_genetic_code'] = self._mito_code.get(taxid)
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
