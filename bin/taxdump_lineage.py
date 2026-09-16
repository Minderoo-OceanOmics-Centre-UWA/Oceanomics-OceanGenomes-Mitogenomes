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

This file is kept byte-identical between the mitogenome and draft-genomes
pipelines. The two behaviours the draft-genomes samplesheet needs and the
mitogenome one must not have are constructor flags defaulting to off:
`index_all_ranks` and `class_falls_back_to_phylum` -- see their docstrings.
"""
import os
import sys


# Ranks worth indexing by name. A nominal_species_id is usually a binomial, but
# the species table's own fallbacks mean it can legitimately be a genus, family
# or order name, and each of those still pins down a class.
INDEXED_RANKS = frozenset({'species', 'genus', 'family', 'order', 'class'})

# Ranks read back out of a lineage walk.
WANTED_RANKS = ('class', 'order', 'family', 'genus', 'species')

# NCBI taxid for Metazoa. Used as the default bounding clade for index_all_ranks:
# everything either pipeline sequences is an animal, so names outside it are only
# ever a source of homonyms and memory.
METAZOA_TAXID = 33208

# Open-nomenclature qualifiers. NCBI holds species-rank placeholder nodes such
# as 'Exocoetus sp.', one per submitter's unidentified organism. Matching one
# would resolve our sample to that particular record; drop to the genus instead.
OPEN_NOMENCLATURE = frozenset({'sp', 'sp.', 'spp', 'spp.', 'cf', 'cf.',
                               'aff', 'aff.', 'nr', 'nr.'})


class TaxdumpLineage:
    """Name -> lineage lookups over nodes.dmp / names.dmp.

    Parsing is lazy and happens once per instance: the two files are ~400 MB of
    text, so callers should construct one resolver and reuse it.

    index_all_ranks
        False (default) indexes only INDEXED_RANKS, which is every rank a curated
        `species` row can carry. True indexes every scientific name whose node
        descends from `clade_taxid`, which is what it takes to resolve a nominal
        name sitting at a rank NCBI uses but the species table does not -- e.g.
        'Caridea' is an infraorder, so the default index misses it entirely. The
        lineage walk still reports class/order/family regardless of the rank of
        the node that matched, so nothing downstream changes shape.

        Bounding the index to a clade is what keeps this affordable: the full dump
        holds ~2.87M scientific names, ~1.32M of them metazoan. It also resolves
        cross-kingdom homonyms for free (Morus the bird vs Morus the mulberry).

    class_falls_back_to_phylum
        NCBI leaves the class rank empty for a number of invertebrate lineages.
        True reports the phylum in the 'class' key when no class node exists on
        the lineage, matching what scripts/taxonomy/load_taxonomy.py already
        writes into `species.class` so the two agree. 'phylum' stays in the
        returned dict either way, so a caller can tell the two apart.

        Leave it off for any caller that routes on the class being a real class
        name (the mitogenome pipeline's INVERT_CLASSES and MitoGeneticCode both
        do); a phylum string there reads as an unknown class, i.e. vertebrate.
    """

    def __init__(self, taxdump_dir, index_all_ranks=False,
                 clade_taxid=METAZOA_TAXID, class_falls_back_to_phylum=False):
        self.taxdump_dir = taxdump_dir
        self.index_all_ranks = index_all_ranks
        # Only meaningful alongside index_all_ranks; None means "no bound".
        self.clade_taxid = clade_taxid
        self.class_falls_back_to_phylum = class_falls_back_to_phylum
        self._loaded = False
        self._parent = {}          # taxid -> parent taxid
        self._rank = {}            # taxid -> rank
        self._name = {}            # taxid -> scientific name (indexed ranks only)
        self._name_to_taxid = {}   # lowercased scientific name -> taxid
        self._ambiguous = set()    # names shared by more than one taxon
        self._in_clade_memo = {}   # taxid -> bool, for the clade bound

        # A phylum is only reportable if its name was indexed, so the fallback
        # has to widen the index too when it is not already indexing everything.
        self._indexed_ranks = INDEXED_RANKS
        self._wanted_ranks = WANTED_RANKS
        if class_falls_back_to_phylum:
            self._indexed_ranks = INDEXED_RANKS | {'phylum'}
            self._wanted_ranks = WANTED_RANKS + ('phylum',)

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

    def _in_clade(self, taxid):
        """Is taxid at or below self.clade_taxid? Memoised over the whole path.

        Walks parent links once per lineage rather than once per node: names.dmp
        asks this ~2.9M times, and without the memo each call would climb the
        tree from scratch.
        """
        target = self.clade_taxid
        if target is None:
            return True
        memo = self._in_clade_memo
        parent = self._parent
        path = []
        current = taxid
        while True:
            if current in memo:
                verdict = memo[current]
                break
            if current == target:
                verdict = True
                break
            nxt = parent.get(current)
            # Root, self-parented root, or a dangling node: not in the clade.
            if nxt is None or nxt == current or nxt == 1:
                verdict = False
                break
            path.append(current)
            current = nxt
        for node in path:
            memo[node] = verdict
        memo[taxid] = verdict
        return verdict

    def _parse_names(self, names_file):
        # nodes.dmp is parsed first, so ranks are known here and only the ranks
        # we actually look up need indexing.
        rank = self._rank
        index_all = self.index_all_ranks
        indexed_ranks = self._indexed_ranks
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
                if index_all:
                    if not self._in_clade(taxid):
                        continue
                elif rank.get(taxid) not in indexed_ranks:
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
        wanted = set(self._wanted_ranks)
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
        if (self.class_falls_back_to_phylum and
                'class' not in found and found.get('phylum')):
            found['class'] = found['phylum']
        return found

    def lineage_for_name(self, name):
        """
        Resolve a taxon name to {'class': ..., 'order': ..., 'family': ...}.

        Tries the name as given, then its genus (first token). Returns a dict
        with a 'matched_name' / 'matched_rank' / 'matched_taxid' triple describing
        what was hit, or {} when the name is absent from NCBI or is an
        unresolvable homonym.
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
            # The taxid of the node that actually matched, at whatever rank that
            # is. FCS-GX takes a --tax-id at any rank, so a class- or genus-rank
            # id is a usable answer where the species table had none at all.
            lineage['matched_taxid'] = taxid
            return lineage
        return {}

    def suggest_near_misses(self, name, limit=3, max_distance=2):
        """Indexed names within a small edit distance of `name`.

        A one-character typo in a specimen label is silently catastrophic: in
        NOVA_260909_LA, 'Actinaria sp.' (for Actiniaria, the sea anemone order)
        failed to resolve, which left invertebrates=false and an empty genetic
        code -- indistinguishable downstream from a resolved vertebrate, so a sea
        anemone was seeded with the animal_mt database, annotated with EMMA and
        judged against the vertebrate 37-gene profile.

        The abort message already names the sample and tells the operator to fix
        the label. Offering the candidate turns a hunt through NCBI into a one
        character edit.

        Returns [] when the taxdump is unavailable, so callers need no guard.
        """
        if not self.available:
            return []
        try:
            self.load()
        except Exception:
            return []

        target = (name or "").strip().split()
        if not target:
            return []
        target = target[0].lower()
        if not target:
            return []

        # Cheap prefilter before the O(n*m) distance: a name within edit distance
        # 2 must share a length band, and in practice the first two characters.
        near = []
        for candidate in self._name_to_taxid:
            cand = candidate          # keys are already lowercased
            if abs(len(cand) - len(target)) > max_distance:
                continue
            if cand[:2] != target[:2]:
                continue
            if cand == target:
                continue
            if _edit_distance_within(cand, target, max_distance):
                near.append(candidate)
                if len(near) >= limit * 4:
                    break
        near.sort(key=lambda c: (_edit_distance_within(c.lower(), target, max_distance) or 99, c))
        return near[:limit]

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


def _edit_distance_within(a, b, limit):
    """Levenshtein distance between a and b if it is <= limit, else None.

    Banded: only the diagonal +/- limit is computed, so a long name pair costs
    O(len * limit) rather than O(len^2).
    """
    la, lb = len(a), len(b)
    if abs(la - lb) > limit:
        return None
    prev = list(range(lb + 1))
    for i in range(1, la + 1):
        cur = [i] + [0] * lb
        lo = max(1, i - limit)
        hi = min(lb, i + limit)
        if lo > 1:
            cur[lo - 1] = limit + 1
        for j in range(lo, hi + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
        for j in range(hi + 1, lb + 1):
            cur[j] = limit + 1
        if min(cur[lo:hi + 1] or [limit + 1]) > limit:
            return None
        prev = cur
    return prev[lb] if prev[lb] <= limit else None
