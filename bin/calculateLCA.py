#!/usr/bin/env python3
"""
BLAST LCA Analysis Tool

Parses BLAST-tabular output and produces Lowest Common Ancestor (LCA) assignments
by querying Fishbase, WoRMS, and NCBI Taxonomy databases for taxonomic lineages.
Fishbase is queried first, if not in Fishbase, it goes to WoRMS, then it goes to NCBI Taxonomy.

Expected BLAST format:
-outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"
"""

import logging
import sys
import tarfile
from argparse import ArgumentParser
from collections import Counter, OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.error import URLError
from urllib.request import urlretrieve
import pandas as pd


@dataclass
class Config:
    """Configuration settings for the LCA analysis."""

    DEFAULT_CUTOFF: float = 1.0  # for the basepair identity
    DEFAULT_COVER_MIN: float = 90.0
    DEFAULT_PIDENT_CUTOFF: float = 90.0
    BLAST_COLUMNS = [
        "qseqid",
        "sseqid",
        "staxids",
        "sscinames",
        "scomnames",
        "sskingdoms",
        "pident",
        "length",
        "qlen",
        "slen",
        "mismatch",
        "gapopen",
        "gaps",
        "qstart",
        "qend",
        "sstart",
        "send",
        "stitle",
        "evalue",
        "bitscore",
        "qcovs",
        "qcovhsp",
    ]
    PIDENT_COLUMN_INDEX: int = 6
    COV_COLUMN_INDEX: int = 20
    BITSCORE_COLUMN_INDEX: int = -3
    NCBI_TAXDUMP_URL: str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


@dataclass
class TaxonomicLineage:
    """Represents a complete taxonomic lineage."""

    class_name: str
    order: str
    family: str
    genus: str
    species: str

    def to_list(self) -> List[Tuple[str, str]]:
        """Convert to list of (rank, name) tuples."""
        return [
            ("C", self.class_name),
            ("O", self.order),
            ("F", self.family),
            ("G", self.genus),
            ("S", self.species),
        ]


@dataclass
class LCAResult:
    """Result of LCA calculation."""

    percentage: float
    coverage: float
    assignment: str
    included_taxa: Set[str]


class NCBITaxdumpParser:
    """Parses NCBI taxdump files for taxonomic information."""

    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(__name__)
        self.taxid_to_lineage = {}
        self.taxid_to_parent = {}
        self.taxid_to_rank = {}
        self.taxid_to_name = {}

    def download_and_extract_taxdump(self) -> Path:
        """Download and extract NCBI taxdump if not already cached."""
        taxdump_dir = self.cache_dir / "taxdump"
        nodes_file = taxdump_dir / "nodes.dmp"
        names_file = taxdump_dir / "names.dmp"

        if nodes_file.exists() and names_file.exists():
            self.logger.info("Using cached NCBI taxdump files")
            return taxdump_dir

        taxdump_tar = self.cache_dir / "taxdump.tar.gz"

        if not taxdump_tar.exists():
            self.logger.info("Downloading NCBI taxdump...")
            try:
                urlretrieve(Config.NCBI_TAXDUMP_URL, taxdump_tar)
                self.logger.info(f"Downloaded taxdump to: {taxdump_tar}")
            except Exception as e:
                self.logger.error(f"Failed to download taxdump: {e}")
                raise

        self.logger.info("Extracting taxdump...")
        taxdump_dir.mkdir(exist_ok=True)

        with tarfile.open(taxdump_tar, "r:gz") as tar:
            tar.extractall(taxdump_dir)

        self.logger.info(f"Extracted taxdump to: {taxdump_dir}")
        return taxdump_dir

    def parse_nodes_file(self, nodes_file: Path):
        """Parse nodes.dmp file to build taxonomy tree structure."""
        self.logger.info("Parsing nodes.dmp...")

        with open(nodes_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("\t|\t")]
                if len(parts) >= 3:
                    taxid = parts[0]
                    parent_taxid = parts[1]
                    rank = parts[2]

                    self.taxid_to_parent[taxid] = parent_taxid
                    self.taxid_to_rank[taxid] = rank

        self.logger.info(f"Parsed {len(self.taxid_to_parent)} taxonomy nodes")

    def parse_names_file(self, names_file: Path):
        """Parse names.dmp file to get scientific names."""
        self.logger.info("Parsing names.dmp...")

        with open(names_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("\t|\t")]
                if len(parts) >= 4:
                    taxid = parts[0]
                    name = parts[1]
                    name_class = parts[3].rstrip("\t|")

                    if name_class == "scientific name":
                        self.taxid_to_name[taxid] = name

        self.logger.info(f"Parsed {len(self.taxid_to_name)} scientific names")

    def build_lineage(self, taxid: str) -> Optional[TaxonomicLineage]:
        """Build complete taxonomic lineage for a given taxid."""
        if not taxid or taxid == "N/A" or taxid not in self.taxid_to_parent:
            return None

        # Check cache first
        if taxid in self.taxid_to_lineage:
            return self.taxid_to_lineage[taxid]

        lineage_data = {
            "superkingdom": None,
            "kingdom": None,
            "phylum": None,
            "class": None,
            "order": None,
            "family": None,
            "genus": None,
            "species": None,
        }

        current_taxid = taxid
        visited = set()  # Prevent infinite loops

        while current_taxid and current_taxid != "1" and current_taxid not in visited:
            visited.add(current_taxid)

            if (
                current_taxid in self.taxid_to_rank
                and current_taxid in self.taxid_to_name
            ):
                rank = self.taxid_to_rank[current_taxid]
                name = self.taxid_to_name[current_taxid]

                if rank in lineage_data:
                    lineage_data[rank] = name

            # Move to parent
            if current_taxid in self.taxid_to_parent:
                current_taxid = self.taxid_to_parent[current_taxid]
            else:
                break

        lineage = TaxonomicLineage(
            class_name=lineage_data["class"] or "Unknown",
            order=lineage_data["order"] or "Unknown",
            family=lineage_data["family"] or "Unknown",
            genus=lineage_data["genus"] or "Unknown",
            species=lineage_data["species"] or self.taxid_to_name.get(taxid, "Unknown"),
        )

        self.taxid_to_lineage[taxid] = lineage
        return lineage

    def load_taxdump(self):
        """Load and parse the complete NCBI taxdump."""
        taxdump_dir = self.download_and_extract_taxdump()

        nodes_file = taxdump_dir / "nodes.dmp"
        names_file = taxdump_dir / "names.dmp"

        self.parse_nodes_file(nodes_file)
        self.parse_names_file(names_file)

        self.logger.info("NCBI taxdump loaded successfully")


class DatabaseManager:
    """Manages downloading and caching of taxonomic databases."""

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or Path.cwd() / "cache"
        self.cache_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(__name__)
        self.ncbi_parser = NCBITaxdumpParser(self.cache_dir)

    def _download_with_cache(self, url: str, filename: str) -> pd.DataFrame:
        """Download file with caching support."""
        cache_path = self.cache_dir / filename

        if cache_path.exists():
            self.logger.info(f"Loading cached file: {cache_path}")
            return pd.read_parquet(cache_path)

        try:
            self.logger.info(f"Downloading: {url}")
            df = pd.read_parquet(url)
            df.to_parquet(cache_path)
            self.logger.info(f"Cached to: {cache_path}")
            return df
        except (URLError, Exception) as e:
            self.logger.error(f"Failed to download {url}: {e}")
            raise

    def load_fishbase_data(
        self,
    ) -> Tuple[Dict[str, List[str]], Dict[str, str], Dict[str, str]]:
        """Load and process Fishbase data."""
        self.logger.info("Loading Fishbase data...")

        species_df = self._download_with_cache(
            #"https://fishbase.ropensci.org/fishbase/species.parquet",
            "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/species.parquet?download=true",
            "fishbase_species.parquet",
        )

        families_df = self._download_with_cache(
            #"https://fishbase.ropensci.org/fishbase/families.parquet",
            "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/families.parquet?download=true",
            "fishbase_families.parquet",
        )

        synonyms_df = self._download_with_cache(
            #"https://fishbase.ropensci.org/fishbase/synonyms.parquet",
            "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/synonyms.parquet?download=true",
            "fishbase_synonyms.parquet",
        )

        self.logger.debug(f"Species columns: {species_df.columns.tolist()}")
        self.logger.debug(f"Families columns: {families_df.columns.tolist()}")

        merged = species_df.merge(families_df, on="FamCode")
        merged = merged.rename(columns={"Species_x": "Species"})

        merged = merged[
            ["SpecCode", "Species", "Genus", "Family", "Order", "Class"]
        ].copy()

        merged["full_species"] = merged["Genus"] + " " + merged["Species"]

        # Key: genus. value: the whole lineage
        genera_to_lineage = (
            merged.groupby("Genus")[["Family", "Order", "Class"]]
            .first()
            .apply(lambda x: x.tolist(), axis=1)
            .to_dict()
        )

        speccode_to_species = merged.set_index("SpecCode")["full_species"].to_dict()

        # Key: Synonym, value: SpecCode of the 'real' species
        synonyms_to_speccode = dict(
            zip(
                synonyms_df["SynGenus"] + " " + synonyms_df["SynSpecies"],
                synonyms_df["SpecCode"],
            )
        )

        self.logger.info(f"Loaded {len(genera_to_lineage)} Fishbase genera")
        return genera_to_lineage, speccode_to_species, synonyms_to_speccode

    def load_worms_data(
        self, worms_file: Path
    ) -> Tuple[Dict[str, List[str]], Set[str]]:
        """Load and process WoRMS data."""
        if not worms_file.exists():
            self.logger.warning(f"WoRMS file not found: {worms_file}")
            return {}, set()

        self.logger.info(f"Loading WoRMS data from {worms_file}")

        try:
            worms_df = pd.read_csv(
                worms_file,
                sep="\t",
                header=None,
                names=[
                    "Species",
                    "Genus",
                    "Kingdom",
                    "Phylum",
                    "Class",
                    "Order",
                    "Family",
                    "g",
                    "trash",
                    "sp",
                ],
            )

            genera_to_lineage = (
                worms_df.groupby("Genus")[["Family", "Order", "Class"]]
                .first()
                .apply(lambda x: x.tolist(), axis=1)
                .to_dict()
            )

            species_set = set(worms_df["Species"])

            self.logger.info(f"Loaded {len(genera_to_lineage)} WoRMS genera")
            return genera_to_lineage, species_set

        except Exception as e:
            self.logger.error(f"Failed to load WoRMS data: {e}")
            return {}, set()

    def load_ncbi_taxdump(self):
        """Load NCBI taxdump data."""
        self.ncbi_parser.load_taxdump()

    def query_ncbi_taxonomy(self, taxid: str) -> Optional[TaxonomicLineage]:
        """
        Query NCBI Taxonomy database for a given taxid using local taxdump.

        Args:
            taxid: NCBI Taxonomy ID

        Returns:
            TaxonomicLineage object or None if not found
        """
        return self.ncbi_parser.build_lineage(taxid)


class SpeciesNameCorrector:
    """Handles species name corrections and standardization."""

    def __init__(self, corrections_file: Optional[Path] = None):
        self.corrections = {}
        if corrections_file and corrections_file.exists():
            self._load_corrections(corrections_file)
        else:
            # Default corrections - should be moved to config file
            self.corrections = {
                "Petroschmidtia albonotatus": "Petroschmidtia albonotata"
            }

    def _load_corrections(self, corrections_file: Path):
        """Load corrections from file."""
        # TODO - currently I have found only one typo...
        # PROBLEM: If it's a typo, it will always go to the NCBI Taxonomy because
        # that's where the typo is. There may be a correct version on Fishbase/WoRMS
        # that does not have the typo in the name, but we don't find that.
        # Takes alot of manual work!
        pass

    def correct_line(self, line: str) -> str:
        """Apply corrections to a line."""
        for wrong, correct in self.corrections.items():
            line = line.replace(wrong, correct)
        return line


class TaxonomicAssigner:
    """Handles taxonomic assignment logic."""

    def __init__(
        self,
        fishbase_genera: Dict,
        fishbase_speccode: Dict,
        fishbase_synonyms: Dict,
        worms_genera: Dict,
        worms_species: Set,
        db_manager: DatabaseManager,
    ):
        self.fishbase_genera = fishbase_genera
        self.fishbase_speccode = fishbase_speccode
        self.fishbase_synonyms = fishbase_synonyms
        self.worms_genera = worms_genera
        self.worms_species = worms_species
        self.db_manager = db_manager
        self.logger = logging.getLogger(__name__)

    def find_species_info(
        self, line_elements: List[str], taxid: Optional[str] = None
    ) -> Optional[Tuple[str, str, str, TaxonomicLineage]]:
        """
        Find species information from BLAST line elements.

        Args:
            line_elements: Split BLAST line elements
            taxid: NCBI taxonomy ID if available

        Returns:
            Tuple of (genus, species, source, lineage) or None if not found
        """

        # Fishbase comes first,
        genus, species, lineage = self._search_fishbase(line_elements)
        if genus:
            return genus, species, "fishbase", lineage

        # then WoRMs,
        genus, species, lineage = self._search_worms(line_elements)
        if genus:
            return genus, species, "worms", lineage

        # and lastly, NCBI Taxonomy via the taxid
        if taxid:
            lineage = self._search_ncbi(taxid)
            if lineage:
                genus = lineage.genus
                species = lineage.species
                return genus, species, "ncbi", lineage

        return None

    def _search_fishbase(
        self, elements: List[str]
    ) -> Tuple[Optional[str], Optional[str], Optional[TaxonomicLineage]]:
        """Search in Fishbase data."""
        # Direct genus match
        for i, element in enumerate(elements[:-1]):
            if element in self.fishbase_genera:
                genus = element
                species_part = elements[i + 1]
                family, order, class_name = self.fishbase_genera[genus]

                lineage = TaxonomicLineage(
                    class_name=class_name,
                    order=order,
                    family=family,
                    genus=genus,
                    species=f"{genus} {species_part}",
                )
                return genus, species_part, lineage

        # Synonym match
        # IF we have a synonym, we need the correct species' lineage
        for i, element in enumerate(elements[:-1]):
            potential_species = f"{element} {elements[i + 1]}"
            if potential_species in self.fishbase_synonyms:
                speccode = self.fishbase_synonyms[potential_species]
                if speccode in self.fishbase_speccode:
                    correct_species = self.fishbase_speccode[speccode]
                    genus, species_part = correct_species.split(" ", 1)
                    family, order, class_name = self.fishbase_genera[genus]

                    lineage = TaxonomicLineage(
                        class_name=class_name,
                        order=order,
                        family=family,
                        genus=genus,
                        species=correct_species,
                    )
                    return genus, species_part, lineage

        return None, None, None

    def _search_ncbi(self, taxid: str) -> Optional[TaxonomicLineage]:
        """Search in NCBI Taxonomy database."""
        return self.db_manager.query_ncbi_taxonomy(taxid)

    def _search_worms(
        self, elements: List[str]
    ) -> Tuple[Optional[str], Optional[str], Optional[TaxonomicLineage]]:
        """Search in WoRMS data."""
        # I BELIEVE that WoRMS automatically has the correct lineage
        # even for synonyms. I may be wrong about that.
        for i, element in enumerate(elements[:-1]):
            if element in self.worms_genera:
                genus = element
                species_part = elements[i + 1]
                family, order, class_name = self.worms_genera[genus]

                lineage = TaxonomicLineage(
                    class_name=class_name,
                    order=order,
                    family=family,
                    genus=genus,
                    species=f"{genus} {species_part}",
                )
                return genus, species_part, lineage

        return None, None, None


class LCACalculator:
    """Calculates Lowest Common Ancestor from taxonomic hits."""

    def __init__(self, cutoff: float, normalise_identity: bool = False):
        self.cutoff = cutoff
        self.normalise_identity = normalise_identity

    def calculate_lca(self, entries: List[Tuple[float, str, float]]) -> LCAResult:
        """
        Calculate LCA from a list of (percentage, taxon, query_cov, bitscore) tuples.

        Args:
            entries: List of (percentage_identity, taxon_name, query_cov, bitscore) tuples

        Returns:
            LCAResult with percentage, assignment, and included taxa
        """
        if not entries:
            return LCAResult(0.0, 0.0, "no_hits", set())

        # Optionally adjust the bp identity by the coverage
        if self.normalise_identity:
            entries = [
                ((entry[0] / 100) * (entry[2] / 100) * 100, entry[1], entry[2], entry[3])
                for entry in entries
            ]

        #sorted_entries = sorted(entries) # used to score entries by percentage identity alone
        sorted_entries = sorted(entries, key=lambda x: x[-1]) # sort by bitscore, the last number in every tuple
        top_percentage = sorted_entries[-1][0]
        percentage_threshold = top_percentage - self.cutoff
        # top_coverage = sorted_entries_by_coverage[-1][2]
        # coverage_threshold = top_coverage - self.cov_cutoff

        filtered_taxa = set()
        valid_percentages = []
        valid_coverages = []

        for percentage, taxon, coverage, bitscore in sorted_entries:
            if (
                percentage >= percentage_threshold
            ):
                filtered_taxa.add(taxon)
                valid_percentages.append(percentage)
                valid_coverages.append(coverage)

        # mean_percentage = statistics.mean(valid_percentages)
        max_percentage = max(valid_percentages)
        max_coverage = max(valid_coverages)

        if len(filtered_taxa) == 1:
            assignment = list(filtered_taxa)[0]
        else:
            assignment = "dropped"

        # CHANGE: used to return mean_percentage here
        return LCAResult(max_percentage, max_coverage, assignment, filtered_taxa)


class BLASTLCAAnalyzer:
    """Main analyzer class."""

    def __init__(self, config: Config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.db_manager = DatabaseManager()
        self.corrector = SpeciesNameCorrector()
        self.lca_calculator = LCACalculator(config.DEFAULT_CUTOFF)

    def setup_logging(self, log_level: str = "INFO"):
        """Setup logging configuration."""
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler("blast_lca.log"),
            ],
        )

    def load_databases(self, worms_file: Optional[Path] = None):
        """Load all required databases."""
        # Load Fishbase data
        fishbase_genera, fishbase_speccode, fishbase_synonyms = (
            self.db_manager.load_fishbase_data()
        )

        # Load WoRMS data
        worms_genera, worms_species = self.db_manager.load_worms_data(
            worms_file or Path("data/worms_species.txt.gz")
        )
        # Load NCBI
        self.db_manager.load_ncbi_taxdump()

        self.assigner = TaxonomicAssigner(
            fishbase_genera,
            fishbase_speccode,
            fishbase_synonyms,
            worms_genera,
            worms_species,
            self.db_manager,
        )

    def process_blast_file(
        self,
        input_file: Path,
        pident_cutoff: float,
        coverage_cutoff: float,
        missing_file: Path,
    ) -> Dict[str, List]:
        """Process BLAST results file."""
        asv_hits = OrderedDict()
        missing_count = Counter()

        self.logger.info(f"Processing BLAST file: {input_file}")

        try:
            with (
                open(input_file, "r") as infile,
                open(missing_file, "w") as missing_out,
            ):
                for line_num, line in enumerate(infile, 1):
                    try:
                        # fix the bad species names
                        corrected_line = self.corrector.correct_line(line.rstrip())
                        elements = corrected_line.split("\t")

                        #  do we have the right number of columns?
                        # If not, that's a user error. Quit
                        if len(elements) < len(self.config.BLAST_COLUMNS):
                            self.logger.warning(
                                f"Line {line_num}: Insufficient columns"
                            )
                            continue

                        try:
                            pident = float(elements[self.config.PIDENT_COLUMN_INDEX])
                        except (ValueError, IndexError):
                            self.logger.warning(
                                f"Line {line_num}: Invalid pident value"
                            )
                            continue

                        if pident < pident_cutoff:
                            continue

                        try:
                            query_coverage = float(
                                elements[self.config.COV_COLUMN_INDEX]
                            )
                        except (ValueError, IndexError):
                            self.logger.warning(
                                f"Line {line_num}: Invalid query coverage value {elements[self.config.COV_COLUMN_INDEX]}"
                            )
                            continue

                        if query_coverage < coverage_cutoff:
                            continue


                        try:
                            bitscore = float(
                                elements[self.config.BITSCORE_COLUMN_INDEX]
                            )
                        except (ValueError, IndexError):
                            self.logger.warning(
                                f"Line {line_num}: Invalid bitscore value {elements[self.config.BITSCORE_COLUMN_INDEX]}"
                            )
                            continue

                        # Split for species name parsing
                        # BOLD uses | as delimiter so let's also split by that
                        line_elements = corrected_line.replace("|", " ").split()
                        asv_name = line_elements[0]

                        taxid = None
                        if len(elements) > 2 and elements[2] != "N/A":
                            # Sometimes you get multiple taxids
                            taxids = elements[2].split(";")
                            taxid = taxids[0] if taxids else None

                        species_info = self.assigner.find_species_info(
                            line_elements, taxid
                        )

                        if species_info is None:
                            missing_out.write(line)
                            missing_count[asv_name] += 1
                            continue

                        genus, species_part, source, lineage = species_info

                        hit_info = (source, pident, lineage, query_coverage, bitscore)
                        if asv_name not in asv_hits:
                            asv_hits[asv_name] = [hit_info]
                        else:
                            asv_hits[asv_name].append(hit_info)

                    except Exception as e:
                        self.logger.error(f"Error processing line {line_num}: {e}")
                        continue

        except FileNotFoundError:
            self.logger.error(f"Input file not found: {input_file}")
            raise
        except Exception as e:
            self.logger.error(f"Error reading input file: {e}")
            raise

        self.logger.info(f"Processed {len(asv_hits)} ASVs with hits")
        if missing_count:
            self.logger.info(f"Missing species for {len(missing_count)} ASVs")
        else:
            self.logger.info(
                "Found species names or taxonomy IDs in *all* lines. The missing file is empty."
            )

        return asv_hits

    def calculate_lca_assignments(self, asv_hits: Dict) -> List[Dict]:
        """Calculate LCA assignments for all ASVs."""
        results = []

        for asv_name, hits in asv_hits.items():
            species_hits = []
            genus_hits = []
            family_hits = []
            order_hits = []
            class_hits = []
            sources = set()

            for source, pident, lineage, query_cov, bitscore in hits:
                sources.add(source)
                lineage_list = lineage.to_list()

                for rank, name in lineage_list:
                    if not name:
                        # Some entries in worms/fishbase have no Order or Class - ignore
                        continue
                    if rank == "S":
                        species_hits.append((pident, name, query_cov, bitscore))
                    elif rank == "G":
                        genus_hits.append((pident, name, query_cov, bitscore))
                    elif rank == "F":
                        family_hits.append((pident, name, query_cov, bitscore))
                    elif rank == "O":
                        order_hits.append((pident, name, query_cov, bitscore))
                    elif rank == "C":
                        class_hits.append((pident, name, query_cov, bitscore))

            # Calculate LCA at each level
            species_lca = self.lca_calculator.calculate_lca(species_hits)
            genus_lca = self.lca_calculator.calculate_lca(genus_hits)
            family_lca = self.lca_calculator.calculate_lca(family_hits)
            order_lca = self.lca_calculator.calculate_lca(order_hits)
            class_lca = self.lca_calculator.calculate_lca(class_hits)
            results.append(
                {
                    "ASV_name": asv_name,
                    "Class": class_lca.assignment,
                    "Order": order_lca.assignment,
                    "Family": family_lca.assignment,
                    "Genus": genus_lca.assignment,
                    "Species": species_lca.assignment,
                    "PercentageID": f"{species_lca.percentage:.2f}",
                    "Coverage": f"{species_lca.coverage:.2f}",
                    "Species_In_LCA": ", ".join(
                        sorted(species_lca.included_taxa)
                    ),  # BUGFIX: sometimes the order changes.
                    "Sources": ", ".join(sorted(sources)),
                }
            )  # TODO: I suspect that sometimes, it lists too many sources? Might need to filter down

        return results

    def write_results(self, results: List[Dict], output_file: Path):
        """Write results to output file."""
        try:
            with open(output_file, "w") as out:
                header = [
                    "ASV_name",
                    "Class",
                    "Order",
                    "Family",
                    "Genus",
                    "Species",
                    "PercentageID",
                    "Coverage",
                    "Species_In_LCA",
                    "Sources",
                    "LCA_run_date",
                    "gene_region"
                ]
                out.write("\t".join(header) + "\n")

                for result in results:
                    row = [str(result[col]) for col in header]
                    out.write("\t".join(row) + "\n")

            self.logger.info(f"Results written to: {output_file}")

        except Exception as e:
            self.logger.error(f"Error writing results: {e}")
            raise

    def run_analysis(
        self,
        input_file: Path,
        output_file: Path,
        cutoff: float,
        pident_cutoff: float,
        cover_minimum: float,
        missing_file: Path,
        worms_file: Optional[Path] = None,
        normalise_identity: bool = False,
    ):
        """Run the complete analysis pipeline."""
        self.logger.info("Starting BLAST LCA analysis")

        self.lca_calculator.cutoff = cutoff
        self.lca_calculator.normalise_identity = normalise_identity

        self.load_databases(worms_file)

        asv_hits = self.process_blast_file(
            input_file, pident_cutoff, cover_minimum, missing_file
        )

        if not asv_hits:
            self.logger.warning("No valid hits found in input file")
            return

        results = self.calculate_lca_assignments(asv_hits)

        self.write_results(results, output_file)

        self.logger.info("Analysis complete")


def main():
    """Main entry point."""
    config = Config()

    parser = ArgumentParser(
        description="Parses BLAST-tabular output and produces LCAs using Fishbase and WoRMS APIs"
    )
    parser.add_argument(
        "-f", "--file", type=Path, help="Input file of BLAST results", required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file of LCAs (tab-delimited)",
        required=True,
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=config.DEFAULT_CUTOFF,
        help=f"Basepair identity percentage cutoff for LCA calculation (default: {config.DEFAULT_CUTOFF})",
    )
    parser.add_argument(
        "--pident",
        type=float,
        default=config.DEFAULT_PIDENT_CUTOFF,
        help=f"Minimum percentage identity for BLAST hits (default: {config.DEFAULT_PIDENT_CUTOFF})",
    )
    parser.add_argument(
        "--min_coverage",
        type=float,
        default=config.DEFAULT_PIDENT_CUTOFF,
        help=f"Minimum query coverage identity for BLAST hits (default: {config.DEFAULT_COVER_MIN})",
    )

    parser.add_argument(
        "--missing_out",
        type=Path,
        default=Path("missing.csv"),
        help="File to write missing species to (default: missing.csv)",
    )
    parser.add_argument(
        "--worms_file",
        type=Path,
        help="Path to WoRMS species file (optional). Default is worms_species.txt.gz, included in the Github repository.",
    )
    parser.add_argument(
        "--log_level",
        choices=["ERROR", "WARNING", "INFO", "DEBUG"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--normalise_identity",
        action="store_true",
        help="Disable identity normalisation by coverage (default: normalisation disabled). Otherwise bp identity is multiplied by coverage.",
    )

    args = parser.parse_args()

    if not args.file.exists():
        print(f"Error: Input file does not exist: {args.file}")
        sys.exit(1)

    if not (0 <= args.pident <= 100):
        print("Error: --pident must be between 0 and 100")
        sys.exit(1)

    if args.cutoff < 0:
        print("Error: --cutoff must be non-negative")
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.missing_out.parent.mkdir(parents=True, exist_ok=True)

    analyzer = BLASTLCAAnalyzer(config)
    analyzer.setup_logging(args.log_level)

    analyzer.run_analysis(
        input_file=args.file,
        output_file=args.output,
        cutoff=args.cutoff,
        cover_minimum=args.min_coverage,
        pident_cutoff=args.pident,
        missing_file=args.missing_out,
        worms_file=args.worms_file,
        normalise_identity=args.normalise_identity,
    )


if __name__ == "__main__":
    main()