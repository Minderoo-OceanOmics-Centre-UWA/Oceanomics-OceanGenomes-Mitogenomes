#!/usr/bin/env python3
"""
Rank-keyed expectations for genes outside the vertebrate set.

Loads assets/lineage_gene_expectations.json and answers one question: given a
sample's taxonomy, which non-vertebrate genes should be present, and is a missing
one blocking or merely worth reporting?

Deliberately small and deliberately weak. Most invertebrates in NOVA_260909_LA
have a blank order and several have wrong family-level labels, so a rule keyed on
a taxonomic rank is inert or unreliable for them. The rule that actually catches
the batch's losses is unconditional and lives in annotation_stats.py: if MITOS
called a gene and conversion lost it, that is an error regardless of taxonomy.
This module only covers the other direction -- MITOS did not call it, and for
this lineage it should have.
"""

import json
import os

_DEFAULT_ASSET = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "assets", "lineage_gene_expectations.json")

# Ranks from most to least specific: the most specific matching rule wins, so an
# order-level rule is not diluted by a phylum-level one for the same gene.
_RANK_ORDER = ("order", "class", "phylum")


def load_rules(path=None):
    path = path or os.environ.get("LINEAGE_GENE_EXPECTATIONS") or _DEFAULT_ASSET
    if not os.path.exists(path):
        # A missing asset must not fail the run: the unconditional
        # "MITOS called it and we lost it" rule does the load-bearing work.
        return []
    with open(path) as fh:
        return json.load(fh).get("rules", [])


def expected_genes(taxonomy, rules=None):
    """Return {GENE: severity} for this sample.

    `taxonomy` is a mapping with any of 'order', 'class', 'phylum' (missing or
    blank values simply match nothing). severity is 'blocking' or 'advisory'.
    """
    rules = load_rules() if rules is None else rules
    norm = {r: str((taxonomy or {}).get(r) or "").strip().lower() for r in _RANK_ORDER}

    best = {}
    for rule in rules:
        rank = str(rule.get("rank", "")).lower()
        taxon = str(rule.get("taxon", "")).strip().lower()
        gene = str(rule.get("gene", "")).strip().upper()
        if rank not in _RANK_ORDER or not gene or not taxon:
            continue
        if norm.get(rank) != taxon:
            continue
        # Most specific rank wins for a given gene.
        prev = best.get(gene)
        if prev is None or _RANK_ORDER.index(rank) < _RANK_ORDER.index(prev[0]):
            best[gene] = (rank, str(rule.get("severity", "advisory")).lower())
    return {g: sev for g, (_rank, sev) in best.items()}


def classify_missing(taxonomy, present_genes, rules=None):
    """Split the expected-but-absent genes into (blocking, advisory) lists."""
    present = {str(g).strip().upper() for g in (present_genes or [])}
    blocking, advisory = [], []
    for gene, severity in sorted(expected_genes(taxonomy, rules).items()):
        if gene in present:
            continue
        (blocking if severity == "blocking" else advisory).append(gene)
    return blocking, advisory
