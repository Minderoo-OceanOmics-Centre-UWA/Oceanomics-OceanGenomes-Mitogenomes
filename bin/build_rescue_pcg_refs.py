#!/usr/bin/env python3
"""Build the reference protein set for the post-EMMA CDS rescue.

EMMA's ``rationalise_matches!`` step drops a short protein-coding gene when its
computed circular overlap with a longer neighbour exceeds half the shorter
feature's length. In practice this loses ND4L (vs ND4) and ATP8 (vs ATP6) from
otherwise-complete vertebrate mitogenomes, which then fail the annotation QC gate
and never reach ENA. ``modules/local/emma_gene_rescue`` recovers the dropped gene
from the assembly using the flanking-gene coordinates EMMA already produced,
guarded by BLAST identity/coverage against a small reference set. This script
builds that reference set.

Source: NCBI RefSeq mitochondrion CDS translations for Actinopterygii and
Chondrichthyes (the OceanOmics target space, plus a few tetrapod outgroups for
breadth). Stdlib only (urllib) so it runs on a bare login node; NCBI E-utilities
is queried anonymously. Pass --api-key / --email to raise the rate limit.

Outputs (into --out-dir, default assets/):
    rescue_pcg_refs.faa           FASTA, headers '>{GENE}_{Genus_species}__{acc}'
    rescue_pcg_refs.manifest.tsv  gene, organism, taxon_group, accession, aa_len

Re-run to refresh; the committed copy is the artifact the pipeline ships.
"""

import argparse
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Genes to collect and the header tokens that identify them in a RefSeq
# fasta_cds_aa record ([gene=...] first, then a [protein=...] fallback).
GENES = {
    "ND4L": {"gene": {"ND4L", "NAD4L"}, "protein_sub": "subunit 4l"},
    "ATP8": {"gene": {"ATP8", "ATPASE8"}, "protein_sub": "subunit 8"},
}

# NCBI taxon ids -> label kept in the manifest. Actinopterygii + Chondrichthyes
# cover the pipeline's fish; a small tetrapod outgroup keeps the guard honest for
# the rare divergent sample.
TAXON_GROUPS = [
    ("7898", "Actinopterygii", 26),
    ("7777", "Chondrichthyes", 10),
    ("32523", "Tetrapoda", 4),
]

PER_GENE_CAP = 34  # keep the set small; BLAST only needs frame + a sanity bound

# One reference per genus keeps the set from filling up with congeners from a few
# recently-sequenced groups (loaches, gudgeons) and spreads it across the tree.
DEDUP_RANK = "genus"


def eutils(endpoint, params, api_key=None, email=None, retries=4):
    q = dict(params)
    q.setdefault("tool", "oceanomics-mitogenomes")
    if api_key:
        q["api_key"] = api_key
    if email:
        q["email"] = email
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(q)}"
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - network flake, just retry
            last = exc
            print(f"[build_refs] {endpoint} attempt {attempt + 1} failed: {exc}",
                  file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"[build_refs] {endpoint} failed after {retries} tries: {last}")


def esearch_mito_ids(txid, retmax, api_key, email):
    term = (f"txid{txid}[Organism:exp] AND refseq[filter] AND mitochondrion[filter] "
            f'AND "complete genome"[Title] AND biomol_genomic[PROP]')
    xml = eutils("esearch.fcgi",
                 {"db": "nuccore", "term": term, "retmax": retmax, "retmode": "json"},
                 api_key, email)
    import json
    ids = json.loads(xml).get("esearchresult", {}).get("idlist", [])
    return ids


def esummary_organisms(uids, api_key, email):
    """uid -> (accession.version, organism)."""
    out = {}
    import json
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        data = json.loads(eutils("esummary.fcgi",
                                 {"db": "nuccore", "id": ",".join(chunk),
                                  "retmode": "json"}, api_key, email))
        res = data.get("result", {})
        for uid in res.get("uids", []):
            rec = res.get(uid, {})
            out[uid] = (rec.get("accessionversion") or rec.get("caption", ""),
                        (rec.get("organism") or "").strip())
        time.sleep(0.34)
    return out


def efetch_cds_aa(uids, api_key, email):
    text = []
    for i in range(0, len(uids), 25):
        chunk = uids[i:i + 25]
        text.append(eutils("efetch.fcgi",
                           {"db": "nuccore", "id": ",".join(chunk),
                            "rettype": "fasta_cds_aa", "retmode": "text"},
                           api_key, email))
        time.sleep(0.34)
    return "".join(text)


def parse_fasta(blob):
    name, seq = None, []
    for line in blob.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def header_tokens(header):
    """'[gene=ND4L] [protein=NADH ...]' -> {'gene': 'ND4L', 'protein': 'nadh ...'}."""
    tokens = {}
    depth, key, buf = 0, None, []
    i = 0
    while i < len(header):
        c = header[i]
        if c == "[":
            depth += 1
            if depth == 1:
                key, buf = None, []
                i += 1
                continue
        if c == "]" and depth == 1:
            frag = "".join(buf)
            if "=" in frag:
                k, v = frag.split("=", 1)
                tokens[k.strip().lower()] = v.strip()
            depth -= 1
            i += 1
            continue
        if depth == 1:
            buf.append(c)
        i += 1
    return tokens


def accession_base(header):
    # 'lcl|NC_002333.2_prot_NP_059331.1_11 ...' -> 'NC_002333.2'
    first = header.split()[0]
    if first.startswith("lcl|"):
        first = first[4:]
    marker = "_prot_"
    return first.split(marker, 1)[0] if marker in first else first


def classify_gene(tokens):
    gene = (tokens.get("gene") or "").upper().replace("-", "").replace("_", "")
    protein = (tokens.get("protein") or "").lower()
    for label, spec in GENES.items():
        if gene in {g.replace("-", "").replace("_", "") for g in spec["gene"]}:
            return label
    for label, spec in GENES.items():
        if spec["protein_sub"] in protein and _protein_family_matches(label, protein):
            return label
    return None


def _protein_family_matches(label, protein):
    if label == "ND4L":
        return "nadh" in protein or "dehydrogenase" in protein
    if label == "ATP8":
        return "atp" in protein
    return False


def genus_species(organism):
    parts = organism.split()
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}".replace("'", "").replace(".", "")
    return organism.replace(" ", "_") or "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="assets/panels", type=Path)
    ap.add_argument("--api-key", default=None, help="NCBI E-utilities API key")
    ap.add_argument("--email", default=None, help="contact email for NCBI E-utilities")
    ap.add_argument("--seed", type=int, default=1348, help="RNG seed for sub-sampling")
    ap.add_argument("--search-retmax", type=int, default=500,
                    help="candidate mitogenomes to pull per taxon before sub-sampling")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    per_gene = defaultdict(dict)  # gene -> {genus_species: (organism, group, acc, seq)}

    for txid, group, want in TAXON_GROUPS:
        uids = esearch_mito_ids(txid, args.search_retmax, args.api_key, args.email)
        if not uids:
            print(f"[build_refs] no mitogenomes for {group} (txid{txid})", file=sys.stderr)
            continue
        rng.shuffle(uids)
        # Over-sample: not every mitogenome parses cleanly for both genes, and
        # genus-level dedup discards a lot of the recently-sequenced clusters.
        picked = uids[: max(want * 6, want + 20)]
        organisms = esummary_organisms(picked, args.api_key, args.email)
        blob = efetch_cds_aa(picked, args.api_key, args.email)

        added = 0
        for header, seq in parse_fasta(blob):
            if not seq or "X" * 5 in seq:
                continue
            label = classify_gene(header_tokens(header))
            if not label:
                continue
            acc = accession_base(header)
            uid = next((u for u, (a, _o) in organisms.items() if a == acc), None)
            organism = organisms.get(uid, (acc, ""))[1] if uid else ""
            if not organism:
                continue
            gs = genus_species(organism)
            key = gs.split("_", 1)[0] if DEDUP_RANK == "genus" else gs
            if key in per_gene[label]:
                continue
            per_gene[label][key] = (organism, group, acc, seq.rstrip("*"))
            added += 1
        print(f"[build_refs] {group}: {added} records", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    faa = args.out_dir / "rescue_pcg_refs.faa"
    manifest = args.out_dir / "rescue_pcg_refs.manifest.tsv"

    n_written = 0
    with faa.open("w") as fa, manifest.open("w") as mf:
        mf.write("gene\torganism\ttaxon_group\taccession\taa_len\n")
        for label in GENES:
            recs = list(per_gene[label].values())
            recs.sort(key=lambda r: (r[1], r[0]))
            if len(recs) > PER_GENE_CAP:
                recs = rng.sample(recs, PER_GENE_CAP)
                recs.sort(key=lambda r: (r[1], r[0]))
            if len(recs) < 5:
                raise SystemExit(f"[build_refs] only {len(recs)} {label} refs; "
                                 "check the NCBI queries before committing")
            for organism, group, acc, seq in recs:
                gs = genus_species(organism)
                fa.write(f">{label}_{gs}__{acc}\n")
                for j in range(0, len(seq), 60):
                    fa.write(seq[j:j + 60] + "\n")
                mf.write(f"{label}\t{organism}\t{group}\t{acc}\t{len(seq)}\n")
                n_written += 1

    print(f"[build_refs] wrote {n_written} sequences -> {faa}", file=sys.stderr)
    for label in GENES:
        print(f"[build_refs]   {label}: {len(per_gene[label])} organisms", file=sys.stderr)


if __name__ == "__main__":
    main()
