#!/usr/bin/env bash
set -euo pipefail

species=''
attempt=1
while [ "$#" -gt 0 ]; do
    case "$1" in
        --species) species="$2"; shift 2 ;;
        --test-attempt) attempt="$2"; shift 2 ;;
        --outfolder) shift 2 ;;
        *) shift ;;
    esac
done

case "$species" in
    success)
        printf '>reference\nACGT\n' > accession.fasta
        printf 'LOCUS       reference 4 bp DNA\n' > accession.gb
        ;;
    no_reference)
        ;;
    transient_then_success)
        if [ "$attempt" -eq 1 ]; then
            exit 137
        fi
        printf '>reference\nACGT\n' > accession.fasta
        printf 'LOCUS       reference 4 bp DNA\n' > accession.gb
        ;;
    exhausted)
        exit 137
        ;;
    *)
        exit 1
        ;;
esac
