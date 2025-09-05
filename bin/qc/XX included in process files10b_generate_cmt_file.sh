#!/bin/bash

    # ── inside your  while read -r dir; do  block ──────────────────────────
    fa=$(find "$dir" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) -print -quit)
    seqid=$(basename "${fa%.*}")   # full FASTA ID, e.g. OG228.ilmn.230324.v177getorg.emma102

    ############################################################################
    # 1 · Derive ASSEMBLY METHOD  (from 4th dot-field, e.g. v177getorg)
    ############################################################################
    field4=$(echo "$seqid" | cut -d'.' -f4)      # -> v177getorg
    clean=${field4#v}                            # -> 177getorg

    digits=$(echo "$clean"   | tr -cd '0-9')     # -> 177
    letters=$(echo "$clean"  | tr -cd '[:alpha:]' | tr '[:upper:]' '[:lower:]')  # -> getorg

    # digits 177  → 1.7.7   (insert dots between every digit)
    version=$(echo "$digits" | sed 's/./&./g; s/\.$//')

    case "$letters" in
    getorg*)   asm_prog="GetOrganelle" ;;
    mitohifi*) asm_prog="MitoHifi"     ;;
    *)         echo "❌  Unknown assembler code in '$seqid' (field '$field4')" >&2; exit 1 ;;
    esac

    assembly_method="$asm_prog v.$version"

    ############################################################################
    # 2 · Derive SEQUENCING TECHNOLOGY  (from seqid pattern)
    ############################################################################
    case "$seqid" in
    *hifi*|*HiFi*)   seq_tech="PacBio HiFi" ;;
    *ilmn*|*Ilmn*)   seq_tech="Illumina"    ;;
    *hic*|*HiC*)     seq_tech="Hi-C"        ;;
    *)               echo "❌  Unknown sequencing tech in '$seqid'" >&2; exit 1 ;;
    esac

    cmt="$dir/$seqid.cmt"
    # --- write the TAB-delimited file in one printf (no external tools) ----
    printf '%s\n' \
    "SeqID	StructuredCommentPrefix	Assembly Method	Sequencing Technology	StructuredCommentSuffix" \
    "$seqid	Assembly-Data	$assembly_method	$seq_tech	Assembly-Data" \
    > "$cmt"

    echo "✅  Wrote $(basename "$cmt")   (assembly = '$assembly_method' ; tech = '$seq_tech')"
    # ----------------------------------------------------------------------
done