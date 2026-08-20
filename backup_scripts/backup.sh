#!/bin/bash --login

#SBATCH --account=pawsey1348
#SBATCH --job-name=mito-backup
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x-%j.out

# -----------------------------------------------------------------------------
# Mitogenome pipeline backup
#
# Two destinations, run in this order:
#
#   1. AWS  : $RUNDIR/mitogenomes/$OG/$ASSEMBLY/ena/package/*
#             -> s3://ocom-oceangenomes/analysed-data/draft-genomes/$OG/$ASSEMBLY/
#             (the final ENA submission files, flattened into the assembly dir,
#             alongside the draft genome files already held there)
#
#   2. Acacia: $RUNDIR/mitogenomes/*
#             -> pawsey0964:oceanomics-mitochondrial-genomes/
#             (full pipeline output, existing OG/ASSEMBLY structure preserved)
#
# Stage 1 must succeed before stage 2 runs in --move mode, otherwise the
# package files would be gone from scratch before they reach AWS.
#
# Usage:
#   sbatch backup.sh                                   # copy in the pipeline
#                                                      # outdir: run dir baked in
#   sbatch backup.sh -r /scratch/pawsey1348/$USER/batch-01
#   bash   backup.sh -r /scratch/pawsey1348/$USER/batch-01 --dry-run
#
# Options:
#   -r, --rundir DIR   Pipeline --outdir (the dir containing mitogenomes/)
#   -m, --move         Move to Acacia instead of copy (frees scratch)
#   -n, --dry-run      Show what would transfer, change nothing
#       --skip-aws     Skip stage 1
#       --skip-acacia  Skip stage 2
#       --audit-only   Run the audit against both remotes, transfer nothing
#   -h, --help         This message
# -----------------------------------------------------------------------------

set -uo pipefail

module load rclone/1.68.1 2>/dev/null

# --- Defaults ----------------------------------------------------------------
# The pipeline run script rewrites the RUNDIR_DEFAULT line below when it copies
# backup_scripts/ into the outdir, so the copied script needs no arguments:
#   sbatch backup_scripts/backup.sh
# It stays empty in the repo copy, where -r/--rundir (or $RUNDIR) is required.
RUNDIR_DEFAULT=""

RUNDIR="${RUNDIR:-$RUNDIR_DEFAULT}"
AWS_DEST="${AWS_DEST:-s3://ocom-oceangenomes/analysed-data/mitogenomes/curated}"
ACACIA_DEST="${ACACIA_DEST:-pawsey0964:oceanomics-mitochondrial-genomes}"

MOVE=false
DRY_RUN=false
SKIP_AWS=false
SKIP_ACACIA=false
AUDIT_ONLY=false

# Print the header comment block (between the two dashed comment rules).
usage() {
    awk '/^# -{10,}/ { n++; next } n >= 2 { exit } n == 1 { sub(/^# ?/, ""); print }' \
        "${BASH_SOURCE[0]}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--rundir)    RUNDIR="$2"; shift 2 ;;
        -m|--move)      MOVE=true; shift ;;
        -n|--dry-run)   DRY_RUN=true; shift ;;
        --skip-aws)     SKIP_AWS=true; shift ;;
        --skip-acacia)  SKIP_ACACIA=true; shift ;;
        --audit-only)   AUDIT_ONLY=true; shift ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$RUNDIR" ]]; then
    echo "ERROR: no run directory given. Use -r/--rundir or export RUNDIR." >&2
    exit 2
fi

RUNDIR_ABS="$(cd "$RUNDIR" 2>/dev/null && pwd -P)" || {
    echo "ERROR: run directory does not exist: $RUNDIR" >&2; exit 2; }
RUNDIR="$RUNDIR_ABS"

MITODIR="$RUNDIR/mitogenomes"
[[ -d "$MITODIR" ]] || { echo "ERROR: not found: $MITODIR" >&2; exit 2; }

STAMP=$(date +%y%m%d_%H%M%S)
LOGDIR="$RUNDIR/backup_logs"
mkdir -p "$LOGDIR"

RCLONE_FLAGS=(--checksum --transfers 16 --checkers 16 --retries 5)

# --progress is unreadable in a slurm log (one rclone call per assembly), so
# only use it interactively.
if [[ -t 1 ]]; then
    RCLONE_FLAGS+=(--progress)
else
    RCLONE_FLAGS+=(--stats 60s --stats-one-line)
fi

$DRY_RUN && RCLONE_FLAGS+=(--dry-run)

log() { echo "[$(date '+%F %T')] $*"; }

log "run dir     : $RUNDIR"
log "mitogenomes : $MITODIR"
log "aws dest    : $AWS_DEST"
log "acacia dest : $ACACIA_DEST"
log "mode        : $($MOVE && echo move || echo copy)$($DRY_RUN && echo ' (dry run)')"

AWS_FAILED=0
ACACIA_FAILED=0

# --- Stage 1: ENA packages to AWS --------------------------------------------
# One rclone call per assembly: the package contents are flattened into
# draft-genomes/$OG/$ASSEMBLY/ so nothing is nested under an ena/package path.
if ! $SKIP_AWS && ! $AUDIT_ONLY; then
    log "=== stage 1: ena packages -> $AWS_DEST ==="

    for ogdir in "$MITODIR"/*/; do
        [[ -d "$ogdir" ]] || continue
        OG=$(basename "$ogdir")

        for asmdir in "$ogdir"*/; do
            [[ -d "$asmdir" ]] || continue
            ASSEMBLY=$(basename "$asmdir")
            PKG="$asmdir/ena/package"

            [[ -d "$PKG" ]] || continue

            if [[ -z "$(find "$PKG" -mindepth 1 -maxdepth 1 -type f -print -quit)" ]]; then
                log "WARN  $OG/$ASSEMBLY: ena/package is empty, skipping"
                continue
            fi

            log "  $OG/$ASSEMBLY -> $AWS_DEST/$OG/$ASSEMBLY/"

            if ! rclone copy "$PKG" "$AWS_DEST/$OG/$ASSEMBLY/" "${RCLONE_FLAGS[@]}"; then
                log "ERROR $OG/$ASSEMBLY: rclone copy to AWS failed"
                AWS_FAILED=$((AWS_FAILED + 1))
            fi
        done
    done

    log "stage 1 done, $AWS_FAILED failure(s)"
fi

# --- Stage 2: full mitogenomes tree to Acacia --------------------------------
if ! $SKIP_ACACIA && ! $AUDIT_ONLY; then
    if $MOVE && [[ "$AWS_FAILED" -ne 0 ]]; then
        log "ERROR: --move requested but stage 1 had $AWS_FAILED failure(s)."
        log "       Refusing to move: the ENA packages would be lost from scratch."
        exit 1
    fi

    if $MOVE && $SKIP_AWS; then
        log "ERROR: --move with --skip-aws would delete un-backed-up ENA packages."
        exit 1
    fi

    RCLONE_OP=$($MOVE && echo move || echo copy)
    log "=== stage 2: $RCLONE_OP $MITODIR -> $ACACIA_DEST ==="

    if ! rclone "$RCLONE_OP" "$MITODIR" "$ACACIA_DEST" "${RCLONE_FLAGS[@]}"; then
        log "ERROR: rclone $RCLONE_OP to Acacia failed"
        ACACIA_FAILED=1
    fi

    log "stage 2 done"
fi

# --- Audit -------------------------------------------------------------------
# Nothing to compare against locally once a move has run, so the audit is
# skipped in that case unless it was asked for on its own.
if $DRY_RUN || { $MOVE && ! $AUDIT_ONLY; }; then
    log "audit skipped ($($DRY_RUN && echo 'dry run' || echo 'source moved'))"
else
    log "=== audit ==="

    AWS_TSV="$LOGDIR/mito_backup_aws.$STAMP.tsv"
    ACACIA_TSV="$LOGDIR/mito_backup_acacia.$STAMP.tsv"

    if ! $SKIP_AWS; then
        printf 'OG\tASSEMBLY\tLOCAL_N\tAWS_N\tLOCAL_BYTES\tAWS_BYTES\tSTATUS\n' > "$AWS_TSV"

        for ogdir in "$MITODIR"/*/; do
            [[ -d "$ogdir" ]] || continue
            OG=$(basename "$ogdir")

            for asmdir in "$ogdir"*/; do
                [[ -d "$asmdir" ]] || continue
                ASSEMBLY=$(basename "$asmdir")
                PKG="$asmdir/ena/package"

                [[ -d "$PKG" ]] || continue

                read -r L_N L_B < <(rclone size --json "$PKG" 2>/dev/null |
                    sed -E 's/.*"count":([0-9]+).*"bytes":([0-9]+).*/\1 \2/')
                read -r R_N R_B < <(rclone size --json "$AWS_DEST/$OG/$ASSEMBLY" 2>/dev/null |
                    sed -E 's/.*"count":([0-9]+).*"bytes":([0-9]+).*/\1 \2/')

                L_N=${L_N:-0}; L_B=${L_B:-0}; R_N=${R_N:-0}; R_B=${R_B:-0}

                # The AWS assembly dir also holds draft genome files, so remote
                # counts are >= local. Verify by hash instead of by count.
                if rclone check "$PKG" "$AWS_DEST/$OG/$ASSEMBLY" --checksum --one-way >/dev/null 2>&1; then
                    STATUS=OK
                else
                    STATUS=MISMATCH
                    AWS_FAILED=$((AWS_FAILED + 1))
                fi

                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$OG" "$ASSEMBLY" "$L_N" "$R_N" "$L_B" "$R_B" "$STATUS" >> "$AWS_TSV"
            done
        done

        log "aws audit    : $AWS_TSV ($(grep -c 'MISMATCH$' "$AWS_TSV") mismatch(es))"
    fi

    if ! $SKIP_ACACIA; then
        printf 'OG\tLOCAL_N\tACACIA_N\tLOCAL_BYTES\tACACIA_BYTES\tSTATUS\n' > "$ACACIA_TSV"

        for ogdir in "$MITODIR"/*/; do
            [[ -d "$ogdir" ]] || continue
            OG=$(basename "$ogdir")

            read -r L_N L_B < <(rclone size --json "$ogdir" 2>/dev/null |
                sed -E 's/.*"count":([0-9]+).*"bytes":([0-9]+).*/\1 \2/')
            read -r R_N R_B < <(rclone size --json "$ACACIA_DEST/$OG" 2>/dev/null |
                sed -E 's/.*"count":([0-9]+).*"bytes":([0-9]+).*/\1 \2/')

            L_N=${L_N:-0}; L_B=${L_B:-0}; R_N=${R_N:-0}; R_B=${R_B:-0}

            if rclone check "$ogdir" "$ACACIA_DEST/$OG" --checksum --one-way >/dev/null 2>&1; then
                STATUS=OK
            else
                STATUS=MISMATCH
                ACACIA_FAILED=1
            fi

            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$OG" "$L_N" "$R_N" "$L_B" "$R_B" "$STATUS" >> "$ACACIA_TSV"
        done

        log "acacia audit : $ACACIA_TSV ($(grep -c 'MISMATCH$' "$ACACIA_TSV") mismatch(es))"
    fi
fi

if [[ "$AWS_FAILED" -ne 0 || "$ACACIA_FAILED" -ne 0 ]]; then
    log "BACKUP INCOMPLETE: aws failures=$AWS_FAILED acacia failures=$ACACIA_FAILED"
    exit 1
fi

log "backup complete"
