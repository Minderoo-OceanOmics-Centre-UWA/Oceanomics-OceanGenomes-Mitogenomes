#!/bin/bash -l
# Calculate total SU usage for a Nextflow run, sourcing the job list from
# Nextflow's own execution database (`nextflow log <session-id>`) rather
# than the rotating .nextflow.log text files. Log rotation caps at 10
# files by default, so on a long, many-times-resumed run the text logs
# can silently lose the record of tasks that ran once early and have
# been CACHED (not resubmitted) on every resume since - the execution
# database has no such limit and is the authoritative source.
# Usage: ./nf_workflow_cost.sh [--tag-filter <file>] <nextflow-run-dir> [output_csv]
#
# <nextflow-run-dir>    The directory Nextflow was launched from for this
#                       run - must contain a .nextflow/ subdirectory with
#                       history and cache/<session-id>/ intact.
# --tag-filter <file>   Only include jobs whose Nextflow task tag (the
#                       text in parentheses, e.g. the sample ID "OG5" in
#                       "GETORGANELLE_RESEED (OG5)") exactly matches one
#                       line in <file>. Useful for scoping a cost/
#                       efficiency report to a subset of samples from a
#                       shared pipeline run without re-running anything.
#
# Outputs (alongside <output_csv>, default metrics_summary.csv):
#   <output_csv>              per-job usage, now including a `tag` column
#   <output_csv basename>.report.txt   human-readable summary
#   cost_per_sample.csv       SU grouped per OceanGenomes sample (cost per
#                             genome): direct SU + amortised shared overhead

set -euo pipefail

TAG_FILTER=""
POSITIONAL=()
while [ $# -gt 0 ]; do
    case "$1" in
        --tag-filter)
            TAG_FILTER="${2:?--tag-filter requires a file argument}"
            shift 2
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL[@]}"

RUN_DIR="${1:?Usage: $0 [--tag-filter <file>] <nextflow-run-dir> [output_csv]}"
OUTFILE="${2:-metrics_summary.csv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -n "$TAG_FILTER" ] && [ ! -f "$TAG_FILTER" ]; then
    echo "Error: --tag-filter file '$TAG_FILTER' not found." >&2
    exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Error: 'nextflow' not found in PATH (module load it first)." >&2
    exit 1
fi

if [ ! -f "$RUN_DIR/.nextflow/history" ]; then
    echo "Error: '$RUN_DIR/.nextflow/history' not found - is this a Nextflow run directory?" >&2
    exit 1
fi

# nextflow log needs to create a local plugins dir, which fails if
# $RUN_DIR isn't writable (common - it's often someone else's project
# dir). Copy just the .nextflow metadata (history + cache, both small -
# task results/work dirs are NOT copied) to a writable temp dir instead.
NXF_WORKDIR=$(mktemp -d)
trap 'rm -rf "$NXF_WORKDIR"' EXIT
cp -r "$RUN_DIR/.nextflow" "$NXF_WORKDIR/"

# A run directory can contain more than one independent session (not
# just resumes of one another); query every distinct session ID on
# record so nothing is missed.
mapfile -t SESSION_IDS < <(awk -F'\t' '{print $6}' "$NXF_WORKDIR/.nextflow/history" | sort -u)
if [ "${#SESSION_IDS[@]}" -eq 0 ]; then
    echo "Error: no sessions found in $RUN_DIR/.nextflow/history" >&2
    exit 1
fi
echo "Found ${#SESSION_IDS[@]} session(s) in .nextflow/history: ${SESSION_IDS[*]}"

# Pull process, tag and native (SLURM) job ID for every task Nextflow has
# ever recorded for these sessions - including CACHED tasks, whose
# original submission may no longer exist in any rotated text log.
JOBID_META_MAP=$(mktemp)
(
    cd "$NXF_WORKDIR"
    for sid in "${SESSION_IDS[@]}"; do
        nextflow log "$sid" -f native_id,process,tag 2>/dev/null
    done
) | awk -F'\t' '$1 != "" && $1 != "-" { print $1","$2","$3 }' \
  | sort -t',' -k1,1 -u > "$JOBID_META_MAP"

# Build the job ID list, optionally restricted to a tag subset.
if [ -n "$TAG_FILTER" ]; then
    mapfile -t JOB_IDS < <(awk -F',' -v OFS=',' '
        NR==FNR { tags[$1]=1; next }
        ($3 in tags) { print $1 }
    ' "$TAG_FILTER" "$JOBID_META_MAP" | sort -un)
    echo "Applied --tag-filter '$TAG_FILTER' ($(wc -l < "$TAG_FILTER") tag(s))"
else
    mapfile -t JOB_IDS < <(cut -d',' -f1 "$JOBID_META_MAP" | sort -un)
fi

if [ "${#JOB_IDS[@]}" -eq 0 ]; then
    echo "No matching submitted SLURM jobs found for session(s): ${SESSION_IDS[*]}" >&2
    exit 1
fi

echo "Found ${#JOB_IDS[@]} unique submitted SLURM job(s) to process"

# Load jobId -> fully-qualified Nextflow process name for inline joining
# (kept in memory so each row can be written straight to $OUTFILE as it's
# computed, instead of buffering the whole run in a temp file — lets you
# tail -f the CSV to watch progress on long runs).
declare -A PROC_MAP
declare -A TAG_MAP
while IFS=',' read -r jid proc tag; do
    PROC_MAP["$jid"]="$proc"
    # tag is the last field read, so `read` keeps any internal commas intact
    # (multi-sample tags like the samplesheet job legitimately contain them).
    TAG_MAP["$jid"]="$tag"
done < "$JOBID_META_MAP"
rm -f "$JOBID_META_MAP"

# sacct occasionally hiccups under rapid repeated calls (connection reset,
# transient empty result) - worth retrying. But job_usage.sh also exits 1
# deterministically for jobs in a state it doesn't recognise as finished
# (e.g. NODE_FAIL isn't in its COMPLETED|FAILED|CANCELLED|OUT_OF_MEM|TIMEOUT
# check) - retrying that is pointless since the state won't change, so skip
# straight away and just record it.
readonly JOB_USAGE_MAX_ATTEMPTS=3
readonly JOB_USAGE_RETRY_DELAY=3
run_job_usage() {
    local jobid="$1" attempt=1 out
    while [ "$attempt" -le "$JOB_USAGE_MAX_ATTEMPTS" ]; do
        if out=$(bash "$SCRIPT_DIR/job_usage.sh" "$jobid" --format csv --no-csv-header 2>&1); then
            printf '%s' "$out"
            return 0
        fi
        if [[ "$out" == *"not yet in a finished state"* ]]; then
            echo "  job_usage.sh cannot process job $jobid (unrecognised terminal state), skipping: $out" >&2
            return 1
        fi
        echo "  Warning: job_usage.sh failed for job $jobid (attempt $attempt/$JOB_USAGE_MAX_ATTEMPTS): $out" >&2
        sleep "$JOB_USAGE_RETRY_DELAY"
        attempt=$((attempt + 1))
    done
    return 1
}

# Header is fixed (schema of vendored job_usage.sh, see README) rather than
# taken from the first job's output, so a failed/skipped first job can't
# leave the CSV headerless.
echo "generated_at,job_id,project,partition,exit_status,job_state,nodes_requested,gcds_requested,ncpus_requested,ncpus_allocated,ncpus_allocated_raw,cpu_time_available,cpu_time_available_s,cpu_time_used,cpu_time_used_s,memory_requested,memory_requested_gb,memory_used,memory_used_gb,walltime_requested,walltime_used,walltime_requested_h,walltime_used_h,walltime_efficiency_pct,cpu_efficiency_pct,memory_efficiency_pct,service_units,job_submitted,job_started,job_ended,process,tag" > "$OUTFILE"

FAILED_JOBS_FILE=$(mktemp)
for jobid in "${JOB_IDS[@]}"; do
    echo "Processing job $jobid..."
    proc="${PROC_MAP[$jobid]:-UNKNOWN}"
    tag="${TAG_MAP[$jobid]:-}"
    if ! data=$(run_job_usage "$jobid"); then
        echo "  Skipping job $jobid" >&2
        echo "$jobid" >> "$FAILED_JOBS_FILE"
        continue
    fi
    # Double any embedded quotes so the quoted CSV cells stay well-formed.
    printf '%s,"%s","%s"\n' "$data" "${proc//\"/\"\"}" "${tag//\"/\"\"}" >> "$OUTFILE"
done

n_failed=$(wc -l < "$FAILED_JOBS_FILE")
if [ "$n_failed" -gt 0 ]; then
    echo ""
    echo "Warning: $n_failed job(s) could not be looked up after retries and were skipped:"
    cat "$FAILED_JOBS_FILE" >&2
fi
rm -f "$FAILED_JOBS_FILE"

# Locate columns by header name (robust to future column reordering)
REPORT_FILE="${OUTFILE%.csv}.report.txt"

SUMMARY=$(awk -F',' '
    function unquote(s) { gsub(/^"|"$/,"",s); return s }
    # Sorts arr in place (1..n) and returns it; used so min/median/max/IQR
    # share one sort instead of re-sorting per stat.
    function sort_arr(arr, n,    i, j, tmp) {
        for (i=1;i<=n;i++) for (j=i+1;j<=n;j++) if (arr[j]<arr[i]) { tmp=arr[i]; arr[i]=arr[j]; arr[j]=tmp }
    }
    # Linear-interpolation percentile (like numpy default), arr must be pre-sorted 1..n
    function pct(arr, n, p,    rank, lo, hi, frac) {
        if (n==0) return 0
        if (n==1) return arr[1]
        rank = p/100.0 * (n-1) + 1
        lo = int(rank); hi = lo+1
        frac = rank - lo
        if (hi > n) return arr[n]
        return arr[lo] + frac * (arr[hi]-arr[lo])
    }
    # Copies a per-process 2D-array slice arr[process,1..n] into a flat 1..n array
    function copy_slice(src, process, n, dest,    i) {
        for (i=1;i<=n;i++) dest[i] = src[process SUBSEP i]
    }
    NR==1 {
        for (i=1;i<=NF;i++) {
            h=unquote($i)
            if (h=="service_units") su_col=i
            else if (h=="cpu_efficiency_pct") cpu_col=i
            else if (h=="memory_efficiency_pct") mem_col=i
            else if (h=="walltime_efficiency_pct") wt_col=i
            else if (h=="job_state") state_col=i
            else if (h=="process") proc_col=i
        }
        next
    }
    {
        n_jobs++
        su = unquote($su_col) + 0; total_su += su; su_arr[n_jobs]=su
        state = unquote($state_col)
        state_count[state]++
        if (state != "COMPLETED") failed_su += su

        cpu = unquote($cpu_col); if (cpu != "") { n_cpu++; cpu_arr[n_cpu]=cpu+0 }
        mem = unquote($mem_col); if (mem != "") { n_mem++; mem_arr[n_mem]=mem+0 }
        wt  = unquote($wt_col);  if (wt  != "") { n_wt++;  wt_arr[n_wt]=wt+0 }

        # Per-process breakdown, keyed by short process name (last :-segment)
        full_proc = unquote($proc_col)
        nseg = split(full_proc, segs, ":")
        proc = segs[nseg]
        if (!(proc in proc_seen)) { proc_seen[proc]=1; proc_order[++n_procs]=proc }
        proc_count[proc]++
        proc_total_su[proc] += su
        if (cpu != "") { proc_cpu_n[proc]++; proc_cpu_arr[proc, proc_cpu_n[proc]] = cpu+0 }
        if (mem != "") { proc_mem_n[proc]++; proc_mem_arr[proc, proc_mem_n[proc]] = mem+0 }
        if (wt  != "") { proc_wt_n[proc]++;  proc_wt_arr[proc, proc_wt_n[proc]] = wt+0 }
    }
    END {
        sort_arr(su_arr, n_jobs); sort_arr(cpu_arr, n_cpu); sort_arr(mem_arr, n_mem); sort_arr(wt_arr, n_wt)

        printf "n_jobs=%d\n", n_jobs
        printf "total_su=%.2f\n", total_su
        printf "failed_su=%.2f\n", failed_su
        printf "median_su=%.4f\n", pct(su_arr, n_jobs, 50)

        printf "cpu_min=%.2f\n", (n_cpu>0 ? cpu_arr[1] : 0)
        printf "cpu_q1=%.2f\n", pct(cpu_arr, n_cpu, 25)
        printf "cpu_median=%.2f\n", pct(cpu_arr, n_cpu, 50)
        printf "cpu_q3=%.2f\n", pct(cpu_arr, n_cpu, 75)
        printf "cpu_max=%.2f\n", (n_cpu>0 ? cpu_arr[n_cpu] : 0)

        printf "mem_min=%.2f\n", (n_mem>0 ? mem_arr[1] : 0)
        printf "mem_q1=%.2f\n", pct(mem_arr, n_mem, 25)
        printf "mem_median=%.2f\n", pct(mem_arr, n_mem, 50)
        printf "mem_q3=%.2f\n", pct(mem_arr, n_mem, 75)
        printf "mem_max=%.2f\n", (n_mem>0 ? mem_arr[n_mem] : 0)

        printf "wt_min=%.2f\n", (n_wt>0 ? wt_arr[1] : 0)
        printf "wt_q1=%.2f\n", pct(wt_arr, n_wt, 25)
        printf "wt_median=%.2f\n", pct(wt_arr, n_wt, 50)
        printf "wt_q3=%.2f\n", pct(wt_arr, n_wt, 75)
        printf "wt_max=%.2f\n", (n_wt>0 ? wt_arr[n_wt] : 0)
        for (s in state_count) printf "state:%s=%d\n", s, state_count[s]

        # Per-process rows, sorted by total SU descending (biggest cost first)
        for (i=1;i<=n_procs;i++) proc_sorted[i]=proc_order[i]
        for (i=1;i<=n_procs;i++) for (j=i+1;j<=n_procs;j++)
            if (proc_total_su[proc_sorted[j]] > proc_total_su[proc_sorted[i]]) {
                tmp=proc_sorted[i]; proc_sorted[i]=proc_sorted[j]; proc_sorted[j]=tmp
            }
        proc_order_line = ""
        for (i=1;i<=n_procs;i++) {
            p = proc_sorted[i]
            proc_order_line = (i==1) ? p : proc_order_line "\x1f" p
            copy_slice(proc_cpu_arr, p, proc_cpu_n[p], tmp_arr); sort_arr(tmp_arr, proc_cpu_n[p])
            proc_cpu_med = pct(tmp_arr, proc_cpu_n[p], 50); delete tmp_arr
            copy_slice(proc_mem_arr, p, proc_mem_n[p], tmp_arr); sort_arr(tmp_arr, proc_mem_n[p])
            proc_mem_med = pct(tmp_arr, proc_mem_n[p], 50); delete tmp_arr
            copy_slice(proc_wt_arr, p, proc_wt_n[p], tmp_arr); sort_arr(tmp_arr, proc_wt_n[p])
            proc_wt_med = pct(tmp_arr, proc_wt_n[p], 50); delete tmp_arr
            printf "proc:%s=%d|%.2f|%.2f|%.2f|%.2f\n", p, proc_count[p], proc_total_su[p], proc_cpu_med, proc_mem_med, proc_wt_med
        }
        printf "proc_order=%s\n", proc_order_line
    }
' "$OUTFILE")

# Parse the awk summary into shell vars
declare -A STATS
while IFS='=' read -r key val; do
    [ -z "$key" ] && continue
    STATS["$key"]="$val"
done <<< "$SUMMARY"

eff_row() {
    # eff_row <label> <prefix>  -> one formatted efficiency-stat table row
    local label="$1" prefix="$2"
    printf "   %-12s %8s%%  %8s%%  %8s%%  %8s%%  %8s%%\n" \
        "$label" "${STATS[${prefix}_min]}" "${STATS[${prefix}_q1]}" \
        "${STATS[${prefix}_median]}" "${STATS[${prefix}_q3]}" "${STATS[${prefix}_max]}"
}

{
    echo "======================================================================================"
    echo "Nextflow workflow compute usage report"
    echo "======================================================================================"
    printf "Generated:    %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    printf "Run source:   %s\n" "$RUN_DIR"
    printf "Per-job CSV:  %s\n" "$OUTFILE"
    echo ""
    echo "JOB SUMMARY"
    echo "--------------------------------------------------------------------------------------"
    printf "   Total Jobs:              %s\n" "${STATS[n_jobs]}"
    printf "   Total Service Units:     %s\n" "${STATS[total_su]}"
    printf "   Median SU per job:       %s\n" "${STATS[median_su]}"
    printf "   SU from non-COMPLETED jobs: %s\n" "${STATS[failed_su]}"
    echo ""
    echo "   Job states:"
    for key in "${!STATS[@]}"; do
        [[ "$key" == state:* ]] && printf "     %-20s %s\n" "${key#state:}" "${STATS[$key]}"
    done
    echo ""
    echo "EFFICIENCY METRICS (across all ${STATS[n_jobs]} jobs)"
    echo "--------------------------------------------------------------------------------------"
    printf "   %-12s %9s   %9s   %9s   %9s   %9s\n" "" "Min" "Q1" "Median" "Q3" "Max"
    eff_row "CPU"      "cpu"
    eff_row "Memory"   "mem"
    eff_row "Walltime" "wt"
    echo ""
    echo "PER-PROCESS BREAKDOWN (sorted by total SU, descending; efficiencies are per-process medians)"
    echo "--------------------------------------------------------------------------------------"
    printf "   %-30s %6s   %10s   %8s   %8s   %8s\n" "Process" "Jobs" "Total SU" "CPU%" "Mem%" "Wall%"
    IFS=$'\x1f' read -r -a PROC_ORDER <<< "${STATS[proc_order]:-}"
    for p in "${PROC_ORDER[@]}"; do
        [ -z "$p" ] && continue
        IFS='|' read -r p_count p_su p_cpu p_mem p_wt <<< "${STATS[proc:$p]}"
        printf "   %-30s %6s   %10s   %7s%%   %7s%%   %7s%%\n" "$p" "$p_count" "$p_su" "$p_cpu" "$p_mem" "$p_wt"
    done
    echo "======================================================================================"
} | tee "$REPORT_FILE"

# ---------------------------------------------------------------------------
# Cost per genome (sample) rollup
# ---------------------------------------------------------------------------
# Group per-job SU by OceanGenomes sample ID. A sample tag is "OG" followed by
# digits, optionally with a suffix (e.g. OG1422.hifi.260227.v10oatk -> OG1422).
# Jobs whose tag does not start with OG (DB downloads, the samplesheet job that
# spans many samples, MultiQC, untagged) are shared overhead, amortised evenly
# across the distinct samples in the run so a genome's cost is honest whether it
# ran solo or in a batch.
COST_FILE="$(dirname "$OUTFILE")/cost_per_sample.csv"

COST_REPORT=$(awk -F',' -v costfile="$COST_FILE" '
    function unquote(s) { gsub(/^"|"$/,"",s); return s }
    NR==1 {
        for (i=1;i<=NF;i++) {
            h=unquote($i)
            if (h=="service_units") su_col=i
            else if (h=="job_state") state_col=i
            else if (h=="tag") tag_col=i
        }
        next
    }
    {
        su = unquote($su_col) + 0
        state = unquote($state_col)
        # tag is the final column; a multi-sample tag may contain commas and
        # spill into extra fields, but its leading chunk sits at tag_col and such
        # tags never start with OG, so they still fall through to overhead.
        tag = unquote($tag_col)
        if (match(tag, /^OG[0-9]+/)) {
            og = substr(tag, RSTART, RLENGTH)
            if (!(og in seen)) { seen[og]=1; order[++n_samples]=og }
            njobs[og]++
            direct[og] += su
            if (state != "COMPLETED") failed[og]++
        } else {
            overhead_su += su
            overhead_jobs++
        }
    }
    END {
        share = (n_samples > 0) ? overhead_su / n_samples : 0

        # Sort samples by direct SU descending (biggest cost first)
        for (i=1;i<=n_samples;i++) srt[i]=order[i]
        for (i=1;i<=n_samples;i++) for (j=i+1;j<=n_samples;j++)
            if (direct[srt[j]] > direct[srt[i]]) { t=srt[i]; srt[i]=srt[j]; srt[j]=t }

        print "sample_id,n_jobs,direct_su,overhead_su_share,total_su,n_failed_jobs" > costfile
        for (i=1;i<=n_samples;i++) {
            og = srt[i]
            printf "%s,%d,%.4f,%.4f,%.4f,%d\n", og, njobs[og], direct[og], share, direct[og]+share, (og in failed ? failed[og] : 0) > costfile
        }
        # If nothing was sample-tagged, keep the overhead SU visible rather than
        # silently dropping it.
        if (n_samples == 0 && overhead_jobs > 0)
            printf "_UNATTRIBUTED,%d,%.4f,%.4f,%.4f,0\n", overhead_jobs, overhead_su, 0, overhead_su > costfile

        # Machine-readable lines echoed back to the shell to build the report.
        printf "n_samples=%d\n", n_samples
        printf "overhead_su=%.4f\n", overhead_su
        printf "overhead_jobs=%d\n", overhead_jobs
        printf "overhead_share=%.4f\n", share
        for (i=1;i<=n_samples;i++) {
            og = srt[i]
            printf "row=%s|%d|%.4f|%.4f|%.4f|%d\n", og, njobs[og], direct[og], share, direct[og]+share, (og in failed ? failed[og] : 0)
        }
    }
' "$OUTFILE")

# Render the cost-per-genome table into the report and stdout.
{
    echo ""
    echo "======================================================================================"
    echo "COST PER GENOME (Service Units; shared overhead amortised across samples)"
    echo "--------------------------------------------------------------------------------------"
    cps_n=$(sed -n 's/^n_samples=//p' <<< "$COST_REPORT")
    cps_oh_su=$(sed -n 's/^overhead_su=//p' <<< "$COST_REPORT")
    cps_oh_jobs=$(sed -n 's/^overhead_jobs=//p' <<< "$COST_REPORT")
    cps_oh_share=$(sed -n 's/^overhead_share=//p' <<< "$COST_REPORT")
    printf "   Samples: %s    Shared overhead: %s SU across %s job(s) => %s SU/sample\n" \
        "${cps_n:-0}" "${cps_oh_su:-0}" "${cps_oh_jobs:-0}" "${cps_oh_share:-0}"
    echo ""
    printf "   %-24s %6s   %12s   %12s   %12s   %8s\n" "Sample" "Jobs" "Direct SU" "Overhead SU" "Total SU" "Failed"
    while IFS='|' read -r c_og c_njobs c_direct c_share c_tot c_failed; do
        [ -z "$c_og" ] && continue
        printf "   %-24s %6s   %12s   %12s   %12s   %8s\n" "$c_og" "$c_njobs" "$c_direct" "$c_share" "$c_tot" "$c_failed"
    done < <(sed -n 's/^row=//p' <<< "$COST_REPORT")
    echo "======================================================================================"
} | tee -a "$REPORT_FILE"

echo ""
echo "Wrote per-job usage to: $OUTFILE"
echo "Wrote summary report to: $REPORT_FILE"
echo "Wrote cost per sample to: $COST_FILE"
