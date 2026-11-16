#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly CONFIG_HELPER="${SCRIPT_DIR}/lib/load_dataset_config.sh"
[[ -r "$CONFIG_HELPER" ]] || { printf 'ERROR: config helper not found: %s\n' "$CONFIG_HELPER" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_HELPER"
dataset_selector="$(dataset_selector_from_args "$@")" || exit 1
load_dataset_config "$dataset_selector" "$SCRIPT_DIR" || exit 1
readonly PROJECT="$DATASET_ID"
readonly RUN_MANIFEST
readonly WORKER="${SCRIPT_DIR}/script_workflow_03_process_one_library.sh"

stop_on_error=0
worker_args=()

usage() {
    cat <<EOF
Usage: $(basename "$0") [--dataset ID] [--stop-on-error] [worker options...]

Run script_workflow_03_process_one_library.sh over every run in
$(basename "$RUN_MANIFEST"), in manifest order.

Libraries that already have a complete final BAM are skipped by the worker, so
re-running this script resumes an interrupted batch without redoing work.

Options:
  --dataset ID      Dataset ID or dataset.conf path (default: PRJNA548952)
  --stop-on-error   Abort at the first failing library
                    (default: process the rest, then report failures)
  -h, --help        Show this help

Any other option is passed through to the worker unchanged, for example:
  $(basename "$0") --threads 8
  $(basename "$0") --dry-run
  $(basename "$0") --dataset-root /path/to/dataset

Exit status is non-zero if any library failed.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --dataset)
            (($# >= 2)) || die "--dataset requires an ID or config file"
            shift 2
            ;;
        --stop-on-error)
            stop_on_error=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --run)
            die "--run is a single-library option; use $(basename "$WORKER") directly"
            ;;
        *)
            worker_args+=("$1")
            shift
            ;;
    esac
done

[[ -r "$RUN_MANIFEST" ]] || die "run manifest not found: $RUN_MANIFEST"
[[ -x "$WORKER" ]] || die "worker script is not executable: $WORKER"

runs=()
while IFS=$'\t' read -r run _rest; do
    [[ "$run" == "Run" || -z "$run" ]] && continue
    runs+=("$run")
done < "$RUN_MANIFEST"

((${#runs[@]})) || die "no runs found in $RUN_MANIFEST"

printf 'Project:   %s\n' "$PROJECT"
printf 'Libraries: %s\n' "${#runs[@]}"
printf 'Worker:    %s\n' "$WORKER"
if ((${#worker_args[@]})); then
    printf 'Passing:   %s\n' "${worker_args[*]}"
fi

succeeded=()
failed=()
for index in "${!runs[@]}"; do
    run="${runs[$index]}"
    printf '\n===== [%s/%s] %s =====\n' "$((index + 1))" "${#runs[@]}" "$run"

    # The worker gets /dev/null on stdin so nothing it spawns can consume the
    # remainder of this loop's input.
    if "$WORKER" --dataset "$DATASET_CONFIG" --run "$run" ${worker_args[@]+"${worker_args[@]}"} </dev/null; then
        succeeded+=("$run")
    else
        status=$?
        failed+=("$run")
        printf '\n[%s] FAILED with status %s\n' "$run" "$status" >&2
        if ((stop_on_error)); then
            printf 'Stopping because --stop-on-error was given.\n' >&2
            break
        fi
        printf 'Continuing with the remaining libraries.\n' >&2
    fi
done

printf '\n===== batch summary =====\n'
printf 'Succeeded: %s of %s\n' "${#succeeded[@]}" "${#runs[@]}"
if ((${#failed[@]})); then
    printf 'Failed:    %s\n' "${failed[*]}" >&2
    printf 'Re-run this script to retry; completed libraries are skipped.\n' >&2
    exit 1
fi

printf 'All libraries processed.\n'
printf 'Next step: script_workflow_04_quantify_PRDX4t_isoform_usage.sh --dataset %s\n' "$PROJECT"
