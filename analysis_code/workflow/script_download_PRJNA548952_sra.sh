#!/usr/bin/env bash
set -Eeuo pipefail

readonly PROJECT="PRJNA548952"
readonly DEFAULT_OUTPUT_DIR="${SRA_ROOT:-${PWD}/data}/PRJNA548952"
readonly DEFAULT_TOOLKIT_BIN="${PRDX4T_ENV_BIN:-${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin}"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly MANIFEST="${SCRIPT_DIR}/datasets/${PROJECT}/runs.tsv"

output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
selected_run=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [--output-dir DIR] [--run SRR_ACCESSION]

Download the 12 SRA archives in ${PROJECT}. Re-running is safe: prefetch
resumes partial downloads, and archives that pass vdb-validate are skipped.

Options:
  --output-dir DIR  Dataset root (default: ${DEFAULT_OUTPUT_DIR})
  --run SRR...      Download only one run from the project manifest
  -h, --help        Show this help

Environment:
  OUTPUT_DIR        Alternative dataset root
  SRA_TOOLKIT_BIN   Directory containing prefetch and vdb-validate
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --output-dir)
            (($# >= 2)) || die "--output-dir requires a value"
            output_dir="$2"
            shift 2
            ;;
        --run)
            (($# >= 2)) || die "--run requires an accession"
            selected_run="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1 (use --help)"
            ;;
    esac
done

[[ -r "$MANIFEST" ]] || die "manifest not found: $MANIFEST"

toolkit_bin="${SRA_TOOLKIT_BIN:-}"
if [[ -n "$toolkit_bin" ]]; then
    :
elif command -v prefetch >/dev/null 2>&1 && command -v vdb-validate >/dev/null 2>&1; then
    toolkit_bin="$(dirname -- "$(command -v prefetch)")"
elif [[ -x "${DEFAULT_TOOLKIT_BIN}/prefetch" && -x "${DEFAULT_TOOLKIT_BIN}/vdb-validate" ]]; then
    toolkit_bin="$DEFAULT_TOOLKIT_BIN"
else
    die "SRA Toolkit was not found. Set SRA_TOOLKIT_BIN or add its bin directory to PATH."
fi

readonly PREFETCH="${toolkit_bin}/prefetch"
readonly VDB_VALIDATE="${toolkit_bin}/vdb-validate"
[[ -x "$PREFETCH" ]] || die "prefetch is not executable: $PREFETCH"
[[ -x "$VDB_VALIDATE" ]] || die "vdb-validate is not executable: $VDB_VALIDATE"

if [[ -n "$selected_run" ]] && ! awk -F '\t' -v run="$selected_run" 'NR > 1 && $1 == run { found=1 } END { exit !found }' "$MANIFEST"; then
    die "$selected_run is not listed in $MANIFEST"
fi

readonly SRA_DIR="${output_dir}/sra"
readonly METADATA_DIR="${output_dir}/metadata"
mkdir -p "$SRA_DIR" "$METADATA_DIR"
cp -f "$MANIFEST" "${METADATA_DIR}/${PROJECT}_runs.tsv"

printf 'Project:     %s\n' "$PROJECT"
printf 'Toolkit:     %s\n' "$("$PREFETCH" --version | tail -n 1)"
printf 'SRA output:  %s\n' "$SRA_DIR"
df -h "$output_dir" | tail -n 1

failures=()
while IFS=$'\t' read -r run sample library layout spots bases size_mb _metadata; do
    [[ "$run" == "Run" ]] && continue
    [[ -z "$selected_run" || "$run" == "$selected_run" ]] || continue

    archive="${SRA_DIR}/${run}/${run}.sra"
    printf '\n[%s] %s / %s (reported SRA size: %s MB)\n' "$run" "$sample" "$library" "$size_mb"

    # This loop reads the manifest on stdin, so every child gets /dev/null
    # instead; a tool that consumed stdin would otherwise swallow run records.
    if [[ -s "$archive" ]] && "$VDB_VALIDATE" "$archive" >/dev/null 2>&1 </dev/null; then
        printf '[%s] Existing archive is valid; skipping.\n' "$run"
        continue
    fi

    if ! "$PREFETCH" --max-size u --progress --output-directory "$SRA_DIR" "$run" </dev/null; then
        printf '[%s] Download failed; continuing with the remaining runs.\n' "$run" >&2
        failures+=("$run")
        continue
    fi

    if ! "$VDB_VALIDATE" "${SRA_DIR}/${run}" </dev/null; then
        printf '[%s] Validation failed. Leave the files in place so prefetch can resume/repair them.\n' "$run" >&2
        failures+=("$run")
        continue
    fi
done < "$MANIFEST"

if ((${#failures[@]})); then
    printf '\nFailed or invalid runs: %s\n' "${failures[*]}" >&2
    printf 'Re-run this script; successful runs will be skipped.\n' >&2
    exit 1
fi

printf '\nAll requested SRA archives downloaded and validated.\n'
