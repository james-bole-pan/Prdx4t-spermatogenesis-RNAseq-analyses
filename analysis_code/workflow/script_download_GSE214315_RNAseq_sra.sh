#!/usr/bin/env bash
set -Eeuo pipefail

readonly GEO_ACCESSION="GSE214315"
readonly RNA_SUBSERIES="GSE214315"
readonly BIOPROJECT="PRJNA884732"
readonly RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${BIOPROJECT}"
readonly DEFAULT_OUTPUT_DIR="${SRA_ROOT:-${PWD}/data}/GSE214315"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly WORKFLOW_MANIFEST="${SCRIPT_DIR}/datasets/${GEO_ACCESSION}/runs.tsv"

output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
selected_run=""
readonly STRATEGY="RNA-Seq"
metadata_only=0
refresh_metadata=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Download the 32 raw RNA-seq SRA archives in ${GEO_ACCESSION}. These runs belong
to the RNA-seq SubSeries ${RNA_SUBSERIES} and its SRA BioProject ${BIOPROJECT};
the three H3K4me3 ChIP-seq runs in the SuperSeries are deliberately excluded.

The 32 runs represent 17 GEO RNA-seq samples. Fifteen samples have two
technical runs, and the final two samples are single-end. Files and metadata
are therefore retained by SRR accession.

Options:
  --output-dir DIR    Dataset root (default: ${DEFAULT_OUTPUT_DIR})
  --run SRR...        Download one run only
  --metadata-only     Prepare the SRA run manifest without downloading archives
  --refresh-metadata  Replace the cached NCBI SRA RunInfo
  -h, --help          Show this help

Environment:
  OUTPUT_DIR          Alternative dataset root
  SRA_TOOLKIT_BIN     Directory containing prefetch and vdb-validate
  PRDX4T_ENV_BIN      Alternative prdx4t-analysis executable directory

Outputs:
  OUTPUT_DIR/metadata/SraRunInfo.csv
  OUTPUT_DIR/metadata/${GEO_ACCESSION}_runs.tsv
  OUTPUT_DIR/sra/SRR.../SRR....sra

Re-running is safe. Archives that pass vdb-validate are skipped, and prefetch
can resume interrupted downloads.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --output-dir)
            (($# >= 2)) || die "--output-dir requires a directory"
            output_dir="$2"
            shift 2
            ;;
        --run)
            (($# >= 2)) || die "--run requires an accession"
            selected_run="$2"
            shift 2
            ;;
        --metadata-only)
            metadata_only=1
            shift
            ;;
        --refresh-metadata)
            refresh_metadata=1
            shift
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

if [[ -n "$selected_run" ]]; then
    [[ "$selected_run" =~ ^SRR[0-9]+$ ]] || die "invalid SRA run accession: $selected_run"
fi

command -v curl >/dev/null 2>&1 || die "curl is required"
[[ -r "$WORKFLOW_MANIFEST" ]] || die "workflow manifest not found: $WORKFLOW_MANIFEST"
python_bin=""
if [[ -n "${PRDX4T_ENV_BIN:-}" && -x "${PRDX4T_ENV_BIN}/python" ]]; then
    python_bin="${PRDX4T_ENV_BIN}/python"
elif command -v python3 >/dev/null 2>&1; then
    python_bin="$(command -v python3)"
elif [[ -x "${DEFAULT_ENV_BIN}/python" ]]; then
    python_bin="${DEFAULT_ENV_BIN}/python"
else
    die "Python 3 was not found; run script_environment_create_prdx4t_analysis.sh first"
fi
readonly PYTHON="$python_bin"

readonly METADATA_DIR="${output_dir}/metadata"
readonly SRA_DIR="${output_dir}/sra"
readonly RUNINFO="${METADATA_DIR}/SraRunInfo.csv"
readonly RUN_MANIFEST="${METADATA_DIR}/${GEO_ACCESSION}_runs.tsv"
mkdir -p "$METADATA_DIR"

download_metadata() {
    local url="$1"
    local destination="$2"
    local label="$3"
    local partial="${destination}.part"

    if ((refresh_metadata == 0)) && [[ -s "$destination" ]]; then
        printf 'Using cached %s: %s\n' "$label" "$destination"
        return
    fi

    rm -f -- "$partial"
    printf 'Downloading %s\n' "$label"
    curl --fail --location --retry 5 --retry-delay 5 \
        --output "$partial" "$url"
    [[ -s "$partial" ]] || die "$label download was empty: $url"
    mv -f "$partial" "$destination"
}

download_metadata "$RUNINFO_URL" "$RUNINFO" "NCBI SRA RunInfo"

manifest_tmp="${RUN_MANIFEST}.tmp"
"$PYTHON" - "$RUNINFO" "$manifest_tmp" "$BIOPROJECT" "$GEO_ACCESSION" \
    "$STRATEGY" "$selected_run" <<'PY'
import csv
import sys
from pathlib import Path

(
    runinfo_path,
    output_path,
    bioproject,
    geo_accession,
    wanted_strategy,
    wanted_run,
) = sys.argv[1:]


def clean(value):
    return (value or "").replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


required = {
    "Run",
    "spots",
    "bases",
    "size_MB",
    "Experiment",
    "LibraryStrategy",
    "LibraryLayout",
    "BioProject",
    "BioSample",
    "ScientificName",
    "SampleName",
}
with Path(runinfo_path).open(newline="") as handle:
    reader = csv.DictReader(handle)
    missing = sorted(required.difference(reader.fieldnames or []))
    if missing:
        raise SystemExit(f"ERROR: RunInfo is missing columns: {', '.join(missing)}")
    run_rows = list(reader)

if not run_rows:
    raise SystemExit(f"ERROR: NCBI returned no runs for {bioproject}")

seen_runs = set()
rows = []
for record in run_rows:
    run = clean(record["Run"])
    if not run.startswith("SRR") or not run[3:].isdigit():
        raise SystemExit(f"ERROR: invalid run accession in RunInfo: {run!r}")
    if run in seen_runs:
        raise SystemExit(f"ERROR: duplicate run in RunInfo: {run}")
    seen_runs.add(run)

    strategy = clean(record["LibraryStrategy"])
    if strategy not in {"RNA-Seq", "ChIP-Seq"}:
        raise SystemExit(f"ERROR: unexpected library strategy for {run}: {strategy}")
    if strategy != wanted_strategy:
        continue
    if clean(record["BioProject"]) != bioproject:
        raise SystemExit(f"ERROR: {run} is not assigned to {bioproject}")
    if wanted_run and run != wanted_run:
        continue

    gsm = clean(record["SampleName"])
    if not gsm.startswith("GSM") or not gsm[3:].isdigit():
        raise SystemExit(f"ERROR: invalid GEO sample accession for {run}: {gsm!r}")

    for integer_field in ("spots", "bases"):
        value = clean(record[integer_field])
        if not value.isdigit() or int(value) <= 0:
            raise SystemExit(f"ERROR: invalid {integer_field} for {run}: {value!r}")
    try:
        size_mb = float(clean(record["size_MB"]))
    except ValueError as exc:
        raise SystemExit(f"ERROR: invalid size_MB for {run}") from exc
    if size_mb <= 0:
        raise SystemExit(f"ERROR: non-positive size_MB for {run}")

    rows.append(
        {
            "Run": run,
            "GEO_sample": gsm,
            "LibraryStrategy": strategy,
            "LibraryLayout": clean(record["LibraryLayout"]),
            "spots": clean(record["spots"]),
            "bases": clean(record["bases"]),
            "size_MB": clean(record["size_MB"]),
            "Experiment": clean(record["Experiment"]),
            "BioSample": clean(record["BioSample"]),
            "ScientificName": clean(record["ScientificName"]),
        }
    )

if wanted_run and not rows:
    if wanted_run not in seen_runs:
        raise SystemExit(f"ERROR: {wanted_run} is not part of {geo_accession}/{bioproject}")
    raise SystemExit(f"ERROR: {wanted_run} is excluded by the requested sequencing subset")
if not rows:
    raise SystemExit("ERROR: no RNA-seq runs were found")

fieldnames = [
    "Run",
    "GEO_sample",
    "LibraryStrategy",
    "LibraryLayout",
    "spots",
    "bases",
    "size_MB",
    "Experiment",
    "BioSample",
    "ScientificName",
]
with Path(output_path).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY
mv -f "$manifest_tmp" "$RUN_MANIFEST"

# The downloaded NCBI metadata is authoritative for archive properties; the
# version-controlled workflow manifest adds biological group/replicate labels.
# Require exact agreement on run, layout, and spot count before downloading.
"$PYTHON" - "$RUN_MANIFEST" "$WORKFLOW_MANIFEST" "$selected_run" <<'PY'
import csv
import sys

downloaded_path, workflow_path, selected_run = sys.argv[1:]

def load(path):
    with open(path, newline="") as handle:
        rows = {row["Run"]: row for row in csv.DictReader(handle, delimiter="\t")}
    return rows

downloaded = load(downloaded_path)
workflow = load(workflow_path)
if selected_run:
    workflow = {selected_run: workflow[selected_run]} if selected_run in workflow else {}
if set(downloaded) != set(workflow):
    missing = sorted(set(downloaded) - set(workflow))
    extra = sorted(set(workflow) - set(downloaded))
    raise SystemExit(
        f"ERROR: workflow/NCBI run mismatch; absent from workflow={missing}, absent from NCBI={extra}"
    )
for run in sorted(downloaded):
    for column in ("LibraryLayout", "spots"):
        if downloaded[run][column].strip() != workflow[run][column].strip():
            raise SystemExit(
                f"ERROR: {run} {column} differs: NCBI={downloaded[run][column]!r}, "
                f"workflow={workflow[run][column]!r}"
            )
PY

run_count="$(awk 'END {print NR-1}' "$RUN_MANIFEST")"
sample_count="$(awk -F '\t' 'NR > 1 {seen[$2]=1} END {print length(seen)}' "$RUN_MANIFEST")"
reported_gib="$(awk -F '\t' 'NR > 1 {mb += $7} END {printf "%.1f", mb/1024}' "$RUN_MANIFEST")"

printf '\nAccession:       %s\n' "$GEO_ACCESSION"
printf 'BioProject:      %s\n' "$BIOPROJECT"
printf 'RNA SubSeries:   %s\n' "$RNA_SUBSERIES"
printf 'GEO samples:     %s\n' "$sample_count"
printf 'RNA-seq runs:    %s\n' "$run_count"
printf 'Reported size:   %s GiB compressed\n' "$reported_gib"
printf 'Run manifest:    %s\n' "$RUN_MANIFEST"

if ((metadata_only)); then
    printf '\nMetadata preparation complete; no SRA archives were downloaded.\n'
    exit 0
fi

toolkit_bin="${SRA_TOOLKIT_BIN:-${PRDX4T_ENV_BIN:-}}"
if [[ -n "$toolkit_bin" ]]; then
    :
elif command -v prefetch >/dev/null 2>&1 && command -v vdb-validate >/dev/null 2>&1; then
    toolkit_bin="$(dirname -- "$(command -v prefetch)")"
elif [[ -x "${DEFAULT_ENV_BIN}/prefetch" && -x "${DEFAULT_ENV_BIN}/vdb-validate" ]]; then
    toolkit_bin="$DEFAULT_ENV_BIN"
else
    die "SRA Toolkit was not found; run script_environment_create_prdx4t_analysis.sh or set SRA_TOOLKIT_BIN"
fi

readonly PREFETCH="${toolkit_bin}/prefetch"
readonly VDB_VALIDATE="${toolkit_bin}/vdb-validate"
[[ -x "$PREFETCH" ]] || die "prefetch is not executable: $PREFETCH"
[[ -x "$VDB_VALIDATE" ]] || die "vdb-validate is not executable: $VDB_VALIDATE"

mkdir -p "$SRA_DIR"
available_gib="$(df -Pk "$output_dir" | awk 'NR == 2 {printf "%d", $4/1048576}')"
recommended_gib="$(awk -v size="$reported_gib" 'BEGIN {printf "%d", size * 1.10 + 2.999}')"
printf 'SRA Toolkit:     %s\n' "$("$PREFETCH" --version | awk 'NF {line=$0} END {print line}')"
printf 'SRA directory:   %s\n' "$SRA_DIR"
printf 'Available disk:  %s GiB\n' "$available_gib"
if ((available_gib < recommended_gib)); then
    die "insufficient free space: about ${recommended_gib} GiB is recommended for this selection"
fi

failures=()
while IFS=$'\t' read -r run gsm library_strategy layout spots bases size_mb experiment biosample species; do
    [[ "$run" == "Run" ]] && continue
    archive="${SRA_DIR}/${run}/${run}.sra"
    printf '\n[%s] %s | %s | %s (reported %.1f GiB)\n' \
        "$run" "$gsm" "$library_strategy" "$layout" \
        "$(awk -v mb="$size_mb" 'BEGIN {printf "%.1f", mb/1024}')"

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
        printf '[%s] Validation failed; retained files can be repaired by re-running.\n' "$run" >&2
        failures+=("$run")
        continue
    fi
done < "$RUN_MANIFEST"

if ((${#failures[@]})); then
    printf '\nFailed or invalid runs: %s\n' "${failures[*]}" >&2
    printf 'Re-run the same command; successful runs will be skipped.\n' >&2
    exit 1
fi

printf '\nAll %s selected SRA runs downloaded and validated.\n' "$run_count"
