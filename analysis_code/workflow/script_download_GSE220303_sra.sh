#!/usr/bin/env bash
set -Eeuo pipefail

readonly GEO_ACCESSION="GSE220303"
readonly BIOPROJECT="PRJNA909605"
readonly SRA_STUDY="SRP411852"
readonly RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${BIOPROJECT}"
readonly DEFAULT_OUTPUT_DIR="${SRA_ROOT:-${PWD}/data}/GSE220303"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly STRATEGY="RNA-Seq"

output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
selected_run=""
selected_sex=""
metadata_only=0
refresh_metadata=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Download the 20 raw RNA-seq SRA archives in ${GEO_ACCESSION} (3D genome
remodeling during meiotic prophase of mouse oogenesis and spermatogenesis).
The runs form SRA study ${SRA_STUDY} under BioProject ${BIOPROJECT}.

${BIOPROJECT} also holds 20 sisHi-C and 4 ATAC-seq runs of the same series,
which are deliberately excluded: they add roughly 770 GiB and carry no
transcript-level information. Only the RNA-seq runs are downloaded.

The 20 runs are 20 GEO samples, one run each, all paired-end HiSeq X Ten
libraries from FACS-sorted B6xPWK F1 germ cells, 2 replicates per stage:
  male    SGA spermatogonia A, then leptotene, zygotene, pachytene and
          diplotene spermatocytes
  female  PGC, then leptotene, zygotene, pachytene and diplotene oocytes
Use --sex male to fetch only the 10 spermatogenic runs.

Options:
  --output-dir DIR    Dataset root (default: ${DEFAULT_OUTPUT_DIR})
  --run SRR...        Download one run only
  --sex male|female   Restrict to one germ line
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
        --sex)
            (($# >= 2)) || die "--sex requires male or female"
            selected_sex="$2"
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
if [[ -n "$selected_sex" ]]; then
    case "$selected_sex" in
        male|female) ;;
        *) die "unknown sex: $selected_sex (male or female)" ;;
    esac
fi

command -v curl >/dev/null 2>&1 || die "curl is required"
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

# NCBI RunInfo carries no stage, sex or replicate for this project, so the GEO
# sample decoding is version controlled here. The Hi-C and ATAC-seq runs of the
# same BioProject are dropped on LibraryStrategy before this table is consulted.
manifest_tmp="${RUN_MANIFEST}.tmp"
"$PYTHON" - "$RUNINFO" "$manifest_tmp" "$BIOPROJECT" "$SRA_STUDY" "$GEO_ACCESSION" \
    "$STRATEGY" "$selected_run" "$selected_sex" <<'PY'
import csv
import sys
from pathlib import Path

(
    runinfo_path,
    output_path,
    bioproject,
    sra_study,
    geo_accession,
    wanted_strategy,
    wanted_run,
    wanted_sex,
) = sys.argv[1:]

# stage code -> (sample name, label, ordering, sex)
STAGES = {
    "SGA": ("SGA_spermatogoniaA", "Spermatogonia A", 0, "male"),
    "MLEP": ("MLEP_leptotenespermatocytes", "Leptotene spermatocytes", 1, "male"),
    "MZYG": ("MZYG_zygotenespermatocytes", "Zygotene spermatocytes", 2, "male"),
    "MPAC": ("MPAC_pachytenespermatocytes", "Pachytene spermatocytes", 3, "male"),
    "MDIP": ("MDIP_diplotenespermatocytes", "Diplotene spermatocytes", 4, "male"),
    "FPGC": ("FPGC_primordialgermcells", "Female primordial germ cells", 5, "female"),
    "FLEP": ("FLEP_leptoteneoocytes", "Leptotene oocytes", 6, "female"),
    "FZYG": ("FZYG_zygoteneoocytes", "Zygotene oocytes", 7, "female"),
    "FPAC": ("FPAC_pachyteneoocytes", "Pachytene oocytes", 8, "female"),
    "FDIP": ("FDIP_diploteneoocytes", "Diplotene oocytes", 9, "female"),
}
# GEO sample -> (stage, replicate), decoded from the sample titles
# "B6xPWK_<stage>_rna_rep<n>" and the GEO tissue/Stage characteristics.
SAMPLES = {
    "GSM6798053": ("FPGC", 1),
    "GSM6798054": ("FPGC", 2),
    "GSM6798055": ("FLEP", 1),
    "GSM6798056": ("FLEP", 2),
    "GSM6798057": ("FZYG", 1),
    "GSM6798058": ("FZYG", 2),
    "GSM6798059": ("FPAC", 1),
    "GSM6798060": ("FPAC", 2),
    "GSM6798061": ("FDIP", 1),
    "GSM6798062": ("FDIP", 2),
    "GSM6798063": ("SGA", 1),
    "GSM6798064": ("SGA", 2),
    "GSM6798065": ("MLEP", 1),
    "GSM6798066": ("MLEP", 2),
    "GSM6798067": ("MZYG", 1),
    "GSM6798068": ("MZYG", 2),
    "GSM6798069": ("MPAC", 1),
    "GSM6798070": ("MPAC", 2),
    "GSM6798071": ("MDIP", 1),
    "GSM6798072": ("MDIP", 2),
}
OTHER_STRATEGIES = {"Hi-C", "ATAC-seq"}


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
    "SRAStudy",
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
by_sample = {}
for record in run_rows:
    run = clean(record["Run"])
    if not run.startswith("SRR") or not run[3:].isdigit():
        raise SystemExit(f"ERROR: invalid run accession in RunInfo: {run!r}")
    if run in seen_runs:
        raise SystemExit(f"ERROR: duplicate run in RunInfo: {run}")
    seen_runs.add(run)

    strategy = clean(record["LibraryStrategy"])
    if strategy not in OTHER_STRATEGIES | {wanted_strategy}:
        raise SystemExit(f"ERROR: unexpected library strategy for {run}: {strategy}")
    if strategy != wanted_strategy:
        continue
    if clean(record["BioProject"]) != bioproject:
        raise SystemExit(f"ERROR: {run} is not assigned to {bioproject}")
    if clean(record["SRAStudy"]) != sra_study:
        raise SystemExit(f"ERROR: {run} is not part of study {sra_study}")
    if clean(record["ScientificName"]) != "Mus musculus":
        raise SystemExit(f"ERROR: {run} is not a mouse library")

    layout = clean(record["LibraryLayout"])
    if layout != "PAIRED":
        raise SystemExit(f"ERROR: {run} is not paired-end as expected: {layout}")

    gsm = clean(record["SampleName"])
    if gsm not in SAMPLES:
        raise SystemExit(f"ERROR: {run} belongs to unrecognised GEO sample {gsm!r}")

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

    if gsm in by_sample:
        raise SystemExit(f"ERROR: {gsm} has more than one run; one run per sample was expected")
    by_sample[gsm] = record

if set(by_sample) != set(SAMPLES):
    absent = sorted(set(SAMPLES) - set(by_sample))
    raise SystemExit(f"ERROR: {geo_accession} RNA-seq samples absent from NCBI RunInfo: {absent}")

rows = []
for gsm in sorted(
    SAMPLES, key=lambda key: (STAGES[SAMPLES[key][0]][2], SAMPLES[key][1])
):
    stage, replicate = SAMPLES[gsm]
    sample_name, stage_label, stage_order, sex = STAGES[stage]
    record = by_sample[gsm]
    rows.append(
        {
            "Run": clean(record["Run"]),
            "SampleName": sample_name,
            "LibraryName": f"{stage}_rep{replicate}",
            "LibraryLayout": clean(record["LibraryLayout"]),
            "spots": clean(record["spots"]),
            "bases": clean(record["bases"]),
            "size_MB": clean(record["size_MB"]),
            "SampleID": gsm,
            "Group": stage,
            "GroupLabel": stage_label,
            "GroupOrder": str(stage_order),
            "Replicate": str(replicate),
            "TechnicalRun": "1",
            "Sex": sex,
            "Experiment": clean(record["Experiment"]),
            "BioSample": clean(record["BioSample"]),
        }
    )

if wanted_sex:
    rows = [row for row in rows if row["Sex"] == wanted_sex]
    if not rows:
        raise SystemExit(f"ERROR: no runs remain for sex {wanted_sex}")
if wanted_run:
    if wanted_run not in seen_runs:
        raise SystemExit(f"ERROR: {wanted_run} is not part of {geo_accession}/{bioproject}")
    rows = [row for row in rows if row["Run"] == wanted_run]
    if not rows:
        raise SystemExit(
            f"ERROR: {wanted_run} is excluded by the requested subset "
            "(it may be a Hi-C or ATAC-seq run)"
        )

fieldnames = [
    "Run",
    "SampleName",
    "LibraryName",
    "LibraryLayout",
    "spots",
    "bases",
    "size_MB",
    "SampleID",
    "Group",
    "GroupLabel",
    "GroupOrder",
    "Replicate",
    "TechnicalRun",
    "Sex",
    "Experiment",
    "BioSample",
]
with Path(output_path).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY
mv -f "$manifest_tmp" "$RUN_MANIFEST"

run_count="$(awk 'END {print NR-1}' "$RUN_MANIFEST")"
sample_count="$(awk -F '\t' 'NR > 1 {seen[$8]=1} END {print length(seen)}' "$RUN_MANIFEST")"
stage_count="$(awk -F '\t' 'NR > 1 {seen[$9]=1} END {print length(seen)}' "$RUN_MANIFEST")"
reported_gib="$(awk -F '\t' 'NR > 1 {mb += $7} END {printf "%.1f", mb/1024}' "$RUN_MANIFEST")"

printf '\nAccession:       %s\n' "$GEO_ACCESSION"
printf 'BioProject:      %s\n' "$BIOPROJECT"
printf 'SRA study:       %s\n' "$SRA_STUDY"
printf 'GEO samples:     %s\n' "$sample_count"
printf 'Germ stages:     %s\n' "$stage_count"
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
# This loop reads the manifest on stdin, so every child gets /dev/null instead;
# a tool that consumed stdin would otherwise swallow run records.
while IFS=$'\t' read -r run sample_name library layout spots bases size_mb gsm group group_label _rest; do
    [[ "$run" == "Run" ]] && continue
    archive="${SRA_DIR}/${run}/${run}.sra"
    printf '\n[%s] %s | %s | %s (reported %.1f GiB)\n' \
        "$run" "$gsm" "$library" "$layout" \
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
