#!/usr/bin/env bash
set -Eeuo pipefail

readonly GEO_ACCESSION="GSE281037"
readonly BIOPROJECT="PRJNA1181871"
readonly SRA_STUDY="SRP543150"
readonly RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${BIOPROJECT}"
readonly DEFAULT_OUTPUT_DIR="${SRA_ROOT:-${PWD}/data}/GSE281037"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly STRATEGY="RNA-Seq"
# This dataset is downloaded for the wild-type arm only; the Zfy1 KO, Zfy2 KO
# and Zfy DKO samples of the same study are deliberately excluded.
readonly GENOTYPE="WT"

output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
selected_run=""
selected_stage=""
metadata_only=0
refresh_metadata=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Download the 15 wild-type RNA-seq SRA archives in ${GEO_ACCESSION} (ZFY1/ZFY2
target genes in male germ cells). The runs form SRA study ${SRA_STUDY} under
BioProject ${BIOPROJECT}; every run in the project is RNA-seq.

${GEO_ACCESSION} holds 48 runs across 4 genotypes. Only the 15 ${GENOTYPE} runs are
downloaded here; the Zfy1 KO, Zfy2 KO and Zfy DKO runs are excluded.

Each run is one GEO sample, paired-end NextSeq 2000. The wild-type arm covers
3 sorted germ-cell stages - SC1 primary spermatocytes (6), SC2 secondary
spermatocytes (3) and RS round spermatids (6) - in 3 replicates per batch.

The manifest records a Batch column, A or B, taken from the two GEO title
conventions: batch A is titled "sampleN_SN" and batch B "N_genotype_stage".
GEO documents no difference between them - one submission, one protocol, one
processing pipeline, and all 48 are distinct libraries - but the deposited
batch A reads are untrimmed (uniform 301 bp per pair) while batch B reads are
trimmed (276-293 bp, varying per sample), and batch A yields roughly half the
counted depth of batch B at similar spot counts. Same-cell-type expression
correlates at r~0.99 within a batch and r~0.98 across, so the batch effect is
small but systematic: model it rather than pooling A and B as six replicates.

Complete download is about 55 GiB compressed; --stage subsets it further.

Options:
  --output-dir DIR    Dataset root (default: ${DEFAULT_OUTPUT_DIR})
  --run SRR...        Download one run only (must be a ${GENOTYPE} run)
  --stage CODE        Restrict to SC1, SC2 or RS
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
        --stage)
            (($# >= 2)) || die "--stage requires a code"
            selected_stage="$2"
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
if [[ -n "$selected_stage" ]]; then
    case "$selected_stage" in
        SC1|SC2|RS) ;;
        *) die "unknown stage: $selected_stage (SC1, SC2 or RS)" ;;
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

# NCBI RunInfo carries no genotype, stage, batch or replicate for this project,
# and the GEO sample titles of batch A carry none either, so the cell type and
# genotype read from the GEO sample characteristics are version controlled here.
# Downloading fails rather than guesses if NCBI reports a sample outside it.
manifest_tmp="${RUN_MANIFEST}.tmp"
"$PYTHON" - "$RUNINFO" "$manifest_tmp" "$BIOPROJECT" "$SRA_STUDY" "$GEO_ACCESSION" \
    "$STRATEGY" "$selected_run" "$GENOTYPE" "$selected_stage" <<'PY'
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
    wanted_genotype,
    wanted_stage,
) = sys.argv[1:]

GENOTYPES = {
    "WT": ("Wild type", 0),
    "Zfy1KO": ("Zfy1 knockout", 1),
    "Zfy2KO": ("Zfy2 knockout", 2),
    "ZfyDKO": ("Zfy1/Zfy2 double knockout", 3),
}
GENOTYPES_ALLOWED = set(GENOTYPES)
if wanted_genotype not in GENOTYPES_ALLOWED:
    raise SystemExit(f"ERROR: unknown genotype: {wanted_genotype}")
STAGES = {
    "SC1": ("SC1_primaryspermatocytes", "Primary spermatocytes", 0),
    "SC2": ("SC2_secondaryspermatocytes", "Secondary spermatocytes", 1),
    "RS": ("RS_roundspermatids", "Round spermatids", 2),
}
# GEO sample -> (genotype, stage, submission batch, biological replicate),
# read from the GEO sample characteristics "cell type" and "genotype" fields.
SAMPLES = {
    "GSM8609620": ("WT", "SC1", "A", 1),
    "GSM8609621": ("WT", "SC1", "A", 2),
    "GSM8609622": ("WT", "SC1", "A", 3),
    "GSM8609623": ("WT", "RS", "A", 1),
    "GSM8609624": ("WT", "RS", "A", 2),
    "GSM8609625": ("WT", "RS", "A", 3),
    "GSM8609626": ("ZfyDKO", "SC1", "A", 1),
    "GSM8609627": ("ZfyDKO", "SC1", "A", 2),
    "GSM8609628": ("ZfyDKO", "SC1", "A", 3),
    "GSM8609629": ("ZfyDKO", "RS", "A", 1),
    "GSM8609630": ("ZfyDKO", "RS", "A", 2),
    "GSM8609631": ("ZfyDKO", "RS", "A", 3),
    "GSM8609632": ("WT", "SC2", "B", 1),
    "GSM8609633": ("WT", "SC2", "B", 2),
    "GSM8609634": ("WT", "SC2", "B", 3),
    "GSM8609635": ("ZfyDKO", "SC2", "B", 1),
    "GSM8609636": ("ZfyDKO", "SC2", "B", 2),
    "GSM8609637": ("ZfyDKO", "SC2", "B", 3),
    "GSM8609638": ("WT", "RS", "B", 1),
    "GSM8609639": ("WT", "RS", "B", 2),
    "GSM8609640": ("WT", "RS", "B", 3),
    "GSM8609641": ("WT", "SC1", "B", 1),
    "GSM8609642": ("WT", "SC1", "B", 2),
    "GSM8609643": ("WT", "SC1", "B", 3),
    "GSM8609644": ("Zfy2KO", "SC1", "B", 1),
    "GSM8609645": ("Zfy2KO", "SC1", "B", 2),
    "GSM8609646": ("Zfy2KO", "SC1", "B", 3),
    "GSM8609647": ("Zfy2KO", "SC2", "B", 1),
    "GSM8609648": ("Zfy2KO", "SC2", "B", 2),
    "GSM8609649": ("Zfy2KO", "SC2", "B", 3),
    "GSM8609650": ("Zfy2KO", "RS", "B", 1),
    "GSM8609651": ("Zfy2KO", "RS", "B", 2),
    "GSM8609652": ("Zfy2KO", "RS", "B", 3),
    "GSM8609653": ("Zfy1KO", "SC1", "B", 1),
    "GSM8609654": ("Zfy1KO", "SC1", "B", 2),
    "GSM8609655": ("Zfy1KO", "SC1", "B", 3),
    "GSM8609656": ("Zfy1KO", "SC2", "B", 1),
    "GSM8609657": ("Zfy1KO", "SC2", "B", 2),
    "GSM8609658": ("Zfy1KO", "SC2", "B", 3),
    "GSM8609659": ("Zfy1KO", "RS", "B", 1),
    "GSM8609660": ("Zfy1KO", "RS", "B", 2),
    "GSM8609661": ("Zfy1KO", "RS", "B", 3),
    "GSM8609662": ("ZfyDKO", "SC1", "B", 1),
    "GSM8609663": ("ZfyDKO", "SC1", "B", 2),
    "GSM8609664": ("ZfyDKO", "SC1", "B", 3),
    "GSM8609665": ("ZfyDKO", "RS", "B", 1),
    "GSM8609666": ("ZfyDKO", "RS", "B", 2),
    "GSM8609667": ("ZfyDKO", "RS", "B", 3),
}


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
    if strategy != wanted_strategy:
        raise SystemExit(f"ERROR: unexpected library strategy for {run}: {strategy}")
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
    raise SystemExit(f"ERROR: {geo_accession} samples absent from NCBI RunInfo: {absent}")

rows = []
for gsm in sorted(
    SAMPLES,
    key=lambda key: (
        GENOTYPES[SAMPLES[key][0]][1],
        STAGES[SAMPLES[key][1]][2],
        SAMPLES[key][2],
        SAMPLES[key][3],
    ),
):
    genotype, stage, batch, replicate = SAMPLES[gsm]
    genotype_label, _genotype_order = GENOTYPES[genotype]
    sample_name, stage_label, stage_order = STAGES[stage]
    record = by_sample[gsm]
    rows.append(
        {
            "Run": clean(record["Run"]),
            "SampleName": sample_name,
            "LibraryName": f"{genotype}_{stage}_{batch}{replicate}",
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
            "Genotype": genotype,
            "GenotypeLabel": genotype_label,
            "Batch": batch,
            "Experiment": clean(record["Experiment"]),
            "BioSample": clean(record["BioSample"]),
        }
    )

rows = [row for row in rows if row["Genotype"] == wanted_genotype]
if not rows:
    raise SystemExit(f"ERROR: no runs remain for genotype {wanted_genotype}")
if wanted_stage:
    rows = [row for row in rows if row["Group"] == wanted_stage]
    if not rows:
        raise SystemExit(f"ERROR: no runs remain for stage {wanted_stage}")
if wanted_run:
    if wanted_run not in seen_runs:
        raise SystemExit(f"ERROR: {wanted_run} is not part of {geo_accession}/{bioproject}")
    rows = [row for row in rows if row["Run"] == wanted_run]
    if not rows:
        raise SystemExit(
            f"ERROR: {wanted_run} is excluded by the requested subset "
            f"(only {wanted_genotype} runs are downloaded by this script)"
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
    "Genotype",
    "GenotypeLabel",
    "Batch",
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
printf 'Genotype:        %s\n' "$GENOTYPE"
printf 'Cell stages:     %s\n' "$stage_count"
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
