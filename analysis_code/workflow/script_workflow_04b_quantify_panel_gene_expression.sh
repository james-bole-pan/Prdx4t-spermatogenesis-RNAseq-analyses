#!/usr/bin/env bash
set -Eeuo pipefail

readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly CONFIG_HELPER="${SCRIPT_DIR}/lib/load_dataset_config.sh"
[[ -r "$CONFIG_HELPER" ]] || { printf 'ERROR: config helper not found: %s\n' "$CONFIG_HELPER" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_HELPER"
dataset_selector="$(dataset_selector_from_args "$@")" || exit 1
load_dataset_config "$dataset_selector" "$SCRIPT_DIR" || exit 1
readonly RUN_MANIFEST

dataset_root="${DATASET_ROOT:-$DEFAULT_DATASET_ROOT}"
reference_dir="${REFERENCE_DIR:-$DEFAULT_REFERENCE_DIR}"
output_dir=""
feature_file=""
min_mapq="${MIN_MAPQ:-255}"
strandedness="${STRANDEDNESS:-auto}"
dry_run=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Quantify gene-level expression for every gene of the workflow step 01b panel
from the panel BAMs written by workflow step 03, then group technical runs by
SampleID. Biological replicates are never combined.

Counts are exon-union fragment counts: a fragment is counted once for a gene if
any aligned block touches any annotated exon of that gene. This is the same
union model the PRDX4t shared-body counting uses, so PRDX4 gene-level values are
directly comparable with PCMT1 and CRBN.

Three levels are written for every gene and sample, per decision on raw versus
comparable values:
  fragments                            raw, no normalization at all
  fragments_per_million_spots          library depth only
  fragments_per_kb_per_million_spots   depth and exonic length, cross-gene

Retained exonic lengths differ among genes. Use the per-kb values rather than
raw fragment counts when comparing the expression profiles of different genes.

Options:
  --dataset ID         Dataset ID or dataset.conf path (default: PRJNA548952)
  --dataset-root DIR   Dataset root (default: ${DEFAULT_DATASET_ROOT})
  --reference-dir DIR  Shared reference root (default: ${DEFAULT_REFERENCE_DIR})
  --output-dir DIR     Result directory
                       (default: DATASET_ROOT/analysis/quantification)
  --features FILE      Panel feature table from workflow step 01b
                       (default: REFERENCE_DIR/gene_panel_features.tsv)
  --min-mapq N         Minimum alignment MAPQ (default: 255)
  --strandedness MODE  auto, unstranded, forward, or reverse (default: auto)
  --dry-run            Validate every input without writing results
  -h, --help           Show this help

Outputs:
  ${DATASET_ID}_panel_gene_raw_counts_per_library.tsv
  ${DATASET_ID}_panel_gene_expression_per_biological_sample.tsv
  ${DATASET_ID}_panel_gene_quantification_manifest.tsv
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --dataset) (($# >= 2)) || die "--dataset requires a value"; shift 2 ;;
        --dataset-root) (($# >= 2)) || die "--dataset-root requires a directory"; dataset_root="$2"; shift 2 ;;
        --reference-dir) (($# >= 2)) || die "--reference-dir requires a directory"; reference_dir="$2"; shift 2 ;;
        --output-dir) (($# >= 2)) || die "--output-dir requires a directory"; output_dir="$2"; shift 2 ;;
        --features) (($# >= 2)) || die "--features requires a file"; feature_file="$2"; shift 2 ;;
        --min-mapq) (($# >= 2)) || die "--min-mapq requires an integer"; min_mapq="$2"; shift 2 ;;
        --strandedness) (($# >= 2)) || die "--strandedness requires a mode"; strandedness="$2"; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1 (use --help)" ;;
    esac
done

case "$strandedness" in
    auto|unstranded|forward|reverse) ;;
    *) die "--strandedness must be auto, unstranded, forward, or reverse" ;;
esac
[[ "$min_mapq" =~ ^[0-9]+$ ]] || die "--min-mapq must be a non-negative integer"

[[ -n "$feature_file" ]] || feature_file="${reference_dir}/gene_panel_features.tsv"
[[ -r "$feature_file" ]] || die "panel feature table not found; run script_workflow_01b_prepare_gene_panel.sh: $feature_file"
[[ -n "$output_dir" ]] || output_dir="${dataset_root}/analysis/quantification"

python_bin=""
if [[ -n "${PRDX4T_ENV_BIN:-}" && -x "${PRDX4T_ENV_BIN}/python" ]]; then
    python_bin="${PRDX4T_ENV_BIN}/python"
elif [[ -x "${DEFAULT_ENV_BIN}/python" ]]; then
    python_bin="$DEFAULT_ENV_BIN/python"
else
    die "the prdx4t-analysis Python was not found; run script_environment_create_prdx4t_analysis.sh first"
fi

mkdir -p "$output_dir"
"$python_bin" "${SCRIPT_DIR}/lib/quantify_panel_genes.py" \
    --dataset-id "$DATASET_ID" \
    --dataset-label "$DATASET_LABEL" \
    --main-stage-groups "$MAIN_STAGE_GROUPS" \
    --wt-sample-filter "${WT_SAMPLE_FILTER:-}" \
    --manifest "$RUN_MANIFEST" \
    --bam-dir "${dataset_root}/analysis/bam" \
    --features "$feature_file" \
    --output-dir "$output_dir" \
    --min-mapq "$min_mapq" \
    --strandedness "$strandedness" \
    $( ((dry_run)) && printf -- '--dry-run' )
