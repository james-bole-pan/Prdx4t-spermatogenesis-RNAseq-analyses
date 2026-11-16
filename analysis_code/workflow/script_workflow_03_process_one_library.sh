#!/usr/bin/env bash
set -Eeuo pipefail

# Processing parameters that define the output BAM. These are declared once and
# used both in the tool invocations and in the compatibility fingerprint, so a
# change to how reads are processed cannot silently leave the fingerprint
# certifying a stale checkpoint as current.
readonly FASTP_FLAGS=(--detect_adapter_for_pe)
readonly BOWTIE2_FLAGS=(--local --very-sensitive-local --no-unal)
readonly STAR_OUT_FILTER_MULTIMAP_NMAX=20
readonly STAR_OUT_SAM_ATTRIBUTES=(NH HI AS nM)

# Bumped whenever the fingerprint's field set changes. A checkpoint written by an
# older format is not rejected outright: it is re-proven with the same structural
# check used for pre-fingerprint BAMs, then rewritten in the current format.
# Format 3 is the single-PRDX4-locus checkpoint; format 4 adds the gene panel.
readonly FINGERPRINT_FORMAT_VERSION_PRDX4=3
readonly FINGERPRINT_FORMAT_VERSION_PANEL=4

readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly CONFIG_HELPER="${SCRIPT_DIR}/lib/load_dataset_config.sh"
[[ -r "$CONFIG_HELPER" ]] || { printf 'ERROR: config helper not found: %s\n' "$CONFIG_HELPER" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG_HELPER"
dataset_selector="$(dataset_selector_from_args "$@")" || exit 1
load_dataset_config "$dataset_selector" "$SCRIPT_DIR" || exit 1
readonly PROJECT="$DATASET_ID"
readonly RUN_MANIFEST

dataset_root="${DATASET_ROOT:-$DEFAULT_DATASET_ROOT}"
reference_dir="${REFERENCE_DIR:-$DEFAULT_REFERENCE_DIR}"
run=""
threads="${THREADS:-12}"
output_flank="${OUTPUT_FLANK:-5000}"
reference_mode="${PRDX4T_REFERENCE_MODE:-panel}"
local_tmp_root="${PRDX4T_LOCAL_TMP:-/private/tmp}"
resume_work=""
keep_work=0
dry_run=0

usage() {
    cat <<EOF
Usage: $(basename "$0") --dataset ID --run SRR_ACCESSION [options]

Process one paired- or single-end ${PROJECT} run:
  SRA -> FASTQ -> fastp -> Bowtie2 gene-panel prefilter -> candidate fragments
      -> full-GRCm39 STAR alignment -> panel-region BAM

In the default panel mode every gene of the workflow step 01b panel is
prefiltered and retained, so PRDX4t can be quantified alongside co-expression
partners such as PCMT1 and CRBN. Mode prdx4 reproduces the original
single-locus behaviour and its file names byte for byte.

Options:
  --dataset ID      Dataset ID or dataset.conf path (default: PRJNA548952)
  --run SRR...       One run accession listed in $(basename "$RUN_MANIFEST")
  --dataset-root DIR Dataset root (default: ${DEFAULT_DATASET_ROOT})
  --reference-dir DIR Shared reference root (default: ${DEFAULT_REFERENCE_DIR})
  --threads N        Worker threads (default: 12)
  --output-flank N   Bases retained on each side of the GTF-derived gene
                     interval in the final BAM (default: 5000)
  --reference-mode M panel (default) or prdx4; panel uses the step 01b gene
                     panel and writes RUN.panel.bam, prdx4 uses the single
                     PRDX4 locus and writes RUN.prdx4.bam
  --local-tmp DIR    Local filesystem for STAR FIFO files (default: /private/tmp)
  --resume-work DIR  Resume at STAR from a retained failed work directory
  --keep-work        Retain large intermediate files after success
  --dry-run          Validate inputs and print the plan without processing reads
  -h, --help         Show this help

Environment:
  DATASET_ROOT       Alternative dataset root
  REFERENCE_DIR       Alternative shared reference root
  THREADS            Alternative thread count
  OUTPUT_FLANK       Alternative final-BAM flank size
  PRDX4T_ENV_BIN     Directory containing the prdx4t-analysis executables
  PRDX4T_LOCAL_TMP   Alternative local filesystem for STAR temporary files
  PRDX4T_REFERENCE_MODE  Alternative default reference mode
  PRDX4T_STAR_BIN    Explicit STAR executable (overrides automatic selection)

Final outputs:
  DATASET_ROOT/analysis/bam/SRR....prdx4.bam[.bai]
  DATASET_ROOT/analysis/logs/SRR.../
  DATASET_ROOT/analysis/metrics/SRR....processing_summary.tsv

Successful runs remove temporary FASTQ and whole-genome BAM files by default.
If a stage fails, its work directory is retained for troubleshooting.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

sha256_file() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        die "neither shasum nor sha256sum is available"
    fi
}

while (($#)); do
    case "$1" in
        --dataset)
            (($# >= 2)) || die "--dataset requires an ID or config file"
            shift 2
            ;;
        --run)
            (($# >= 2)) || die "--run requires an accession"
            run="$2"
            shift 2
            ;;
        --dataset-root)
            (($# >= 2)) || die "--dataset-root requires a directory"
            dataset_root="$2"
            shift 2
            ;;
        --reference-dir)
            (($# >= 2)) || die "--reference-dir requires a directory"
            reference_dir="$2"
            shift 2
            ;;
        --threads)
            (($# >= 2)) || die "--threads requires a value"
            threads="$2"
            shift 2
            ;;
        --output-flank)
            (($# >= 2)) || die "--output-flank requires a value"
            output_flank="$2"
            shift 2
            ;;
        --reference-mode)
            (($# >= 2)) || die "--reference-mode requires panel or prdx4"
            reference_mode="$2"
            shift 2
            ;;
        --local-tmp)
            (($# >= 2)) || die "--local-tmp requires a directory"
            local_tmp_root="$2"
            shift 2
            ;;
        --resume-work)
            (($# >= 2)) || die "--resume-work requires a directory"
            resume_work="$2"
            shift 2
            ;;
        --keep-work)
            keep_work=1
            shift
            ;;
        --dry-run)
            dry_run=1
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

[[ -n "$run" ]] || die "--run is required"
[[ "$run" =~ ^SRR[0-9]+$ ]] || die "invalid run accession: $run"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$output_flank" =~ ^[0-9]+$ ]] || die "--output-flank must be a non-negative integer"
[[ -r "$RUN_MANIFEST" ]] || die "run manifest not found: $RUN_MANIFEST"
command -v cmp >/dev/null 2>&1 || die "cmp is required"
[[ -d "$local_tmp_root" && -w "$local_tmp_root" ]] || \
    die "STAR local temporary directory is not writable: $local_tmp_root"

# STAR uses named pipes to decompress input reads. ExFAT and some network
# filesystems do not support them, so verify the selected local filesystem
# before creating any large FASTQ files.
fifo_test="${local_tmp_root%/}/.prdx4t_fifo_test.$$.${RANDOM}"
if ! mkfifo "$fifo_test" 2>/dev/null; then
    die "$local_tmp_root does not support FIFO files; choose a local APFS directory with --local-tmp"
fi
rm -f -- "$fifo_test"

metadata_line="$(awk -F '\t' -v wanted="$run" 'BEGIN {OFS="\t"} NR > 1 && $1 == wanted {print $1,$2,$3,$4,$5,$6,$7; found=1; exit} END {if (!found) exit 1}' "$RUN_MANIFEST")" || \
    die "$run is not listed in $RUN_MANIFEST"
IFS=$'\t' read -r manifest_run sample_name library_name library_layout spots bases size_mb <<< "$metadata_line"
[[ "$manifest_run" == "$run" ]] || die "manifest lookup failed for $run"
case "$library_layout" in
    PAIRED) layout_label="paired-end" ;;
    SINGLE) layout_label="single-end" ;;
    *) die "$run has unsupported LibraryLayout '$library_layout' in $RUN_MANIFEST" ;;
esac

env_bin="${PRDX4T_ENV_BIN:-}"
if [[ -n "$env_bin" ]]; then
    :
elif command -v STAR >/dev/null 2>&1 && command -v bowtie2 >/dev/null 2>&1; then
    env_bin="$(dirname -- "$(command -v STAR)")"
elif [[ -x "${DEFAULT_ENV_BIN}/STAR" ]]; then
    env_bin="$DEFAULT_ENV_BIN"
else
    die "prdx4t-analysis tools were not found; run script_environment_create_prdx4t_analysis.sh first"
fi

readonly FASTERQ="${env_bin}/fasterq-dump"
readonly FASTP="${env_bin}/fastp"
readonly BOWTIE2="${env_bin}/bowtie2"
readonly SAMTOOLS="${env_bin}/samtools"
readonly PYTHON="${env_bin}/python"
readonly PIGZ="${env_bin}/pigz"
readonly BUNDLED_ARM64_STAR="${SCRIPT_DIR}/bin/STAR-2.7.11b-macos-arm64"

star_bin="${PRDX4T_STAR_BIN:-}"
if [[ -z "$star_bin" ]]; then
    if [[ "$(uname -s)-$(uname -m)" == "Darwin-arm64" && -x "$BUNDLED_ARM64_STAR" ]]; then
        star_bin="$BUNDLED_ARM64_STAR"
    else
        star_bin="${env_bin}/STAR"
    fi
fi
readonly STAR_BIN="$star_bin"

for executable in "$FASTERQ" "$FASTP" "$BOWTIE2" "$STAR_BIN" "$SAMTOOLS" "$PYTHON" "$PIGZ"; do
    [[ -x "$executable" ]] || die "required executable is missing: $executable"
done
[[ "$("$STAR_BIN" --version)" == "2.7.11b" ]] || die "STAR 2.7.11b is required: $STAR_BIN"

readonly SRA_ARCHIVE="${dataset_root}/sra/${run}/${run}.sra"
readonly REFERENCE_DIR="$reference_dir"
case "$reference_mode" in
    panel)
        readonly FEATURE_TABLE="${REFERENCE_DIR}/gene_panel_features.tsv"
        readonly INTERVALS_TABLE="${REFERENCE_DIR}/gene_panel_intervals.tsv"
        readonly BT2_PREFIX="${REFERENCE_DIR}/bt2_panel/panel"
        readonly BT2_SOURCE_SHA="${REFERENCE_DIR}/bt2_panel/source_fasta.sha256"
        readonly BAM_SUFFIX="panel"
        readonly METRICS_PREFIX="panel_"
        readonly FINGERPRINT_FORMAT_VERSION="$FINGERPRINT_FORMAT_VERSION_PANEL"
        readonly PREPARE_STEP="script_workflow_01b_prepare_gene_panel.sh"
        ;;
    prdx4)
        readonly FEATURE_TABLE="${REFERENCE_DIR}/PRDX4_locus_features.tsv"
        readonly INTERVALS_TABLE=""
        readonly BT2_PREFIX="${REFERENCE_DIR}/bt2_locus/locus"
        readonly BT2_SOURCE_SHA="${REFERENCE_DIR}/bt2_locus/source_fasta.sha256"
        readonly BAM_SUFFIX="prdx4"
        readonly METRICS_PREFIX=""
        readonly FINGERPRINT_FORMAT_VERSION="$FINGERPRINT_FORMAT_VERSION_PRDX4"
        readonly PREPARE_STEP="script_workflow_01_prepare_PRDX4_locus.sh"
        ;;
    *)
        die "--reference-mode must be panel or prdx4, not: $reference_mode"
        ;;
esac
readonly STAR_INDEX="${REFERENCE_DIR}/star_GRCm39_Ensembl110_full"
readonly STAR_MANIFEST="${REFERENCE_DIR}/GRCm39_STAR_reference_manifest.tsv"
readonly ANALYSIS_DIR="${dataset_root}/analysis"
readonly BAM_DIR="${ANALYSIS_DIR}/bam"
readonly LOG_ROOT="${ANALYSIS_DIR}/logs"
readonly METRICS_DIR="${ANALYSIS_DIR}/metrics"
readonly WORK_ROOT="${ANALYSIS_DIR}/work"
readonly FINAL_BAM="${BAM_DIR}/${run}.${BAM_SUFFIX}.bam"
readonly FINAL_BAI="${FINAL_BAM}.bai"
readonly FINAL_SUMMARY="${METRICS_DIR}/${run}.${METRICS_PREFIX}processing_summary.tsv"
readonly FINAL_FINGERPRINT="${METRICS_DIR}/${run}.${METRICS_PREFIX}processing_compatibility.tsv"

[[ -r "$FEATURE_TABLE" ]] || die "GTF-derived feature table is missing; rerun ${PREPARE_STEP}: $FEATURE_TABLE"
gene_definition="$(awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) column[$i] = i
        for (name in required) if (!(name in column)) exit 1
        next
    }
    $1 == "gene" {
        print $column["chrom"], $column["start"], $column["end"], $column["strand"]
        found++
    }
    BEGIN {required["chrom"]; required["start"]; required["end"]; required["strand"]}
    END {if (found != 1) exit 1}
' "$FEATURE_TABLE")" || die "could not derive a unique gene interval from $FEATURE_TABLE"
read -r output_chrom gene_start gene_end gene_strand <<< "$gene_definition"
[[ "$gene_start" =~ ^[0-9]+$ && "$gene_end" =~ ^[0-9]+$ && "$gene_start" -le "$gene_end" ]] || \
    die "invalid gene interval in $FEATURE_TABLE: $gene_definition"
[[ "$gene_strand" == "+" || "$gene_strand" == "-" ]] || die "invalid gene strand in $FEATURE_TABLE: $gene_strand"
output_start=$((gene_start - output_flank))
((output_start >= 1)) || output_start=1
readonly OUTPUT_LOCUS="${output_chrom}:${output_start}-$((gene_end + output_flank))"

# Regions retained in the final BAM. In prdx4 mode that is the PRDX4 locus
# alone; in panel mode it is every gene of the step 01b panel, so one STAR
# alignment serves PRDX4t and its co-expression partners.
retained_genes=()
retained_regions=()
if [[ "$reference_mode" == "panel" ]]; then
    [[ -r "$INTERVALS_TABLE" ]] || die "panel interval table is missing; rerun ${PREPARE_STEP}: $INTERVALS_TABLE"
    while IFS=$'\t' read -r panel_gene _panel_id _role _quantification _contig panel_chrom \
        _strand panel_gene_start panel_gene_end _locus_start _locus_end _locus_length _rationale; do
        [[ "$panel_gene" == "gene_name" ]] && continue
        [[ "$panel_gene_start" =~ ^[0-9]+$ && "$panel_gene_end" =~ ^[0-9]+$ ]] || \
            die "invalid interval for $panel_gene in $INTERVALS_TABLE"
        ((panel_gene_start <= panel_gene_end)) || die "inverted interval for $panel_gene"
        panel_start=$((panel_gene_start - output_flank))
        ((panel_start >= 1)) || panel_start=1
        retained_genes+=("$panel_gene")
        retained_regions+=("${panel_chrom}:${panel_start}-$((panel_gene_end + output_flank))")
    done < "$INTERVALS_TABLE"
    ((${#retained_regions[@]} > 0)) || die "panel interval table lists no genes: $INTERVALS_TABLE"
    # The panel and feature tables must agree about PRDX4 itself.
    prdx4_index=-1
    for index in "${!retained_genes[@]}"; do
        [[ "${retained_genes[index]}" == "Prdx4" ]] && prdx4_index="$index"
    done
    ((prdx4_index >= 0)) || die "the panel does not contain Prdx4: $INTERVALS_TABLE"
    [[ "${retained_regions[prdx4_index]}" == "$OUTPUT_LOCUS" ]] || \
        die "panel and feature tables disagree about the PRDX4 region: ${retained_regions[prdx4_index]} vs $OUTPUT_LOCUS"
else
    retained_genes=("Prdx4")
    retained_regions=("$OUTPUT_LOCUS")
fi

if [[ -n "$resume_work" ]]; then
    case "$resume_work" in
        "${WORK_ROOT}/${run}."*) ;;
        *) die "--resume-work must be a retained ${run} directory under ${WORK_ROOT}" ;;
    esac
    [[ -d "$resume_work" ]] || die "retained work directory does not exist: $resume_work"
fi

[[ -s "$SRA_ARCHIVE" ]] || die "SRA archive is missing: $SRA_ARCHIVE"
[[ -s "$BT2_SOURCE_SHA" ]] || die "Bowtie2 source fingerprint is missing: $BT2_SOURCE_SHA"
[[ -s "$STAR_MANIFEST" ]] || die "STAR reference manifest is missing: $STAR_MANIFEST"
# bowtie2-build emits .bt2l instead of .bt2 for large references, so accept
# either; workflow step 01 already validates the index the same way.
for suffix in 1 2 3 4 rev.1 rev.2; do
    [[ -s "${BT2_PREFIX}.${suffix}.bt2" || -s "${BT2_PREFIX}.${suffix}.bt2l" ]] || \
        die "Bowtie2 locus index is incomplete (missing ${BT2_PREFIX}.${suffix}.bt2[l])"
done
for required in Genome SA SAindex chrLength.txt chrName.txt chrStart.txt genomeParameters.txt sjdbInfo.txt; do
    [[ -s "${STAR_INDEX}/${required}" ]] || die "full-genome STAR index is incomplete; missing $required"
done
for region in "${retained_regions[@]}"; do
    grep -Fxq "${region%%:*}" "${STAR_INDEX}/chrName.txt" || \
        die "chromosome ${region%%:*} from the reference tables is absent from the STAR index"
done

# samtools view emits records in the order its regions are given, so a panel
# listed in gene order would produce an unsorted, unindexable BAM. Re-order the
# regions into the alignment's own coordinate order: chromosome as the STAR
# index declares it, then start position. Panel loci never overlap, which
# workflow step 01b enforces, so this ordering is total.
sorted_genes=()
sorted_regions=()
while IFS=$'\t' read -r _chrom_order _region_start sorted_gene sorted_region; do
    [[ -n "$sorted_region" ]] || continue
    sorted_genes+=("$sorted_gene")
    sorted_regions+=("$sorted_region")
done < <(
    for index in "${!retained_regions[@]}"; do
        region="${retained_regions[index]}"
        region_chrom="${region%%:*}"
        region_start="${region#*:}"
        region_start="${region_start%%-*}"
        chrom_order="$(awk -v chrom="$region_chrom" '$0 == chrom {print NR; exit}' \
            "${STAR_INDEX}/chrName.txt")"
        printf '%s\t%s\t%s\t%s\n' "$chrom_order" "$region_start" \
            "${retained_genes[index]}" "$region"
    done | LC_ALL=C sort -k1,1n -k2,2n
)
((${#sorted_regions[@]} == ${#retained_regions[@]})) || \
    die "could not order the retained regions by alignment coordinate"
retained_genes=("${sorted_genes[@]}")
retained_regions=("${sorted_regions[@]}")
readonly retained_genes retained_regions

# Only inputs that can change the output BAM belong here. The feature table is
# deliberately NOT hashed as a whole: this script reads a single row from it (the
# gene interval, pinned below), so hashing the file would invalidate every
# checkpoint whenever workflow step 01 revises an exon annotation that only
# workflow step 04 consumes, forcing a full reprocess to reproduce byte-identical
# BAMs.
write_compatibility_fingerprint() {
    printf 'field\tvalue\n'
    printf 'format_version\t%s\n' "$FINGERPRINT_FORMAT_VERSION"
    printf 'project\t%s\n' "$PROJECT"
    printf 'run\t%s\n' "$run"
    printf 'library_layout\t%s\n' "$library_layout"
    printf 'gene_interval\t%s:%s-%s(%s)\n' \
        "$output_chrom" "$gene_start" "$gene_end" "$gene_strand"
    printf 'output_flank_bp\t%s\n' "$output_flank"
    printf 'reference_mode\t%s\n' "$reference_mode"
    printf 'required_output_region\t%s\n' "$OUTPUT_LOCUS"
    printf 'retained_regions\t%s\n' "${retained_regions[*]}"
    printf 'bowtie_locus_fasta_sha256\t%s\n' "$(tr -d '[:space:]' < "$BT2_SOURCE_SHA")"
    printf 'star_reference_manifest_sha256\t%s\n' "$(sha256_file "$STAR_MANIFEST")"
    if [[ "$library_layout" == "PAIRED" ]]; then
        printf 'fastp_flags\t%s\n' "${FASTP_FLAGS[*]}"
        printf 'bowtie2_input_mode\tpaired (-1/-2)\n'
    else
        printf 'fastp_flags\t(single-end defaults)\n'
        printf 'bowtie2_input_mode\tsingle (-U)\n'
    fi
    printf 'bowtie2_flags\t%s\n' "${BOWTIE2_FLAGS[*]}"
    printf 'STAR_outFilterMultimapNmax\t%s\n' "$STAR_OUT_FILTER_MULTIMAP_NMAX"
    printf 'STAR_outSAMattributes\t%s\n' "${STAR_OUT_SAM_ATTRIBUTES[*]}"
    printf 'STAR_version\t%s\n' "$("$STAR_BIN" --version)"
    printf 'Bowtie2_version\t%s\n' "$("$BOWTIE2" --version 2>&1 | sed -n '1p')"
    printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
    printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1 | sed -n '1p')"
}

legacy_checkpoint_is_compatible() {
    # Panel BAMs are a new output under new file names, so no earlier format of
    # them exists and a saved single-region summary must never satisfy them.
    [[ "$reference_mode" == "prdx4" ]] || return 1
    [[ -s "$FINAL_SUMMARY" ]] || return 1
    local saved_region saved_chrom saved_start saved_end
    saved_region="$(awk -F '\t' '$1 == "final_region" {print $2; found=1; exit} END {if (!found) exit 1}' "$FINAL_SUMMARY")" || return 1
    if [[ "$saved_region" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        saved_chrom="${BASH_REMATCH[1]}"
        saved_start="${BASH_REMATCH[2]}"
        saved_end="${BASH_REMATCH[3]}"
    else
        return 1
    fi
    [[ "$saved_chrom" == "$output_chrom" ]] || return 1
    ((saved_start <= output_start && saved_end >= gene_end + output_flank)) || return 1

    local bam_length star_length
    bam_length="$("$SAMTOOLS" view -H "$FINAL_BAM" | awk -F '\t' -v chrom="$output_chrom" '
        $1 == "@SQ" {
            name=""; seq_length=""
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SN:/) name=substr($i,4)
                if ($i ~ /^LN:/) seq_length=substr($i,4)
            }
            if (name == chrom) {print seq_length; exit}
        }
    ')"
    star_length="$(paste "${STAR_INDEX}/chrName.txt" "${STAR_INDEX}/chrLength.txt" | awk -v chrom="$output_chrom" '$1 == chrom {print $2; exit}')"
    [[ -n "$bam_length" && "$bam_length" == "$star_length" ]]
}

if [[ -s "$FINAL_BAM" && -s "$FINAL_BAI" && -s "$FINAL_SUMMARY" ]]; then
    "$SAMTOOLS" quickcheck "$FINAL_BAM" || die "existing final BAM is corrupt: $FINAL_BAM"
    "$SAMTOOLS" idxstats "$FINAL_BAM" >/dev/null || die "existing BAM index cannot be read: $FINAL_BAI"
    recorded_fingerprint_version=""
    if [[ -s "$FINAL_FINGERPRINT" ]]; then
        recorded_fingerprint_version="$(awk -F '\t' '$1 == "format_version" {print $2; exit}' "$FINAL_FINGERPRINT")"
    fi

    if [[ "$recorded_fingerprint_version" == "$FINGERPRINT_FORMAT_VERSION" ]]; then
        cmp -s <(write_compatibility_fingerprint) "$FINAL_FINGERPRINT" || \
            die "existing BAM checkpoint does not match the current references or processing parameters: $FINAL_BAM"
    elif legacy_checkpoint_is_compatible; then
        if [[ -n "$recorded_fingerprint_version" ]]; then
            checkpoint_kind="format-${recorded_fingerprint_version} checkpoint"
        else
            checkpoint_kind="legacy checkpoint"
        fi
        if ((dry_run)); then
            printf '%s is a compatible %s; the real run will record its compatibility fingerprint.\n' \
                "$run" "$checkpoint_kind"
        else
            fingerprint_tmp="${FINAL_FINGERPRINT}.tmp"
            write_compatibility_fingerprint > "$fingerprint_tmp"
            mv -f "$fingerprint_tmp" "$FINAL_FINGERPRINT"
            printf '%s is a compatible %s; recorded a format-%s fingerprint.\n' \
                "$run" "$checkpoint_kind" "$FINGERPRINT_FORMAT_VERSION"
        fi
    else
        die "existing BAM cannot be proven compatible with the current reference interval: $FINAL_BAM"
    fi
    printf '%s is already complete; no work was performed.\n' "$run"
    printf 'Final BAM: %s\n' "$FINAL_BAM"
    exit 0
fi
if [[ -e "$FINAL_BAM" || -e "$FINAL_BAI" || -e "$FINAL_SUMMARY" || -e "$FINAL_FINGERPRINT" ]]; then
    die "a partial final checkpoint exists; move it aside before retrying: $FINAL_BAM"
fi

fastp_threads="$threads"
if ((fastp_threads > 16)); then
    fastp_threads=16
fi

required_disk_gib="$(awk -v mb="$size_mb" 'BEGIN {printf "%d", (mb * 6 / 1024) + 21}')"
if [[ -d "$dataset_root" ]]; then
    available_disk_gib="$(df -Pk "$dataset_root" | awk 'NR == 2 {printf "%d", $4/1048576}')"
else
    available_disk_gib="0"
fi

cat <<EOF
Run:                 ${run}
Sample:              ${sample_name}
Library:             ${library_name}
Layout:              ${library_layout}
Reported spots:      ${spots}
SRA archive:         ${SRA_ARCHIVE}
Bowtie2 locus index: ${BT2_PREFIX}
Full STAR index:     ${STAR_INDEX}
STAR executable:    ${STAR_BIN}
STAR local temp:     ${local_tmp_root}
Reference mode:      ${reference_mode}
Retained regions:    ${retained_regions[*]}
Threads:             ${threads}
Available disk:      ${available_disk_gib} GiB
Estimated minimum:   ${required_disk_gib} GiB
EOF
if [[ -n "$resume_work" ]]; then
    printf 'Resume work:         %s\n' "$resume_work"
fi

((available_disk_gib >= required_disk_gib)) || \
    die "insufficient free disk for temporary FASTQ files (need about ${required_disk_gib} GiB)"

if ((dry_run)); then
    cat <<EOF

Dry run passed. No reads were processed and no output directories were made.

The real run will retain only the panel-region BAM, its index, QC/alignment
logs, and a processing summary. Large intermediates are deleted after success.
EOF
    exit 0
fi

mkdir -p "$BAM_DIR" "$LOG_ROOT" "$METRICS_DIR" "$WORK_ROOT"
readonly RUN_LOG_DIR="${LOG_ROOT}/${run}"
mkdir -p "$RUN_LOG_DIR"

if [[ -n "$resume_work" ]]; then
    work_dir="$resume_work"
else
    work_dir="$(mktemp -d "${WORK_ROOT}/${run}.XXXXXX")"
fi
readonly work_dir
readonly RAW_DIR="${work_dir}/raw_fastq"
readonly TRIM_DIR="${work_dir}/trimmed_fastq"
readonly CAND_DIR="${work_dir}/candidate_fastq"
readonly TEMP_DIR="${work_dir}/fasterq_tmp"
readonly STAR_OUT_DIR="${work_dir}/star"
mkdir -p "$RAW_DIR" "$TRIM_DIR" "$CAND_DIR" "$TEMP_DIR" "$STAR_OUT_DIR"
readonly RAW_R1="${RAW_DIR}/${run}_1.fastq"
readonly RAW_R2="${RAW_DIR}/${run}_2.fastq"
readonly RAW_SINGLE="${RAW_DIR}/${run}.fastq"
readonly TRIM_R1="${TRIM_DIR}/${run}_1.trim.fq.gz"
readonly TRIM_R2="${TRIM_DIR}/${run}_2.trim.fq.gz"
readonly TRIM_SINGLE="${TRIM_DIR}/${run}.trim.fq.gz"
readonly CANDIDATE_NAMES="${work_dir}/${run}.candidate_names.txt"
readonly CAND_R1="${CAND_DIR}/${run}_candidate_1.fq.gz"
readonly CAND_R2="${CAND_DIR}/${run}_candidate_2.fq.gz"
readonly CAND_SINGLE="${CAND_DIR}/${run}_candidate.fq.gz"

success=0
stage="initialization"
star_tmp_parent=""
cleanup_on_exit() {
    local rc=$?
    if [[ -n "$star_tmp_parent" ]]; then
        case "$star_tmp_parent" in
            "${local_tmp_root%/}/prdx4t_${run}."*) rm -rf -- "$star_tmp_parent" ;;
            *) printf 'WARNING: refused to remove unexpected STAR temp path: %s\n' "$star_tmp_parent" >&2 ;;
        esac
    fi
    if ((success == 1 && keep_work == 0)); then
        case "$work_dir" in
            "${WORK_ROOT}/${run}."*) rm -rf -- "$work_dir" ;;
            *) printf 'WARNING: refused to remove unexpected work path: %s\n' "$work_dir" >&2 ;;
        esac
    elif ((rc != 0)); then
        printf '\nFAILED during: %s\n' "$stage" >&2
        printf 'Diagnostic work files retained at: %s\n' "$work_dir" >&2
    fi
}
trap cleanup_on_exit EXIT

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

if [[ -z "$resume_work" ]]; then
    stage="SRA-to-FASTQ conversion"
    printf '\n[%s] %s: converting the SRA archive to %s FASTQ\n' "$(timestamp)" "$run" "$layout_label"
    fasterq_args=(
        "$SRA_ARCHIVE"
        --skip-technical
        --outdir "$RAW_DIR"
        --temp "$TEMP_DIR"
        --threads "$threads"
    )
    if [[ "$library_layout" == "PAIRED" ]]; then
        fasterq_args+=(--split-files)
    fi
    "$FASTERQ" "${fasterq_args[@]}" \
        >"${RUN_LOG_DIR}/${run}.fasterq.stdout.log" \
        2>"${RUN_LOG_DIR}/${run}.fasterq.stderr.log"
    if [[ "$library_layout" == "PAIRED" ]]; then
        [[ -s "$RAW_R1" && -s "$RAW_R2" ]] || die "fasterq-dump did not produce both paired FASTQ files"
    else
        [[ -s "$RAW_SINGLE" ]] || die "fasterq-dump did not produce the single-end FASTQ file"
    fi

    stage="fastp trimming and QC"
    printf '[%s] %s: trimming and quality filtering\n' "$(timestamp)" "$run"
    fastp_args=(
        --json "${RUN_LOG_DIR}/${run}.fastp.json"
        --html "${RUN_LOG_DIR}/${run}.fastp.html"
        --thread "$fastp_threads"
    )
    if [[ "$library_layout" == "PAIRED" ]]; then
        fastp_args+=(--in1 "$RAW_R1" --in2 "$RAW_R2" --out1 "$TRIM_R1" --out2 "$TRIM_R2" "${FASTP_FLAGS[@]}")
    else
        fastp_args+=(--in1 "$RAW_SINGLE" --out1 "$TRIM_SINGLE")
    fi
    "$FASTP" "${fastp_args[@]}" \
        2>"${RUN_LOG_DIR}/${run}.fastp.stderr.log"
    if [[ "$library_layout" == "PAIRED" ]]; then
        [[ -s "$TRIM_R1" && -s "$TRIM_R2" ]] || die "fastp did not produce both trimmed FASTQ files"
    else
        [[ -s "$TRIM_SINGLE" ]] || die "fastp did not produce the trimmed single-end FASTQ file"
    fi

    # The raw FASTQs are the largest redundant files and are no longer needed.
    rm -f -- "$RAW_R1" "$RAW_R2" "$RAW_SINGLE"

    stage="Bowtie2 PRDX4-locus prefilter"
    printf '[%s] %s: finding fragments with evidence of PRDX4-locus origin\n' "$(timestamp)" "$run"
    bowtie2_args=("${BOWTIE2_FLAGS[@]}" --threads "$threads" -x "$BT2_PREFIX")
    if [[ "$library_layout" == "PAIRED" ]]; then
        bowtie2_args+=(-1 "$TRIM_R1" -2 "$TRIM_R2")
    else
        bowtie2_args+=(-U "$TRIM_SINGLE")
    fi
    "$BOWTIE2" "${bowtie2_args[@]}" \
    2>"${RUN_LOG_DIR}/${run}.bowtie2_locus.log" \
    | "$SAMTOOLS" view - \
    | cut -f 1 \
    | LC_ALL=C sort -u > "${CANDIDATE_NAMES}.tmp"
    mv -f "${CANDIDATE_NAMES}.tmp" "$CANDIDATE_NAMES"
    candidate_fragments="$(wc -l < "$CANDIDATE_NAMES" | tr -d '[:space:]')"
    [[ "$candidate_fragments" =~ ^[0-9]+$ ]] || die "could not count candidate fragments"
    ((candidate_fragments > 0)) || die "Bowtie2 found no PRDX4-locus candidate fragments"
    printf '[%s] %s: %s candidate fragments\n' "$(timestamp)" "$run" "$candidate_fragments"

    stage="candidate-fragment extraction"
    printf '[%s] %s: extracting candidate %s fragments\n' "$(timestamp)" "$run" "$layout_label"
    "$PYTHON" - "$library_layout" "$CANDIDATE_NAMES" "$TRIM_R1" "$TRIM_R2" "$TRIM_SINGLE" "$CAND_R1" "$CAND_R2" "$CAND_SINGLE" \
    >"${RUN_LOG_DIR}/${run}.candidate_extraction.tsv" <<'PY'
import gzip
import sys

layout, names_file, in_r1, in_r2, in_single, out_r1, out_r2, out_single = sys.argv[1:]
with open(names_file) as handle:
    wanted = {line.rstrip("\n") for line in handle if line.strip()}

def normalized_name(header):
    name = header[1:].split(None, 1)[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name

def extract(source, destination):
    written = 0
    with gzip.open(source, "rt") as src, gzip.open(destination, "wt", compresslevel=3) as dst:
        while True:
            header = src.readline()
            if not header:
                break
            sequence = src.readline()
            plus = src.readline()
            quality = src.readline()
            if not (sequence and plus and quality):
                raise RuntimeError(f"truncated FASTQ record in {source}")
            if normalized_name(header) in wanted:
                dst.write(header)
                dst.write(sequence)
                dst.write(plus)
                dst.write(quality)
                written += 1
    return written

print("metric\tvalue")
print(f"requested_candidate_names\t{len(wanted)}")
if layout == "PAIRED":
    count_r1 = extract(in_r1, out_r1)
    count_r2 = extract(in_r2, out_r2)
    print(f"mate1_records_written\t{count_r1}")
    print(f"mate2_records_written\t{count_r2}")
    if count_r1 != len(wanted) or count_r2 != len(wanted):
        raise RuntimeError(
            f"candidate extraction mismatch: wanted={len(wanted)}, R1={count_r1}, R2={count_r2}"
        )
elif layout == "SINGLE":
    count = extract(in_single, out_single)
    print(f"single_end_records_written\t{count}")
    if count != len(wanted):
        raise RuntimeError(f"candidate extraction mismatch: wanted={len(wanted)}, reads={count}")
else:
    raise RuntimeError(f"unsupported layout: {layout}")
PY
    if [[ "$library_layout" == "PAIRED" ]]; then
        [[ -s "$CAND_R1" && -s "$CAND_R2" ]] || die "candidate FASTQ extraction produced an empty mate file"
    else
        [[ -s "$CAND_SINGLE" ]] || die "candidate FASTQ extraction produced an empty single-end file"
    fi
else
    stage="retained candidate validation"
    [[ -s "$CANDIDATE_NAMES" ]] || die "candidate-name checkpoint is missing: $CANDIDATE_NAMES"
    if [[ "$library_layout" == "PAIRED" ]]; then
        [[ -s "$CAND_R1" && -s "$CAND_R2" ]] || die "paired candidate FASTQ checkpoints are missing from $CAND_DIR"
    else
        [[ -s "$CAND_SINGLE" ]] || die "single-end candidate FASTQ checkpoint is missing from $CAND_DIR"
    fi
    candidate_fragments="$(wc -l < "$CANDIDATE_NAMES" | tr -d '[:space:]')"
    [[ "$candidate_fragments" =~ ^[1-9][0-9]*$ ]] || die "retained candidate count is invalid"
    printf '\n[%s] %s: resuming with %s retained candidate fragments\n' "$(timestamp)" "$run" "$candidate_fragments"
fi

stage="full-genome STAR alignment"
printf '[%s] %s: aligning candidate fragments against full GRCm39\n' "$(timestamp)" "$run"
readonly STAR_PREFIX="${STAR_OUT_DIR}/${run}."
star_tmp_parent="$(mktemp -d "${local_tmp_root%/}/prdx4t_${run}.XXXXXX")"
readonly STAR_TMP_DIR="${star_tmp_parent}/STARtmp"
if [[ "$library_layout" == "PAIRED" ]]; then
    star_reads=("$CAND_R1" "$CAND_R2")
else
    star_reads=("$CAND_SINGLE")
fi
"$STAR_BIN" \
    --genomeDir "$STAR_INDEX" \
    --runThreadN "$threads" \
    --sysShell /bin/bash \
    --outTmpDir "$STAR_TMP_DIR" \
    --readFilesIn "${star_reads[@]}" \
    --readFilesCommand "$PIGZ" -dc \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$STAR_PREFIX" \
    --outFilterMultimapNmax "$STAR_OUT_FILTER_MULTIMAP_NMAX" \
    --outSAMattributes "${STAR_OUT_SAM_ATTRIBUTES[@]}" \
    >"${RUN_LOG_DIR}/${run}.star.stdout.log" \
    2>"${RUN_LOG_DIR}/${run}.star.stderr.log"
readonly STAR_BAM="${STAR_PREFIX}Aligned.out.bam"
[[ -s "$STAR_BAM" ]] || die "STAR did not produce an alignment BAM"
cp -f "${STAR_PREFIX}Log.final.out" "${RUN_LOG_DIR}/${run}.STAR.Log.final.out"

stage="BAM sorting and PRDX4-region extraction"
printf '[%s] %s: sorting and retaining %s\n' "$(timestamp)" "$run" "${retained_regions[*]}"
readonly SORTED_BAM="${work_dir}/${run}.candidate_genome.sorted.bam"
readonly LOCUS_BAM="${work_dir}/${run}.${BAM_SUFFIX}.bam"
"$SAMTOOLS" sort --threads "$threads" -o "$SORTED_BAM" "$STAR_BAM"
"$SAMTOOLS" index --threads "$threads" "$SORTED_BAM"
"$SAMTOOLS" view --threads "$threads" -b -o "$LOCUS_BAM" "$SORTED_BAM" "${retained_regions[@]}"
"$SAMTOOLS" index --threads "$threads" "$LOCUS_BAM"
"$SAMTOOLS" quickcheck "$LOCUS_BAM"
"$SAMTOOLS" flagstat --threads "$threads" "$LOCUS_BAM" > "${RUN_LOG_DIR}/${run}.${BAM_SUFFIX}.flagstat.txt"

# Records count each mate separately; fragments count each read pair once, which
# is the unit workflow step 04 reports. Both are recorded so the two are comparable.
# Per-gene records let a failed prefilter be caught per gene rather than only
# in the panel total. The unqualified locus_* metrics below stay PRDX4-only so
# they remain directly comparable with format-3 summaries.
gene_record_lines=()
for index in "${!retained_genes[@]}"; do
    gene_records="$("$SAMTOOLS" view -c "$LOCUS_BAM" "${retained_regions[index]}")"
    gene_unique="$("$SAMTOOLS" view -c -q 255 "$LOCUS_BAM" "${retained_regions[index]}")"
    gene_record_lines+=("$(printf 'gene_records.%s\t%s (MAPQ255 %s) %s' \
        "${retained_genes[index]}" "$gene_records" "$gene_unique" "${retained_regions[index]}")")
done

locus_alignment_records="$("$SAMTOOLS" view -c "$LOCUS_BAM" "$OUTPUT_LOCUS")"
unique_locus_alignment_records="$("$SAMTOOLS" view -c -q 255 "$LOCUS_BAM" "$OUTPUT_LOCUS")"
locus_fragments="$("$SAMTOOLS" view "$LOCUS_BAM" "$OUTPUT_LOCUS" | cut -f 1 | LC_ALL=C sort -u | wc -l | tr -d '[:space:]')"
unique_locus_fragments="$("$SAMTOOLS" view -q 255 "$LOCUS_BAM" "$OUTPUT_LOCUS" | cut -f 1 | LC_ALL=C sort -u | wc -l | tr -d '[:space:]')"

stage="final output checkpoint"
summary_tmp="${FINAL_SUMMARY}.tmp"
fingerprint_tmp="${FINAL_FINGERPRINT}.tmp"
{
    printf 'metric\tvalue\n'
    printf 'project\t%s\n' "$PROJECT"
    printf 'run\t%s\n' "$run"
    printf 'sample_name\t%s\n' "$sample_name"
    printf 'library_name\t%s\n' "$library_name"
    printf 'manifest_spots\t%s\n' "$spots"
    printf 'library_layout\t%s\n' "$library_layout"
    printf 'bowtie2_candidate_fragments\t%s\n' "$candidate_fragments"
    printf 'reference_mode\t%s\n' "$reference_mode"
    printf 'final_region\t%s\n' "$OUTPUT_LOCUS"
    printf 'retained_regions\t%s\n' "${retained_regions[*]}"
    for summary_line in "${gene_record_lines[@]}"; do
        printf '%s\n' "$summary_line"
    done
    printf 'final_region_flank_bp\t%s\n' "$output_flank"
    printf 'feature_definition\t%s\n' "$FEATURE_TABLE"
    printf 'count_unit_note\tfragment is one read pair for PAIRED or one read for SINGLE; query names are counted once per window\n'
    printf 'locus_alignment_records\t%s\n' "$locus_alignment_records"
    printf 'unique_locus_alignment_records_mapq255\t%s\n' "$unique_locus_alignment_records"
    printf 'locus_fragments\t%s\n' "$locus_fragments"
    printf 'unique_locus_fragments_mapq255\t%s\n' "$unique_locus_fragments"
    printf 'STAR_version\t%s\n' "$("$STAR_BIN" --version)"
    printf 'STAR_binary\t%s\n' "$STAR_BIN"
    printf 'Bowtie2_version\t%s\n' "$("$BOWTIE2" --version 2>&1 | sed -n '1p')"
    printf 'samtools_version\t%s\n' "$("$SAMTOOLS" --version | sed -n '1p')"
    printf 'fastp_version\t%s\n' "$("$FASTP" --version 2>&1 | sed -n '1p')"
} > "$summary_tmp"
write_compatibility_fingerprint > "$fingerprint_tmp"

cp -f "$LOCUS_BAM" "${FINAL_BAM}.partial"
cp -f "${LOCUS_BAM}.bai" "${FINAL_BAI}.partial"
"$SAMTOOLS" quickcheck "${FINAL_BAM}.partial"
mv -f "${FINAL_BAI}.partial" "$FINAL_BAI"
mv -f "${FINAL_BAM}.partial" "$FINAL_BAM"
mv -f "$summary_tmp" "$FINAL_SUMMARY"
mv -f "$fingerprint_tmp" "$FINAL_FINGERPRINT"

success=1
printf '\n[%s] %s COMPLETE\n' "$(timestamp)" "$run"
printf 'Candidate fragments:         %s\n' "$candidate_fragments"
printf 'PRDX4 alignment records:     %s\n' "$locus_alignment_records"
printf 'Unique records (MAPQ 255):   %s\n' "$unique_locus_alignment_records"
printf 'PRDX4 fragments:             %s\n' "$locus_fragments"
printf 'Unique fragments (MAPQ 255): %s\n' "$unique_locus_fragments"
for summary_line in "${gene_record_lines[@]}"; do
    printf '%s\n' "$summary_line"
done
printf 'Final BAM:                   %s\n' "$FINAL_BAM"
printf 'Processing summary:          %s\n' "$FINAL_SUMMARY"
if ((keep_work)); then
    printf 'Intermediate work retained:  %s\n' "$work_dir"
else
    printf 'Large intermediates will now be removed.\n'
fi
