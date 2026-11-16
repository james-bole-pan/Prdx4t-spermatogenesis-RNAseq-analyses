#!/usr/bin/env bash
set -Eeuo pipefail

readonly ENSEMBL_RELEASE="110"
readonly ASSEMBLY="GRCm39"
readonly READ_LENGTH="150"
readonly SJDB_OVERHANG="149"

readonly GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
readonly GENOME_SUM="16996"
readonly GENOME_BLOCKS="787519"
readonly GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"
readonly GTF_SUM="42178"
readonly GTF_BLOCKS="31546"

readonly DEFAULT_REFERENCE_DIR="${PRDX4T_REFERENCE_DIR:-${SRA_ROOT:-${PWD}/data}/reference}"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly MIN_RAM_GIB="32"
readonly MIN_DISK_GIB="50"

reference_dir="${REFERENCE_DIR:-$DEFAULT_REFERENCE_DIR}"
threads="${THREADS:-12}"
dry_run=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [--reference-dir DIR] [--threads N] [--dry-run]

Download and verify the Ensembl release-${ENSEMBL_RELEASE} GRCm39 primary
assembly, reuse/verify the matching GTF prepared by workflow step 01, and build a
full-genome STAR ${ASSEMBLY} index for 150-bp RNA-seq reads.
The index is dataset-independent and is shared by PRJNA548952 and GSE214315.

Options:
  --reference-dir DIR  Reference root
                       (default: ${DEFAULT_REFERENCE_DIR})
  --threads N          STAR genome-generation threads (default: 12)
  --dry-run            Check tools/resources and show planned paths without downloading
  -h, --help           Show this help

Environment:
  REFERENCE_DIR        Alternative reference root
  THREADS              Alternative thread count
  PRDX4T_ENV_BIN       Directory containing STAR and samtools

Expected resources:
  RAM:  at least ${MIN_RAM_GIB} GiB
  Disk: at least ${MIN_DISK_GIB} GiB free in the reference destination
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --reference-dir)
            (($# >= 2)) || die "--reference-dir requires a value"
            reference_dir="$2"
            shift 2
            ;;
        --threads)
            (($# >= 2)) || die "--threads requires a value"
            threads="$2"
            shift 2
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

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"

env_bin="${PRDX4T_ENV_BIN:-}"
if [[ -n "$env_bin" ]]; then
    :
elif command -v STAR >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then
    env_bin="$(dirname -- "$(command -v STAR)")"
elif [[ -x "${DEFAULT_ENV_BIN}/STAR" && -x "${DEFAULT_ENV_BIN}/samtools" ]]; then
    env_bin="$DEFAULT_ENV_BIN"
else
    die "prdx4t-analysis tools were not found; run script_environment_create_prdx4t_analysis.sh first"
fi

readonly STAR_BIN="${env_bin}/STAR"
readonly SAMTOOLS="${env_bin}/samtools"
[[ -x "$STAR_BIN" ]] || die "STAR is not executable: $STAR_BIN"
[[ -x "$SAMTOOLS" ]] || die "samtools is not executable: $SAMTOOLS"
command -v curl >/dev/null 2>&1 || die "curl is required"
command -v gzip >/dev/null 2>&1 || die "gzip is required"
command -v sum >/dev/null 2>&1 || die "sum is required"
command -v cmp >/dev/null 2>&1 || die "cmp is required"

readonly SOURCE_DIR="${reference_dir}/source_ensembl_release_${ENSEMBL_RELEASE}"
readonly GENOME_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
readonly GENOME_FA="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
readonly GTF_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.110.gtf.gz"
readonly GTF_FILE="${SOURCE_DIR}/Mus_musculus.GRCm39.110.gtf"
readonly STAR_DIR="${reference_dir}/star_GRCm39_Ensembl110_full"
readonly MANIFEST="${reference_dir}/GRCm39_STAR_reference_manifest.tsv"

mkdir -p "$SOURCE_DIR" "$reference_dir"

available_disk_gib() {
    df -Pk "$reference_dir" | awk 'NR == 2 {printf "%d", $4/1048576}'
}

physical_ram_gib() {
    local bytes=""
    bytes="$(sysctl -n hw.memsize 2>/dev/null || true)"
    if [[ -n "$bytes" ]]; then
        awk -v b="$bytes" 'BEGIN {printf "%d", b/1073741824}'
        return
    fi
    if command -v hostinfo >/dev/null 2>&1; then
        hostinfo 2>/dev/null | awk '/Primary memory available:/ {printf "%d", $4}'
        return
    fi
    printf '0'
}

ram_gib="$(physical_ram_gib)"
disk_gib="$(available_disk_gib)"
printf 'Detected RAM:       %s GiB\n' "${ram_gib:-unknown}"
printf 'Available disk:     %s GiB\n' "$disk_gib"
printf 'STAR threads:       %s\n' "$threads"
star_version="$("$STAR_BIN" --version)"
printf 'STAR version:       %s\n' "$star_version"
printf 'Reference root:     %s\n' "$reference_dir"

if [[ -n "$ram_gib" && "$ram_gib" != "0" ]] && ((ram_gib < MIN_RAM_GIB)); then
    die "STAR genome generation requires at least ${MIN_RAM_GIB} GiB RAM"
fi
((disk_gib >= MIN_DISK_GIB)) || die "at least ${MIN_DISK_GIB} GiB free disk is required"

if ((dry_run)); then
    cat <<EOF

Dry-run plan:
  Genome download: ${GENOME_URL}
  Reused GTF:      ${GTF_GZ}
  Genome FASTA:    ${GENOME_FA}
  STAR index:      ${STAR_DIR}
  sjdbOverhang:    ${SJDB_OVERHANG}

No downloads or index generation were performed.
EOF
    exit 0
fi

sha256_file() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        die "neither shasum nor sha256sum is available"
    fi
}

verify_ensembl_sum() {
    local file="$1"
    local expected_sum="$2"
    local expected_blocks="$3"
    local got_sum got_blocks ignored
    read -r got_sum got_blocks ignored < <(sum "$file")
    [[ "$got_sum" == "$expected_sum" && "$got_blocks" == "$expected_blocks" ]] || \
        die "Ensembl checksum failed for $file (got: $got_sum $got_blocks; expected: $expected_sum $expected_blocks)"
}

download_reference() {
    local url="$1"
    local destination="$2"
    local expected_sum="$3"
    local expected_blocks="$4"
    local partial="${destination}.part"

    if [[ -s "$destination" ]]; then
        verify_ensembl_sum "$destination" "$expected_sum" "$expected_blocks"
        gzip -t "$destination"
        printf 'Verified existing download: %s\n' "$destination"
        return
    fi

    printf 'Downloading %s\n' "$url"
    curl --fail --location --retry 5 --retry-delay 5 \
        --continue-at - --output "$partial" "$url"
    verify_ensembl_sum "$partial" "$expected_sum" "$expected_blocks"
    gzip -t "$partial"
    mv -f "$partial" "$destination"
    printf 'Downloaded and verified: %s\n' "$destination"
}

decompress_if_needed() {
    local source_gz="$1"
    local destination="$2"
    local source_sha marker marker_tmp temporary
    source_sha="$(sha256_file "$source_gz")"
    marker="${destination}.source_gz.sha256"
    if [[ -s "$destination" ]]; then
        if [[ -s "$marker" && "$(tr -d '[:space:]' < "$marker")" == "$source_sha" ]]; then
            printf 'Verified decompressed reference provenance: %s\n' "$destination"
            return
        fi
        printf 'Validating existing decompressed reference against %s\n' "$source_gz"
        gzip -dc "$source_gz" | cmp -s - "$destination" || \
            die "decompressed reference does not match its verified source: $destination"
        marker_tmp="${marker}.tmp"
        printf '%s\n' "$source_sha" > "$marker_tmp"
        mv -f "$marker_tmp" "$marker"
        return
    fi
    temporary="${destination}.tmp"
    printf 'Decompressing %s\n' "$source_gz"
    gzip -dc "$source_gz" > "$temporary"
    [[ -s "$temporary" ]] || die "decompression produced an empty file: $temporary"
    mv -f "$temporary" "$destination"
    marker_tmp="${marker}.tmp"
    printf '%s\n' "$source_sha" > "$marker_tmp"
    mv -f "$marker_tmp" "$marker"
}

manifest_value() {
    local key="$1"
    [[ -r "$MANIFEST" ]] || return 1
    awk -F '\t' -v wanted="$key" '$1 == wanted {print $2; found=1; exit} END {if (!found) exit 1}' "$MANIFEST"
}

download_reference "$GENOME_URL" "$GENOME_GZ" "$GENOME_SUM" "$GENOME_BLOCKS"

if [[ ! -s "$GTF_GZ" ]]; then
    printf 'The release-110 GTF from workflow step 01 is absent; downloading it now.\n'
    download_reference "$GTF_URL" "$GTF_GZ" "$GTF_SUM" "$GTF_BLOCKS"
else
    verify_ensembl_sum "$GTF_GZ" "$GTF_SUM" "$GTF_BLOCKS"
    gzip -t "$GTF_GZ"
    printf 'Verified GTF from workflow step 01: %s\n' "$GTF_GZ"
fi

decompress_if_needed "$GENOME_GZ" "$GENOME_FA"
decompress_if_needed "$GTF_GZ" "$GTF_FILE"

printf 'Indexing the genome FASTA with samtools...\n'
"$SAMTOOLS" faidx "$GENOME_FA"
[[ -s "${GENOME_FA}.fai" ]] || die "samtools FASTA index was not created"

chromosome_count="$(awk '$1 ~ /^(1[0-9]?|[2-9]|X|Y|MT)$/ {n++} END {print n+0}' "${GENOME_FA}.fai")"
has_x="$(awk '$1 == "X" {print $2}' "${GENOME_FA}.fai")"
[[ -n "$has_x" ]] || die "chromosome X is missing from the primary assembly"
((chromosome_count >= 21)) || die "too few assembled chromosomes were detected: $chromosome_count"
printf 'Detected assembled chromosomes: %s; chromosome X length: %s bp\n' "$chromosome_count" "$has_x"

star_complete=1
for required in Genome SA SAindex chrLength.txt chrName.txt chrStart.txt genomeParameters.txt sjdbInfo.txt; do
    [[ -s "${STAR_DIR}/${required}" ]] || star_complete=0
done

genome_source_sha="$(sha256_file "$GENOME_GZ")"
gtf_source_sha="$(sha256_file "$GTF_GZ")"
if ((star_complete)); then
    grep -Eq "sjdbOverhang[[:space:]]+${SJDB_OVERHANG}" "${STAR_DIR}/genomeParameters.txt" || star_complete=0
    grep -Eq 'genomeSAindexNbases[[:space:]]+14' "${STAR_DIR}/genomeParameters.txt" || star_complete=0
    [[ "$(manifest_value assembly 2>/dev/null || true)" == "$ASSEMBLY" ]] || star_complete=0
    [[ "$(manifest_value ensembl_release 2>/dev/null || true)" == "$ENSEMBL_RELEASE" ]] || star_complete=0
    [[ "$(manifest_value STAR_version 2>/dev/null || true)" == "$star_version" ]] || star_complete=0
    [[ "$(manifest_value read_length 2>/dev/null || true)" == "$READ_LENGTH" ]] || star_complete=0
    [[ "$(manifest_value sjdbOverhang 2>/dev/null || true)" == "$SJDB_OVERHANG" ]] || star_complete=0
    [[ "$(manifest_value genomeSAindexNbases 2>/dev/null || true)" == "14" ]] || star_complete=0
    [[ "$(manifest_value genome_fasta_gz.sha256 2>/dev/null || true)" == "$genome_source_sha" ]] || star_complete=0
    [[ "$(manifest_value annotation_gtf_gz.sha256 2>/dev/null || true)" == "$gtf_source_sha" ]] || star_complete=0
fi

if ((star_complete)); then
    printf 'Existing full-genome STAR index is complete; skipping rebuild.\n'
else
    if [[ -d "$STAR_DIR" ]]; then
        printf 'Existing STAR index does not match the current reference fingerprint; rebuilding.\n'
    fi
    star_build_dir="$(mktemp -d "${reference_dir}/.star_GRCm39_build.XXXXXX")"
    printf '\nBuilding the full-genome STAR index. This can take considerable RAM and time.\n'
    if ! "$STAR_BIN" \
        --runMode genomeGenerate \
        --runThreadN "$threads" \
        --genomeDir "$star_build_dir" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang "$SJDB_OVERHANG" \
        --genomeSAindexNbases 14 \
        2>&1 | tee "${star_build_dir}/genomeGenerate.log"; then
        printf 'STAR genome generation failed; build files retained at %s\n' "$star_build_dir" >&2
        exit 1
    fi

    for required in Genome SA SAindex chrLength.txt chrName.txt chrStart.txt genomeParameters.txt sjdbInfo.txt; do
        [[ -s "${star_build_dir}/${required}" ]] || die "STAR index is incomplete; missing $required"
    done

    if [[ -d "$STAR_DIR" ]]; then
        incomplete_dir="${STAR_DIR}.incomplete_$(date '+%Y%m%d_%H%M%S')"
        mv "$STAR_DIR" "$incomplete_dir"
        printf 'Preserved the previous incomplete index at %s\n' "$incomplete_dir"
    fi
    mv "$star_build_dir" "$STAR_DIR"
fi

manifest_tmp="${MANIFEST}.tmp"
{
    printf 'field\tvalue\n'
    printf 'species\tMus musculus\n'
    printf 'assembly\t%s\n' "$ASSEMBLY"
    printf 'ensembl_release\t%s\n' "$ENSEMBL_RELEASE"
    printf 'reference_scope\tfull primary assembly\n'
    printf 'genome_fasta_url\t%s\n' "$GENOME_URL"
    printf 'annotation_gtf_url\t%s\n' "$GTF_URL"
    printf 'STAR_version\t%s\n' "$star_version"
    printf 'read_length\t%s\n' "$READ_LENGTH"
    printf 'sjdbOverhang\t%s\n' "$SJDB_OVERHANG"
    printf 'genomeSAindexNbases\t14\n'
    printf 'assembled_chromosome_count\t%s\n' "$chromosome_count"
    printf 'genome_fasta_gz.sha256\t%s\n' "$genome_source_sha"
    printf 'annotation_gtf_gz.sha256\t%s\n' "$gtf_source_sha"
    printf 'STAR_index_directory\t%s\n' "$STAR_DIR"
} > "$manifest_tmp"
mv -f "$manifest_tmp" "$MANIFEST"

printf '\nFull-genome STAR reference preparation complete.\n'
printf 'STAR index: %s\n' "$STAR_DIR"
printf 'Manifest:   %s\n' "$MANIFEST"
