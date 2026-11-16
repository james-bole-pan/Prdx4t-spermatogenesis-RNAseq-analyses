#!/usr/bin/env bash
set -Eeuo pipefail

readonly ENSEMBL_RELEASE="110"
readonly ASSEMBLY="GRCm39"
readonly GENE_ID="ENSMUSG00000025289"
readonly GENE_NAME="Prdx4"
readonly PRDX4_TRANSCRIPT_ID="ENSMUST00000026328"
readonly PRDX4T_TRANSCRIPT_ID="ENSMUST00000130349"

readonly CHRX_URL="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.X.fa.gz"
readonly GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"

# Ensembl CHECKSUMS values (BSD sum checksum and 1024-byte block count).
readonly CHRX_SUM="59093"
readonly CHRX_BLOCKS="48414"
readonly GTF_SUM="42178"
readonly GTF_BLOCKS="31546"

readonly DEFAULT_OUTPUT_DIR="${PRDX4T_REFERENCE_DIR:-${SRA_ROOT:-${PWD}/data}/reference}"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"

output_dir="${REFERENCE_DIR:-$DEFAULT_OUTPUT_DIR}"
threads="${THREADS:-8}"
prefilter_flank="${PREFILTER_FLANK:-10000}"
build_index=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [--output-dir DIR] [--threads N] [--prefilter-flank N]
       [--skip-bowtie2-index]

Prepare the shared PRDX4 locus reference from Ensembl release 110:
  - GRCm39 Prdx4 locus FASTA with configurable flanking sequence
  - original chromosome-X locus GTF
  - reference-derived PRDX4/PRDX4t feature-coordinate table
  - FASTA index
  - Bowtie2 locus index (unless explicitly skipped)

The resulting reference is dataset-independent and can be reused by both
PRJNA548952 and GSE214315.

Downloaded source files are checked against Ensembl's published checksums.
Generated locus files are validated for contig name, length, gene identity,
coordinate bounds, and non-empty content.

Options:
  --output-dir DIR        Reference output directory
                          (default: ${DEFAULT_OUTPUT_DIR})
  --threads N             Bowtie2-build threads (default: 8)
  --prefilter-flank N     Bases added on each side of the GTF gene interval
                          (default: 10000)
  --skip-bowtie2-index    Prepare reference files without building the index
  -h, --help              Show this help

Environment:
  REFERENCE_DIR           Alternative reference output directory
  THREADS                 Alternative thread count
  PREFILTER_FLANK         Alternative prefilter flank size
  PRDX4T_ENV_BIN          Directory containing samtools and bowtie2-build
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
        --threads)
            (($# >= 2)) || die "--threads requires a value"
            threads="$2"
            shift 2
            ;;
        --prefilter-flank)
            (($# >= 2)) || die "--prefilter-flank requires a value"
            prefilter_flank="$2"
            shift 2
            ;;
        --skip-bowtie2-index)
            build_index=0
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
[[ "$prefilter_flank" =~ ^[0-9]+$ ]] || die "--prefilter-flank must be a non-negative integer"

env_bin="${PRDX4T_ENV_BIN:-}"
if [[ -n "$env_bin" ]]; then
    :
elif command -v samtools >/dev/null 2>&1 && command -v bowtie2-build >/dev/null 2>&1; then
    env_bin="$(dirname -- "$(command -v samtools)")"
elif [[ -x "${DEFAULT_ENV_BIN}/samtools" && -x "${DEFAULT_ENV_BIN}/bowtie2-build" ]]; then
    env_bin="$DEFAULT_ENV_BIN"
else
    die "prdx4t-analysis tools were not found; run script_environment_create_prdx4t_analysis.sh first"
fi

readonly SAMTOOLS="${env_bin}/samtools"
readonly BOWTIE2_BUILD="${env_bin}/bowtie2-build"
readonly PYTHON="${env_bin}/python"
[[ -x "$SAMTOOLS" ]] || die "samtools is not executable: $SAMTOOLS"
[[ -x "$PYTHON" ]] || die "Python is not executable: $PYTHON"
if ((build_index)); then
    [[ -x "$BOWTIE2_BUILD" ]] || die "bowtie2-build is not executable: $BOWTIE2_BUILD"
fi
command -v curl >/dev/null 2>&1 || die "curl is required"
command -v gzip >/dev/null 2>&1 || die "gzip is required"
command -v sum >/dev/null 2>&1 || die "sum is required"

readonly SOURCE_DIR="${output_dir}/source_ensembl_release_${ENSEMBL_RELEASE}"
readonly CHRX_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.chromosome.X.fa.gz"
readonly CHRX_FA="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.chromosome.X.fa"
readonly FULL_GTF_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.110.gtf.gz"
readonly LOCUS_FASTA="${output_dir}/prdx4_locus.fa"
readonly LOCUS_FAI="${LOCUS_FASTA}.fai"
readonly LOCUS_GTF="${output_dir}/prdx4_locus.gtf"
readonly FEATURE_TABLE="${output_dir}/PRDX4_locus_features.tsv"
readonly BT2_DIR="${output_dir}/bt2_locus"
readonly BT2_PREFIX="${BT2_DIR}/locus"
readonly BT2_SOURCE_SHA="${BT2_DIR}/source_fasta.sha256"
readonly MANIFEST="${output_dir}/PRDX4_locus_reference_manifest.tsv"

mkdir -p "$SOURCE_DIR" "$output_dir"

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

validate_locus_files() {
    [[ -s "$LOCUS_FASTA" && -s "$LOCUS_FAI" && -s "$LOCUS_GTF" ]] || return 1

    local fasta_header fasta_length fai_name fai_length
    fasta_header="$(sed -n '1{s/^>//;p;}' "$LOCUS_FASTA")"
    fasta_length="$(awk '!/^>/ {gsub(/[[:space:]]/, ""); n+=length($0)} END {print n+0}' "$LOCUS_FASTA")"
    read -r fai_name fai_length _ < "$LOCUS_FAI"
    [[ "$fasta_header" == "prdx4locus" ]] || return 1
    [[ "$fasta_length" == "$LOCUS_LENGTH" ]] || return 1
    [[ "$fai_name" == "prdx4locus" && "$fai_length" == "$LOCUS_LENGTH" ]] || return 1

    awk -F '\t' -v gene="$GENE_ID" -v chrom="$CHROM" \
        -v gene_start="$GENE_START" -v gene_end="$GENE_END" '
        BEGIN {found=0}
        $0 !~ /^#/ && NF >= 9 && index($0, "gene_id \"" gene "\"") {
            if ($1 != chrom || $4 < gene_start || $5 > gene_end) exit 1
            found=1
        }
        END {exit !found}
    ' "$LOCUS_GTF" || return 1
}

write_feature_table() {
    local temporary="${FEATURE_TABLE}.tmp"
    "$PYTHON" - "$LOCUS_GTF" "$temporary" "$GENE_ID" \
        "$PRDX4_TRANSCRIPT_ID" "$PRDX4T_TRANSCRIPT_ID" <<'PY'
import csv
import re
import sys
from pathlib import Path


gtf_path, output_path, gene_id, prdx4_id, prdx4t_id = sys.argv[1:]


def attributes(text):
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]*)"', text))


records = []
with Path(gtf_path).open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = attributes(fields[8])
        if attrs.get("gene_id") != gene_id:
            continue
        records.append(
            {
                "chrom": fields[0],
                "type": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "transcript_id": attrs.get("transcript_id", "."),
                "exon_number": attrs.get("exon_number", "."),
            }
        )

gene_records = [record for record in records if record["type"] == "gene"]
if len(gene_records) != 1:
    raise SystemExit(f"ERROR: expected one {gene_id} gene record, found {len(gene_records)}")
gene = gene_records[0]


def exon(transcript_id, exon_number):
    matches = [
        record
        for record in records
        if record["type"] == "exon"
        and record["transcript_id"] == transcript_id
        and record["exon_number"] == str(exon_number)
    ]
    unique = {(record["chrom"], record["start"], record["end"], record["strand"]) for record in matches}
    if len(unique) != 1:
        raise SystemExit(
            f"ERROR: expected one coordinate for exon {exon_number} of {transcript_id}, found {len(unique)}"
        )
    chrom, start, end, strand = unique.pop()
    if chrom != gene["chrom"] or strand != gene["strand"]:
        raise SystemExit(f"ERROR: inconsistent chromosome or strand for {transcript_id} exon {exon_number}")
    return start, end


prdx4_exon1 = exon(prdx4_id, 1)
prdx4t_exon1 = exon(prdx4t_id, 1)
canonical_body_exons = {number: exon(prdx4_id, number) for number in range(2, 8)}

# Ensembl release 110 marks Prdx4-203 as mRNA_end_NF/cds_end_NF and only
# annotates its alternative exon 1 plus downstream exons 2-5. Verify every
# downstream exon that Ensembl does annotate against the canonical transcript.
# The experimentally established PRDX4t model is exon 1A followed by canonical
# exons 2-7 (Yim et al. 2011, doi:10.1074/jbc.M111.257220).
prdx4t_body_numbers = sorted(
    {
        int(record["exon_number"])
        for record in records
        if record["type"] == "exon"
        and record["transcript_id"] == prdx4t_id
        and record["exon_number"].isdigit()
        and int(record["exon_number"]) >= 2
    }
)
if not prdx4t_body_numbers or prdx4t_body_numbers[0] != 2:
    raise SystemExit("ERROR: PRDX4t lacks an annotated exon 2")
for number in prdx4t_body_numbers:
    if number not in canonical_body_exons:
        raise SystemExit(f"ERROR: unexpected PRDX4t downstream exon number: {number}")
    if exon(prdx4t_id, number) != canonical_body_exons[number]:
        raise SystemExit(
            f"ERROR: PRDX4 and PRDX4t do not share the same annotated exon {number}"
        )
ordered_exon1 = sorted((prdx4_exon1, prdx4t_exon1))
gap_start = ordered_exon1[0][1] + 1
gap_end = ordered_exon1[1][0] - 1
if gap_start > gap_end:
    raise SystemExit("ERROR: PRDX4 and PRDX4t exon 1 intervals overlap; no inter-exon gap exists")

rows = [
    ("gene", gene["start"], gene["end"], ".", ".", "GTF gene feature"),
    ("prdx4t_exon1", *prdx4t_exon1, prdx4t_id, "1", "GTF exon feature"),
    ("prdx4_exon1", *prdx4_exon1, prdx4_id, "1", "GTF exon feature"),
    ("gap_between_tss", gap_start, gap_end, ".", ".", "derived gap between exon 1 intervals"),
]

# Exons 2-7 are emitted individually rather than as one span. The span from the
# start of exon 7 to the end of exon 2 is mostly intron, so a span would count
# intronic and unspliced signal as shared-body expression. Downstream code takes
# the union of these intervals. Ensembl's
# incomplete Prdx4-203 model directly confirms exons 2-5; exons 6-7 follow the
# experimentally established PRDX4t transcript structure.
for number, (start, end) in canonical_body_exons.items():
    if number in prdx4t_body_numbers:
        transcript_id = f"{prdx4_id},{prdx4t_id}"
        derivation = "matching GTF exon features"
    else:
        transcript_id = prdx4_id
        derivation = (
            "canonical GTF exon; shared with Prdx4t per Yim et al. 2011 "
            "(doi:10.1074/jbc.M111.257220; Prdx4-203 is mRNA_end_NF)"
        )
    rows.append(
        (f"shared_exon{number}", start, end, transcript_id, str(number), derivation)
    )

fieldnames = [
    "feature_id",
    "chrom",
    "start",
    "end",
    "strand",
    "gene_id",
    "transcript_id",
    "exon_number",
    "derivation",
]
with Path(output_path).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for feature_id, start, end, transcript_id, exon_number, derivation in rows:
        writer.writerow(
            {
                "feature_id": feature_id,
                "chrom": gene["chrom"],
                "start": start,
                "end": end,
                "strand": gene["strand"],
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "exon_number": exon_number,
                "derivation": derivation,
            }
        )
PY
    mv -f "$temporary" "$FEATURE_TABLE"
}

validate_feature_table() {
    [[ -s "$FEATURE_TABLE" ]] || return 1
    awk -F '\t' -v chrom="$CHROM" -v strand="$GENE_STRAND" \
        -v locus_start="$LOCUS_START" -v locus_end="$LOCUS_END" '
        NR == 1 {
            if ($1 != "feature_id" || $2 != "chrom" || $3 != "start" || $4 != "end") exit 1
            next
        }
        {
            if (seen[$1]++) exit 1
            if ($2 != chrom || $3 < locus_start || $4 < $3 || $4 > locus_end || $5 != strand) exit 1
        }
        END {
            required["gene"]=1
            required["prdx4t_exon1"]=1
            required["prdx4_exon1"]=1
            required["gap_between_tss"]=1
            required["shared_exon2"]=1
            required["shared_exon3"]=1
            required["shared_exon4"]=1
            required["shared_exon5"]=1
            required["shared_exon6"]=1
            required["shared_exon7"]=1
            for (feature in required) if (!(feature in seen)) exit 1
        }
    ' "$FEATURE_TABLE"
}

download_reference "$GTF_URL" "$FULL_GTF_GZ" "$GTF_SUM" "$GTF_BLOCKS"

gene_definition="$(gzip -dc "$FULL_GTF_GZ" | awk -F '\t' -v gene="$GENE_ID" '
    $0 !~ /^#/ && $3 == "gene" && index($0, "gene_id \"" gene "\"") {
        print $1, $4, $5, $7
        found++
    }
    END {if (found != 1) exit 1}
')" || die "could not derive a unique $GENE_ID gene interval from the Ensembl GTF"
read -r CHROM GENE_START GENE_END GENE_STRAND <<< "$gene_definition"
[[ "$CHROM" == "X" ]] || die "$GENE_ID is on chromosome $CHROM, but the configured FASTA source is chromosome X"
[[ "$GENE_START" =~ ^[0-9]+$ && "$GENE_END" =~ ^[0-9]+$ && "$GENE_START" -le "$GENE_END" ]] || \
    die "invalid gene interval derived from the Ensembl GTF: $gene_definition"
[[ "$GENE_STRAND" == "+" || "$GENE_STRAND" == "-" ]] || die "invalid gene strand: $GENE_STRAND"

LOCUS_START=$((GENE_START - prefilter_flank))
((LOCUS_START >= 1)) || LOCUS_START=1
readonly CHROM GENE_START GENE_END GENE_STRAND LOCUS_START
readonly LOCUS_END=$((GENE_END + prefilter_flank))
readonly LOCUS_LENGTH=$((LOCUS_END - LOCUS_START + 1))

printf 'GTF-derived Prdx4 interval: %s:%s-%s (%s)\n' "$CHROM" "$GENE_START" "$GENE_END" "$GENE_STRAND"
printf 'Bowtie2 prefilter interval: %s:%s-%s (%s-bp flanks)\n' \
    "$CHROM" "$LOCUS_START" "$LOCUS_END" "$prefilter_flank"

download_reference "$CHRX_URL" "$CHRX_GZ" "$CHRX_SUM" "$CHRX_BLOCKS"

references_complete=0
if validate_locus_files; then
    references_complete=1
    printf 'Existing PRDX4 locus files passed GTF-derived structural validation.\n'
fi

if ((references_complete == 0)); then

    if [[ ! -s "$CHRX_FA" ]]; then
        chr_tmp="${CHRX_FA}.tmp"
        gzip -dc "$CHRX_GZ" > "$chr_tmp"
        mv -f "$chr_tmp" "$CHRX_FA"
    fi
    "$SAMTOOLS" faidx "$CHRX_FA"
    chr_length="$(awk -F '\t' '$1 == "X" {print $2}' "${CHRX_FA}.fai")"
    [[ -n "$chr_length" ]] || die "chromosome X was not found in $CHRX_FA"
    printf 'Chromosome X length: %s bp\n' "$chr_length"

    fasta_tmp="${LOCUS_FASTA}.tmp"
    "$SAMTOOLS" faidx "$CHRX_FA" "${CHROM}:${LOCUS_START}-${LOCUS_END}" \
        | awk 'NR == 1 {$0=">prdx4locus"} {print}' > "$fasta_tmp"
    observed_length="$(awk '!/^>/ {gsub(/[[:space:]]/, ""); n+=length($0)} END {print n+0}' "$fasta_tmp")"
    [[ "$observed_length" == "$LOCUS_LENGTH" ]] || \
        die "unexpected locus FASTA length: $observed_length (expected $LOCUS_LENGTH)"
    mv -f "$fasta_tmp" "$LOCUS_FASTA"
    "$SAMTOOLS" faidx "$LOCUS_FASTA"

    gtf_tmp="${LOCUS_GTF}.tmp"
    gzip -dc "$FULL_GTF_GZ" \
        | awk -v gene="$GENE_ID" '$0 !~ /^#/ && index($0, "gene_id \"" gene "\"") {print}' \
        > "$gtf_tmp"
    [[ -s "$gtf_tmp" ]] || die "$GENE_ID records were not found in the Ensembl GTF"
    mv -f "$gtf_tmp" "$LOCUS_GTF"

    validate_locus_files || die "generated PRDX4 locus files failed structural validation"
    printf 'Generated PRDX4 locus files passed structural validation.\n'
fi

write_feature_table
validate_feature_table || die "GTF-derived PRDX4 feature table failed validation"
printf 'Wrote GTF-derived feature table: %s\n' "$FEATURE_TABLE"

bt2_complete=1
for suffix in 1 2 3 4; do
    [[ -s "${BT2_PREFIX}.${suffix}.bt2" || -s "${BT2_PREFIX}.${suffix}.bt2l" ]] || bt2_complete=0
done
for suffix in 1 2; do
    [[ -s "${BT2_PREFIX}.rev.${suffix}.bt2" || -s "${BT2_PREFIX}.rev.${suffix}.bt2l" ]] || bt2_complete=0
done
current_fasta_sha="$(sha256_file "$LOCUS_FASTA")"
if [[ ! -s "$BT2_SOURCE_SHA" || "$(tr -d '[:space:]' < "$BT2_SOURCE_SHA")" != "$current_fasta_sha" ]]; then
    bt2_complete=0
fi

if ((build_index)); then
    if ((bt2_complete)); then
        printf 'Existing Bowtie2 locus index is complete; skipping rebuild.\n'
    else
        mkdir -p "$BT2_DIR"
        bt2_tmp="$(mktemp -d "${output_dir}/.bt2_locus.XXXXXX")"
        "$BOWTIE2_BUILD" --threads "$threads" "$LOCUS_FASTA" "${bt2_tmp}/locus"
        printf '%s\n' "$current_fasta_sha" > "${bt2_tmp}/source_fasta.sha256"
        find "$bt2_tmp" -type f -name 'locus*.bt2*' -exec mv -f {} "$BT2_DIR"/ \;
        mv -f "${bt2_tmp}/source_fasta.sha256" "$BT2_SOURCE_SHA"
        rmdir "$bt2_tmp"
        printf 'Built Bowtie2 locus index: %s\n' "$BT2_PREFIX"
    fi
fi

manifest_tmp="${MANIFEST}.tmp"
{
    printf 'field\tvalue\n'
    printf 'species\tMus musculus\n'
    printf 'assembly\t%s\n' "$ASSEMBLY"
    printf 'ensembl_release\t%s\n' "$ENSEMBL_RELEASE"
    printf 'gene_id\t%s\n' "$GENE_ID"
    printf 'gene_name\t%s\n' "$GENE_NAME"
    printf 'gene_interval_from_gtf\t%s:%s-%s\n' "$CHROM" "$GENE_START" "$GENE_END"
    printf 'gene_strand_from_gtf\t%s\n' "$GENE_STRAND"
    printf 'prdx4_transcript_id\t%s\n' "$PRDX4_TRANSCRIPT_ID"
    printf 'prdx4t_transcript_id\t%s\n' "$PRDX4T_TRANSCRIPT_ID"
    printf 'prdx4t_biological_model\talternative exon 1A followed by shared exons 2-7\n'
    printf 'prdx4t_model_source\tYim et al. 2011; doi:10.1074/jbc.M111.257220\n'
    printf 'prefilter_flank_bp\t%s\n' "$prefilter_flank"
    printf 'genomic_interval\t%s:%s-%s\n' "$CHROM" "$LOCUS_START" "$LOCUS_END"
    printf 'locus_length_bp\t%s\n' "$LOCUS_LENGTH"
    printf 'chromosome_fasta_url\t%s\n' "$CHRX_URL"
    printf 'annotation_gtf_url\t%s\n' "$GTF_URL"
    printf 'prdx4_locus.fa.sha256\t%s\n' "$(sha256_file "$LOCUS_FASTA")"
    printf 'prdx4_locus.fa.fai.sha256\t%s\n' "$(sha256_file "$LOCUS_FAI")"
    printf 'prdx4_locus.gtf.sha256\t%s\n' "$(sha256_file "$LOCUS_GTF")"
    printf 'PRDX4_locus_features.tsv.sha256\t%s\n' "$(sha256_file "$FEATURE_TABLE")"
} > "$manifest_tmp"
mv -f "$manifest_tmp" "$MANIFEST"

printf '\nPRDX4 locus preparation complete.\n'
printf 'Reference directory: %s\n' "$output_dir"
printf 'Manifest:            %s\n' "$MANIFEST"
