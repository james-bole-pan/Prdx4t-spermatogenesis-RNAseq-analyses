#!/usr/bin/env bash
set -Eeuo pipefail

readonly ENSEMBL_RELEASE="110"
readonly ASSEMBLY="GRCm39"

readonly GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
readonly GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"

# Ensembl CHECKSUMS values (BSD sum checksum and 1024-byte block count).
readonly GENOME_SUM="16996"
readonly GENOME_BLOCKS="787519"
readonly GTF_SUM="42178"
readonly GTF_BLOCKS="31546"

readonly DEFAULT_OUTPUT_DIR="${PRDX4T_REFERENCE_DIR:-${SRA_ROOT:-${PWD}/data}/reference}"
readonly DEFAULT_ENV_BIN="${CONDA_PREFIX:-${HOME}/miniforge3/envs/prdx4t-analysis}/bin"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly DEFAULT_PANEL="${SCRIPT_DIR}/panels/prdx4t_coexpression.tsv"

output_dir="${REFERENCE_DIR:-$DEFAULT_OUTPUT_DIR}"
panel_file="${GENE_PANEL:-$DEFAULT_PANEL}"
threads="${THREADS:-8}"
prefilter_flank="${PREFILTER_FLANK:-10000}"
build_index=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [--panel FILE] [--output-dir DIR] [--threads N]
       [--prefilter-flank N] [--skip-bowtie2-index]

Prepare a multi-gene panel reference from Ensembl release ${ENSEMBL_RELEASE}:
  - one locus FASTA contig per panel gene, with configurable flanking sequence
  - the panel subset of the Ensembl GTF
  - a panel feature table with PRDX4/PRDX4t first-exon features and exon-union
    intervals for every gene-level member
  - FASTA index
  - one Bowtie2 prefilter index spanning every panel locus (unless skipped)

This generalizes script_workflow_01_prepare_PRDX4_locus.sh from the single
PRDX4 locus to a gene panel, so that steps 03-05 can quantify PRDX4t together
with co-expression partners such as PCMT1 and CRBN. It writes gene_panel.*
files and leaves the existing prdx4_locus.* files untouched, so results
produced against the PRDX4-only reference stay reproducible.

The PRDX4 rows of the panel feature table use the same feature identifiers,
coordinates and derivations as the single-locus script, so the validated
PRDX4t quantification logic reads them unchanged.

Downloaded source files are checked against Ensembl's published checksums.

Options:
  --panel FILE            Panel definition TSV
                          (default: ${DEFAULT_PANEL})
  --output-dir DIR        Reference output directory
                          (default: ${DEFAULT_OUTPUT_DIR})
  --threads N             Bowtie2-build threads (default: 8)
  --prefilter-flank N     Bases added on each side of each GTF gene interval
                          (default: 10000)
  --skip-bowtie2-index    Prepare reference files without building the index
  -h, --help              Show this help

Environment:
  REFERENCE_DIR           Alternative reference output directory
  GENE_PANEL              Alternative panel definition TSV
  THREADS                 Alternative thread count
  PREFILTER_FLANK         Alternative prefilter flank size
  PRDX4T_ENV_BIN          Directory containing samtools and bowtie2-build

Panel definition format (tab separated, '#' comments allowed):
  gene_name  gene_id  role  quantification  rationale
where quantification is prdx4_isoform (PRDX4 only) or gene_exonunion.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --panel)
            (($# >= 2)) || die "--panel requires a file"
            panel_file="$2"
            shift 2
            ;;
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
[[ -r "$panel_file" ]] || die "panel definition is not readable: $panel_file"

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
readonly GENOME_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
readonly GENOME_FA="${SOURCE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
readonly FULL_GTF_GZ="${SOURCE_DIR}/Mus_musculus.GRCm39.110.gtf.gz"
readonly PANEL_FASTA="${output_dir}/gene_panel.fa"
readonly PANEL_FAI="${PANEL_FASTA}.fai"
readonly PANEL_GTF="${output_dir}/gene_panel.gtf"
readonly PANEL_INTERVALS="${output_dir}/gene_panel_intervals.tsv"
readonly FEATURE_TABLE="${output_dir}/gene_panel_features.tsv"
readonly BT2_DIR="${output_dir}/bt2_panel"
readonly BT2_PREFIX="${BT2_DIR}/panel"
readonly BT2_SOURCE_SHA="${BT2_DIR}/source_fasta.sha256"
readonly MANIFEST="${output_dir}/gene_panel_reference_manifest.tsv"

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

download_reference "$GTF_URL" "$FULL_GTF_GZ" "$GTF_SUM" "$GTF_BLOCKS"
download_reference "$GENOME_URL" "$GENOME_GZ" "$GENOME_SUM" "$GENOME_BLOCKS"

if [[ ! -s "$GENOME_FA" ]]; then
    printf 'Decompressing the primary assembly (this takes a few minutes)...\n'
    genome_tmp="${GENOME_FA}.tmp"
    gzip -dc "$GENOME_GZ" > "$genome_tmp"
    mv -f "$genome_tmp" "$GENOME_FA"
fi
[[ -s "${GENOME_FA}.fai" ]] || "$SAMTOOLS" faidx "$GENOME_FA"

# Panel members and their GTF-derived intervals. Emitted as a table so the
# FASTA extraction, the GTF subset and the feature table all read one source.
printf 'Deriving panel gene intervals from the Ensembl GTF...\n'
intervals_tmp="${PANEL_INTERVALS}.tmp"
gzip -dc "$FULL_GTF_GZ" > "${output_dir}/.panel_gtf.tmp"
"$PYTHON" - "$panel_file" "${output_dir}/.panel_gtf.tmp" "${GENOME_FA}.fai" \
    "$prefilter_flank" "$intervals_tmp" <<'PY'
import csv
import re
import sys
from pathlib import Path

panel_path, gtf_path, fai_path, flank_text, output_path = sys.argv[1:]
flank = int(flank_text)

QUANTIFICATIONS = {"prdx4_isoform", "gene_exonunion"}


def read_panel(path):
    members = []
    seen_ids = set()
    seen_names = set()
    with Path(path).open() as handle:
        rows = [line for line in handle if not line.lstrip().startswith("#")]
    reader = csv.DictReader(rows, delimiter="\t")
    required = {"gene_name", "gene_id", "role", "quantification"}
    missing = sorted(required.difference(reader.fieldnames or []))
    if missing:
        raise SystemExit(f"ERROR: panel definition is missing columns: {', '.join(missing)}")
    for row in reader:
        gene_id = row["gene_id"].strip()
        gene_name = row["gene_name"].strip()
        quantification = row["quantification"].strip()
        if not re.fullmatch(r"ENSMUSG[0-9]{11}", gene_id):
            raise SystemExit(f"ERROR: invalid Ensembl gene id in panel: {gene_id!r}")
        if not re.fullmatch(r"[A-Za-z0-9._-]+", gene_name):
            raise SystemExit(f"ERROR: invalid gene name in panel: {gene_name!r}")
        if quantification not in QUANTIFICATIONS:
            raise SystemExit(
                f"ERROR: {gene_name} has unknown quantification {quantification!r}; "
                f"expected one of {sorted(QUANTIFICATIONS)}"
            )
        if gene_id in seen_ids or gene_name.lower() in seen_names:
            raise SystemExit(f"ERROR: duplicate panel member: {gene_name} / {gene_id}")
        seen_ids.add(gene_id)
        seen_names.add(gene_name.lower())
        members.append(
            {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "role": row["role"].strip(),
                "quantification": quantification,
                "rationale": (row.get("rationale") or "").strip(),
            }
        )
    if not members:
        raise SystemExit("ERROR: panel definition contains no genes")
    isoform_members = [m for m in members if m["quantification"] == "prdx4_isoform"]
    if len(isoform_members) != 1:
        raise SystemExit(
            "ERROR: exactly one panel member must use quantification prdx4_isoform, "
            f"found {len(isoform_members)}"
        )
    return members


members = read_panel(panel_path)
wanted = {member["gene_id"]: member for member in members}

gene_records = {}
with Path(gtf_path).open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        match = re.search(r'gene_id "([^"]*)"', fields[8])
        if match is None or match.group(1) not in wanted:
            continue
        gene_id = match.group(1)
        if gene_id in gene_records:
            raise SystemExit(f"ERROR: {gene_id} has more than one gene record in the GTF")
        name_match = re.search(r'gene_name "([^"]*)"', fields[8])
        gtf_name = name_match.group(1) if name_match else ""
        if gtf_name != wanted[gene_id]["gene_name"]:
            raise SystemExit(
                f"ERROR: {gene_id} is named {gtf_name!r} in the GTF but "
                f"{wanted[gene_id]['gene_name']!r} in the panel"
            )
        gene_records[gene_id] = fields

absent = sorted(set(wanted) - set(gene_records))
if absent:
    raise SystemExit(f"ERROR: panel genes absent from the Ensembl GTF: {', '.join(absent)}")

contig_lengths = {}
with Path(fai_path).open() as handle:
    for line in handle:
        name, length = line.split("\t")[:2]
        contig_lengths[name] = int(length)

rows = []
for member in members:
    fields = gene_records[member["gene_id"]]
    chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
    if chrom not in contig_lengths:
        raise SystemExit(f"ERROR: {member['gene_name']} is on contig {chrom}, absent from the assembly")
    if strand not in {"+", "-"} or start > end:
        raise SystemExit(f"ERROR: invalid GTF gene record for {member['gene_name']}")
    locus_start = max(1, start - flank)
    locus_end = min(contig_lengths[chrom], end + flank)
    rows.append(
        {
            "gene_name": member["gene_name"],
            "gene_id": member["gene_id"],
            "role": member["role"],
            "quantification": member["quantification"],
            "rationale": member["rationale"],
            "contig": f"{member['gene_name'].lower()}locus",
            "chrom": chrom,
            "strand": strand,
            "gene_start": start,
            "gene_end": end,
            "locus_start": locus_start,
            "locus_end": locus_end,
            "locus_length": locus_end - locus_start + 1,
        }
    )

# Overlapping panel loci would put the same genomic sequence into the prefilter
# index twice and make locus-BAM extraction ambiguous.
by_chrom = {}
for row in rows:
    by_chrom.setdefault(row["chrom"], []).append(row)
for chrom, group in by_chrom.items():
    group.sort(key=lambda item: item["locus_start"])
    for earlier, later in zip(group, group[1:]):
        if later["locus_start"] <= earlier["locus_end"]:
            raise SystemExit(
                f"ERROR: panel loci overlap on {chrom}: "
                f"{earlier['gene_name']} and {later['gene_name']}"
            )

fieldnames = [
    "gene_name", "gene_id", "role", "quantification", "contig", "chrom", "strand",
    "gene_start", "gene_end", "locus_start", "locus_end", "locus_length", "rationale",
]
with Path(output_path).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY
mv -f "$intervals_tmp" "$PANEL_INTERVALS"

printf '\n%-8s %-20s %-22s %s\n' "GENE" "GENE ID" "GENE INTERVAL" "PREFILTER LOCUS"
awk -F '\t' 'NR > 1 {printf "%-8s %-20s %-22s %s\n", $1, $2, $6":"$8"-"$9"("$7")", $6":"$10"-"$11}' \
    "$PANEL_INTERVALS"
printf '\n'

# Panel FASTA: one renamed contig per gene, in panel order.
fasta_tmp="${PANEL_FASTA}.tmp"
: > "$fasta_tmp"
while IFS=$'\t' read -r gene_name gene_id role quantification contig chrom strand \
    gene_start gene_end locus_start locus_end locus_length rationale; do
    [[ "$gene_name" == "gene_name" ]] && continue
    "$SAMTOOLS" faidx "$GENOME_FA" "${chrom}:${locus_start}-${locus_end}" \
        | awk -v name="$contig" 'NR == 1 {$0=">" name} {print}' >> "$fasta_tmp"
done < "$PANEL_INTERVALS"
[[ -s "$fasta_tmp" ]] || die "panel FASTA extraction produced no sequence"
mv -f "$fasta_tmp" "$PANEL_FASTA"
rm -f -- "$PANEL_FAI"
"$SAMTOOLS" faidx "$PANEL_FASTA"

# Every contig must have exactly the length the interval table promised.
while IFS=$'\t' read -r gene_name gene_id role quantification contig chrom strand \
    gene_start gene_end locus_start locus_end locus_length rationale; do
    [[ "$gene_name" == "gene_name" ]] && continue
    observed="$(awk -F '\t' -v name="$contig" '$1 == name {print $2}' "$PANEL_FAI")"
    [[ "$observed" == "$locus_length" ]] || \
        die "contig $contig has length ${observed:-none}, expected $locus_length"
done < "$PANEL_INTERVALS"
printf 'Panel FASTA contigs passed length validation: %s\n' "$PANEL_FASTA"

# GTF subset restricted to the panel genes.
gtf_tmp="${PANEL_GTF}.tmp"
awk -F '\t' 'NR > 1 {print $2}' "$PANEL_INTERVALS" > "${output_dir}/.panel_ids.tmp"
awk -F '\t' '
    NR == FNR {wanted[$1]=1; next}
    /^#/ {next}
    {
        if (match($9, /gene_id "[^"]*"/)) {
            id = substr($9, RSTART + 9, RLENGTH - 10)
            if (id in wanted) print
        }
    }
' "${output_dir}/.panel_ids.tmp" "${output_dir}/.panel_gtf.tmp" > "$gtf_tmp"
[[ -s "$gtf_tmp" ]] || die "panel gene records were not found in the Ensembl GTF"
mv -f "$gtf_tmp" "$PANEL_GTF"
rm -f -- "${output_dir}/.panel_ids.tmp"
printf 'Wrote panel GTF subset: %s\n' "$PANEL_GTF"

# Feature table. The PRDX4 rows reproduce script_workflow_01's feature ids,
# coordinates and derivations exactly; the remaining genes contribute merged
# exon-union intervals plus the exonic length needed for length normalization.
printf 'Building the panel feature table...\n'
feature_tmp="${FEATURE_TABLE}.tmp"
"$PYTHON" - "$PANEL_INTERVALS" "$PANEL_GTF" "${output_dir}/.panel_gtf.tmp" "$feature_tmp" <<'PY'
import csv
import re
import sys
from pathlib import Path

intervals_path, gtf_path, full_gtf_path, output_path = sys.argv[1:]

PRDX4_TRANSCRIPT_ID = "ENSMUST00000026328"
PRDX4T_TRANSCRIPT_ID = "ENSMUST00000130349"


def attributes(text):
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]*)"', text))


with Path(intervals_path).open(newline="") as handle:
    panel = list(csv.DictReader(handle, delimiter="\t"))

records_by_gene = {row["gene_id"]: [] for row in panel}
with Path(gtf_path).open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = attributes(fields[8])
        gene_id = attrs.get("gene_id")
        if gene_id not in records_by_gene:
            continue
        records_by_gene[gene_id].append(
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


def opposite_strand_neighbours(full_gtf_path, panel):
    """Genes on the strand opposite a panel gene, overlapping its interval.

    An unstranded library cannot tell such a neighbour's transcripts from the
    panel gene's own, so the exonic bases they cover are unusable for a
    comparison that must hold across stranded and unstranded datasets alike.
    """
    wanted = {
        row["gene_id"]: (row["chrom"], int(row["gene_start"]), int(row["gene_end"]), row["strand"])
        for row in panel
    }
    masks = {gene_id: [] for gene_id in wanted}
    with Path(full_gtf_path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            match = re.search(r'gene_id "([^"]*)"', fields[8])
            if match is None:
                continue
            other_id = match.group(1)
            name_match = re.search(r'gene_name "([^"]*)"', fields[8])
            other_name = name_match.group(1) if name_match else other_id
            for gene_id, (gene_chrom, gene_start, gene_end, gene_strand) in wanted.items():
                if other_id == gene_id or chrom != gene_chrom or strand == gene_strand:
                    continue
                if start <= gene_end and end >= gene_start:
                    masks[gene_id].append((start, end, other_name))
    return masks


def subtract(intervals, masks):
    """Remove every masked base from a list of inclusive intervals."""
    result = []
    for start, end in intervals:
        pieces = [(start, end)]
        for mask_start, mask_end, _ in masks:
            remaining = []
            for piece_start, piece_end in pieces:
                if mask_end < piece_start or mask_start > piece_end:
                    remaining.append((piece_start, piece_end))
                    continue
                if mask_start > piece_start:
                    remaining.append((piece_start, mask_start - 1))
                if mask_end < piece_end:
                    remaining.append((mask_end + 1, piece_end))
            pieces = remaining
        result.extend(pieces)
    return sorted(piece for piece in result if piece[1] >= piece[0])


def merge_exons(records, gene_name):
    exons = sorted(
        {(record["start"], record["end"]) for record in records if record["type"] == "exon"}
    )
    if not exons:
        raise SystemExit(f"ERROR: {gene_name} has no annotated exons")
    merged = [list(exons[0])]
    for start, end in exons[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def prdx4_features(records, gene_name):
    """Reproduce script_workflow_01_prepare_PRDX4_locus.sh feature derivation."""

    def exon(transcript_id, exon_number):
        matches = [
            record
            for record in records
            if record["type"] == "exon"
            and record["transcript_id"] == transcript_id
            and record["exon_number"] == str(exon_number)
        ]
        unique = {
            (record["chrom"], record["start"], record["end"], record["strand"])
            for record in matches
        }
        if len(unique) != 1:
            raise SystemExit(
                f"ERROR: expected one coordinate for exon {exon_number} of "
                f"{transcript_id}, found {len(unique)}"
            )
        chrom, start, end, strand = unique.pop()
        return start, end

    gene_records = [record for record in records if record["type"] == "gene"]
    if len(gene_records) != 1:
        raise SystemExit(f"ERROR: expected one {gene_name} gene record")
    gene = gene_records[0]

    prdx4_exon1 = exon(PRDX4_TRANSCRIPT_ID, 1)
    prdx4t_exon1 = exon(PRDX4T_TRANSCRIPT_ID, 1)
    canonical_body_exons = {number: exon(PRDX4_TRANSCRIPT_ID, number) for number in range(2, 8)}

    # Ensembl release 110 marks Prdx4-203 as mRNA_end_NF/cds_end_NF and only
    # annotates its alternative exon 1 plus downstream exons 2-5. Verify every
    # downstream exon that Ensembl does annotate against the canonical
    # transcript. The experimentally established PRDX4t model is exon 1A
    # followed by canonical exons 2-7 (Yim et al. 2011,
    # doi:10.1074/jbc.M111.257220).
    prdx4t_body_numbers = sorted(
        {
            int(record["exon_number"])
            for record in records
            if record["type"] == "exon"
            and record["transcript_id"] == PRDX4T_TRANSCRIPT_ID
            and record["exon_number"].isdigit()
            and int(record["exon_number"]) >= 2
        }
    )
    if not prdx4t_body_numbers or prdx4t_body_numbers[0] != 2:
        raise SystemExit("ERROR: PRDX4t lacks an annotated exon 2")
    for number in prdx4t_body_numbers:
        if number not in canonical_body_exons:
            raise SystemExit(f"ERROR: unexpected PRDX4t downstream exon number: {number}")
        if exon(PRDX4T_TRANSCRIPT_ID, number) != canonical_body_exons[number]:
            raise SystemExit(
                f"ERROR: PRDX4 and PRDX4t do not share the same annotated exon {number}"
            )

    ordered_exon1 = sorted((prdx4_exon1, prdx4t_exon1))
    gap_start = ordered_exon1[0][1] + 1
    gap_end = ordered_exon1[1][0] - 1
    if gap_start > gap_end:
        raise SystemExit("ERROR: PRDX4 and PRDX4t exon 1 intervals overlap")

    rows = [
        ("gene", gene["start"], gene["end"], ".", ".", "GTF gene feature"),
        ("prdx4t_exon1", *prdx4t_exon1, PRDX4T_TRANSCRIPT_ID, "1", "GTF exon feature"),
        ("prdx4_exon1", *prdx4_exon1, PRDX4_TRANSCRIPT_ID, "1", "GTF exon feature"),
        ("gap_between_tss", gap_start, gap_end, ".", ".", "derived gap between exon 1 intervals"),
    ]
    # Exons 2-7 are emitted individually rather than as one span. The span from
    # the start of exon 7 to the end of exon 2 is mostly intron, so a span would
    # count intronic and unspliced signal as shared-body expression.
    for number, (start, end) in canonical_body_exons.items():
        if number in prdx4t_body_numbers:
            transcript_id = f"{PRDX4_TRANSCRIPT_ID},{PRDX4T_TRANSCRIPT_ID}"
            derivation = "matching GTF exon features"
        else:
            transcript_id = PRDX4_TRANSCRIPT_ID
            derivation = (
                "canonical GTF exon; shared with Prdx4t per Yim et al. 2011 "
                "(doi:10.1074/jbc.M111.257220; Prdx4-203 is mRNA_end_NF)"
            )
        rows.append((f"shared_exon{number}", start, end, transcript_id, str(number), derivation))
    return rows


neighbour_masks = opposite_strand_neighbours(full_gtf_path, panel)

output_rows = []
for entry in panel:
    gene_id = entry["gene_id"]
    gene_name = entry["gene_name"]
    records = records_by_gene[gene_id]
    if not records:
        raise SystemExit(f"ERROR: no GTF records for {gene_name}")

    chroms = {record["chrom"] for record in records}
    strands = {record["strand"] for record in records}
    if chroms != {entry["chrom"]} or strands != {entry["strand"]}:
        raise SystemExit(f"ERROR: inconsistent chromosome or strand for {gene_name}")

    exon_union = merge_exons(records, gene_name)
    # Drop exonic bases an unstranded library could not attribute to this gene.
    masks = neighbour_masks[gene_id]
    masked_by = ""
    if masks:
        before = sum(end - start + 1 for start, end in exon_union)
        exon_union = subtract(exon_union, masks)
        if not exon_union:
            raise SystemExit(
                f"ERROR: every exon of {gene_name} is overlapped by an opposite-strand gene"
            )
        masked_by = ",".join(sorted({name for _, _, name in masks}))
        after = sum(end - start + 1 for start, end in exon_union)
        print(
            f"  {gene_name}: masked {before - after} of {before} exonic bp overlapped "
            f"by opposite-strand {masked_by}",
            file=sys.stderr,
        )
    exonic_length = sum(end - start + 1 for start, end in exon_union)

    prefix = gene_name.lower()
    if entry["quantification"] == "prdx4_isoform":
        feature_rows = prdx4_features(records, gene_name)
    else:
        feature_rows = [(f"{prefix}_gene", int(entry["gene_start"]), int(entry["gene_end"]),
                         ".", ".", "GTF gene feature")]
    # Every gene, PRDX4 included, also gets exon-union intervals. Gene-level
    # PRDX4 is what puts it on the same measurement basis as PCMT1 and CRBN;
    # the isoform features above stay the basis for PRDX4t-versus-PRDX4.
    for index, (start, end) in enumerate(exon_union, start=1):
        feature_rows.append(
            (
                f"{prefix}_exonunion{index}",
                start,
                end,
                ".",
                ".",
                "union of all annotated exons"
                + (", minus opposite-strand overlap" if masked_by else
                   " of every transcript"),
            )
        )

    locus_start = int(entry["locus_start"])
    locus_end = int(entry["locus_end"])
    for feature_id, start, end, transcript_id, exon_number, derivation in feature_rows:
        if start < locus_start or end > locus_end or end < start:
            raise SystemExit(
                f"ERROR: {feature_id} ({start}-{end}) falls outside the "
                f"{gene_name} prefilter locus {locus_start}-{locus_end}"
            )
        output_rows.append(
            {
                "feature_id": feature_id,
                "gene_name": gene_name,
                "role": entry["role"],
                "quantification": entry["quantification"],
                "chrom": entry["chrom"],
                "start": start,
                "end": end,
                "strand": entry["strand"],
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "exon_number": exon_number,
                "gene_exonic_length_bp": exonic_length,
                "exonic_bp_masked_by": masked_by,
                "derivation": derivation,
            }
        )

seen = set()
for row in output_rows:
    if row["feature_id"] in seen:
        raise SystemExit(f"ERROR: duplicate feature id: {row['feature_id']}")
    seen.add(row["feature_id"])

# The PRDX4t quantification in workflow step 04 requires these exact ids.
required = {"gene", "prdx4t_exon1", "prdx4_exon1", "gap_between_tss"} | {
    f"shared_exon{number}" for number in range(2, 8)
}
missing = sorted(required.difference(seen))
if missing:
    raise SystemExit(f"ERROR: feature table is missing PRDX4 features: {', '.join(missing)}")

fieldnames = [
    "feature_id", "gene_name", "role", "quantification", "chrom", "start", "end",
    "strand", "gene_id", "transcript_id", "exon_number", "gene_exonic_length_bp",
    "exonic_bp_masked_by", "derivation",
]
with Path(output_path).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
PY
mv -f "$feature_tmp" "$FEATURE_TABLE"
rm -f -- "${output_dir}/.panel_gtf.tmp"
printf 'Wrote panel feature table: %s\n' "$FEATURE_TABLE"

# The PRDX4 rows must be byte-identical to the single-locus feature table when
# that table is present, so existing PRDX4t results stay reproducible.
readonly LEGACY_FEATURE_TABLE="${output_dir}/PRDX4_locus_features.tsv"
if [[ -s "$LEGACY_FEATURE_TABLE" ]]; then
    "$PYTHON" - "$LEGACY_FEATURE_TABLE" "$FEATURE_TABLE" <<'PY'
import csv
import sys
from pathlib import Path

legacy_path, panel_path = sys.argv[1:]
compared = ("chrom", "start", "end", "strand", "gene_id", "transcript_id",
            "exon_number", "derivation")

def load(path):
    with Path(path).open(newline="") as handle:
        return {row["feature_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

legacy = load(legacy_path)
panel = load(panel_path)
for feature_id, row in legacy.items():
    if feature_id not in panel:
        raise SystemExit(f"ERROR: panel feature table lost legacy feature {feature_id}")
    for column in compared:
        if row[column] != panel[feature_id][column]:
            raise SystemExit(
                f"ERROR: {feature_id} {column} changed: "
                f"legacy={row[column]!r} panel={panel[feature_id][column]!r}"
            )
print(f"PRDX4 features match the single-locus table: {len(legacy)} features")
PY
fi

column -t -s $'\t' "$FEATURE_TABLE" | cut -c1-120

bt2_complete=1
for suffix in 1 2 3 4; do
    [[ -s "${BT2_PREFIX}.${suffix}.bt2" || -s "${BT2_PREFIX}.${suffix}.bt2l" ]] || bt2_complete=0
done
for suffix in 1 2; do
    [[ -s "${BT2_PREFIX}.rev.${suffix}.bt2" || -s "${BT2_PREFIX}.rev.${suffix}.bt2l" ]] || bt2_complete=0
done
current_fasta_sha="$(sha256_file "$PANEL_FASTA")"
if [[ ! -s "$BT2_SOURCE_SHA" || "$(tr -d '[:space:]' < "$BT2_SOURCE_SHA")" != "$current_fasta_sha" ]]; then
    bt2_complete=0
fi

if ((build_index)); then
    if ((bt2_complete)); then
        printf '\nExisting Bowtie2 panel index is complete; skipping rebuild.\n'
    else
        mkdir -p "$BT2_DIR"
        bt2_tmp="$(mktemp -d "${output_dir}/.bt2_panel.XXXXXX")"
        "$BOWTIE2_BUILD" --threads "$threads" "$PANEL_FASTA" "${bt2_tmp}/panel"
        printf '%s\n' "$current_fasta_sha" > "${bt2_tmp}/source_fasta.sha256"
        find "$bt2_tmp" -type f -name 'panel*.bt2*' -exec mv -f {} "$BT2_DIR"/ \;
        mv -f "${bt2_tmp}/source_fasta.sha256" "$BT2_SOURCE_SHA"
        rmdir "$bt2_tmp"
        printf 'Built Bowtie2 panel index: %s\n' "$BT2_PREFIX"
    fi
fi

manifest_tmp="${MANIFEST}.tmp"
{
    printf 'field\tvalue\n'
    printf 'species\tMus musculus\n'
    printf 'assembly\t%s\n' "$ASSEMBLY"
    printf 'ensembl_release\t%s\n' "$ENSEMBL_RELEASE"
    printf 'panel_definition\t%s\n' "$panel_file"
    printf 'panel_definition.sha256\t%s\n' "$(sha256_file "$panel_file")"
    printf 'panel_gene_count\t%s\n' "$(awk 'END {print NR-1}' "$PANEL_INTERVALS")"
    while IFS=$'\t' read -r gene_name gene_id role quantification contig chrom strand \
        gene_start gene_end locus_start locus_end locus_length rationale; do
        [[ "$gene_name" == "gene_name" ]] && continue
        printf 'gene.%s\t%s %s:%s-%s(%s) prefilter %s:%s-%s role=%s quantification=%s\n' \
            "$gene_name" "$gene_id" "$chrom" "$gene_start" "$gene_end" "$strand" \
            "$chrom" "$locus_start" "$locus_end" "$role" "$quantification"
    done < "$PANEL_INTERVALS"
    printf 'prdx4_transcript_id\t%s\n' "ENSMUST00000026328"
    printf 'prdx4t_transcript_id\t%s\n' "ENSMUST00000130349"
    printf 'prdx4t_biological_model\talternative exon 1A followed by shared exons 2-7\n'
    printf 'prdx4t_model_source\tYim et al. 2011; doi:10.1074/jbc.M111.257220\n'
    printf 'prefilter_flank_bp\t%s\n' "$prefilter_flank"
    printf 'genome_fasta_url\t%s\n' "$GENOME_URL"
    printf 'annotation_gtf_url\t%s\n' "$GTF_URL"
    printf 'gene_panel.fa.sha256\t%s\n' "$(sha256_file "$PANEL_FASTA")"
    printf 'gene_panel.fa.fai.sha256\t%s\n' "$(sha256_file "$PANEL_FAI")"
    printf 'gene_panel.gtf.sha256\t%s\n' "$(sha256_file "$PANEL_GTF")"
    printf 'gene_panel_intervals.tsv.sha256\t%s\n' "$(sha256_file "$PANEL_INTERVALS")"
    printf 'gene_panel_features.tsv.sha256\t%s\n' "$(sha256_file "$FEATURE_TABLE")"
} > "$manifest_tmp"
mv -f "$manifest_tmp" "$MANIFEST"

printf '\nGene panel preparation complete.\n'
printf 'Reference directory: %s\n' "$output_dir"
printf 'Manifest:            %s\n' "$MANIFEST"
