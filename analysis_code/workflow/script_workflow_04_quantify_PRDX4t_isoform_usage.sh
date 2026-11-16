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
feature_file="${PRDX4_FEATURE_TABLE:-}"
min_mapq="${MIN_MAPQ:-255}"
junction_slack="${JUNCTION_SLACK:-0}"
scan_flank="${SCAN_FLANK:-3000}"
strandedness="${STRANDEDNESS:-auto}"
dry_run=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Write run-level PRDX4t and canonical PRDX4 counts, then group technical runs by
SampleID to write biological-sample counts and per-million-spot metrics.
Biological replicates are never combined.

Options:
  --dataset ID        Dataset ID or dataset.conf path (default: PRJNA548952)
  --dataset-root DIR   Dataset root (default: ${DEFAULT_DATASET_ROOT})
  --reference-dir DIR Shared reference root (default: ${DEFAULT_REFERENCE_DIR})
  --output-dir DIR     Result directory
                       (default: DATASET_ROOT/analysis/quantification)
  --features FILE      GTF-derived feature table from workflow step 01
                       (default: REFERENCE_DIR/PRDX4_locus_features.tsv)
  --min-mapq N         Minimum alignment MAPQ (default: 255)
  --junction-slack N   Allowed splice-boundary offset in bp (default: 0)
  --scan-flank N       Bases scanned on each side of the GTF gene interval
                       (default: 3000)
  --strandedness MODE  auto, unstranded, forward, or reverse (default: auto)
                       forward  = mate 1 (or single read) lies on the transcribed strand
                                  (fr-secondstrand)
                       reverse  = mate 1 lies opposite it (fr-firststrand/dUTP)
                       auto     = call the protocol from the shared body exons
                                  of every library and require agreement
  --dry-run            Validate every input without writing results
  -h, --help           Show this help

Environment:
  DATASET_ROOT         Alternative dataset root
  REFERENCE_DIR         Alternative shared reference root
  PRDX4T_ENV_BIN       Directory containing the prdx4t-analysis Python
  PRDX4_FEATURE_TABLE  Alternative GTF-derived feature table
  MIN_MAPQ             Alternative minimum alignment MAPQ
  JUNCTION_SLACK       Alternative splice-boundary offset
  SCAN_FLANK           Alternative gene scan flank size
  STRANDEDNESS         Alternative library strandedness mode

Outputs:
  ${DATASET_ID}_PRDX4t_raw_counts_per_library.tsv
  ${DATASET_ID}_PRDX4t_normalized_metrics.tsv
  ${DATASET_ID}_PRDX4t_raw_counts_per_biological_sample.tsv
  ${DATASET_ID}_PRDX4t_normalized_metrics_per_biological_sample.tsv
  ${DATASET_ID}_PRDX4t_quantification_manifest.tsv
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
        --output-dir)
            (($# >= 2)) || die "--output-dir requires a directory"
            output_dir="$2"
            shift 2
            ;;
        --features)
            (($# >= 2)) || die "--features requires a file"
            feature_file="$2"
            shift 2
            ;;
        --min-mapq)
            (($# >= 2)) || die "--min-mapq requires an integer"
            min_mapq="$2"
            shift 2
            ;;
        --junction-slack)
            (($# >= 2)) || die "--junction-slack requires an integer"
            junction_slack="$2"
            shift 2
            ;;
        --scan-flank)
            (($# >= 2)) || die "--scan-flank requires an integer"
            scan_flank="$2"
            shift 2
            ;;
        --strandedness)
            (($# >= 2)) || die "--strandedness requires a mode"
            strandedness="$2"
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

[[ "$min_mapq" =~ ^[0-9]+$ ]] || die "--min-mapq must be an integer from 0 to 255"
((min_mapq <= 255)) || die "--min-mapq must be an integer from 0 to 255"
[[ "$junction_slack" =~ ^[0-9]+$ ]] || die "--junction-slack must be a non-negative integer"
[[ "$scan_flank" =~ ^[0-9]+$ ]] || die "--scan-flank must be a non-negative integer"
case "$strandedness" in
    auto|unstranded|forward|reverse) ;;
    *) die "--strandedness must be auto, unstranded, forward, or reverse" ;;
esac
[[ -r "$RUN_MANIFEST" ]] || die "run manifest not found: $RUN_MANIFEST"
[[ -n "$output_dir" ]] || output_dir="${dataset_root}/analysis/quantification"
[[ -n "$feature_file" ]] || feature_file="${reference_dir}/PRDX4_locus_features.tsv"
[[ -r "$feature_file" ]] || die "GTF-derived feature table is missing; rerun script_workflow_01_prepare_PRDX4_locus.sh: $feature_file"

env_bin="${PRDX4T_ENV_BIN:-}"
if [[ -n "$env_bin" ]]; then
    :
elif [[ -x "${DEFAULT_ENV_BIN}/python" ]]; then
    env_bin="$DEFAULT_ENV_BIN"
elif command -v python >/dev/null 2>&1 && python -c 'import pysam' >/dev/null 2>&1; then
    env_bin="$(dirname -- "$(command -v python)")"
else
    die "prdx4t-analysis Python was not found; run script_environment_create_prdx4t_analysis.sh first"
fi

readonly PYTHON="${env_bin}/python"
[[ -x "$PYTHON" ]] || die "Python is not executable: $PYTHON"

args=(
    --dataset-id "$DATASET_ID"
    --main-stage-groups "$MAIN_STAGE_GROUPS"
    --manifest "$RUN_MANIFEST"
    --bam-dir "${dataset_root}/analysis/bam"
    --features "$feature_file"
    --output-dir "$output_dir"
    --min-mapq "$min_mapq"
    --junction-slack "$junction_slack"
    --scan-flank "$scan_flank"
    --strandedness "$strandedness"
)
if ((dry_run)); then
    args+=(--dry-run)
fi

"$PYTHON" - "${args[@]}" <<'PY'
import argparse
import csv
import importlib.metadata
import os
from pathlib import Path
import platform

import pysam


CHROM = None
GENE_STRAND = None
GENE_START = None
GENE_END = None
SCAN_START = None
SCAN_END = None
ALT_EXON1 = None
CANONICAL_EXON1 = None
SHARED_EXON2 = None
SHARED_BODY = None
REGIONS = None
BAM_SUFFIX = None

# Exons 2-7 are shared by PRDX4 and PRDX4t. They are used as a union of exon
# intervals, never as one span: the span between them is 94% intron, so a span
# would fold unspliced and intronic signal into the shared-body tally and into
# the strandedness call.
SHARED_BODY_FEATURES = (
    "shared_exon2",
    "shared_exon3",
    "shared_exon4",
    "shared_exon5",
    "shared_exon6",
    "shared_exon7",
)

# Strandedness is called against the unstranded null of 0.5, not against an
# expectation of perfection: 0.75 is a 3:1 strand enrichment, which no unstranded
# protocol produces. Libraries between the stranded and unstranded bands are
# treated as undiagnosable rather than quietly assumed.
STRAND_FORWARD_MIN = 0.75
STRAND_REVERSE_MAX = 0.25
STRAND_UNSTRANDED_BAND = (0.40, 0.60)

# Minimum shared-exon fragments needed to call one library on its own. The unit
# is fragments, not alignment records: a pair contributes once, so this is about
# half the corresponding record count. At 50 fragments the binomial 95% interval
# around the 0.75 cutoff excludes 0.5, which is the discrimination this needs.
STRAND_MIN_EVIDENCE = 50

# Stranded libraries below this are still called, but their antisense background
# is high enough to matter for near-zero windows and is reported.
STRAND_LEAKY_WARN = 0.90

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--main-stage-groups", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-mapq", required=True, type=int)
    parser.add_argument("--junction-slack", required=True, type=int)
    parser.add_argument("--scan-flank", required=True, type=int)
    parser.add_argument(
        "--strandedness",
        required=True,
        choices=("auto", "unstranded", "forward", "reverse"),
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def load_feature_table(path):
    required_columns = {
        "feature_id",
        "chrom",
        "start",
        "end",
        "strand",
        "gene_id",
        "transcript_id",
        "exon_number",
        "derivation",
    }
    required_features = {
        "gene",
        "prdx4t_exon1",
        "prdx4_exon1",
        "gap_between_tss",
        *SHARED_BODY_FEATURES,
    }
    features = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = sorted(required_columns.difference(reader.fieldnames or []))
        if missing_columns:
            raise SystemExit(
                f"ERROR: feature table is missing columns: {', '.join(missing_columns)}"
            )
        for row in reader:
            feature_id = row["feature_id"].strip()
            if feature_id in features:
                raise SystemExit(f"ERROR: duplicate feature in {path}: {feature_id}")
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as exc:
                raise SystemExit(f"ERROR: non-integer coordinates for {feature_id} in {path}") from exc
            if start < 1 or end < start:
                raise SystemExit(f"ERROR: invalid coordinates for {feature_id} in {path}: {start}-{end}")
            row["start"] = start
            row["end"] = end
            features[feature_id] = row

    missing_features = sorted(required_features.difference(features))
    if missing_features:
        raise SystemExit(f"ERROR: feature table is missing features: {', '.join(missing_features)}")

    gene = features["gene"]
    for feature_id in required_features:
        feature = features[feature_id]
        if feature["chrom"] != gene["chrom"] or feature["strand"] != gene["strand"]:
            raise SystemExit(f"ERROR: chromosome/strand mismatch for {feature_id} in {path}")
        if feature["start"] < gene["start"] or feature["end"] > gene["end"]:
            raise SystemExit(f"ERROR: {feature_id} lies outside the GTF gene interval in {path}")
    return features


def use_alignment(read, min_mapq):
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.mapping_quality < min_mapq
    )


def overlaps_region(blocks, intervals):
    """True if any aligned block touches any interval of the region.

    Regions are lists of 1-based inclusive intervals so that a multi-exon
    feature is the union of its exons rather than the span across them.
    """
    for start, end in intervals:
        start0 = start - 1
        if any(block_end > start0 and block_start < end for block_start, block_end in blocks):
            return True
    return False


def supports_junction(blocks, exon_a, exon_b, slack):
    low_end = min(exon_a[1], exon_b[1])
    high_start = max(exon_a[0], exon_b[0])

    for (_, upstream_end), (downstream_start, _) in zip(blocks[:-1], blocks[1:]):
        downstream_start_1based = downstream_start + 1
        if (
            abs(upstream_end - low_end) <= slack
            and abs(downstream_start_1based - high_start) <= slack
        ):
            return True
    return False


def mate1_is_reverse(read):
    """Orientation of mate 1, inferred from whichever mate this record is."""
    if not read.is_paired or read.is_read1:
        return read.is_reverse
    if not read.mate_is_unmapped:
        return read.mate_is_reverse
    # A mapped read 2 with an unmapped mate has no mate-orientation flag that can
    # be trusted. Fall back to the expected opposite orientation of the library.
    return not read.is_reverse


def is_sense(read, strandedness):
    """True if the fragment came from the strand this gene is transcribed from.

    Both PRDX4 and PRDX4t lie on the same strand of the same gene, so this
    separates real transcription from antisense background; it does not
    discriminate between the two isoforms.
    """
    if strandedness == "unstranded":
        return True
    gene_on_minus = GENE_STRAND == "-"
    if strandedness == "forward":
        return mate1_is_reverse(read) == gene_on_minus
    return mate1_is_reverse(read) != gene_on_minus


def detect_strandedness(bam_path, min_mapq):
    """Call the protocol from the shared body exons.

    Exons 2-7 are common to PRDX4 and PRDX4t, so they carry the most alignments
    of any window and their strand composition is independent of which isoform
    a library expresses.
    """
    fetch_start = min(start for start, _ in SHARED_BODY)
    fetch_end = max(end for _, end in SHARED_BODY)
    gene_on_minus = GENE_STRAND == "-"
    fragment_orientations = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(CHROM, fetch_start - 1, fetch_end):
            if not use_alignment(read, min_mapq):
                continue
            blocks = read.get_blocks()
            if not blocks or not overlaps_region(blocks, SHARED_BODY):
                continue
            fragment_orientations.setdefault(read.query_name, mate1_is_reverse(read))

    total = len(fragment_orientations)
    if total < STRAND_MIN_EVIDENCE:
        return None, total, None

    mate1_on_gene_strand = sum(
        orientation == gene_on_minus for orientation in fragment_orientations.values()
    )
    fraction = mate1_on_gene_strand / total
    if fraction >= STRAND_FORWARD_MIN:
        return "forward", total, fraction
    if fraction <= STRAND_REVERSE_MAX:
        return "reverse", total, fraction
    if STRAND_UNSTRANDED_BAND[0] <= fraction <= STRAND_UNSTRANDED_BAND[1]:
        return "unstranded", total, fraction
    return "ambiguous", total, fraction


def resolve_strandedness(records, min_mapq, requested):
    """Detect the protocol for every library and require them to agree.

    Detection runs even when the mode is given explicitly, so a mismatch between
    what was asked for and what the data show is reported rather than applied
    silently.
    """
    observations = []
    for record in records:
        call, total, fraction = detect_strandedness(Path(record["bam_path"]), min_mapq)
        record["strandedness_detected"] = call or "undetermined"
        record["strand_evidence_fragments"] = str(total)
        observations.append((record["Run"].strip(), call, total, fraction))
        shown = "n/a" if fraction is None else f"{fraction:.3f}"
        print(
            f"  {record['Run'].strip()}: mate-1 sense fraction {shown} "
            f"over {total} shared-exon fragments -> {call or 'undetermined'}"
        )

    leaky = [
        (run, fraction)
        for run, call, _, fraction in observations
        if call in ("forward", "reverse") and fraction is not None
        and min(fraction, 1 - fraction) > 1 - STRAND_LEAKY_WARN
    ]
    for run, fraction in leaky:
        print(
            f"  WARNING: {run} is stranded but leaky "
            f"({100 * min(fraction, 1 - fraction):.1f}% antisense at the shared exons); "
            "its antisense columns are the background estimate for near-zero windows."
        )

    clear_calls = {"forward", "reverse", "unstranded"}
    called = {call for _, call, _, _ in observations if call in clear_calls}
    ambiguous = [run for run, call, _, _ in observations if call == "ambiguous"]
    undetermined = [run for run, call, _, _ in observations if call is None]

    if requested != "auto":
        if requested == "unstranded":
            detected = next(iter(called)) if len(called) == 1 else "mixed_or_uncertain"
            print("  NOTE: explicit unstranded counting requested; automatic calls are diagnostic only.")
            return requested, detected

        conflicts = [
            f"{run}={call}"
            for run, call, _, _ in observations
            if call in clear_calls and call != requested
        ]
        if conflicts:
            raise SystemExit(
                f"ERROR: --strandedness {requested} conflicts with clear automatic calls: "
                + ", ".join(conflicts)
            )
        uncertain = ambiguous + undetermined
        if uncertain:
            print(
                f"  NOTE: applying explicit '{requested}' to uncertain libraries: "
                + ", ".join(uncertain)
            )
        detected = next(iter(called)) if len(called) == 1 else "ambiguous_or_undetermined"
        return requested, detected

    if ambiguous:
        raise SystemExit(
            "ERROR: strandedness is ambiguous for: "
            f"{', '.join(ambiguous)}. Re-run with an explicit --strandedness."
        )
    if len(called) > 1:
        detail = ", ".join(f"{run}={call}" for run, call, _, _ in observations if call)
        raise SystemExit(f"ERROR: libraries disagree on strandedness: {detail}")
    if not called:
        raise SystemExit(
            "ERROR: no library had enough shared-exon coverage to call strandedness; "
            "re-run with an explicit --strandedness."
        )

    detected = called.pop()
    if undetermined:
        print(
            f"  NOTE: {', '.join(undetermined)} had under {STRAND_MIN_EVIDENCE} "
            f"shared-exon fragments; applying the consensus call '{detected}'."
        )

    return detected, detected


COUNTED_FEATURES = (
    "unique_primary_locus",
    "alt_exon1_to_exon2_junction",
    "canonical_exon1_to_exon2_junction",
)


def analyze_bam(bam_path, min_mapq, junction_slack, strandedness):
    """Count fragments rather than alignment records.

    Fragments are keyed on query name: a paired-end spot contributes once to a
    window even when both mates overlap it, and a single-end spot contributes
    once through its one read. Orientation is resolved consistently at fragment
    level, so a fragment falls entirely on the sense or antisense side.
    """
    groups = {}
    for name in list(REGIONS) + list(COUNTED_FEATURES):
        groups[f"{name}_fragments"] = set()
        groups[f"{name}_antisense_fragments"] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(CHROM, SCAN_START - 1, SCAN_END):
            if not use_alignment(read, min_mapq):
                continue
            blocks = read.get_blocks()
            if not blocks:
                continue

            # Every alignment is tallied; the strand call routes its fragment to
            # the sense set or the antisense set, so nothing is discarded.
            suffix = "_fragments" if is_sense(read, strandedness) else "_antisense_fragments"
            fragment = read.query_name

            groups[f"unique_primary_locus{suffix}"].add(fragment)
            for region, intervals in REGIONS.items():
                if overlaps_region(blocks, intervals):
                    groups[f"{region}{suffix}"].add(fragment)

            if supports_junction(blocks, ALT_EXON1, SHARED_EXON2, junction_slack):
                groups[f"alt_exon1_to_exon2_junction{suffix}"].add(fragment)
            if supports_junction(blocks, CANONICAL_EXON1, SHARED_EXON2, junction_slack):
                groups[f"canonical_exon1_to_exon2_junction{suffix}"].add(fragment)

    return {key: len(value) for key, value in groups.items()}


def resolve_bam_suffix(manifest_path, bam_dir):
    """Pick the BAM flavour this dataset has, and require it to be uniform.

    Workflow step 03 writes RUN.panel.bam in panel mode and RUN.prdx4.bam in
    prdx4 mode. A dataset half-processed in each mode would silently mix two
    reference definitions, so a partial set is an error rather than a fallback.
    """
    runs = []
    with Path(manifest_path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = row["Run"].strip()
            if run:
                runs.append(run)
    if not runs:
        raise SystemExit(f"ERROR: run manifest contains no runs: {manifest_path}")
    present = {
        suffix: sum(1 for run in runs if (Path(bam_dir) / f"{run}.{suffix}.bam").is_file())
        for suffix in ("panel", "prdx4")
    }
    complete = [suffix for suffix, count in present.items() if count == len(runs)]
    if len(complete) == 2:
        # Both exist; the panel BAM is the current output of workflow step 03.
        return "panel"
    if len(complete) == 1:
        return complete[0]
    raise SystemExit(
        "ERROR: no complete set of locus BAMs for this dataset in "
        f"{bam_dir}\n  panel BAMs: {present['panel']}/{len(runs)}"
        f"\n  prdx4 BAMs: {present['prdx4']}/{len(runs)}\n"
        "Run script_workflow_03_process_all_libraries.sh for this dataset first."
    )


def load_and_validate_inputs(manifest_path, bam_dir):
    required = {
        "Run", "SampleName", "LibraryName", "LibraryLayout", "spots",
        "SampleID", "Group", "GroupLabel", "GroupOrder", "Replicate", "TechnicalRun",
    }
    records = []
    seen_runs = set()

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise SystemExit(f"ERROR: run manifest is missing columns: {', '.join(missing)}")

        for record in reader:
            run = record["Run"].strip()
            if not run:
                raise SystemExit("ERROR: run manifest contains an empty Run value")
            if run in seen_runs:
                raise SystemExit(f"ERROR: duplicate run in manifest: {run}")
            seen_runs.add(run)

            if record["LibraryLayout"].strip() not in {"PAIRED", "SINGLE"}:
                raise SystemExit(
                    f"ERROR: {run} has unsupported LibraryLayout: {record['LibraryLayout']!r}"
                )
            for column in ("spots", "GroupOrder", "Replicate", "TechnicalRun"):
                value = record[column].strip()
                if not value.isdigit() or int(value) < 0:
                    raise SystemExit(f"ERROR: {run} has invalid {column}: {value!r}")
            if any(int(record[column]) <= 0 for column in ("spots", "Replicate", "TechnicalRun")):
                raise SystemExit(f"ERROR: {run} has non-positive spots, Replicate, or TechnicalRun")

            bam_path = bam_dir / f"{run}.{BAM_SUFFIX}.bam"
            bai_path = Path(f"{bam_path}.bai")
            if not bam_path.is_file():
                raise SystemExit(f"ERROR: BAM not found: {bam_path}")
            if not bai_path.is_file():
                raise SystemExit(f"ERROR: BAM index not found: {bai_path}")

            try:
                with pysam.AlignmentFile(bam_path, "rb") as bam:
                    if CHROM not in bam.references:
                        raise SystemExit(f"ERROR: chromosome {CHROM} is absent from {bam_path}")
                    bam.check_index()
            except (OSError, ValueError) as exc:
                raise SystemExit(f"ERROR: cannot use {bam_path}: {exc}") from exc

            record["bam_path"] = str(bam_path)
            records.append(record)

    if not records:
        raise SystemExit(f"ERROR: run manifest contains no libraries: {manifest_path}")
    return records


def atomic_write_tsv(path, rows, fieldnames):
    temporary_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


FRAGMENT_COUNT_FIELDNAMES = (
    "unique_primary_locus_fragments",
    "alt_exon1_prdx4t_fragments",
    "canonical_exon1_fragments",
    "gap_between_tss_fragments",
    "shared_body_exons2to7_fragments",
    "alt_exon1_to_exon2_junction_fragments",
    "canonical_exon1_to_exon2_junction_fragments",
    "unique_primary_locus_antisense_fragments",
    "alt_exon1_prdx4t_antisense_fragments",
    "canonical_exon1_antisense_fragments",
    "gap_between_tss_antisense_fragments",
    "shared_body_exons2to7_antisense_fragments",
    "alt_exon1_to_exon2_junction_antisense_fragments",
    "canonical_exon1_to_exon2_junction_antisense_fragments",
)


def derive_metrics(rows):
    """Normalize the four plotted raw counts by sequenced SRA-spot depth."""
    metrics = []
    for row in rows:
        spots = int(row["manifest_spots"])
        if spots <= 0:
            raise SystemExit(f"ERROR: {row['run']} has a non-positive manifest spots value: {spots}")
        million_spots = spots / 1e6
        prdx4t_exon1 = row["alt_exon1_prdx4t_fragments"]
        canonical_exon1 = row["canonical_exon1_fragments"]
        prdx4t_junction = row["alt_exon1_to_exon2_junction_fragments"]
        canonical_junction = row["canonical_exon1_to_exon2_junction_fragments"]
        def per_million(count):
            return f"{count / million_spots:.4f}"

        metrics.append(
            {
                "run": row["run"],
                "sample_id": row["sample_id"],
                "library_name": row["library_name"],
                "library_layout": row["library_layout"],
                "stage": row["stage"],
                "stage_label": row["stage_label"],
                "stage_order": row["stage_order"],
                "replicate": row["replicate"],
                "technical_run": row["technical_run"],
                "manifest_spots": row["manifest_spots"],
                "million_spots": f"{million_spots:.4f}",
                "prdx4t_exon1_per_million_spots": per_million(prdx4t_exon1),
                "canonical_exon1_per_million_spots": per_million(canonical_exon1),
                "prdx4t_junction_per_million_spots": per_million(prdx4t_junction),
                "canonical_junction_per_million_spots": per_million(canonical_junction),
            }
        )
    return metrics


def aggregate_biological_samples(rows):
    """Combine technical SRA runs belonging to the same biological sample."""
    samples = {}
    invariant_fields = (
        "sample_name",
        "sample_id",
        "library_layout",
        "stage",
        "stage_label",
        "stage_order",
        "replicate",
        "include_in_main_plot",
    )

    for row in rows:
        sample_id = row["sample_id"]
        if sample_id not in samples:
            sample = {field: row[field] for field in invariant_fields}
            sample.update(
                {
                    "library_name": sample_id,
                    "technical_run_count": 0,
                    "source_runs": [],
                    "manifest_spots": 0,
                    "strandedness_detected": set(),
                    "strand_evidence_fragments": 0,
                    "bam_paths": [],
                }
            )
            sample.update({field: 0 for field in FRAGMENT_COUNT_FIELDNAMES})
            samples[sample_id] = sample

        sample = samples[sample_id]
        for field in invariant_fields:
            if sample[field] != row[field]:
                raise SystemExit(
                    f"ERROR: technical runs for {sample_id} disagree on {field}: "
                    f"{sample[field]!r} versus {row[field]!r}"
                )
        sample["technical_run_count"] += 1
        sample["source_runs"].append(row["run"])
        sample["manifest_spots"] += int(row["manifest_spots"])
        sample["strandedness_detected"].add(row["strandedness_detected"])
        sample["strand_evidence_fragments"] += int(row["strand_evidence_fragments"])
        sample["bam_paths"].append(row["bam_path"])
        for field in FRAGMENT_COUNT_FIELDNAMES:
            sample[field] += int(row[field])

    aggregated = []
    for sample in samples.values():
        sample["source_runs"] = ",".join(sample["source_runs"])
        sample["strandedness_detected"] = ",".join(sorted(sample["strandedness_detected"]))
        sample["bam_paths"] = ",".join(sample["bam_paths"])
        aggregated.append(sample)
    aggregated.sort(key=lambda row: (row["stage_order"], row["replicate"], row["sample_id"]))
    return aggregated


def derive_sample_metrics(rows):
    metrics = []
    for row in rows:
        spots = int(row["manifest_spots"])
        million_spots = spots / 1e6

        def per_million(field):
            return f"{int(row[field]) / million_spots:.4f}"

        metrics.append(
            {
                "sample_id": row["sample_id"],
                "library_name": row["library_name"],
                "library_layout": row["library_layout"],
                "stage": row["stage"],
                "stage_label": row["stage_label"],
                "stage_order": row["stage_order"],
                "replicate": row["replicate"],
                "technical_run_count": row["technical_run_count"],
                "source_runs": row["source_runs"],
                "include_in_main_plot": row["include_in_main_plot"],
                "manifest_spots": spots,
                "million_spots": f"{million_spots:.4f}",
                "prdx4t_exon1_per_million_spots": per_million("alt_exon1_prdx4t_fragments"),
                "canonical_exon1_per_million_spots": per_million("canonical_exon1_fragments"),
                "prdx4t_junction_per_million_spots": per_million("alt_exon1_to_exon2_junction_fragments"),
                "canonical_junction_per_million_spots": per_million("canonical_exon1_to_exon2_junction_fragments"),
            }
        )
    return metrics


METRIC_FIELDNAMES = [
    "run",
    "sample_id",
    "library_name",
    "library_layout",
    "stage",
    "stage_label",
    "stage_order",
    "replicate",
    "technical_run",
    "manifest_spots",
    "million_spots",
    "prdx4t_exon1_per_million_spots",
    "canonical_exon1_per_million_spots",
    "prdx4t_junction_per_million_spots",
    "canonical_junction_per_million_spots",
]

SAMPLE_METRIC_FIELDNAMES = [
    "sample_id",
    "library_name",
    "library_layout",
    "stage",
    "stage_label",
    "stage_order",
    "replicate",
    "technical_run_count",
    "source_runs",
    "include_in_main_plot",
    "manifest_spots",
    "million_spots",
    "prdx4t_exon1_per_million_spots",
    "canonical_exon1_per_million_spots",
    "prdx4t_junction_per_million_spots",
    "canonical_junction_per_million_spots",
]

SAMPLE_COUNT_FIELDNAMES = [
    "sample_id",
    "sample_name",
    "library_name",
    "library_layout",
    "stage",
    "stage_label",
    "stage_order",
    "replicate",
    "technical_run_count",
    "source_runs",
    "include_in_main_plot",
    "manifest_spots",
    "strandedness_detected",
    "strand_evidence_fragments",
    *FRAGMENT_COUNT_FIELDNAMES,
    "bam_paths",
]


def write_settings(
    path, args, manifest_path, bam_dir, feature_path, library_count, biological_sample_count,
    main_sample_count, strandedness, detected
):
    settings = [
        ("dataset_id", args.dataset_id),
        ("reference", "GRCm39"),
        ("annotation_coordinates", "Ensembl release 110"),
        ("feature_definition", str(feature_path)),
        ("chromosome", CHROM),
        ("gene_strand", GENE_STRAND),
        ("gene_interval_1based_inclusive", f"{CHROM}:{GENE_START}-{GENE_END}"),
        ("scan_flank_bp", str(args.scan_flank)),
        ("scan_interval_1based_inclusive", f"{CHROM}:{SCAN_START}-{SCAN_END}"),
        ("alt_exon1_prdx4t_interval", f"{CHROM}:{ALT_EXON1[0]}-{ALT_EXON1[1]}"),
        ("canonical_exon1_interval", f"{CHROM}:{CANONICAL_EXON1[0]}-{CANONICAL_EXON1[1]}"),
        ("shared_exon2_interval", f"{CHROM}:{SHARED_EXON2[0]}-{SHARED_EXON2[1]}"),
        ("minimum_mapq", str(args.min_mapq)),
        ("junction_boundary_slack_bp", str(args.junction_slack)),
        ("excluded_alignment_types", "unmapped,secondary,supplementary"),
        ("strandedness_requested", args.strandedness),
        ("strandedness_detected", detected),
        ("strandedness_applied", strandedness),
        (
            "strand_assignment_rule",
            "forward: mate 1 (or single read) on the transcribed strand; reverse: opposite it",
        ),
        (
            "shared_body_intervals",
            "union of " + ",".join(f"{CHROM}:{start}-{end}" for start, end in SHARED_BODY),
        ),
        (
            "strand_detection_window",
            "shared exons 2-7, exon union only",
        ),
        (
            "antisense_handling",
            "reported in separate *_antisense_fragments columns; never discarded",
        ),
        ("count_unit", "fragment (one pair for PAIRED; one read for SINGLE; counted once per window)"),
        ("runs_quantified", str(library_count)),
        ("biological_samples_quantified", str(biological_sample_count)),
        ("main_stage_groups", args.main_stage_groups),
        ("main_stage_samples_selected_for_plotting", str(main_sample_count)),
        (
            "technical_run_aggregation",
            "raw fragment counts and manifest spots summed by SampleID before sample-level normalization",
        ),
        ("replicate_aggregation", "none"),
        ("counts_table_normalization", "none (raw counts)"),
        (
            "derived_metrics_table",
            "four plotted counts per million sequenced SRA spots",
        ),
        (
            "depth_denominator",
            "manifest spots (read pairs for PAIRED; reads for SINGLE); mapped-read totals are "
            "unavailable because only the PRDX4 locus is retained",
        ),
        (
            "depth_denominator_caveat",
            "assumes comparable filtering and mapping rates across runs; consult per-run fastp and STAR logs",
        ),
        ("plots_created", "no"),
        ("run_manifest", str(manifest_path)),
        ("bam_directory", str(bam_dir)),
        ("python_version", platform.python_version()),
        ("pysam_version", importlib.metadata.version("pysam")),
    ]
    rows = [{"parameter": key, "value": value} for key, value in settings]
    atomic_write_tsv(path, rows, ["parameter", "value"])


def main():
    global CHROM, GENE_STRAND, GENE_START, GENE_END, SCAN_START, SCAN_END
    global ALT_EXON1, CANONICAL_EXON1, SHARED_EXON2, SHARED_BODY, REGIONS

    args = parse_args()
    main_stage_groups = {value.strip() for value in args.main_stage_groups.split(",") if value.strip()}
    if not main_stage_groups:
        raise SystemExit("ERROR: --main-stage-groups contains no stage codes")
    manifest_path = Path(args.manifest).resolve()
    bam_dir = Path(args.bam_dir).resolve()
    feature_path = Path(args.features).resolve()
    output_dir = Path(args.output_dir).resolve()
    features = load_feature_table(feature_path)

    gene = features["gene"]
    CHROM = gene["chrom"]
    GENE_STRAND = gene["strand"]
    GENE_START = gene["start"]
    GENE_END = gene["end"]
    SCAN_START = max(1, GENE_START - args.scan_flank)
    SCAN_END = GENE_END + args.scan_flank
    ALT_EXON1 = (features["prdx4t_exon1"]["start"], features["prdx4t_exon1"]["end"])
    CANONICAL_EXON1 = (features["prdx4_exon1"]["start"], features["prdx4_exon1"]["end"])
    SHARED_EXON2 = (features["shared_exon2"]["start"], features["shared_exon2"]["end"])
    SHARED_BODY = [
        (features[name]["start"], features[name]["end"]) for name in SHARED_BODY_FEATURES
    ]
    REGIONS = {
        "alt_exon1_prdx4t": [ALT_EXON1],
        "canonical_exon1": [CANONICAL_EXON1],
        "gap_between_tss": [
            (features["gap_between_tss"]["start"], features["gap_between_tss"]["end"])
        ],
        "shared_body_exons2to7": SHARED_BODY,
    }

    global BAM_SUFFIX
    BAM_SUFFIX = resolve_bam_suffix(manifest_path, bam_dir)
    all_records = load_and_validate_inputs(manifest_path, bam_dir)
    print(f"Locus BAMs: RUN.{BAM_SUFFIX}.bam")
    records = [
        record for record in all_records if record["Group"].strip() in main_stage_groups
    ]
    excluded_records = [
        record for record in all_records if record["Group"].strip() not in main_stage_groups
    ]
    if not records:
        raise SystemExit("ERROR: no runs belong to the configured main developmental stages")

    print(f"Validated {len(all_records)} indexed BAMs listed in {manifest_path}")
    if excluded_records:
        excluded = ", ".join(record["Run"].strip() for record in excluded_records)
        print(
            f"Ignoring {len(excluded_records)} non-main-stage runs before strandedness detection: "
            f"{excluded}"
        )
    print(f"Quantifying {len(records)} runs from main stages: {', '.join(sorted(main_stage_groups))}")
    print(f"GTF-derived features: {feature_path}")
    print(f"Minimum MAPQ: {args.min_mapq}; splice-boundary slack: {args.junction_slack} bp")

    print(f"Calling library strandedness (requested: {args.strandedness}):")
    strandedness, detected = resolve_strandedness(records, args.min_mapq, args.strandedness)
    print(f"Strandedness in use: {strandedness} (detected: {detected})")
    if strandedness == "unstranded":
        print("Counting both strands together; antisense columns will read zero.")
    else:
        print("Sense and antisense alignments are reported in separate columns.")

    if args.dry_run:
        print("Dry run complete; no result files were written.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for record in records:
        run = record["Run"].strip()
        library_name = record["LibraryName"].strip()
        counts = analyze_bam(
            Path(record["bam_path"]), args.min_mapq, args.junction_slack, strandedness
        )
        row = {
            "run": run,
            "sample_name": record["SampleName"].strip(),
            "sample_id": record["SampleID"].strip(),
            "library_name": library_name,
            "library_layout": record["LibraryLayout"].strip(),
            "stage": record["Group"].strip(),
            "stage_label": record["GroupLabel"].strip(),
            "stage_order": int(record["GroupOrder"]),
            "replicate": int(record["Replicate"]),
            "technical_run": int(record["TechnicalRun"]),
            "include_in_main_plot": int(record["Group"].strip() in main_stage_groups),
            "manifest_spots": record["spots"].strip(),
            "strandedness_detected": record["strandedness_detected"],
            "strand_evidence_fragments": record["strand_evidence_fragments"],
            **counts,
            "bam_path": record["bam_path"],
        }
        rows.append(row)
        print(
            f"{run}: PRDX4t exon 1={counts['alt_exon1_prdx4t_fragments']} "
            f"(antisense {counts['alt_exon1_prdx4t_antisense_fragments']}), "
            f"canonical exon 1={counts['canonical_exon1_fragments']} "
            f"(antisense {counts['canonical_exon1_antisense_fragments']}), "
            f"PRDX4t junction={counts['alt_exon1_to_exon2_junction_fragments']}, "
            f"canonical junction={counts['canonical_exon1_to_exon2_junction_fragments']}"
        )

    rows.sort(key=lambda row: (row["stage_order"], row["replicate"], row["technical_run"], row["run"]))
    fieldnames = [
        "run",
        "sample_name",
        "sample_id",
        "library_name",
        "library_layout",
        "stage",
        "stage_label",
        "stage_order",
        "replicate",
        "technical_run",
        "include_in_main_plot",
        "manifest_spots",
        "strandedness_detected",
        "strand_evidence_fragments",
        "unique_primary_locus_fragments",
        "alt_exon1_prdx4t_fragments",
        "canonical_exon1_fragments",
        "gap_between_tss_fragments",
        "shared_body_exons2to7_fragments",
        "alt_exon1_to_exon2_junction_fragments",
        "canonical_exon1_to_exon2_junction_fragments",
        "unique_primary_locus_antisense_fragments",
        "alt_exon1_prdx4t_antisense_fragments",
        "canonical_exon1_antisense_fragments",
        "gap_between_tss_antisense_fragments",
        "shared_body_exons2to7_antisense_fragments",
        "alt_exon1_to_exon2_junction_antisense_fragments",
        "canonical_exon1_to_exon2_junction_antisense_fragments",
        "bam_path",
    ]

    sample_rows = aggregate_biological_samples(rows)
    main_sample_count = sum(row["include_in_main_plot"] for row in sample_rows)

    counts_path = output_dir / f"{args.dataset_id}_PRDX4t_raw_counts_per_library.tsv"
    metrics_path = output_dir / f"{args.dataset_id}_PRDX4t_normalized_metrics.tsv"
    sample_counts_path = output_dir / f"{args.dataset_id}_PRDX4t_raw_counts_per_biological_sample.tsv"
    sample_metrics_path = output_dir / f"{args.dataset_id}_PRDX4t_normalized_metrics_per_biological_sample.tsv"
    settings_path = output_dir / f"{args.dataset_id}_PRDX4t_quantification_manifest.tsv"
    atomic_write_tsv(counts_path, rows, fieldnames)
    atomic_write_tsv(metrics_path, derive_metrics(rows), METRIC_FIELDNAMES)
    atomic_write_tsv(sample_counts_path, sample_rows, SAMPLE_COUNT_FIELDNAMES)
    atomic_write_tsv(sample_metrics_path, derive_sample_metrics(sample_rows), SAMPLE_METRIC_FIELDNAMES)
    write_settings(
        settings_path,
        args,
        manifest_path,
        bam_dir,
        feature_path,
        len(rows),
        len(sample_rows),
        main_sample_count,
        strandedness,
        detected,
    )

    print(f"Raw per-library counts: {counts_path}")
    print(f"Normalized metrics: {metrics_path}")
    print(f"Biological-sample counts: {sample_counts_path}")
    print(f"Biological-sample normalized metrics: {sample_metrics_path}")
    print(f"Main-stage biological samples selected for plotting: {main_sample_count}")
    print(f"Quantification settings: {settings_path}")


if __name__ == "__main__":
    main()
PY
