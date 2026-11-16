#!/usr/bin/env python3
"""Gene-level exon-union quantification for the PRDX4t co-expression panel.

Reads the panel BAMs written by workflow step 03 and counts, per library and
per gene, the fragments whose aligned blocks touch any annotated exon of that
gene. Technical runs are then grouped by SampleID; biological replicates are
never combined.

The alignment-filtering, strand-inference and exon-union conventions are the
same ones script_workflow_04_quantify_PRDX4t_isoform_usage.sh applies to the
PRDX4 locus, so gene-level values here are comparable with the PRDX4t values
written there.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# Same thresholds as the PRDX4t quantifier: 0.75 requires 3:1 strand enrichment,
# and the unstranded band brackets equal representation of both orientations.
# Libraries between these bands are left undetermined.
STRAND_FORWARD_MIN = 0.75
STRAND_REVERSE_MAX = 0.25
STRAND_UNSTRANDED_BAND = (0.40, 0.60)
STRAND_MIN_EVIDENCE = 50
STRAND_LEAKY_WARN = 0.90


def parse_args():
    parser = argparse.ArgumentParser()
    for name in ("dataset-id", "dataset-label", "main-stage-groups", "manifest",
                 "bam-dir", "features", "output-dir", "strandedness"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--wt-sample-filter", default="")
    parser.add_argument("--min-mapq", required=True, type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def load_panel(path):
    """Genes with exon-union intervals, keyed by gene name."""
    required = {"feature_id", "gene_name", "role", "chrom", "start", "end",
                "strand", "gene_id", "gene_exonic_length_bp"}
    genes = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"ERROR: feature table is missing columns: {', '.join(missing)}")
        for row in reader:
            if "_exonunion" not in row["feature_id"]:
                continue
            gene = row["gene_name"]
            entry = genes.setdefault(gene, {
                "gene_name": gene,
                "gene_id": row["gene_id"],
                "role": row["role"],
                "chrom": row["chrom"],
                "strand": row["strand"],
                "exonic_length_bp": int(row["gene_exonic_length_bp"]),
                "exons": [],
            })
            if entry["chrom"] != row["chrom"] or entry["strand"] != row["strand"]:
                raise SystemExit(f"ERROR: inconsistent chrom/strand for {gene}")
            entry["exons"].append((int(row["start"]), int(row["end"])))
    if not genes:
        raise SystemExit(f"ERROR: no exon-union features in {path}")
    for gene, entry in genes.items():
        entry["exons"].sort()
        covered = sum(end - start + 1 for start, end in entry["exons"])
        if covered != entry["exonic_length_bp"]:
            raise SystemExit(
                f"ERROR: {gene} exon-union intervals cover {covered} bp but the table "
                f"records {entry['exonic_length_bp']} bp"
            )
    return genes


def use_alignment(read, min_mapq):
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.mapping_quality < min_mapq
    )


def touches_exons(blocks, exons):
    for start, end in exons:
        start0 = start - 1
        if any(block_end > start0 and block_start < end for block_start, block_end in blocks):
            return True
    return False


def mate1_is_reverse(read):
    """Orientation of mate 1, inferred from whichever mate this record is."""
    if not read.is_paired or read.is_read1:
        return read.is_reverse
    if not read.mate_is_unmapped:
        return read.mate_is_reverse
    return not read.is_reverse


def is_sense(read, gene_strand, strandedness):
    if strandedness == "unstranded":
        return True
    gene_on_minus = gene_strand == "-"
    if strandedness == "forward":
        return mate1_is_reverse(read) == gene_on_minus
    return mate1_is_reverse(read) != gene_on_minus


def fragment_orientations(bam_path, genes, min_mapq):
    """Per gene, the mate-1 orientation of every exon-touching fragment."""
    per_gene = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for gene, entry in genes.items():
            orientations = {}
            start = max(0, entry["exons"][0][0] - 1)
            end = entry["exons"][-1][1]
            for read in bam.fetch(entry["chrom"], start, end):
                if not use_alignment(read, min_mapq):
                    continue
                blocks = read.get_blocks()
                if not blocks or not touches_exons(blocks, entry["exons"]):
                    continue
                # One vote per fragment: mates of a pair agree by construction.
                orientations.setdefault(read.query_name, mate1_is_reverse(read))
            per_gene[gene] = orientations
    return per_gene


def detect_strandedness(per_gene, genes):
    """Call the protocol from every panel gene at once.

    Genes on opposite strands are folded together by asking whether mate 1 lies
    on each gene's own transcribed strand, so a mixed-strand panel still gives
    one coherent fraction.
    """
    on_gene_strand = 0
    total = 0
    for gene, orientations in per_gene.items():
        gene_on_minus = genes[gene]["strand"] == "-"
        for reverse in orientations.values():
            total += 1
            if reverse == gene_on_minus:
                on_gene_strand += 1
    if total < STRAND_MIN_EVIDENCE:
        return "undiagnosable", total, None
    fraction = on_gene_strand / total
    if fraction >= STRAND_FORWARD_MIN:
        return "forward", total, fraction
    if fraction <= STRAND_REVERSE_MAX:
        return "reverse", total, fraction
    if STRAND_UNSTRANDED_BAND[0] <= fraction <= STRAND_UNSTRANDED_BAND[1]:
        return "unstranded", total, fraction
    return "undiagnosable", total, fraction


def count_gene_fragments(per_gene, genes, strandedness):
    counts = {}
    for gene, orientations in per_gene.items():
        gene_on_minus = genes[gene]["strand"] == "-"
        sense = antisense = 0
        for reverse in orientations.values():
            if strandedness == "unstranded":
                sense += 1
            elif (reverse == gene_on_minus) == (strandedness == "forward"):
                sense += 1
            else:
                antisense += 1
        counts[gene] = (sense, antisense)
    return counts


def load_manifest(path, bam_dir):
    required = {"Run", "SampleName", "LibraryName", "LibraryLayout", "spots",
                "SampleID", "Group", "GroupLabel", "GroupOrder", "Replicate", "TechnicalRun"}
    records = []
    seen = set()
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"ERROR: manifest is missing columns: {', '.join(missing)}")
        for row in reader:
            run = row["Run"].strip()
            if not run or run in seen:
                raise SystemExit(f"ERROR: empty or duplicate Run in manifest: {run!r}")
            seen.add(run)
            if row["LibraryLayout"].strip() not in {"PAIRED", "SINGLE"}:
                raise SystemExit(f"ERROR: {run} has unsupported LibraryLayout")
            for column in ("spots", "GroupOrder", "Replicate", "TechnicalRun"):
                if not row[column].strip().isdigit():
                    raise SystemExit(f"ERROR: {run} has invalid {column}: {row[column]!r}")
            if int(row["spots"]) <= 0:
                raise SystemExit(f"ERROR: {run} has non-positive spots")
            bam = Path(bam_dir) / f"{run}.panel.bam"
            if not bam.is_file():
                raise SystemExit(
                    f"ERROR: panel BAM not found: {bam}\n"
                    "Run script_workflow_03_process_all_libraries.sh for this dataset first."
                )
            if not Path(f"{bam}.bai").is_file():
                raise SystemExit(f"ERROR: panel BAM index not found: {bam}.bai")
            row["bam_path"] = str(bam)
            records.append(row)
    if not records:
        raise SystemExit(f"ERROR: manifest lists no libraries: {path}")
    return records


def parse_wt_filter(text):
    """'Column=value[,value...]' -> (column, {values}); empty text keeps all."""
    text = (text or "").strip()
    if not text:
        return None, None
    if text.count("=") != 1:
        raise SystemExit(f"ERROR: WT_SAMPLE_FILTER must be Column=value[,value...], got {text!r}")
    column, values = text.split("=", 1)
    wanted = {value.strip() for value in values.split(",") if value.strip()}
    if not column.strip() or not wanted:
        raise SystemExit(f"ERROR: malformed WT_SAMPLE_FILTER: {text!r}")
    return column.strip(), wanted


def atomic_write_tsv(path, rows, fieldnames):
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main():
    args = parse_args()
    genes = load_panel(args.features)
    records = load_manifest(args.manifest, args.bam_dir)
    main_groups = {value.strip() for value in args.main_stage_groups.split(",") if value.strip()}
    filter_column, filter_values = parse_wt_filter(args.wt_sample_filter)

    if filter_column and filter_column not in records[0]:
        raise SystemExit(
            f"ERROR: WT_SAMPLE_FILTER names column {filter_column!r}, which the manifest "
            f"does not have: {args.manifest}"
        )

    gene_order = sorted(genes, key=lambda g: (genes[g]["role"] != "primary", g))
    print(f"Dataset: {args.dataset_id} ({args.dataset_label})")
    print(f"Panel genes: {', '.join(gene_order)}")
    print(f"Libraries: {len(records)}; minimum MAPQ {args.min_mapq}")
    if filter_column:
        print(f"Wild-type sample filter: {filter_column} in {sorted(filter_values)}")
    if args.dry_run:
        print("Dry run complete; no result files were written.")
        return

    # --- per library ------------------------------------------------------
    library_rows = []
    calls = {}
    for record in records:
        run = record["Run"].strip()
        per_gene = fragment_orientations(record["bam_path"], genes, args.min_mapq)
        call, evidence, fraction = detect_strandedness(per_gene, genes)
        calls[run] = (call, evidence, fraction)
        mode = args.strandedness if args.strandedness != "auto" else call
        if mode == "undiagnosable":
            raise SystemExit(
                f"ERROR: {run} strandedness could not be called from {evidence} fragments"
                + (f" (mate-1 sense fraction {fraction:.3f})" if fraction is not None else "")
                + "; pass --strandedness explicitly"
            )
        counts = count_gene_fragments(per_gene, genes, mode)
        shown = f"{fraction:.3f}" if fraction is not None else "n/a"
        print(f"  {run}: sense fraction {shown} over {evidence} fragments -> {mode}")
        for gene in gene_order:
            sense, antisense = counts[gene]
            library_rows.append({
                "run": run, "sample_id": record["SampleID"].strip(),
                "library_name": record["LibraryName"].strip(),
                "library_layout": record["LibraryLayout"].strip(),
                "technical_run": record["TechnicalRun"].strip(),
                "stage": record["Group"].strip(), "gene_name": gene,
                "gene_id": genes[gene]["gene_id"], "role": genes[gene]["role"],
                "strandedness": mode, "manifest_spots": record["spots"].strip(),
                "exonic_length_bp": genes[gene]["exonic_length_bp"],
                "sense_fragments": sense, "antisense_fragments": antisense,
            })

    if args.strandedness == "auto":
        # Each library is counted with its own call, which is what auto means.
        # Disagreement still matters, but only fatally when it reaches the
        # samples the analysis actually uses: a differently prepared library
        # sitting in an excluded group cannot corrupt a main-stage result.
        by_run = {record["Run"].strip(): record for record in records}
        distinct = {call for call, _, _ in calls.values()}
        if len(distinct) > 1:
            analysed = {
                call
                for run, (call, _, _) in calls.items()
                if by_run[run]["Group"].strip() in main_groups
            }
            odd = sorted(
                run for run, (call, _, _) in calls.items()
                if call != max(analysed or distinct, key=lambda c: sum(
                    1 for r, (x, _, _) in calls.items() if x == c))
            )
            if len(analysed) > 1:
                raise SystemExit(
                    "ERROR: libraries inside the main stage groups disagree about "
                    f"strandedness: {sorted(analysed)}; pass --strandedness explicitly"
                )
            print(
                f"NOTE: {len(odd)} librar(ies) outside the main stage groups were "
                f"prepared differently and are counted with their own call: "
                + ", ".join(f"{run} ({calls[run][0]})" for run in odd)
            )
        leaky = [run for run, (_, _, fraction) in calls.items()
                 if fraction is not None and 0.5 < fraction < STRAND_LEAKY_WARN]
        if leaky and distinct != {"unstranded"}:
            print(f"NOTE: {len(leaky)} stranded librar(ies) below "
                  f"{STRAND_LEAKY_WARN:.2f} sense fraction; antisense background is reported.")

    # --- grouped by biological sample -------------------------------------
    samples = {}
    for row in library_rows:
        key = (row["sample_id"], row["gene_name"])
        entry = samples.setdefault(key, {
            "sample_id": row["sample_id"], "gene_name": row["gene_name"],
            "gene_id": row["gene_id"], "role": row["role"],
            "exonic_length_bp": row["exonic_length_bp"],
            "library_layout": row["library_layout"], "strandedness": row["strandedness"],
            "sense_fragments": 0, "antisense_fragments": 0, "manifest_spots": 0,
            "technical_runs": 0, "source_runs": [],
        })
        entry["sense_fragments"] += row["sense_fragments"]
        entry["antisense_fragments"] += row["antisense_fragments"]
        entry["manifest_spots"] += int(row["manifest_spots"])
        entry["technical_runs"] += 1
        entry["source_runs"].append(row["run"])

    by_sample = {record["SampleID"].strip(): record for record in records}
    sample_rows = []
    for (sample_id, gene), entry in samples.items():
        record = by_sample[sample_id]
        million_spots = entry["manifest_spots"] / 1e6
        exonic_kb = entry["exonic_length_bp"] / 1000
        fragments = entry["sense_fragments"]
        in_main = record["Group"].strip() in main_groups
        in_wt = in_main and (
            filter_column is None or record.get(filter_column, "").strip() in filter_values
        )
        sample_rows.append({
            "sample_id": sample_id,
            "library_name": record["LibraryName"].strip(),
            "library_layout": entry["library_layout"],
            "stage": record["Group"].strip(),
            "stage_label": record["GroupLabel"].strip(),
            "stage_order": record["GroupOrder"].strip(),
            "replicate": record["Replicate"].strip(),
            "genotype": record.get("Genotype", "").strip(),
            "batch": record.get("Batch", "").strip(),
            "gene_name": gene,
            "gene_id": entry["gene_id"],
            "role": entry["role"],
            "exonic_length_bp": entry["exonic_length_bp"],
            "strandedness": entry["strandedness"],
            "technical_run_count": entry["technical_runs"],
            "source_runs": ",".join(sorted(entry["source_runs"])),
            "include_in_main_plot": int(in_main),
            "include_in_wt_analysis": int(in_wt),
            "manifest_spots": entry["manifest_spots"],
            "million_spots": f"{million_spots:.4f}",
            "fragments": fragments,
            "antisense_fragments": entry["antisense_fragments"],
            "fragments_per_million_spots": f"{fragments / million_spots:.4f}",
            "fragments_per_kb_per_million_spots":
                f"{fragments / exonic_kb / million_spots:.4f}",
        })
    sample_rows.sort(key=lambda row: (int(row["stage_order"]), row["sample_id"],
                                      gene_order.index(row["gene_name"])))

    output_dir = Path(args.output_dir)
    library_path = output_dir / f"{args.dataset_id}_panel_gene_raw_counts_per_library.tsv"
    sample_path = output_dir / f"{args.dataset_id}_panel_gene_expression_per_biological_sample.tsv"
    manifest_path = output_dir / f"{args.dataset_id}_panel_gene_quantification_manifest.tsv"

    atomic_write_tsv(library_path, library_rows, list(library_rows[0]))
    atomic_write_tsv(sample_path, sample_rows, list(sample_rows[0]))

    wt_samples = len({row["sample_id"] for row in sample_rows if row["include_in_wt_analysis"]})
    provenance = [
        ("dataset_id", args.dataset_id), ("dataset_label", args.dataset_label),
        ("panel_feature_table", args.features),
        ("panel_genes", ",".join(gene_order)),
        ("run_manifest", args.manifest), ("bam_directory", args.bam_dir),
        ("bam_name_pattern", "RUN.panel.bam (workflow step 03, reference mode panel)"),
        ("libraries_quantified", str(len(records))),
        ("biological_samples", str(len({row["sample_id"] for row in sample_rows}))),
        ("main_stage_groups", args.main_stage_groups),
        ("wt_sample_filter", args.wt_sample_filter or "(none)"),
        ("wt_samples_selected", str(wt_samples)),
        ("minimum_mapq", str(args.min_mapq)),
        ("strandedness_requested", args.strandedness),
        ("strandedness_applied", sorted({row["strandedness"] for row in sample_rows})[0]),
        ("counting_model", "exon-union fragments; one count per fragment per gene"),
        ("fragments", "raw sense fragments; NOT comparable between genes"),
        ("fragments_per_million_spots", "normalized for library depth only"),
        ("fragments_per_kb_per_million_spots",
         "normalized for depth and exonic length; use for cross-gene comparison"),
        ("normalization_denominator",
         "manifest spots (read pairs for PAIRED, reads for SINGLE)"),
    ]
    atomic_write_tsv(manifest_path, [{"field": f, "value": v} for f, v in provenance],
                     ["field", "value"])

    print(f"\nLibraries quantified: {len(records)}")
    print(f"Biological samples:   {len({row['sample_id'] for row in sample_rows})}"
          f" ({wt_samples} in the wild-type analysis)")
    print(f"Per-library counts:   {library_path}")
    print(f"Per-sample expression:{sample_path}")
    print(f"Manifest:             {manifest_path}")


if __name__ == "__main__":
    main()
