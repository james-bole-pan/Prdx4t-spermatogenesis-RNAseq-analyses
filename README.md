# Prdx4t spermatogenesis RNA-seq analyses

Source data, analysis scripts and figures for *Prdx4t*, canonical *Prdx4*, *Pcmt1* and *Crbn* across five mouse spermatogenesis studies: Soumillon (GSE43717), He (GSE220303), Trovero (PRJNA548952), Holmlund (GSE281037) and Gill (GSE214315). The analysis includes 51 biological samples or cell pools and 76 sequencing runs. *Prdx4t* refers to the alternative-first-exon transcript of *Prdx4*.

## Contents

Illustrator-compatible `.ai` versions of the figures are not tracked in this repository; the plotting scripts below regenerate them.


- `Source_Data_Supplementary_Figure_4.xlsx`: counts, expression estimates, stage means and SEMs, study selection, strandedness results and column definitions.
- `figures/RNAseq_nonstacked_relative.pdf`: the original multi-panel figure with relative expression values.
- `figures/RNAseq_stacked_relative.pdf`: relative expression, with the five studies stacked on a common stage axis.
- `figures/RNAseq_stacked_abundance_linear.pdf`: abundance before peak scaling, on a shared linear axis.
- `figures/RNAseq_stacked_abundance_asinh.pdf`: the same abundance estimates on a shared inverse-hyperbolic-sine axis.
- `analysis_code/`: plotting scripts, the raw-read processing workflow and the RSeQC strandedness assessment.
- `tables/`: stage mappings, abundance means and SEMs, and checks for the stacked figures.
- `reference/`: Ensembl-derived first-exon and retained exon-union coordinates used for counting.

## Reproduce the figures

From this folder, using Python 3.11:

```bash
python -m pip install -r analysis_code/requirements.txt
python analysis_code/reproduce_figure.py
python analysis_code/reproduce_stacked_figures.py
```

The first script creates the non-stacked figure; the second creates all three stacked figures. Both recalculate expression from the workbook and check the results against its sample and stage summaries. Outputs are saved in `figures/` as PDF and Illustrator-compatible EPS files with an `.ai` extension. Use Arial to match the supplied figures.

To check the calculations without writing files:

```bash
python analysis_code/reproduce_figure.py --check-only
python analysis_code/reproduce_stacked_figures.py --check-only
```

Both scripts accept `--workbook FILE`, `--output-dir DIR` and `--tables-dir DIR`. Default paths are relative to the scripts, so they can also be run from another directory.

The stacked script also writes:

- `tables/stacked_stage_alignment.tsv`: original stage labels, their spans on the common axis and marker positions.
- `tables/stacked_abundance_stage_summary.tsv`: abundance means, SEMs, sample sizes and units.
- `tables/stacked_figure_checks.json`: validation results, source-file hashes and plotting settings.

## Reading the stacked figures

The studies run from top to bottom in the order Soumillon, He, Trovero, Holmlund and Gill. Each heading gives the number of biological samples or cell pools per annotated stage. Technical sequencing runs from the same sample are combined before analysis.

The horizontal axis shows ordered stage categories. Its spacing does not represent developmental time. Circles mark single-stage or control populations. Diamonds mark combined-stage populations at a schematic midpoint; grey brackets show the stages covered by the original annotation. Open markers indicate *Prdx4t* and filled blue markers indicate canonical *Prdx4*. Lines connect the measured populations across gaps; no values are imputed for missing stages.

Grey shading marks Sertoli cells and the mixed 2C fraction. Red shading marks diplotene spermatocytes, secondary spermatocytes and spermatids. This interval was selected from the stage annotations and the observed expression patterns, without an automated cutoff. The shading describes the proposed expression window; a highlighted stage need not have a measurement in every study.

The linear and asinh figures show abundance before peak scaling. The asinh axis uses `linear_width=0.1`, includes zero and adds no pseudocount. Its tick labels give abundance values on the original scale.

## Analysis definitions

Abundance is measured in fragments per kilobase of retained exonic sequence per million sequenced SRA spots. One spot is one read pair for paired-end data or one read for single-end data. These are length- and depth-normalized estimates, not raw counts, TPM, conventional mapped-read FPKM or molecules per cell. Comparisons are strongest within studies.

For each biological sample, the *Prdx4t* fraction is the exon 1t-exon 2 junction count divided by the sum of alternative and canonical first-exon junction counts. All included samples meet the minimum of five total junction fragments. This fraction and its complement split total *Prdx4* abundance into estimated *Prdx4t* and canonical *Prdx4* components.

In the relative-expression figures, both isoforms are divided by the highest stage-mean total *Prdx4* abundance in their dataset. *Pcmt1* and *Crbn* each use their own highest stage mean. Error bars show SEM across biological samples. Soumillon has one pool per stage, so no SEM is shown.

The Trovero 2C fraction contains spermatogonia and testicular somatic cells. Holmlund is restricted to wild-type batch B. Study selection and limitations are recorded in the workbook.

## Reprocess raw sequencing data

The read-processing workflow uses fastp 0.23.4, Bowtie2 2.5.4, STAR 2.7.11b, samtools 1.20, SRA Toolkit 3.1.1 and pysam 0.22.1. Plotting uses NumPy 1.26.4 and Matplotlib 3.8.4. To prepare the Mamba environment:

```bash
bash analysis_code/workflow/script_environment_create_prdx4t_analysis.sh
mamba activate prdx4t-analysis
export SRA_ROOT="$PWD/data"
export PRDX4T_ENV_BIN="$CONDA_PREFIX/bin"
export PRDX4T_REFERENCE_DIR="$SRA_ROOT/reference"
```

Dataset configurations and SRA run manifests are in `analysis_code/workflow/datasets/`. Each accession has its own folder under `SRA_ROOT`. Process the reads in this order:

1. Download the five datasets with their `script_download_*` scripts.
2. Prepare the GRCm39/Ensembl release 110 reference and six-gene panel with workflows 01 and 01b, then build the STAR index with workflow 02.
3. Run workflow 03 in panel mode for each dataset.
4. Count isoform-specific junction fragments with workflow 04 and gene-level fragments with workflow 04b.

Each script provides `--help` and accepts explicit input/output paths. Reference-building scripts use `--output-dir "$PRDX4T_REFERENCE_DIR"`; the STAR-index script uses `--reference-dir "$PRDX4T_REFERENCE_DIR"`. For each accession, run:

```bash
bash analysis_code/workflow/script_workflow_03_process_all_libraries.sh --dataset ACCESSION --reference-mode panel
bash analysis_code/workflow/script_workflow_04_quantify_PRDX4t_isoform_usage.sh --dataset ACCESSION
bash analysis_code/workflow/script_workflow_04b_quantify_panel_gene_expression.sh --dataset ACCESSION
```

To plot from the resulting biological-sample counts:

```bash
python analysis_code/reproduce_figure.py --counts-root "$SRA_ROOT"
```

The panel contains *Prdx4*, *Pcmt1*, *Crbn*, *Pgk1*, *Pgk2* and *Actb* with 10-kb flanks. Processing enables fastp `--detect_adapter_for_pe` for paired-end libraries, uses Bowtie2 `--local --very-sensitive-local --no-unal`, and realigns retained candidate fragments to the full genome with STAR `--outFilterMultimapNmax 20`. The STAR index uses `--sjdbOverhang 149`. Counting retains primary MAPQ-255 alignments, excludes opposite-strand gene-span overlaps from all exon unions, and counts each fragment once per gene. Other analytical settings are defaults. Raw FASTQs, BAMs and genome indexes are obtained or generated by these scripts and are not duplicated in this package.

## Reassess library strandedness

The workbook's `Strandedness` sheet contains results for all 76 sequencing runs. RSeQC `infer_experiment.py` was run with the required BAM and BED12 paths and default analytical settings. The results agree with the strand assignments used for counting.

Use a separate Python 3.12 environment for the RSeQC dependencies:

```bash
python -m pip install -r analysis_code/requirements_strandedness.txt
python analysis_code/prepare_transcript_bed.py --gtf /path/to/Mus_musculus.GRCm39.110.gtf --output /path/to/transcripts.bed12
python analysis_code/infer_strandedness.py --bam-root /path/to/data --reference-bed /path/to/transcripts.bed12 --output /path/to/strandedness.tsv
```

The BED12 conversion includes all annotated transcripts. The assessment uses the retained panel-region BAMs, not whole-transcriptome alignments. Forward/reverse fractions are renormalized after excluding annotation-ambiguous records. A fraction of at least 0.8 identifies a stranded orientation; an absolute forward/reverse difference below 0.1 identifies an unstranded library. Intermediate results are undetermined. These interpretation rules follow [nf-core/rnaseq](https://nf-co.re/rnaseq/3.15.0/docs/usage#strandedness-prediction) and are separate from the default parameters of [RSeQC](https://rseqc.sourceforge.net/#infer-experiment-py).

## Citation

Manuscript in preparation. Citation details will be added on publication.

## License

Code in this repository is released under the [MIT License](LICENSE).
