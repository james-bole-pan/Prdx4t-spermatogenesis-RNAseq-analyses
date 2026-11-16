#!/usr/bin/env python3
"""Recalculate the published source data and draw Supplementary Figure 4."""

from pathlib import Path
from collections import defaultdict
import argparse
import csv
import math
import statistics
import textwrap
import sys

sys.dont_write_bytecode = True

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.offsetbox import TextArea, HPacker, VPacker, AnnotationBbox

from source_data import DEFAULT_WORKBOOK, DATASET_ORDER, read_workbook, calculate_samples, stage_summary, validate_summary

BLUE, GREEN, MAGENTA, PINK = '#0072B2', '#009E73', '#CC79A7', '#F3E2E8'
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42, 'ps.fonttype': 42})



def format_gene_text(figure, text, lines, *, bold=False, center=False):
    """Draw mixed text using real Arial italic glyphs, including numerals.

    Keep the reserved text position, but use ordinary font runs instead of
    mathematical text, whose treatment of digits differs from gene typography.
    Each line contains (text, italic) pairs.
    """
    renderer = figure.canvas.get_renderer()
    bounds = text.get_window_extent(renderer)
    rows = []
    for line in lines:
        runs = [TextArea(value, textprops={
            'fontfamily': 'Arial', 'fontsize': text.get_fontsize(),
            'fontweight': 'bold' if bold else 'normal',
            'fontstyle': 'italic' if italic else 'normal',
        }) for value, italic in line]
        rows.append(HPacker(children=runs, align='baseline', pad=0, sep=0))
    box = VPacker(children=rows, align='left', pad=0, sep=0)
    x = (bounds.x0 + bounds.x1) / 2 if center else bounds.x0
    y = (bounds.y0 + bounds.y1) / 2
    anchor = figure.transFigure.inverted().transform((x, y))
    replacement = AnnotationBbox(
        box, anchor, xycoords=figure.transFigure,
        box_alignment=(0.5 if center else 0, 0.5),
        frameon=False, pad=0, annotation_clip=False,
    )
    replacement.set_in_layout(False)
    figure.add_artist(replacement)
    text.set_visible(False)

def write_table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'figures')
    parser.add_argument('--counts-root', type=Path, help='Optional raw-workflow output root containing the five accession-named dataset folders.')
    parser.add_argument('--tables-dir', type=Path, help='Optionally export calculated sample and stage tables to this separate folder.')
    parser.add_argument('--check-only', action='store_true', help='Check calculations without writing figure files.')
    args = parser.parse_args()
    book = read_workbook(args.workbook)
    samples = calculate_samples(book, args.counts_root)
    summary = stage_summary(samples)
    validate_summary(summary, book['Stage_summary'])
    if args.tables_dir:
        args.tables_dir.mkdir(parents=True, exist_ok=True)
        write_table(args.tables_dir / 'samples.tsv', samples)
        write_table(args.tables_dir / 'stage_summary.tsv', summary)
    if args.check_only:
        print('Checked 51 biological samples or pools and 88 plotted stage means.')
        return
    figure, axes = plt.subplots(2, 3, figsize=(17, 10.5), constrained_layout=True)
    title = figure.suptitle('Prdx4t, Prdx4, Pcmt1 and Crbn across mouse spermatogenesis', fontsize=17, fontweight='bold', y=1.06)
    figure.get_layout_engine().set(hspace=0.05, wspace=0.04, h_pad=0.08)
    studies = {r['dataset_id']: r['study'] for r in book['Datasets']}
    series = [('PRDX4t_scaled', BLUE, 'Prdx4t', '-', 'o', True),
              ('canonical_PRDX4_scaled', BLUE, 'Prdx4 (canonical)', (0, (4, 2)), 's', False),
              ('Pcmt1_scaled', GREEN, 'Pcmt1', '-', 'o', True),
              ('Crbn_scaled', MAGENTA, 'Crbn', '-', 'o', True)]
    for axis, dataset in zip(axes.flat, DATASET_ORDER):
        selected = [r for r in samples if r['dataset_id'] == dataset]
        stages = sorted({(r['stage_order'], r['stage'], r['stage_label']) for r in selected})
        for key, color, label, line, marker, filled in series:
            means, errors = [], []
            for _, stage, _ in stages:
                values = [r[key] for r in selected if r['stage'] == stage]
                means.append(sum(values) / len(values))
                errors.append(statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0)
            axis.errorbar(range(len(stages)), means, yerr=errors, marker=marker, markersize=6.5,
                          linewidth=2.2, linestyle=line, color=color, label=label,
                          capsize=3.2, elinewidth=1.1, markerfacecolor=color if filled else 'white',
                          markeredgecolor='white' if filled else color, markeredgewidth=0.8 if filled else 1.5)
        flags = [next(r['in_shaded_window'] for r in selected if r['stage'] == stage) for _, stage, _ in stages]
        pos = 0
        while pos < len(flags):
            if not flags[pos]: pos += 1; continue
            start = pos
            while pos + 1 < len(flags) and flags[pos + 1]: pos += 1
            axis.axvspan(start - 0.5, pos + 0.5, color=PINK, zorder=0, linewidth=0)
            pos += 1
        labels = []
        for _, stage, label in stages:
            n = sum(r['stage'] == stage for r in selected)
            wrapped = textwrap.wrap(label, width=22, break_long_words=False, break_on_hyphens=False)
            labels.append('\n'.join(wrapped + [f'(n = {n})']))
        axis.set_title(f'{studies[dataset]} et al. ({dataset})', fontsize=15, fontweight='bold', pad=14)
        axis.set_xticks(range(len(stages)))
        axis.set_xticklabels(labels, fontsize=10.5, rotation=35, ha='right', rotation_mode='anchor', multialignment='right')
        axis.set_xlim(-0.45, len(stages) - 0.55); axis.set_ylim(-0.08, 1.15)
        axis.tick_params(axis='y', labelsize=10)
        # Reserve the published label space before fixing the plot positions.
        if axis in axes[:, 0]:
            axis.set_ylabel('Relative expression\n(fraction of each gene\u2019s dataset maximum)', fontsize=12, labelpad=12)
        else:
            axis.tick_params(axis='y', labelleft=False)
        axis.grid(color='#E4E4E4', linewidth=0.7); axis.set_axisbelow(True)
        for side in ('top', 'right'): axis.spines[side].set_visible(False)
    axes.flat[5].set_axis_off()
    handles, labels = axes.flat[0].get_legend_handles_labels()
    handles.insert(0, Patch(facecolor=PINK, label='Prdx4t window (diplotene, secondary' + '\n' + 'spermatocyte, spermatid)'))
    legend = axes.flat[5].legend(handles=handles, fontsize=14, frameon=False, loc='upper left', bbox_to_anchor=(0.0, 0.88), labelspacing=1.2, handlelength=2.6)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    crop = figure.get_tightbbox(renderer).padded(plt.rcParams['savefig.pad_inches'])
    centers = []
    for axis in axes[:, 0]:
        bounds = axis.yaxis.label.get_window_extent(renderer)
        centers.append(figure.transFigure.inverted().transform(((bounds.x0 + bounds.x1)/2, (bounds.y0 + bounds.y1)/2)))
    positions = [axis.get_position().frozen() for axis in axes.flat]
    figure.set_layout_engine(None)
    for axis, position in zip(axes.flat, positions): axis.set_position(position)
    for axis, center in zip(axes[:, 0], centers):
        axis.yaxis.label.set_text('Scaled relative expression')
        axis.yaxis.label.set_horizontalalignment('center')
        axis.yaxis.label.set_verticalalignment('center')
        axis.yaxis.set_label_coords(center[0], center[1], transform=figure.transFigure)
    # Apply true italic font runs after fixing the published panel positions.
    figure.canvas.draw()
    format_gene_text(figure, title, [[
        ('Prdx4t', True), (', ', False), ('Prdx4', True), (', ', False),
        ('Pcmt1', True), (' and ', False), ('Crbn', True),
        (' across mouse spermatogenesis', False),
    ]], bold=True, center=True)
    legend_lines = [
        [[('Prdx4t', True), (' window (diplotene, secondary', False)],
         [('spermatocyte, spermatid)', False)]],
        [[('Prdx4t', True)]],
        [[('Prdx4', True), (' (canonical)', False)]],
        [[('Pcmt1', True)]],
        [[('Crbn', True)]],
    ]
    for text, lines in zip(legend.get_texts(), legend_lines):
        format_gene_text(figure, text, lines)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_dir / 'RNAseq_nonstacked_relative.pdf', bbox_inches=crop, dpi=300)
    # Illustrator opens the same EPS-compatible vector format as the supplied AI.
    figure.savefig(args.output_dir / 'RNAseq_nonstacked_relative.ai', format='eps', bbox_inches=crop, dpi=300)
    plt.close(figure)
    print(f'Checked 51 samples and 88 plotted means; saved PDF and AI in {args.output_dir}.')


if __name__ == '__main__': main()
