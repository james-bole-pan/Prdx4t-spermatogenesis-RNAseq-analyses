#!/usr/bin/env python3
"""Validate source data and reproduce all three stacked RNA-seq figures."""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import statistics
import sys

sys.dont_write_bytecode = True

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import transforms
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as path_effects

from source_data import (
    DEFAULT_WORKBOOK, DATASET_ORDER, read_workbook, calculate_samples,
    stage_summary, validate_summary,
)
from reproduce_figure import format_gene_text

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
UNITS = 'fragments per kilobase of retained exonic sequence per million sequenced SRA spots'
PLOT_DATASET_ORDER = ['GSE43717', 'GSE220303', 'PRJNA548952', 'GSE281037', 'GSE214315']
STAGES = [
    'Sertoli\ncells', 'Mixed 2C\nfraction', 'Spermatogonia', 'Leptotene',
    'Zygotene', 'Pachytene', 'Diplotene', 'Secondary\nspermatocytes',
    'Round\nspermatids', 'Early elongating\nspermatids',
    'Late elongating\nspermatids', 'Mature\nspermatozoa',
]
# Bounds reflect annotation breadth, not measured developmental durations or
# a known mixture composition. A broad population contributes only one point.
MAPPING = {
    'PRJNA548952': {'2C': (1, 1), 'LZ': (3, 4), 'PS': (5, 5), 'RS': (8, 8)},
    'GSE220303': {'SGA': (2, 2), 'MLEP': (3, 3), 'MZYG': (4, 4), 'MPAC': (5, 5), 'MDIP': (6, 6)},
    'GSE281037': {'SC1': (3, 6), 'SC2': (7, 7), 'RS': (8, 8)},
    'GSE214315': {'LZ': (3, 4), 'PD': (5, 6), 'RS': (8, 8), 'EES': (9, 9), 'LES': (10, 10)},
    'GSE43717': {'SL': (0, 0), 'SG': (2, 2), 'SC': (3, 7), 'SD': (8, 10), 'SZ': (11, 11)},
}
RANGE_LABELS = {
    'LZ': 'Leptotene/\nzygotene', 'PD': 'Pachytene/\ndiplotene',
    'SC1': 'Primary spermatocytes', 'SC': 'Spermatocytes', 'SD': 'Spermatids',
}
SERIES = [
    ('Prdx4t', '#0072B2', '-', False),
    ('Prdx4 canonical', '#0072B2', (0, (4, 2)), True),
    ('Pcmt1', '#009E73', '-', True),
    ('Crbn', '#CC79A7', '-', True),
]
ABUNDANCE_KEYS = {
    'Prdx4t': ('PRDX4t_estimated_component', 'Prdx4_peak_stage_mean'),
    'Prdx4 canonical': ('canonical_PRDX4_estimated_component', 'Prdx4_peak_stage_mean'),
    'Pcmt1': ('Pcmt1_fragments_per_kb_per_million_spots', 'Pcmt1_peak_stage_mean'),
    'Crbn': ('Crbn_fragments_per_kb_per_million_spots', 'Crbn_peak_stage_mean'),
}
GREY_REGION, RED_REGION = '#F3F4F6', '#FCE2E2'
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 14.5, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'axes.linewidth': 0.7,
})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def summarize_abundance(samples, relative):
    reference = {(r['dataset_id'], r['stage'], r['gene_or_isoform']): r for r in relative}
    rows = []
    for dataset in DATASET_ORDER:
        selected = [r for r in samples if r['dataset_id'] == dataset]
        for order, stage, label in sorted({(r['stage_order'], r['stage'], r['stage_label']) for r in selected}):
            group = [r for r in selected if r['stage'] == stage]
            for gene, (value_key, peak_key) in ABUNDANCE_KEYS.items():
                values = [r[value_key] for r in group]
                mean = statistics.mean(values)
                sem = statistics.stdev(values)/math.sqrt(len(values)) if len(values) > 1 else None
                peaks = {r[peak_key] for r in group}
                assert len(peaks) == 1
                peak = peaks.pop()
                old = reference[(dataset, stage, gene)]
                assert math.isclose(mean/peak, old['mean_scaled_expression'], rel_tol=1e-10, abs_tol=1e-10)
                if sem is None:
                    assert old['SEM_scaled_expression'] is None
                else:
                    assert math.isclose(sem/peak, old['SEM_scaled_expression'], rel_tol=1e-10, abs_tol=1e-10)
                rows.append({
                    'dataset_id': dataset, 'stage': stage, 'stage_label': label,
                    'stage_order': order, 'gene_or_isoform': gene, 'n': len(values),
                    'mean_abundance': mean, 'SEM_abundance': sem, 'units': UNITS,
                })
    assert len(rows) == 88
    for sample in samples:
        assert math.isclose(sample['PRDX4t_estimated_component'] + sample['canonical_PRDX4_estimated_component'],
                            sample['Prdx4_fragments_per_kb_per_million_spots'], rel_tol=1e-12, abs_tol=1e-12)
    return rows


def stage_alignment(relative):
    rows = []
    for dataset in PLOT_DATASET_ORDER:
        populations = [r for r in relative if r['dataset_id'] == dataset and r['gene_or_isoform'] == 'Prdx4t']
        assert {r['stage'] for r in populations} == set(MAPPING[dataset])
        for row in sorted(populations, key=lambda r: r['stage_order']):
            low, high = MAPPING[dataset][row['stage']]
            rows.append({
                'dataset_id': dataset, 'original_stage': row['stage'],
                'original_label': row['stage_label'], 'n': row['n'],
                'axis_start': low, 'axis_end': high, 'marker_x': (low+high)/2,
                'display': 'single stage/control' if low == high else 'one diamond at a schematic midpoint',
            })
    assert len(rows) == 22
    return rows


def add_bottom_legends(fig, center):
    gene_handles = [Line2D([], [], color=color, marker='o', linestyle=line,
                          markersize=7.3, markerfacecolor=color if filled else 'white',
                          markeredgecolor=color, lw=1.7, label=name)
                    for name, color, line, filled in SERIES]
    gene_legend = fig.legend(handles=gene_handles, loc='upper center',
                             bbox_to_anchor=(center, 0.145), ncol=4,
                             frameon=False, fontsize=16, columnspacing=1.3, handlelength=1.8)
    single_handles = (
        Line2D([], [], marker='o', linestyle='none', color='#4B5563', markersize=7.2),
        Line2D([], [], marker='o', linestyle='none', color='#4B5563', markerfacecolor='white', markersize=7.2),
    )
    broad_handle = Line2D([], [], marker='D', linestyle='none', color='#303B48',
                          markerfacecolor='#7A8795', markersize=8.2)
    fig.legend(handles=[single_handles, broad_handle],
               labels=['Single-stage', 'Combined stages (schematic midpoint)'],
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.35)},
               loc='upper center', bbox_to_anchor=(center, 0.105), ncol=2,
               frameon=False, fontsize=15, columnspacing=1.5, handlelength=1.8)
    fig.legend(handles=[
        Patch(facecolor=GREY_REGION, edgecolor='#AFB6BF', label='Somatic / mixed\npopulations'),
        Patch(facecolor=RED_REGION, edgecolor='#EAA4A4',
              label='Highlighted stages: diplotene / secondary\nspermatocytes and spermatids'),
    ], loc='upper center', bbox_to_anchor=(center, 0.066), ncol=2,
       frameon=False, fontsize=15, columnspacing=1.8, handlelength=1.8)
    return gene_legend


def draw_series(ax, rows, dataset, series, mean_key, sem_key):
    gene, color, linestyle, filled = series
    points = []
    for row in sorted(rows, key=lambda r: r['stage_order']):
        low, high = MAPPING[dataset][row['stage']]
        points.append({'dataset': dataset, 'gene': gene, 'stage': row['stage'],
                       'x': (low+high)/2, 'mean': row[mean_key], 'sem': row[sem_key], 'broad': low != high})
    ax.plot([p['x'] for p in points], [p['mean'] for p in points], color=color, lw=1.9, linestyle=linestyle, zorder=2)
    for broad in [False, True]:
        selected = [p for p in points if p['broad'] == broad]
        if not selected:
            continue
        marker = ax.scatter([p['x'] for p in selected], [p['mean'] for p in selected],
                            s=72 if broad else 40, marker='D' if broad else 'o',
                            facecolors=color if filled else 'white',
                            edgecolors=('#303B48' if filled else color) if broad else color,
                            linewidths=1.1 if broad else 1, zorder=6 if broad else 5)
        if broad:
            marker.set_path_effects([path_effects.Stroke(linewidth=2.8, foreground='white'), path_effects.Normal()])
    for point in points:
        if point['sem'] is not None:
            ax.errorbar(point['x'], point['mean'], yerr=point['sem'], fmt='none',
                        ecolor=color, capsize=2.7, elinewidth=0.95, zorder=4)
    return points


def verify_layout(fig, axes):
    renderer = fig.canvas.get_renderer()
    title_box = fig._suptitle.get_window_extent(renderer)
    assert fig._suptitle.get_ha() == 'center'
    assert abs((title_box.x0+title_box.x1)/2-fig.bbox.width/2) < 0.5, 'Overall title is not centered on the figure'
    assert title_box.x0 >= 0 and title_box.x1 <= fig.bbox.width, 'Clipped overall title'
    min_gap = 6*fig.dpi/72
    for ax in axes:
        ticks = sorted([t.get_window_extent(renderer) for t in ax.get_yticklabels() if t.get_visible()], key=lambda b: b.y0)
        assert all(b.y0 >= a.y1+1 for a, b in zip(ticks, ticks[1:])), 'Overlapping y-axis labels'
        assert ax.title.get_ha() == 'center'
        title_box = ax.title.get_window_extent(renderer)
        for annotation in ax.texts:
            annotation_box = annotation.get_window_extent(renderer)
            assert title_box.y0-annotation_box.y1 >= min_gap, 'Insufficient gap between title and combined-stage annotation'
        for text in [ax.title, *ax.get_xticklabels(), *ax.get_yticklabels(), *ax.texts]:
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                assert box.x0 >= 0 and box.x1 <= fig.bbox.width+1, f'Horizontal text clipping: {text.get_text()}'
                assert box.y0 >= 0 and box.y1 <= fig.bbox.height+1, f'Vertical text clipping: {text.get_text()}'
    for previous, current in zip(axes, axes[1:]):
        assert previous.bbox.y0-current.title.get_window_extent(renderer).y1 >= 18*fig.dpi/72, 'Insufficient gap above the next panel title'
    boxes = [legend.get_window_extent(renderer) for legend in fig.legends]
    assert all(b.x0 >= 0 and b.x1 <= fig.bbox.width and b.y0 >= 0 for b in boxes), 'Clipped legend'
    assert all(a.y0 > b.y1 for a, b in zip(boxes, boxes[1:])), 'Overlapping legend rows'


def draw_figure(rows, studies, mode, output_dir):
    relative = mode == 'relative'
    mean_key = 'mean_scaled_expression' if relative else 'mean_abundance'
    sem_key = 'SEM_scaled_expression' if relative else 'SEM_abundance'
    left = 0.120 if relative else 0.150
    center = (left+0.980)/2
    fig, axes = plt.subplots(5, 1, sharex=True, sharey=True, figsize=(10.8, 17.4))
    fig.subplots_adjust(left=left, right=0.980, top=0.915, bottom=0.300, hspace=0.75)
    max_upper = max(row[mean_key]+(row[sem_key] or 0) for row in rows)
    ceiling = None if relative else math.ceil(max_upper/50)*50
    points = []
    for ax, dataset in zip(axes, PLOT_DATASET_ORDER):
        selected = [row for row in rows if row['dataset_id'] == dataset]
        populations = [row for row in selected if row['gene_or_isoform'] == 'Prdx4t']
        sizes = {row['n'] for row in populations}
        assert len(sizes) == 1
        n = sizes.pop()
        n_label = f'n = {n} per stage'
        ax.set_title(f'{studies[dataset]} et al.  ({dataset}; {n_label})',
                     loc='center', fontsize=16.2, fontweight='bold', pad=24)
        if relative:
            ticks, tick_labels = [0, 0.5, 1], ['0', '0.5', '1.0']
            ax.set_ylim(-0.065, 1.36)
        elif mode == 'asinh':
            ax.set_yscale('asinh', linear_width=0.1, base=10)
            forward, inverse = ax.yaxis.get_transform().transform, ax.yaxis.get_transform().inverted().transform
            ax.set_ylim(float(inverse(-0.045*forward(ceiling))), float(inverse(1.36*forward(ceiling))))
            ticks, tick_labels = [0, 1, 10, 100], ['0', '1', '10', '100']
        else:
            ticks = [0, 100, 200, ceiling]
            tick_labels = [str(int(tick)) for tick in ticks]
            ax.set_ylim(-0.045*ceiling, 1.36*ceiling)
        ax.set_yticks(ticks); ax.set_yticklabels(tick_labels, fontsize=14); ax.minorticks_off()
        ax.axvspan(-0.5, 1.5, color=GREY_REGION, zorder=-5)
        ax.axvspan(5.5, 10.5, color=RED_REGION, zorder=-4, linewidth=0)
        ax.axvline(1.5, color='#C4C8CF', lw=0.9, ls=(0, (3, 3)), zorder=-3)
        ax.set_axisbelow(True); ax.set_xticks(range(len(STAGES)))
        ax.grid(axis='x', color='#E7E9EC', linewidth=0.7)
        for y in ticks:
            ax.axhline(y, color='#E4E7EB', lw=0.65, zorder=-3)
        ax.set_xlim(-0.52, 11.52)
        ax.tick_params(axis='x', length=0, pad=8)
        ax.tick_params(axis='y', length=3, color='#9DA3AB')
        for edge in ['top', 'right']:
            ax.spines[edge].set_visible(False)
        for edge in ['left', 'bottom']:
            ax.spines[edge].set_color('#B0B6BE')
        for population in populations:
            low, high = MAPPING[dataset][population['stage']]
            if low == high:
                continue
            x1, x2 = low-0.27, high+0.27
            if relative:
                bracket_y, label_y, transform = [1.11, 1.15, 1.15, 1.11], 1.18, ax.transData
            else:
                bracket_y, label_y = [0.755, 0.790, 0.790, 0.755], 0.810
                transform = transforms.blended_transform_factory(ax.transData, ax.transAxes)
            ax.plot([x1, x1, x2, x2], bracket_y, color='#89919C', lw=0.8, transform=transform, clip_on=False)
            ax.text((low+high)/2, label_y, RANGE_LABELS[population['stage']], ha='center', va='bottom',
                    fontsize=12, color='#626B77', transform=transform)
        for series in SERIES:
            gene_rows = [row for row in selected if row['gene_or_isoform'] == series[0]]
            points.extend(draw_series(ax, gene_rows, dataset, series, mean_key, sem_key))
    keys = {(p['dataset'], p['gene'], p['stage']) for p in points}
    markers = sum(len(collection.get_offsets()) for ax in axes for collection in ax.collections if isinstance(collection, PathCollection))
    assert len(points) == len(keys) == markers == 88
    assert sum(point['broad'] for point in points) == 24
    if not relative:
        assert max_upper < ceiling
    axes[-1].set_xticklabels(STAGES, rotation=65, ha='right', rotation_mode='anchor', fontsize=14.5)
    y_label = 'Scaled relative expression' if relative else 'Abundance estimate\nFragments / kb / million sequenced spots'
    fig.text(0.025 if relative else 0.034, 0.608, y_label, rotation=90,
             va='center', ha='center', fontsize=17.5 if relative else 16.5, linespacing=1.55)
    title = 'Gene expression across mouse spermatogenesis' if relative else f'RNA abundance before peak scaling - {mode} scale'
    fig.suptitle(title, x=0.5, y=0.982, fontsize=22, fontweight='bold', ha='center', va='top')
    fig.text(center, 0.177, 'Developmental stage', fontsize=17, ha='center')
    gene_legend = add_bottom_legends(fig, center)
    fig.canvas.draw()
    legend_runs = [[[('Prdx4t', True)]], [[('Prdx4', True), (' (canonical)', False)]],
                   [[('Pcmt1', True)]], [[('Crbn', True)]]]
    for label, runs in zip(gene_legend.get_texts(), legend_runs):
        format_gene_text(fig, label, runs)
    fig.canvas.draw()
    verify_layout(fig, axes)
    name = 'RNAseq_stacked_relative' if relative else f'RNAseq_stacked_abundance_{mode}'
    files = []
    for extension, export_format in [('pdf', 'pdf'), ('ai', 'eps')]:
        target = output_dir/f'{name}.{extension}'
        fig.savefig(target, format=export_format, dpi=130, facecolor='white')
        files.append(target.name)
    plt.close(fig)
    return {'mode': mode, 'files': files, 'markers_verified': markers,
            'overall_title_centered': True, 'panel_headings_centered': True, 'legends_at_bottom': True,
            'heading_sample_size_label': 'n = x per stage',
            'figure_size_inches': [10.8, 17.4], 'panel_hspace': 0.75,
            'title_padding_points': 24, 'minimum_checked_label_gap_points': 6,
            'minimum_checked_gap_above_title_points': 18,
            'shared_y_axis': True, 'asinh_linear_width': 0.1 if mode == 'asinh' else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument('--output-dir', type=Path, default=PACKAGE_ROOT/'figures')
    parser.add_argument('--tables-dir', type=Path, default=PACKAGE_ROOT/'tables')
    parser.add_argument('--check-only', action='store_true', help='Validate all values without writing any outputs.')
    args = parser.parse_args()
    book = read_workbook(args.workbook)
    samples = calculate_samples(book)
    relative = stage_summary(samples)
    validate_summary(relative, book['Stage_summary'])
    abundance = summarize_abundance(samples, relative)
    alignment = stage_alignment(relative)
    if args.check_only:
        print('Checked 51 biological samples/pools, 88 relative means/SEMs, 88 abundance means/SEMs, and 22 stage alignments; no files written.')
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    studies = {row['dataset_id']: row['study'] for row in book['Datasets']}
    outputs = [draw_figure(relative, studies, 'relative', args.output_dir)]
    outputs.extend(draw_figure(abundance, studies, mode, args.output_dir) for mode in ['linear', 'asinh'])
    write_table(args.tables_dir/'stacked_stage_alignment.tsv', alignment)
    write_table(args.tables_dir/'stacked_abundance_stage_summary.tsv', abundance)
    report = {
        'biological_samples_or_pools_validated': len(samples), 'relative_means_and_sems_validated': len(relative),
        'abundance_means_and_sems_validated': len(abundance), 'original_populations': len(alignment),
        'panel_order': [studies[dataset] for dataset in PLOT_DATASET_ORDER],
        'abundance_units': UNITS,
        'abundance_interpretation': 'Length/depth-normalized abundance before peak scaling; not raw counts or molecule counts. Comparisons are strongest within studies.',
        'broad_stage_display': 'One diamond plus SEM at a schematic midpoint; grey brackets indicate original annotation breadth. Single-stage/control points use circles. Prdx4t is open; canonical Prdx4 is filled.',
        'grey_region': 'Somatic/mixed populations: Sertoli cells and mixed 2C fraction.',
        'red_region': 'Stage-based highlight: diplotene and secondary spermatocytes and spermatids; not a cutoff-based expression call.',
        'connecting_lines': 'Connect observed populations; missing stages have no measurements or imputed markers.',
        'asinh_interpretation': 'Zero-preserving inverse hyperbolic sine scale; linear_width=0.1; no pseudocount.',
        'ai_format': 'Illustrator-compatible EPS, matching the existing package convention; not a native Illustrator document.',
        'source_hashes': {args.workbook.name: sha(args.workbook), 'source_data.py': sha(Path(__file__).with_name('source_data.py')),
                          'reproduce_figure.py': sha(Path(__file__).with_name('reproduce_figure.py')), Path(__file__).name: sha(__file__)},
        'figures': outputs,
    }
    (args.tables_dir/'stacked_figure_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Created 3 stacked figures, each in PDF and Illustrator-compatible AI; validated all 51 samples and 88 stage estimates.')


if __name__ == '__main__':
    main()
