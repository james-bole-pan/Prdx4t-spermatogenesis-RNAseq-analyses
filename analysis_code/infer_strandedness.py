#!/usr/bin/env python3
"""Assess the included sequencing runs using default-parameter RSeQC."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse
import csv
import os
import re
import subprocess
import tempfile
import sys

sys.dont_write_bytecode = True

from source_data import DEFAULT_WORKBOOK, read_workbook


def parse_result(stdout, stderr):
    values = [float(value) for value in re.findall(r'Fraction of reads explained by [^\n]+: ([0-9.eE+-]+)', stdout)]
    ambiguous = re.search(r'Fraction of reads failed to determine: ([0-9.eE+-]+)', stdout)
    sampled = re.search(r'Total (\d+) usable reads were sampled', stdout+'\n'+stderr)
    if len(values) != 2 or ambiguous is None or sampled is None:
        raise ValueError('Could not parse the RSeQC orientation report.')
    forward, reverse = values
    unclassified = float(ambiguous.group(1))
    if abs(forward+reverse+unclassified-1) > 0.00021:
        raise ValueError('RSeQC orientation fractions do not sum to one.')
    classified = forward+reverse
    normalized = forward/classified if classified else 0
    if classified == 0: call = 'undetermined'
    elif normalized >= 0.8: call = 'forward'
    elif 1-normalized >= 0.8: call = 'reverse'
    elif abs(2*normalized-1) < 0.1: call = 'unstranded'
    else: call = 'undetermined'
    return {'sampled_alignment_records':int(sampled.group(1)),
            'fraction_forward_raw':forward, 'fraction_reverse_raw':reverse,
            'fraction_undetermined':unclassified,
            'fraction_forward_among_classifiable':normalized, 'rseqc_call':call}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument('--bam-root', type=Path, required=True, help='Folder containing the included RUN.panel.bam files, searched recursively.')
    parser.add_argument('--reference-bed', type=Path, required=True)
    parser.add_argument('--rseqc', default='infer_experiment.py')
    parser.add_argument('--workers', type=int, default=6, help='Number of independent sequencing runs processed concurrently.')
    parser.add_argument('--output', type=Path, required=True, help='Tab-separated results file; replaced on a subsequent run.')
    args = parser.parse_args()
    book = read_workbook(args.workbook)
    rows = book['Strandedness']
    if len(rows) != 76 or len({r['run'] for r in rows}) != 76:
        raise ValueError('Expected the 76 sequencing runs included in Figure 4.')
    expected = {run for row in book['Sample_data'] for run in row['source_runs'].split(',')}
    if expected != {row['run'] for row in rows}: raise ValueError('Sample and sequencing-run lists differ.')
    candidates = {}
    for path in args.bam_root.rglob('*.panel.bam'):
        run = path.name.removesuffix('.panel.bam')
        if run not in expected: continue
        if run in candidates: raise ValueError(f'Multiple BAMs found for {run}. Use a more specific --bam-root.')
        candidates[run] = path.resolve()
    missing = sorted(expected-candidates.keys())
    if missing: raise FileNotFoundError('Missing panel BAMs: '+', '.join(missing))
    if not args.reference_bed.is_file(): raise FileNotFoundError(args.reference_bed)
    with tempfile.TemporaryDirectory(prefix='rseqc_') as cache:
        env = {**os.environ, 'MPLCONFIGDIR':cache, 'PYTHONDONTWRITEBYTECODE':'1', 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'}
        def assess(row):
            # Only required paths are supplied; every analytical setting is a RSeQC default.
            command = [args.rseqc, '-i', str(candidates[row['run']]), '-r', str(args.reference_bed.resolve())]
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env, timeout=600)
            if result.returncode:
                raise RuntimeError(f"RSeQC failed for {row['run']}: {result.stderr}")
            return {key:row[key] for key in ['dataset_id','sample_id','run','library_layout']} | parse_result(result.stdout, result.stderr)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(assess, rows))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(results)
    differences = [(row['run'], old['rseqc_call'], row['rseqc_call']) for old,row in zip(rows,results) if old['rseqc_call'] != row['rseqc_call']]
    print(f'Assessed {len(results)} sequencing runs using RSeQC defaults; {len(differences)} calls differ from the source workbook.')
    if differences:
        for run, old, new in differences: print(f'{run}: published {old}; assessed {new}')


if __name__ == '__main__': main()
