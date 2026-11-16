"""Read the source workbook and calculate Supplementary Figure 4 expression."""

from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile
import csv
import math
import posixpath
import statistics
import xml.etree.ElementTree as ET

DATASET_ORDER = ['PRJNA548952', 'GSE220303', 'GSE281037', 'GSE214315', 'GSE43717']
MIN_JUNCTION_FRAGMENTS = 5
DEFAULT_WORKBOOK = Path(__file__).resolve().parents[1] / 'Source_Data_Supplementary_Figure_4.xlsx'
MAIN = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'


def read_workbook(path):
    """Read the flat, value-only worksheets using the Python standard library."""
    with ZipFile(path) as archive:
        strings = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            for item in ET.fromstring(archive.read('xl/sharedStrings.xml')):
                strings.append(''.join(item.itertext()))
        relationships = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        targets = {item.attrib['Id']: item.attrib['Target'] for item in relationships}
        book = ET.fromstring(archive.read('xl/workbook.xml'))
        result = {}
        for sheet in book.find(f'{{{MAIN}}}sheets'):
            target = targets[sheet.attrib[f'{{{REL}}}id']]
            target = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
            root = ET.fromstring(archive.read(target))
            matrix = []
            for row in root.findall(f'{{{MAIN}}}sheetData/{{{MAIN}}}row'):
                values = []
                for cell in row:
                    letters = ''.join(c for c in cell.attrib['r'] if c.isalpha())
                    col = 0
                    for char in letters: col = col * 26 + ord(char) - ord('A') + 1
                    while len(values) < col: values.append(None)
                    value = cell.find(f'{{{MAIN}}}v')
                    kind = cell.attrib.get('t')
                    if kind == 'inlineStr':
                        value = ''.join(cell.find(f'{{{MAIN}}}is').itertext())
                    elif value is None:
                        value = None
                    elif kind == 's':
                        value = strings[int(value.text)]
                    elif kind == 'b':
                        value = value.text == '1'
                    elif kind in ('str', 'e'):
                        value = value.text
                    else:
                        value = float(value.text)
                        if value.is_integer(): value = int(value)
                    values[col - 1] = value
                matrix.append(values)
            header, rows = matrix[0], matrix[1:]
            result[sheet.attrib['name']] = [dict(zip(header, row + [None] * (len(header) - len(row)))) for row in rows]
        return result


def read_tsv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def calculate_samples(workbook, counts_root=None):
    rows = [dict(row) for row in workbook['Sample_data']]
    if len(rows) != 51 or len({(r['dataset_id'], r['sample_id']) for r in rows}) != 51:
        raise ValueError('Expected the 51 distinct biological samples or pools in Figure 4.')
    if counts_root:
        for dataset in DATASET_ORDER:
            quant = Path(counts_root) / dataset / 'analysis/quantification'
            genes = {(r['sample_id'], r['gene_name']): r for r in read_tsv(quant / f'{dataset}_panel_gene_expression_per_biological_sample.tsv')}
            junctions = {r['sample_id']: r for r in read_tsv(quant / f'{dataset}_PRDX4t_raw_counts_per_biological_sample.tsv')}
            for row in (r for r in rows if r['dataset_id'] == dataset):
                sid = row['sample_id']
                for gene in ['Prdx4', 'Pcmt1', 'Crbn']:
                    value = genes[(sid, gene)]
                    row[gene + '_fragments'] = int(value['fragments'])
                    row[gene + '_exon_union_bp'] = int(value['exonic_length_bp'])
                    row['sequenced_SRA_spots'] = int(value['manifest_spots'])
                row['PRDX4t_junction_fragments'] = int(junctions[sid]['alt_exon1_to_exon2_junction_fragments'])
                row['canonical_PRDX4_junction_fragments'] = int(junctions[sid]['canonical_exon1_to_exon2_junction_fragments'])
    for row in rows:
        for gene in ['Prdx4', 'Pcmt1', 'Crbn']:
            depth = row['sequenced_SRA_spots'] / 1e6
            length = row[gene + '_exon_union_bp'] / 1000
            # The gene quantifier writes these values to four decimal places.
            row[gene + '_fragments_per_kb_per_million_spots'] = round(row[gene + '_fragments'] / length / depth, 4)
        alt = row['PRDX4t_junction_fragments']
        canonical = row['canonical_PRDX4_junction_fragments']
        if alt + canonical < MIN_JUNCTION_FRAGMENTS:
            raise ValueError(f"Insufficient junction evidence: {row['sample_id']}")
        row['PRDX4t_fraction'] = alt / (alt + canonical)
        total = row['Prdx4_fragments_per_kb_per_million_spots']
        row['PRDX4t_estimated_component'] = total * row['PRDX4t_fraction']
        row['canonical_PRDX4_estimated_component'] = total * (1 - row['PRDX4t_fraction'])
    for dataset in DATASET_ORDER:
        selected = [r for r in rows if r['dataset_id'] == dataset]
        for gene in ['Prdx4', 'Pcmt1', 'Crbn']:
            stages = defaultdict(list)
            for row in selected:
                stages[row['stage']].append(row[gene + '_fragments_per_kb_per_million_spots'])
            peak = max(statistics.mean(values) for values in stages.values())
            if peak <= 0: raise ValueError(f'No positive expression for {dataset}, {gene}')
            for row in selected:
                row[gene + '_peak_stage_mean'] = peak
                if gene == 'Prdx4':
                    row['PRDX4t_scaled'] = row['PRDX4t_estimated_component'] / peak
                    row['canonical_PRDX4_scaled'] = row['canonical_PRDX4_estimated_component'] / peak
                else:
                    row[gene + '_scaled'] = row[gene + '_fragments_per_kb_per_million_spots'] / peak
    expected = {(r['dataset_id'], r['sample_id']): r for r in workbook['Sample_data']}
    numeric = [k for k in rows[0] if any(term in k for term in ['_fragments', '_exon_union_bp', '_component', '_scaled', '_fraction', '_peak_stage_mean'])]
    numeric.append('sequenced_SRA_spots')
    for row in rows:
        original = expected[(row['dataset_id'], row['sample_id'])]
        for key in numeric:
            if not math.isclose(row[key], original[key], abs_tol=1e-10, rel_tol=1e-10):
                raise ValueError(f"Calculated value differs from published source data: {row['sample_id']}, {key}")
    return rows


def stage_summary(samples):
    output = []
    for dataset in DATASET_ORDER:
        selected = [r for r in samples if r['dataset_id'] == dataset]
        stages = sorted({(r['stage_order'], r['stage'], r['stage_label']) for r in selected})
        for order, stage, label in stages:
            group = [r for r in selected if r['stage'] == stage]
            for name, key in [('Prdx4t', 'PRDX4t_scaled'), ('Prdx4 canonical', 'canonical_PRDX4_scaled'), ('Pcmt1', 'Pcmt1_scaled'), ('Crbn', 'Crbn_scaled')]:
                values = [r[key] for r in group]
                output.append({'dataset_id': dataset, 'stage': stage, 'stage_label': label, 'stage_order': order,
                               'gene_or_isoform': name, 'n': len(values),
                               'mean_scaled_expression': statistics.mean(values),
                               'SEM_scaled_expression': statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else None,
                               'in_shaded_window': group[0]['in_shaded_window']})
    return output


def validate_summary(calculated, published):
    expected = {(r['dataset_id'], r['stage'], r['gene_or_isoform']): r for r in published}
    if len(expected) != len(calculated): raise ValueError('Stage-summary row count differs.')
    for row in calculated:
        old = expected[(row['dataset_id'], row['stage'], row['gene_or_isoform'])]
        if row['n'] != old['n']: raise ValueError('Stage sample-size mismatch.')
        for key in ['mean_scaled_expression', 'SEM_scaled_expression']:
            if row[key] is None and old[key] is None: continue
            if not math.isclose(row[key], old[key], abs_tol=1e-10, rel_tol=1e-10):
                raise ValueError(f'Stage-summary mismatch: {row}, {key}')
