#!/usr/bin/env python3
"""Convert Ensembl exon annotations to transcript BED12 for RSeQC."""

from pathlib import Path
import argparse
import gzip
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf', required=True, type=Path, help='Ensembl release 110 GRCm39 GTF, plain or gzip-compressed.')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    models = {}
    opener = gzip.open if args.gtf.suffix == '.gz' else open
    with opener(args.gtf, 'rt') as stream:
        for line in stream:
            if line.startswith('#'): continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'exon': continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not match: continue
            transcript = match.group(1)
            chrom, start, end, strand = fields[0], int(fields[3])-1, int(fields[4]), fields[6]
            model = models.setdefault(transcript, {'chrom':chrom, 'strand':strand, 'exons':[]})
            if (model['chrom'], model['strand']) != (chrom, strand):
                raise ValueError(f'Inconsistent transcript coordinates: {transcript}')
            model['exons'].append((start, end))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as stream:
        for transcript, model in sorted(models.items(), key=lambda item:(item[1]['chrom'], min(x[0] for x in item[1]['exons']), item[0])):
            blocks = []
            for start, end in sorted(set(model['exons'])):
                if blocks and start <= blocks[-1][1]:
                    blocks[-1] = (blocks[-1][0], max(end, blocks[-1][1]))
                else: blocks.append((start, end))
            start, end = blocks[0][0], blocks[-1][1]
            values = [model['chrom'], start, end, transcript, 0, model['strand'], start, end, 0, len(blocks),
                      ','.join(str(b-a) for a,b in blocks)+',', ','.join(str(a-start) for a,b in blocks)+',']
            stream.write('\t'.join(map(str,values))+'\n')
    print(f'Wrote {len(models)} annotated transcript models to {args.output}.')


if __name__ == '__main__': main()
