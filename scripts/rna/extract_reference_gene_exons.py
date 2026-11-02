#!/usr/bin/env python3
"""Extract merged reference exon intervals for a single GTF gene.

Input GTF coordinates are 1-based inclusive; output BED is 0-based half-open.
The merged exons are the union of exons across reference transcript models,
not a single isoform's exon structure.
"""
import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

GENE_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
GENE_NAME = re.compile(r'(?:^|;\s*)gene_name\s+"([^"]+)"')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference-gtf', required=True)
    ap.add_argument('--gene', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    opener = gzip.open if args.reference_gtf.endswith('.gz') else open
    intervals = defaultdict(list)
    gene_ids = set()
    found = 0
    with opener(args.reference_gtf, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t', 8)
            if len(fields) != 9 or fields[0] != args.chrom or fields[2] != 'exon':
                continue
            attributes = fields[8]
            gene_id = GENE_ID.search(attributes)
            gene_name = GENE_NAME.search(attributes)
            id_value = gene_id.group(1) if gene_id else None
            name_value = gene_name.group(1) if gene_name else None
            if args.gene not in (id_value, name_value):
                continue
            gene_ids.add(id_value)
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            if start0 >= end0:
                raise ValueError(f'Invalid exon coordinates: {line[:100]}')
            intervals[fields[6]].append((start0, end0))
            found += 1

    if not found:
        raise SystemExit(f'No exons found for {args.gene} on {args.chrom}')
    if len(gene_ids) != 1 or len(intervals) != 1:
        raise SystemExit(f'Ambiguous gene identifiers or strands: {gene_ids}, {list(intervals)}')

    strand, exons = next(iter(intervals.items()))
    merged = []
    for start, end in sorted(exons):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as output:
        output.write('chrom\tstart0\tend0\tstrand\n')
        for start, end in merged:
            output.write(f'{args.chrom}\t{start}\t{end}\t{strand}\n')
    print(f'{args.gene}: {found} GTF exon records -> {len(merged)} union-exon blocks, {strand}')
    print(f'BED output: {path}')

if __name__ == '__main__':
    main()
