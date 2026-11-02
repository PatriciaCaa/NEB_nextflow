#!/usr/bin/env python3
"""Compare two read-supported junction-count tables and create annotation look-up files.

Input: junction_counts.tsv files from extract_splice_junctions.py.
Input/output junction coordinates: 0-based, half-open genomic intron intervals.
One case + one control is a descriptive TECHNICAL comparison, not a differential splicing test.

Outputs:
  junction_comparison.tsv      support count per sample; detected in both/one/none
  isopedia_splice_queries.bed  four-column BED, *candidate* batch splice queries
  introverse_lookup.tsv        coordinate lookup table, not an IntroVerse upload API

No isoforms, strand, splice-event PSI or annotation novelty is inferred.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

COLUMNS = {'sample_id', 'chrom', 'intron_start0', 'intron_end0', 'read_support'}


def load_counts(path: Path):
    counts = Counter()
    sample_ids = set()
    with path.open(newline='') as inp:
        reader = csv.DictReader(inp, delimiter='\t')
        if not reader.fieldnames or not COLUMNS.issubset(reader.fieldnames):
            raise ValueError(f'{path}: missing columns {sorted(COLUMNS - set(reader.fieldnames or []))}')
        for row in reader:
            sample_ids.add(row['sample_id'])
            chrom = row['chrom']
            start = int(row['intron_start0'])
            end = int(row['intron_end0'])
            support = int(row['read_support'])
            if not chrom or start < 0 or end <= start or support < 0:
                raise ValueError(f'{path}: invalid junction row: {row}')
            key = (chrom, start, end)
            if key in counts:
                raise ValueError(f'{path}: duplicated junction {key}')
            counts[key] = support
    if len(sample_ids) > 1:
        raise ValueError(f'{path}: expected one sample, found {sorted(sample_ids)}')
    return counts, next(iter(sample_ids), path.stem)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', type=Path, required=True, help='Case junction_counts.tsv')
    p.add_argument('--control', type=Path, required=True, help='Control junction_counts.tsv')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--gene', default='', help='Optional label for lookups, e.g. NEB')
    p.add_argument('--min-support', type=int, default=1,
                   help='Minimum reads in at least one sample to include a junction (default: 1)')
    args = p.parse_args()
    if args.min_support < 1:
        p.error('--min-support must be at least 1')
    case, case_id = load_counts(args.case)
    control, control_id = load_counts(args.control)
    if case_id == control_id:
        p.error(f'Case and control share sample ID {case_id!r}')
    args.outdir.mkdir(parents=True, exist_ok=True)
    keys = sorted((set(case) | set(control)), key=lambda k: (k[0], k[1], k[2]))
    comparison = args.outdir / 'junction_comparison.tsv'
    bed = args.outdir / 'isopedia_splice_queries.bed'
    introverse = args.outdir / 'introverse_lookup.tsv'
    summary = Counter()
    with comparison.open('w', newline='') as a, bed.open('w', newline='') as b, introverse.open('w', newline='') as c:
        w = csv.writer(a, delimiter='\t')
        i = csv.writer(c, delimiter='\t')
        w.writerow(['junction_id', 'gene_label', 'chrom', 'intron_start0', 'intron_end0',
                    'case_sample', 'control_sample', 'case_read_support',
                    'control_read_support', 'detection_class'])
        i.writerow(['junction_id', 'gene_label', 'chrom', 'intron_start0', 'intron_end0',
                    'intron_start1', 'intron_end1', 'coordinate_system',
                    'case_read_support', 'control_read_support', 'detection_class'])
        for chrom, start, end in keys:
            ca, co = case.get((chrom, start, end), 0), control.get((chrom, start, end), 0)
            if max(ca, co) < args.min_support:
                continue
            detection = 'shared' if ca > 0 and co > 0 else ('case_only' if ca > 0 else 'control_only')
            summary[detection] += 1
            j_id = f'{chrom}:{start}-{end}'
            w.writerow([j_id, args.gene, chrom, start, end, case_id, control_id, ca, co, detection])
            # BED 0-based half-open intron interval and stable junction identifier.
            b.write(f'{chrom}\t{start}\t{end}\t{j_id}\n')
            # The first deleted/intronic base is start0+1 in 1-based notation;
            # the last base is end0. Do not confuse these with exon boundary bases.
            i.writerow([j_id, args.gene, chrom, start, end, start+1, end,
                        '0-based half-open intron (also 1-based inclusive intron)', ca, co, detection])
    print(f'Case: {case_id}; control: {control_id}; junctions: {sum(summary.values())}; '
          f'shared={summary["shared"]}, case_only={summary["case_only"]}, control_only={summary["control_only"]}')
    print(f'Outputs: {args.outdir}')


if __name__ == '__main__':
    main()
