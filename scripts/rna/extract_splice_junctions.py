#!/usr/bin/env python3
"""Extract directly observed splice junctions and consecutive internal junction chains from SAM.

Use: samtools view -F 2308 sample.bam chr:START-END | python3 extract_splice_junctions.py --sample S --outdir DIR
Junction coordinates: 0-based half-open intron intervals; strand is not inferred.
A chain requires >=2 consecutive junctions passing the anchor filter. A rejected
junction breaks the chain: we never connect supported junctions across it.
Internal segments are *read-supported partial junction chains*, not complete isoforms.
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
REF_CONSUMING = frozenset('MDN=X')
MATCH = frozenset('M=X')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--min-mapq', type=int, default=20)
    p.add_argument('--min-anchor', type=int, default=8)
    args = p.parse_args()
    if args.min_anchor < 1 or args.min_mapq < 0:
        p.error('--min-anchor must be >=1 and --min-mapq must be >=0')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    chains = Counter()
    complete_chain_support = Counter()
    segment_chain_support = Counter()
    alignments_seen = eligible_alignments = spliced_alignments = 0
    rejected_junction_observations = reads_with_rejected_junctions = 0
    reads_with_multi_junction_chains = 0

    for line in sys.stdin:
        if not line or line[0] == '@':
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            continue
        alignments_seen += 1
        flag = int(fields[1])
        if flag & 0x904 or fields[2] == '*' or int(fields[4]) < args.min_mapq:
            continue
        chrom, cigar = fields[2], fields[5]
        if cigar == '*':
            continue
        pieces = CIGAR.findall(cigar)
        if ''.join(n + op for n, op in pieces) != cigar:
            raise ValueError(f'Unexpected CIGAR {cigar!r} for {fields[0]}')
        eligible_alignments += 1
        position0 = int(fields[3]) - 1
        accepted = []
        segments = []
        consecutive = []
        had_n = False
        had_rejection = False
        for i, (num, op) in enumerate(pieces):
            length = int(num)
            if op == 'N':
                had_n = True
                # Conservative: contiguous matched bases immediately beside each splice.
                left = int(pieces[i - 1][0]) if i > 0 and pieces[i - 1][1] in MATCH else 0
                right = int(pieces[i + 1][0]) if i + 1 < len(pieces) and pieces[i + 1][1] in MATCH else 0
                if min(left, right) >= args.min_anchor:
                    junc = (chrom, position0, position0 + length)
                    accepted.append(junc)
                    consecutive.append((position0, position0 + length))
                else:
                    had_rejection = True
                    rejected_junction_observations += 1
                    if len(consecutive) >= 2:
                        segments.append(tuple(consecutive))
                    consecutive = []
            if op in REF_CONSUMING:
                position0 += length
        if len(consecutive) >= 2:
            segments.append(tuple(consecutive))
        if had_n:
            spliced_alignments += 1
        if had_rejection:
            reads_with_rejected_junctions += 1
        for junc in set(accepted):
            counts[junc] += 1
        if segments:
            reads_with_multi_junction_chains += 1
            for sites in set(segments):  # one observation per chain per read
                key = (chrom, sites)
                chains[key] += 1
                if not had_rejection:
                    complete_chain_support[key] += 1
                else:
                    segment_chain_support[key] += 1

    with (outdir / 'junction_counts.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id', 'chrom', 'intron_start0', 'intron_end0', 'read_support'])
        for (chrom, start, end), support in sorted(counts.items()):
            w.writerow([args.sample, chrom, start, end, support])

    with (outdir / 'junction_chains.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id', 'chrom', 'junction_chain_start0_end0', 'n_junctions',
                    'read_support', 'all_junctions_pass_read_support', 'internal_segment_support'])
        for (chrom, sites), support in sorted(chains.items()):
            chain_str = ','.join(f'{start}-{end}' for start, end in sites)
            key = (chrom, sites)
            w.writerow([args.sample, chrom, chain_str, len(sites), support,
                        complete_chain_support[key], segment_chain_support[key]])

    with (outdir / 'extraction_summary.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample_id', 'alignments_seen', 'eligible_primary_alignments',
                    'spliced_primary_alignments', 'distinct_junctions', 'distinct_chains',
                    'min_mapq', 'min_anchor', 'total_accepted_junction_observations',
                    'rejected_junction_observations', 'reads_with_rejected_junctions',
                    'reads_with_multi_junction_chains'])
        w.writerow([args.sample, alignments_seen, eligible_alignments,
                    spliced_alignments, len(counts), len(chains), args.min_mapq,
                    args.min_anchor, sum(counts.values()), rejected_junction_observations,
                    reads_with_rejected_junctions, reads_with_multi_junction_chains])
    print(f'{args.sample}: {eligible_alignments} eligible alignments, {len(counts)} distinct junctions, '
          f'{sum(counts.values())} junction observations, {len(chains)} multi-junction chains; '
          f'{reads_with_rejected_junctions} reads have low-anchor junction(s)', file=sys.stderr)


if __name__ == '__main__':
    main()
