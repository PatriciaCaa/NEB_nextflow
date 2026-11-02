#!/usr/bin/env python3
"""Reference-guided gene-region selection and structural summary for transcript GTFs.

Selection is evidence-based, NOT a validated transcript-to-gene assignment.
An overlap with the gene locus does not prove transcript identity or novelty.
"""
import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

ATTR_RE = re.compile(r'([\w.:-]+)\s+"([^"]*)"')
FIELDS = [
    "sample_id", "chrom", "gene_id", "gene_name", "transcript_id",
    "transcript_start", "transcript_end", "genomic_span", "spliced_length",
    "exon_count", "total_junctions", "coverage", "FPKM", "TPM",
    "strand", "target_gene_id", "target_gene_name", "assignment_evidence",
    "exonic_overlap_bp", "assembler",
]


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            if fields[2] not in ('gene', 'transcript', 'exon'):
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            yield fields, start, end, dict(ATTR_RE.findall(fields[8]))


def merge_intervals(intervals):
    merged = []
    for a, b in sorted(intervals):
        if merged and a <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def overlapping_bases(query, reference):
    # Both sequences must be sorted, non-overlapping genomic intervals.
    i = j = total = 0
    while i < len(query) and j < len(reference):
        a, b = query[i]
        c, d = reference[j]
        total += max(0, min(b, d) - max(a, c) + 1)
        if b < d:
            i += 1
        else:
            j += 1
    return total


def get_target(reference_gtf, gene):
    candidates = defaultdict(lambda: {'exons': [], 'spans': [], 'names': set()})
    for f, start, end, a in lines(reference_gtf):
        if gene not in (a.get('gene_id'), a.get('gene_name')):
            continue
        gid = a.get('gene_id')
        if not gid:
            continue
        key = (gid, f[0], f[6])
        data = candidates[key]
        data['names'].add(a.get('gene_name') or gene)
        data['spans'].append((start, end))
        if f[2] == 'exon':
            data['exons'].append((start, end))

    if not candidates:
        raise SystemExit(f"Target {gene!r} not found in reference GTF")
    unique_ids = {k[0] for k in candidates}
    if len(unique_ids) != 1:
        raise SystemExit(f"Ambiguous gene name {gene!r}; matches gene IDs: {sorted(unique_ids)}. Use a stable gene ID.")

    targets = {}
    for (gid, chrom, strand), data in candidates.items():
        spans = data['spans']
        targets[(chrom, strand)] = {
            'gene_id': gid,
            'gene_name': sorted(data['names'])[0],
            'start': min(a for a, _ in spans),
            'end': max(b for _, b in spans),
            'exons': merge_intervals(data['exons']),
        }
    return targets


def summarise(args):
    targets = get_target(args.reference_gtf, args.gene)
    records = {}
    for f, start, end, a in lines(args.assembly_gtf):
        if f[2] not in ('transcript', 'exon'):
            continue
        tid = a.get('transcript_id')
        if not tid:
            continue
        sample = a.get('sample_id') or args.sample_id
        if not sample:
            raise SystemExit('GTF lacks sample_id: supply --sample-id for a single-sample GTF')
        key = (sample, tid, f[0], f[6])
        record = records.setdefault(key, {
            'sample_id': sample, 'chrom': f[0], 'strand': f[6],
            'transcript_id': tid, 'gene_id': '', 'gene_name': '',
            'ref_gene_id': '', 'starts': [], 'ends': [], 'exons': [],
            'coverage': '', 'FPKM': '', 'TPM': '',
        })
        record['starts'].append(start)
        record['ends'].append(end)
        if f[2] == 'exon':
            record['exons'].append((start, end))
        for name in ('gene_id', 'coverage', 'FPKM', 'TPM'):
            attribute = 'cov' if name == 'coverage' else name
            if a.get(attribute) not in (None, ''):
                record[name] = a[attribute]
        if a.get('ref_gene_id'):
            record['ref_gene_id'] = a['ref_gene_id']
        if a.get('ref_gene_name') or a.get('gene_name'):
            record['gene_name'] = a.get('ref_gene_name') or a.get('gene_name')

    result = []
    for row in records.values():
        if not row['exons']:
            continue  # Junction and spliced-length metrics need exon features.
        chrom, strand = row['chrom'], row['strand']
        target = targets.get((chrom, strand))
        if target is None:
            continue
        start, end = min(row['starts']), max(row['ends'])
        if end < target['start'] or start > target['end']:
            continue
        exons = merge_intervals(row['exons'])
        bp = overlapping_bases(exons, target['exons'])
        reference_label_match = (
            row['ref_gene_id'] == target['gene_id'] or
            row['gene_name'] == target['gene_name']
        )
        if bp:
            evidence = 'reference_label_and_exon_overlap' if reference_label_match else 'exon_overlap_candidate'
        elif reference_label_match:
            evidence = 'reference_label_only'
        elif args.include_span_only:
            evidence = 'span_only_candidate'
        else:
            continue
        row.update({
            'transcript_start': start, 'transcript_end': end,
            'genomic_span': end - start + 1,
            'spliced_length': sum(b - a + 1 for a, b in exons),
            'exon_count': len(exons), 'total_junctions': max(len(exons) - 1, 0),
            'target_gene_id': target['gene_id'],
            'target_gene_name': target['gene_name'],
            'assignment_evidence': evidence,
            'exonic_overlap_bp': bp,
            'assembler': args.assembler,
        })
        result.append(row)

    result.sort(key=lambda r: (r['sample_id'], r['chrom'], r['transcript_start'], r['transcript_id']))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result)
    print(f'Wrote {len(result)} transcript records to {args.output}')
    print('Evidence categories:', dict((k, sum(r['assignment_evidence'] == k for r in result))
                                       for k in sorted({r['assignment_evidence'] for r in result})))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference-gtf', required=True)
    ap.add_argument('--assembly-gtf', required=True)
    ap.add_argument('--gene', required=True, help='Reference gene symbol or stable gene_id')
    ap.add_argument('--sample-id', help='Required if GTF lacks sample_id')
    ap.add_argument('--assembler', choices=['stringtie', 'flair'], required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--include-span-only', action='store_true',
                    help='Include overlapping genomic spans with no reference-exon overlap, labelled candidates')
    summarise(ap.parse_args())


if __name__ == '__main__':
    main()
