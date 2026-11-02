#!/usr/bin/env python3
"""Describe relative splice-junction support within shared donor/acceptor sites.

Denominators are sums of *observed junction-spanning alignment support* at a
shared splice site, not numbers of all locally informative reads. These
fractions are descriptive, not event-level PSI or statistical evidence of
differential splicing. Stranded gene annotation is required to assign donor
versus acceptor sites.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

REQUIRED = (
    'junction_id', 'gene_label', 'chrom', 'intron_start0', 'intron_end0',
    'case_sample', 'control_sample', 'case_read_support',
    'control_read_support', 'reference_strand',
)
SIDES = ('donor', 'acceptor')
SAMPLES = ('case', 'control')


def read_rows(path):
    with open(path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('Input TSV is empty or lacks a header')
        missing = sorted(set(REQUIRED) - set(reader.fieldnames))
        if missing:
            raise ValueError('Missing required columns: ' + ', '.join(missing))
        rows = list(reader)
    if not rows:
        raise ValueError('Input TSV has no junction rows')
    seen = set()
    pairs = set()
    for i, row in enumerate(rows, 2):
        try:
            start, end = int(row['intron_start0']), int(row['intron_end0'])
            case, control = int(row['case_read_support']), int(row['control_read_support'])
        except (ValueError, TypeError) as exc:
            raise ValueError(f'Line {i}: invalid coordinates or read support') from exc
        if start < 0 or end <= start or case < 0 or control < 0:
            raise ValueError(f'Line {i}: negative counts or invalid half-open intron')
        strand = row['reference_strand'].strip()
        if strand not in ('+', '-', '.', ''):
            raise ValueError(f'Line {i}: unsupported reference strand {strand!r}')
        key = (row['gene_label'], row['chrom'], start, end, strand)
        if key in seen:
            raise ValueError(f'Line {i}: repeated junction {key}; would double-count support')
        seen.add(key)
        pairs.add((row['case_sample'], row['control_sample']))
        row['_start'], row['_end'] = start, end
        row['_case'], row['_control'] = case, control
        row['_strand'] = strand
    if len(pairs) != 1:
        raise ValueError('Input mixes comparison sample pairs: ' + repr(sorted(pairs)))
    return reader.fieldnames, rows, next(iter(pairs))


def site_key(row, side):
    strand = row['_strand']
    if strand not in ('+', '-'):
        return None
    # Splice-site coordinates are *0-based boundaries* of the half-open intron.
    coordinate = (row['_start'] if strand == '+' else row['_end']) if side == 'donor' \
        else (row['_end'] if strand == '+' else row['_start'])
    return (row['gene_label'], row['chrom'], strand, side, coordinate)


def fmt_ratio(numerator, denominator):
    return f'{numerator / denominator:.6f}' if denominator else ''


def site_status(entry, sample, min_support):
    if entry is None:
        return 'unknown_strand'
    if entry['totals'][sample] == 0:
        return 'no_coverage'
    if entry['totals'][sample] < min_support:
        return 'low_site_support'
    if len(entry['observed_union']) < 2:
        return 'no_alternative_observed'
    return 'descriptive_comparable'


def process(input_path, outdir, prefix, min_support):
    fields, rows, (case_name, control_name) = read_rows(input_path)
    site_groups = defaultdict(lambda: {
        'totals': {'case': 0, 'control': 0},
        'observed': {'case': set(), 'control': set()},
        'observed_union': set(),
    })

    for row in rows:
        for side in SIDES:
            key = site_key(row, side)
            if key is None:
                continue
            group = site_groups[key]
            for sample in SAMPLES:
                support = row[f'_{sample}']
                group['totals'][sample] += support
                if support:
                    group['observed'][sample].add(row['junction_id'])
                    group['observed_union'].add(row['junction_id'])

    extra_cols = []
    for side in SIDES:
        extra_cols += [f'{side}_site0', f'{side}_pooled_observed_junctions']
        for sample in SAMPLES:
            extra_cols += [
                f'{sample}_{side}_total_support',
                f'{sample}_{side}_observed_junctions',
                f'{sample}_{side}_usage',
                f'{sample}_{side}_status',
            ]
        extra_cols += [f'{side}_delta_case_minus_control']

    output_rows = []
    for row in rows:
        out = {field: row[field] for field in fields}
        for side in SIDES:
            key = site_key(row, side)
            entry = site_groups.get(key) if key is not None else None
            out[f'{side}_site0'] = key[-1] if key else ''
            out[f'{side}_pooled_observed_junctions'] = len(entry['observed_union']) if entry else ''
            for sample in SAMPLES:
                total = entry['totals'][sample] if entry else 0
                out[f'{sample}_{side}_total_support'] = total if entry else ''
                out[f'{sample}_{side}_observed_junctions'] = len(entry['observed'][sample]) if entry else ''
                out[f'{sample}_{side}_usage'] = fmt_ratio(row[f'_{sample}'], total) if entry else ''
                out[f'{sample}_{side}_status'] = site_status(entry, sample, min_support)
            valid = entry and all(
                site_status(entry, sample, min_support) == 'descriptive_comparable'
                for sample in SAMPLES
            )
            out[f'{side}_delta_case_minus_control'] = (
                f"{row['_case'] / entry['totals']['case'] - row['_control'] / entry['totals']['control']:.6f}"
                if valid else ''
            )
        output_rows.append(out)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    usage_path = outdir / f'{prefix}.junction_relative_usage.tsv'
    with usage_path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fields + extra_cols)
        writer.writeheader()
        writer.writerows(output_rows)

    sites_path = outdir / f'{prefix}.junction_usage_sites.tsv'
    site_fields = [
        'gene_label', 'chrom', 'reference_strand', 'site_type', 'site0',
        'sample_role', 'sample_id', 'total_support',
        'observed_junctions', 'pooled_observed_junctions', 'site_status',
    ]
    with sites_path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=site_fields)
        writer.writeheader()
        for (gene, chrom, strand, side, position), entry in sorted(site_groups.items()):
            for role, sample_id in (('case', case_name), ('control', control_name)):
                writer.writerow({
                    'gene_label': gene, 'chrom': chrom,
                    'reference_strand': strand, 'site_type': side, 'site0': position,
                    'sample_role': role, 'sample_id': sample_id,
                    'total_support': entry['totals'][role],
                    'observed_junctions': len(entry['observed'][role]),
                    'pooled_observed_junctions': len(entry['observed_union']),
                    'site_status': site_status(entry, role, min_support),
                })

    summary_path = outdir / f'{prefix}.junction_usage_summary.tsv'
    summary = {
        'case_sample': case_name,
        'control_sample': control_name,
        'min_site_support': min_support,
        'total_junctions': len(rows),
        'unknown_strand_junctions': sum(row['_strand'] not in ('+', '-') for row in rows),
    }
    for side in SIDES:
        groups = [g for key, g in site_groups.items() if key[3] == side]
        summary[f'{side}_sites'] = len(groups)
        summary[f'{side}_sites_with_alternative'] = sum(
            len(g['observed_union']) >= 2 for g in groups
        )
        summary[f'{side}_sites_comparable_at_threshold'] = sum(
            len(g['observed_union']) >= 2 and
            all(g['totals'][role] >= min_support for role in SAMPLES)
            for g in groups
        )
    with summary_path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)

    print(f'Processed {len(rows)} junctions from {case_name} versus {control_name}')
    print(f'Wrote {usage_path}, {sites_path}, {summary_path}')
    print('Interpretation: descriptive proportions of observed splice-junction support; not PSI.')
    return usage_path, sites_path, summary_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='Ensembl-annotated junction comparison TSV')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--prefix', required=True, help='Output comparison name')
    parser.add_argument('--min-site-support', type=int, default=10)
    args = parser.parse_args()
    if args.min_site_support < 1:
        parser.error('--min-site-support must be >= 1')
    if Path(args.prefix).name != args.prefix or args.prefix in ('.', '..'):
        parser.error('--prefix must be a plain filename component')
    process(args.input, args.outdir, args.prefix, args.min_site_support)


if __name__ == '__main__':
    main()
