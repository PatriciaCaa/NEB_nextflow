#!/usr/bin/env python3
"""Annotate read-supported junctions against adjacent-exon junctions in a gene's reference GTF.

Inputs: junction_comparison.tsv from compare_splice_junctions.py; reference GTF.
Coordinates: chrom, intron_start0, intron_end0 are 0-based half-open.
Outputs: annotated_junction_comparison.tsv, annotation_summary.tsv.

'Unannotated' means absent from THIS GTF, not necessarily biologically novel.
Transcript identity cannot be inferred from a single junction.
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

ATTR_RE = re.compile(r'(?:^|;\s*)([A-Za-z_][\w.]*)\s+"([^"]*)"')

def attrs(value):
    return dict(ATTR_RE.findall(value))

def load_reference(gtf, gene):
    # Key transcript by gene and chromosome too, to avoid accidental merging.
    tx_exons = defaultdict(list)
    gene_loci = []
    matching_exons = 0
    with gtf.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] not in ('exon', 'gene'):
                continue
            if f'gene_name "{gene}"' not in fields[8]:
                continue
            at = attrs(fields[8])
            if at.get('gene_name') != gene:
                continue
            chrom, strand = fields[0], fields[6]
            start, end = int(fields[3]), int(fields[4])
            if fields[2] == 'gene':
                gene_loci.append((chrom, start-1, end, strand))
            else:
                tid = at.get('transcript_id')
                gid = at.get('gene_id', '')
                if tid and strand in ('+', '-'):
                    tx_exons[(gid, tid, chrom, strand)].append((start, end))
                    matching_exons += 1
    if not tx_exons:
        raise ValueError(f'No annotated exons for gene {gene!r} were found in {gtf}')
    if not gene_loci:
        # Reference may lack gene feature; derive a conservative envelope.
        locus_parts = defaultdict(list)
        for (gid, tid, chrom, strand), exon_list in tx_exons.items():
            locus_parts[(chrom, strand)].extend(exon_list)
        gene_loci = [(chrom, min(s for s,e in exons)-1, max(e for s,e in exons), strand)
                     for (chrom, strand), exons in locus_parts.items()]
    known_junctions = defaultdict(set)
    donors = set()
    acceptors = set()
    for (gid, tid, chrom, strand), exons in tx_exons.items():
        exon_list = sorted(set(exons))
        for (_left_s, left_e), (right_s, _right_e) in zip(exon_list, exon_list[1:]):
            start0, end0 = left_e, right_s-1
            if end0 <= start0:
                continue
            j = (chrom, start0, end0)
            known_junctions[j].add(tid)
            donor, acceptor = (start0, end0) if strand == '+' else (end0, start0)
            donors.add((chrom, strand, donor))
            acceptors.add((chrom, strand, acceptor))
    return gene_loci, known_junctions, donors, acceptors, len(tx_exons), matching_exons

def classify(chrom, start, end, loci, known, donors, acceptors):
    key = (chrom, start, end)
    # Reference chromosome and target-locus overlap must agree.
    overlaps = [strand for ch, ls, le, strand in loci if ch == chrom and start < le and end > ls]
    if not overlaps:
        return 'outside_reference_gene_locus', '.', False, False, []
    strands = sorted(set(overlaps))
    strand = strands[0] if len(strands) == 1 else '.'
    matches = sorted(known.get(key, ()))
    if matches:
        return 'annotated_junction', strand, True, True, matches
    # If ambiguous reference strand, don't incorrectly assign donor/acceptor novelty.
    if strand == '.':
        return 'unannotated_ambiguous_strand', strand, False, False, []
    donor_coord, acceptor_coord = (start, end) if strand == '+' else (end, start)
    d = (chrom, strand, donor_coord) in donors
    a = (chrom, strand, acceptor_coord) in acceptors
    if d and a:
        label = 'unannotated_combination_known_splice_sites'
    elif d:
        label = 'unannotated_acceptor_site'
    elif a:
        label = 'unannotated_donor_site'
    else:
        label = 'unannotated_both_splice_sites'
    return label, strand, d, a, []

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison', required=True, type=Path)
    p.add_argument('--reference-gtf', required=True, type=Path)
    p.add_argument('--gene', required=True)
    p.add_argument('--outdir', required=True, type=Path)
    args = p.parse_args()
    loci, known, donors, acceptors, n_tx, n_exons = load_reference(args.reference_gtf, args.gene)
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = Counter()
    out = args.outdir / 'annotated_junction_comparison.tsv'
    with args.comparison.open(newline='') as inp, out.open('w', newline='') as fh:
        reader = csv.DictReader(inp, delimiter='\t')
        required = {'chrom', 'intron_start0', 'intron_end0'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError('Comparison missing required columns: ' + ', '.join(sorted(required)))
        extra = ['reference_annotation_status', 'reference_strand',
                 'known_donor_site', 'known_acceptor_site',
                 'matching_reference_transcripts', 'reference_gtf']
        writer = csv.DictWriter(fh, fieldnames=list(reader.fieldnames)+extra, delimiter='\t')
        writer.writeheader()
        for row in reader:
            chrom = row['chrom']
            start, end = int(row['intron_start0']), int(row['intron_end0'])
            if start < 0 or end <= start:
                raise ValueError(f'Invalid 0-based junction coordinates: {chrom}:{start}-{end}')
            status, strand, d, a, tx = classify(chrom, start, end, loci, known, donors, acceptors)
            summary[status] += 1
            row.update(reference_annotation_status=status, reference_strand=strand,
                       known_donor_site=str(d).lower(), known_acceptor_site=str(a).lower(),
                       matching_reference_transcripts=','.join(tx),
                       reference_gtf=args.reference_gtf.name)
            writer.writerow(row)
    with (args.outdir/'annotation_summary.tsv').open('w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene', 'reference_gtf', 'annotated_transcripts',
                         'reference_exon_records', 'annotation_status', 'junctions'])
        for status, n in sorted(summary.items()):
            writer.writerow([args.gene, args.reference_gtf.name, n_tx, n_exons, status, n])
    print(f'{args.gene}: {sum(summary.values())} compared junctions; '
          f'{summary["annotated_junction"]} match a reference transcript; '
          f'{sum(n for s,n in summary.items() if s.startswith("unannotated"))} absent from this GTF.')
    print(f'Wrote {out} and annotation_summary.tsv')

if __name__ == '__main__':
    main()
