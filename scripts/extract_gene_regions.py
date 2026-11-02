
#!/usr/bin/env python3

import argparse
from pathlib import Path
import gzip
import re
import sys


def parse_attributes(field):
    return dict(re.findall(r'(\w+)\s+"([^"]*)"', field))


def load_reference_lengths(fai_path):
    lengths = {}
    with open(fai_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def main():
    parser = argparse.ArgumentParser(
        description="Extract a gene's genomic span from a GTF."
    )
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    reference_lengths = load_reference_lengths(args.fai)
    matches = []

    opener = gzip.open if args.gtf.endswith(".gz") else open

    with opener(args.gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "gene":
                continue

            attributes = parse_attributes(fields[8])

            if args.gene not in (
                attributes.get("gene_name"),
                attributes.get("gene_id"),
            ):
                continue

            chrom = fields[0]

            if chrom not in reference_lengths:
                prefixed = f"chr{chrom}"
                if prefixed in reference_lengths:
                    chrom = prefixed
                else:
                    raise ValueError(
                        f"Chromosome {fields[0]} is absent "
                        "from the reference FASTA index."
                    )

            start = int(fields[3]) - 1
            end = int(fields[4])

            if start < 0 or end > reference_lengths[chrom] or start >= end:
                raise ValueError(
                    f"Invalid coordinates for {args.gene}: "
                    f"{chrom}:{start}-{end}"
                )

            matches.append(
                (chrom, start, end, attributes.get("gene_id"), fields[6])
            )

    if len(matches) != 1:
        sys.exit(
            f"Expected exactly one gene matching '{args.gene}', "
            f"but found {len(matches)}."
        )

    chrom, start, end, gene_id, strand = matches[0]

    with open(args.output, "w") as handle:
        handle.write(
            f"{chrom}\t{start}\t{end}\t{args.gene}\t0\t{strand}\n"
        )

    print(
        f"Extracted {args.gene} ({gene_id}): "
        f"{chrom}:{start}-{end} [BED coordinates]"
    )

def write_exon_bed(gtf_path, target_gene, output_path):
    """Write merged exon intervals across all transcripts of a gene."""
    import gzip

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    intervals = []

    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            attributes = fields[8]
            gene_name = None
            gene_id = None

            for item in attributes.split(";"):
                item = item.strip()
                if item.startswith('gene_name "'):
                    gene_name = item.split('"')[1]
                elif item.startswith('gene_id "'):
                    gene_id = item.split('"')[1]

            if target_gene not in (gene_name, gene_id):
                continue

            start = int(fields[3]) - 1
            end = int(fields[4])
            intervals.append((fields[0], start, end))

    if not intervals:
        raise ValueError(f"No annotated exons found for {target_gene}")

    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = []

    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            previous = merged[-1]
            merged[-1] = (chrom, previous[1], max(previous[2], end))
        else:
            merged.append((chrom, start, end))

    with open(output_path, "w") as handle:
        for chrom, start, end in merged:
            handle.write(f"{chrom}\t{start}\t{end}\n")

if __name__ == "__main__":
    main()

    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True)
    args, _ = parser.parse_known_args()

    exon_output = args.output[:-4] + ".exons.bed"
    write_exon_bed(args.gtf, args.gene, exon_output)

    # Match exon chromosome names to the reference FASTA.
    with open(args.output) as handle:
        reference_chrom = handle.readline().split("\t")[0]

    exon_path = Path(exon_output)
    lines = exon_path.read_text().splitlines()

    with exon_path.open("w") as handle:
        for line in lines:
            fields = line.split("\t")
            fields[0] = reference_chrom
            handle.write("\t".join(fields) + "\n")


