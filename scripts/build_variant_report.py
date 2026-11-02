from pathlib import Path
import argparse
import gzip
import re
import csv


# Read command-line arguments passed from Nextflow.
# These arguments tell the script which annotated VCF files to use,
# which target gene to report, whether a known variant was supplied,
# and where to write the final CSV.
def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a target-gene DNA variant report from annotated VCF files."
    )

    # VEP-annotated small variant VCF
    parser.add_argument("--small", required=True)

    # VEP-annotated structural variant VCF
    parser.add_argument("--sv", required=True)

    # VEP-annotated CNV VCF from either QDNAseq or Spectre
    parser.add_argument("--cnv", required=True)

    # Gene to retain in the final report, e.g. NEB
    parser.add_argument("--gene", required=True)
    parser.add_argument("--gene_bed", required=True, help="BED file containing the target gene coordinates")

    # Optional VCF containing previously identified WGS variant(s).
    # Used for matching and flagging, not filtering detected variants.
    parser.add_argument(
    	"--known_variant",
    	default="",
    	help="Path to a VCF containing previously identified variants"
    )

    parser.add_argument("--output", required=True)

    return parser.parse_args()


# Open either a compressed or uncompressed VCF.
def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


# Read the VEP CSQ header and return the names of the annotation fields.
def get_csq_fields(vcf_path):
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("##INFO=<ID=CSQ"):

                match = re.search(r'Format: ([^">]+)', line)

                if match:
                    return match.group(1).split("|")

    # If no CSQ header is found, return an empty list.
    return []


# Convert a VCF INFO field into a dictionary.
def parse_info(info_string):
    info = {}

    for item in info_string.split(";"):

        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value

        elif item:
            info[item] = True

    return info


# Parse VEP CSQ annotations into one dictionary per transcript annotation.
def parse_csq(csq_string, csq_fields):
    annotations = []

    if not csq_string or not csq_fields:
        return annotations

    for entry in csq_string.split(","):
        values = entry.split("|")

        annotation = {
            field: values[i] if i < len(values) else ""
            for i, field in enumerate(csq_fields)
        }

        annotations.append(annotation)

    return annotations


# Read VCF records and retain annotations for the target gene.
def read_target_gene_variants(vcf_path, target_gene, variant_type, gene_chrom, gene_start, gene_end):
    rows = []

    csq_fields = get_csq_fields(vcf_path)

    with open_vcf(vcf_path) as vcf:
        for line in vcf:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = parse_info(fields[7])

            # Extract structural variant and copy-number information.
            # These fields may be absent depending on the variant caller.
            variant_end = info.get("END", "")
            sv_type = info.get("SVTYPE", "")
            sv_length = info.get("SVLEN", "")
            copy_number = info.get("CN", "")
            log2_count = info.get("LOG2CNT", "")

            # Additional structural variant evidence.
            support = info.get("SUPPORT", "")
            allele_fraction = info.get("AF", "")
            breakpoint_stdev_pos = info.get("STDEV_POS", "")
            breakpoint_precision = (
                "imprecise" if "IMPRECISE" in info else
                "not_flagged_imprecise" if variant_type == "structural_variant"
                else ""
            )

            sample_data = {}

            # Extract the genotype (GT) from the first sample column,
            # when genotype information is present in the VCF.
            genotype = ""

            if len(fields) > 9:
                format_keys = fields[8].split(":")
                sample_values = fields[9].split(":")

                sample_data = {
                    key: sample_values[i] if i < len(sample_values) else ""
                    for i, key in enumerate(format_keys)
                }

                genotype = sample_data.get("GT", "")

            reference_reads = sample_data.get("DR", "")
            variant_reads = sample_data.get("DV", "")
            genotype_quality = sample_data.get("GQ", "")

            # Convert the raw VCF genotype into a readable allele-state label.
            allele_state = ""

            if genotype in {"0/1", "1/0", "0|1", "1|0"}:
                allele_state = "heterozygous"

            elif genotype in {"1/1", "1|1"}:
                allele_state = "homozygous_alt"

            elif genotype in {"0/0", "0|0"}:
                allele_state = "homozygous_reference"

            elif genotype:
                allele_state = "other"

            # Parse the VEP CSQ annotations for this variant.
            annotations = parse_csq(
                info.get("CSQ", ""),
                csq_fields
            )

            # Keep gene-overlapping variants even without a matching
            # VEP gene annotation, across all three variant classes.
            matching_annotations = [
                annotation for annotation in annotations
                if annotation.get("SYMBOL") == target_gene
            ]

            if not matching_annotations:
                if not overlaps_gene(
                    chrom, pos, variant_end,
                    gene_chrom, gene_start, gene_end
                ):
                    continue
                matching_annotations = [{"SYMBOL": target_gene}]

            # Retain transcript annotations belonging to the requested gene.
            for annotation in matching_annotations:

                rows.append({
                    "variant_type": variant_type,
                    "chromosome": chrom,
                    "position": pos,
                    "id": variant_id,
                    "reference": ref,
                    "alternate": alt,
                    "gene": annotation.get("SYMBOL", ""),
                    "transcript": annotation.get("Feature", ""),
                    "consequence": annotation.get("Consequence", ""),
                    "canonical": annotation.get("CANONICAL", ""),
                    "mane_select": annotation.get("MANE_SELECT", ""),
                    "gnomade_af": annotation.get("gnomADe_AF", ""),
                    "sift": annotation.get("SIFT", ""),
                    "polyphen": annotation.get("PolyPhen", ""),
                    "genotype": genotype,
                    "allele_state": allele_state,
                    "variant_end": variant_end,
                    "sv_type": sv_type,
                    "sv_length": sv_length,
                    "copy_number": copy_number,
                    "log2_count": log2_count,
                    "rsid": ",".join(
                        dict.fromkeys(
                            identifier
                            for identifier in re.split(
                                r"[&,]", annotation.get("Existing_variation", "")
                            )
                            if re.fullmatch(r"rs[0-9]+", identifier)
                        )
                    ),
                    "support": support,
                    "reference_reads": reference_reads,
                    "variant_reads": variant_reads,
                    "allele_fraction": allele_fraction,
                    "genotype_quality": genotype_quality,
                    "breakpoint_stdev_pos": breakpoint_stdev_pos,
                    "breakpoint_precision": breakpoint_precision,
                })

    return rows

# Read previously identified WGS variants from a VCF.
def load_known_variants(vcf_path):
    known = set()

    if not vcf_path:
        return known

    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = (fields[0][3:] if fields[0].startswith("chr") else fields[0])
            pos = fields[1]
            ref = fields[3]

            for alt in fields[4].split(","):
                known.add((chrom, pos, ref, alt))

    return known


def read_gene_bed(bed_path, target_gene):
    with open(bed_path) as bed:
        for line in bed:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split()
            if len(fields) >= 4 and fields[3] == target_gene:
                return (fields[0][3:] if fields[0].startswith("chr") else fields[0]), int(fields[1]), int(fields[2])

    raise ValueError(f"Gene {target_gene} not found in {bed_path}")



def overlaps_gene(chrom, pos, variant_end, gene_chrom, gene_start, gene_end):
    if (chrom[3:] if chrom.startswith("chr") else chrom) != gene_chrom:
        return False

    # Convert VCF 1-based coordinates to BED 0-based coordinates.
    variant_start = int(pos) - 1
    variant_end = int(variant_end) if variant_end else int(pos)

    return variant_start < gene_end and variant_end > gene_start



def main():
    args = parse_args()
    gene_chrom, gene_start, gene_end = read_gene_bed(
        args.gene_bed, args.gene
    )

    all_rows = []
    for vcf, variant_type in [
        (args.small, "small_variant"),
        (args.sv, "structural_variant"),
        (args.cnv, "cnv"),
    ]:
        all_rows.extend(
            read_target_gene_variants(
                vcf, args.gene, variant_type,
                gene_chrom, gene_start, gene_end
            )
        )

    known_variants = load_known_variants(args.known_variant)

    def variant_key(row):
        chrom = str(row["chromosome"]).removeprefix("chr")
        return (
            row["variant_type"],
            chrom,
            str(row["position"]),
            str(row["reference"]),
            str(row["alternate"]),
            str(row["variant_end"]),
        )

    def is_known(row):
        chrom = str(row["chromosome"]).removeprefix("chr")
        return (
            chrom,
            str(row["position"]),
            str(row["reference"]),
            str(row["alternate"]),
        ) in known_variants

    def annotation_rank(row):
        # Prefer MANE Select, then canonical transcripts.
        return (
            0 if row["mane_select"] else 1,
            0 if str(row["canonical"]).upper() == "YES" else 1,
            str(row["transcript"]),
            str(row["consequence"]),
        )

    transcript_fields = [
        "variant_type", "chromosome", "position", "id",
        "reference", "alternate", "gene", "transcript",
        "consequence", "canonical", "mane_select",
        "gnomade_af", "sift", "polyphen",
        "genotype", "allele_state", "variant_end",
        "sv_type", "sv_length", "copy_number", "log2_count",
        "known_variant", "rsid", "support",
        "reference_reads", "variant_reads", "allele_fraction",
        "genotype_quality", "breakpoint_stdev_pos",
        "breakpoint_precision"
    ]

    for row in all_rows:
        row["known_variant"] = "yes" if is_known(row) else "no"

    # Preserve every transcript-level annotation.
    supplementary_path = (
        Path(args.output).with_name(
            Path(args.output).stem + "_transcripts.csv"
        )
    )

    with supplementary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=transcript_fields
        )
        writer.writeheader()
        writer.writerows(all_rows)

    # Select one representative transcript per distinct variant.
    variants = {}
    for row in all_rows:
        key = variant_key(row)
        if key not in variants or (
            annotation_rank(row) < annotation_rank(variants[key])
        ):
            variants[key] = row

    main_fields = [
        field for field in transcript_fields
        if field not in ("transcript", "canonical", "mane_select")
    ] + ["representative_transcript", "annotation_selection"]

    main_rows = []
    for row in variants.values():
        result = {
            field: row.get(field, "")
            for field in main_fields
        }
        result["representative_transcript"] = row["transcript"]
        result["annotation_selection"] = (
            "MANE_Select" if row["mane_select"]
            else "canonical" if str(row["canonical"]).upper() == "YES"
            else "other"
        )
        main_rows.append(result)

    main_rows.sort(
        key=lambda row: (
            row["variant_type"],
            row["chromosome"],
            int(row["position"]),
            row["reference"],
            row["alternate"],
        )
    )

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=main_fields)
        writer.writeheader()
        writer.writerows(main_rows)

    print(
        f"Reported {len(main_rows)} variants and "
        f"{len(all_rows)} transcript annotations."
    )


# Run main() only when this script is executed directly.
if __name__ == "__main__":
    main()
