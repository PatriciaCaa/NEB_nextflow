# QDNAseq CNV analysis for Oxford Nanopore adaptive-sampling DNA data.
# Input: aligned BAM file.

library(QDNAseq)

args <- commandArgs(trailingOnly = TRUE)

bam_file <- args[1]
output_prefix <- args[2]
bin_size_kbp <- as.numeric(args[3])
genome_build <- args[4]

if (!genome_build %in% c("hg19", "hg38")) {
    stop("QDNAseq currently expects genome_build to be 'hg19' or 'hg38'.")
}

# Load bin annotations for the selected genome build and bin size
bins <- getBinAnnotations(
    binSize = bin_size_kbp,
    genome = genome_build
)

# Count aligned reads in each genomic bin
read_counts <- binReadCounts(
    bins = bins,
    bamfiles = bam_file
)