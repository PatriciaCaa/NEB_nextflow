process PYCHOPPER {

    tag "${sample}"
    label 'rna_pychopper'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.oriented.fastq"), emit: trimmed_reads
    tuple val(sample), path("${sample}.unclassified.fastq"), emit: unclassified
    tuple val(sample), path("${sample}.rescued.fastq"), emit: rescued
    tuple val(sample), path("${sample}.pychopper_stats.tsv"), emit: stats
    tuple val(sample), path("${sample}.pychopper_report.pdf"), emit: report

    script:
    def kit = params.rna_pychopper_kit ?: 'PCS109'

    """
    set -euo pipefail

    pychopper \
        -k ${kit} \
        -t ${task.cpus} \
        -r ${sample}.pychopper_report.pdf \
        -S ${sample}.pychopper_stats.tsv \
        -u ${sample}.unclassified.fastq \
        -w ${sample}.rescued.fastq \
        ${fastq} ${sample}.full_length.fastq

    # Combine classified full-length and rescued reads for downstream alignment.
    # Keep unclassified reads separate from the analysis input.

    cat ${sample}.full_length.fastq \
        ${sample}.rescued.fastq \
        > ${sample}.oriented.fastq
    """
}
