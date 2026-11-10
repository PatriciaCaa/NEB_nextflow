process FLAIR_TRANSCRIPTOME {

    tag "${sample}"
    label 'rna_flair'

    publishDir "${params.results_dir}/rna/flair",
        mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path reference_fasta
    path reference_gtf

    output:
    tuple val(sample),
          path("${sample}.flair.isoforms.gtf"),
          emit: isoform_gtf

    tuple val(sample),
          path("${sample}.flair.isoforms.bed"),
          emit: isoform_bed

    tuple val(sample),
          path("${sample}.flair.isoforms.fa"),
          emit: isoform_fasta

    script:
    """
    set -euo pipefail

    # Assemble transcript isoforms from aligned long reads using FLAIR3.
    # Use the reference genome and GTF to guide transcript reconstruction.
    # Output GTF, BED and FASTA files for downstream characterisation.

    flair transcriptome \
        -b ${bam} \
        -g ${reference_fasta} \
        -f ${reference_gtf} \
        -o ${sample}.flair \
        -t ${task.cpus}
    """
}
