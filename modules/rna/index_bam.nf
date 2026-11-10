process INDEX_BAM {

    tag "${sample}"
    label 'samtools'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai"),
          emit: indexed_bam

    script:
    """
    set -euo pipefail

    # The BAM is already coordinate-sorted by MINIMAP2.
    # Create its index for efficient gene-region queries.

    samtools index -@ ${task.cpus} ${bam}
    """
}
