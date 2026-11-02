process STRINGTIE_REQUANT {

    tag "${sample}"

    label 'rna_tools'

    cpus 4
    memory '8 GB'

    input:
    tuple val(sample), path(sorted_bam)
    path merged_gtf

    output:
    tuple val(sample),
          path("${sample}_mergedQuant.gtf"),
          path("${sample}_mergedQuant.abundance.txt"),
          emit: requantified

    script:
    """
    set -euo pipefail

    stringtie \
        ${sorted_bam} \
        -e \
        -G ${merged_gtf} \
        -o ${sample}_mergedQuant.gtf \
        -A ${sample}_mergedQuant.abundance.txt \
        -p ${task.cpus}
    """
}
