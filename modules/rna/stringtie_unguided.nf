process STRINGTIE_UNGUIDED {

    tag "${sample}"

    label 'rna_tools'

    cpus 4
    memory '8 GB'

    input:
    tuple val(sample), path(sorted_bam)

    output:
    tuple val(sample), path("${sample}_unguided.gtf"), emit: unguided_gtf

    script:
    """
    set -euo pipefail

    stringtie \
        ${sorted_bam} \
        -o ${sample}_unguided.gtf \
        -p ${task.cpus}
    """
}
