process MINIMAP2 {

    tag "${sample}"

    input:
    tuple val(sample), path(trimmed_fastq)
    path reference_fasta

    output:
    tuple val(sample), path("${sample}.sorted.bam"), emit: sorted_bam

    script:
    def strandMode = params.rna_minimap2_strand_mode
    if (!(strandMode in ['auto', 'forward'])) {
        error "Unsupported RNA strand mode: ${strandMode}"
    }
    def strandArg = strandMode == 'forward' ? '-uf' : ''
    """
    set -euo pipefail

    minimap2 \
        -ax splice \
        -t ${task.cpus} \
        ${strandArg} \
        --secondary=no \
        -C5 \
        ${reference_fasta} \
        ${trimmed_fastq} \
    | samtools sort \
        -@ ${task.cpus} \
        -o ${sample}.sorted.bam \
        -
    """
}
