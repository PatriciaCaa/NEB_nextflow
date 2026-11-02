process NANOSTAT {

    tag "${sample}"

    label 'process_low'

    publishDir "${params.results_dir}/qc/rna/nanostat", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_stats"), emit: stats

    script:
    """
    NanoStat \
        --fastq ${fastq} \
        --outdir ${sample}_stats \
        --name ${sample}
    """
}
