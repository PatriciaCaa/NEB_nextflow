process CUTADAPT {

    tag "${sample}"

    input:
    tuple val(sample), path(fastq)
    path fwd_primer
    path rev_primer

    output:
    tuple val(sample), path("${sample}.trimmed.fq"), emit: trimmed_reads

    script:
    def mode = params.rna_trim_mode ?: 'independent'

    if (mode == 'linked') {
        """
        FWD=\$(grep -v '^>' ${fwd_primer} | tr -d '\\r\\n')
        REV=\$(grep -v '^>' ${rev_primer} | tr -d '\\r\\n')

        cutadapt \
            -j ${task.cpus} \
            --revcomp \
            -e 0.20 \
            -O 15 \
            -g "\${FWD}...\${REV}" \
            --discard-untrimmed \
            -o ${sample}.trimmed.fq \
            ${fastq}
        """
    } else if (mode == 'independent') {
        """
        cutadapt \
            -j ${task.cpus} \
            -g file:${fwd_primer} \
            -a file:${rev_primer} \
            --discard-untrimmed \
            -o ${sample}.trimmed.fq \
            ${fastq}
        """
    } else {
        error "Unsupported RNA trimming mode: ${mode}"
    }
}
