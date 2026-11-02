

process STRUCTURAL_VARIANTS {

    cpus 8 // I've hard-coded this for now I need to find ways to make it configurable

    label 'sniffles'

    publishDir "${params.results_dir}/structural_variants", mode: 'copy'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    path reference_fasta

    output:
    path "*.sv.vcf.gz", emit: sv_vcf
    path "*.sv.vcf.gz.tbi", emit: sv_vcf_index

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    sniffles \
        --input ${aligned_bam} \
        --reference ${reference_fasta} \
        --threads ${task.cpus} \
        --vcf \${sample_id}.sv.vcf.gz
    """
}

process STRUCTURAL_VARIANTS_QC{

    publishDir "${params.results_dir}/qc/structural_variants", mode: 'copy'

    label 'bcftools'

    input:
    path sv_vcf

    output:
    path "*.structural_variants.stats.txt", emit: stats

    script:
    """
    sample_id=\$(basename ${sv_vcf} .sv.vcf.gz)

    bcftools stats ${sv_vcf} \
        > \${sample_id}.structural_variants.stats.txt
    """
}
