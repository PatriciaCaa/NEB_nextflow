
process EXTRACT_GENE_REGIONS {

    tag "${target_gene}"
    label 'python'

    publishDir "${params.results_dir}/qc/gene_coverage", mode: 'copy'

    input:
    val target_gene
    path annotation_gtf
    path reference_fai
    path gene_region_script

    output:
    path "${target_gene}.gene.bed", emit: gene_bed
    path "${target_gene}.gene.exons.bed", emit: exon_bed

    script:
    """
    python3 ${gene_region_script} \
        --gtf ${annotation_gtf} \
        --fai ${reference_fai} \
        --gene "${target_gene}" \
        --output ${target_gene}.gene.bed
    """
}


process GENE_COVERAGE_QC {

    tag "Gene coverage"
    label 'samtools'

    publishDir "${params.results_dir}/qc/gene_coverage", mode: 'copy'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    path gene_bed

    output:
    path "*.gene.depth.tsv", emit: gene_depth

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    samtools depth -a -b ${gene_bed} ${aligned_bam} \
        > \${sample_id}.gene.depth.tsv
    """
}


