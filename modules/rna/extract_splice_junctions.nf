process EXTRACT_SPLICE_JUNCTIONS {
    tag "${sample}:${gene}"
    label 'samtools'

    publishDir "${params.outdir ?: 'results'}/rna/junctions/per_sample", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam), path(bai)
    path reference_gtf
    path extractor_script
    val gene

    output:
    tuple val(sample), path("${sample}.${gene}.junction_counts.tsv"), emit: counts
    tuple val(sample), path("${sample}.${gene}.junction_chains.tsv"), emit: chains
    tuple val(sample), path("${sample}.${gene}.extraction_summary.tsv"), emit: summaries

    script:
    """
    set -euo pipefail

    region=\$(awk -F '\t' -v query='gene_name "${gene}";' '\$3=="gene" && index(\$9, query) {print \$1 ":" \$4 "-" \$5; exit}' ${reference_gtf})
    if [[ -z "\$region" ]]; then
        echo 'ERROR: Gene ${gene} was not found in the provided reference GTF' >&2
        exit 1
    fi
    echo "Junction extraction: ${sample}, ${gene}, \$region"

    mkdir -p output

    # Analyse reads overlapping the target gene, excluding unmapped,
    # Retain primary mapped alignments by excluding unmapped reads (4),
    # secondary alignments (256) and supplementary alignments (2048).

    samtools view -F 2308 ${bam} "\$region" | python3 ${extractor_script} --sample "${sample}" --min-mapq ${params.rna_junction_min_mapq ?: 20} --min-anchor ${params.rna_junction_min_anchor ?: 8} --outdir output

    mv output/junction_counts.tsv "${sample}.${gene}.junction_counts.tsv"
    mv output/junction_chains.tsv "${sample}.${gene}.junction_chains.tsv"
    mv output/extraction_summary.tsv "${sample}.${gene}.extraction_summary.tsv"
    """
}
