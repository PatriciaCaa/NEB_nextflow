process EXTRACT_GENE_TRANSCRIPT_SUMMARY {
    tag "${sample}:${assembler}:${target_gene}"
    label 'rna_report'

    publishDir "${params.results_dir}/rna/characterisation/summaries",
        mode: 'copy',
        saveAs: { name -> "${assembler}/${sample}/${name}" }


    input:
    tuple val(sample), val(assembler), path(assembly_gtf)
    path reference_gtf
    path extractor_script
    val target_gene

    output:
    tuple val(sample), val(assembler), path('transcript_summary.tsv'), emit: summary

    // Summarise assembled transcripts associated with the target gene.
    // Gene assignment is based on annotation evidence and genomic overlap, and does not independently validate transcript identity.

    script:
    """
    set -euo pipefail
    python3 ${extractor_script} \
        --reference-gtf ${reference_gtf} \
        --assembly-gtf ${assembly_gtf} \
        --gene '${target_gene}' \
        --sample-id '${sample}' \
        --assembler '${assembler}' \
        --output transcript_summary.tsv
    """
}
