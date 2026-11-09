// Generate gene-focused variant reports from annotated SNPs, indels, SVs and CNVs.
// Include available variant annotations and genotype information.
// Highlight the specified known variant without excluding other findings.
// Visualise detected variants relative to the target gene and exon structure.

process GENERATE_VARIANT_REPORT {

    label 'python'    

    tag "DNA variant report"

    publishDir "${params.results_dir}/variant_report", mode: 'copy'

    input:
    path small_variants_vcf
    path structural_variants_vcf
    path cnv_vcf
    val target_gene
    path gene_bed
    path exon_bed
    val known_variant
    path report_script
    path plot_script

    output:
    path "target_gene_variants.csv", emit: variant_csv
    path "target_gene_variants_transcripts.csv", emit: transcript_csv
    path "target_gene_variants.png", emit: variant_png
    path "target_gene_variants.pdf", emit: variant_pdf

    script:
    """
    python3 ${report_script} \
        --small ${small_variants_vcf} \
        --sv ${structural_variants_vcf} \
        --cnv ${cnv_vcf} \
        --gene ${target_gene} \
        --gene_bed ${gene_bed} \
        --known_variant "${known_variant ?: ''}" \
        --output target_gene_variants.csv

    python3 ${plot_script} \
        --csv target_gene_variants.csv \
        --gene_bed ${gene_bed} --exon_bed ${exon_bed} \
        --gene "${target_gene}" \
        --output_prefix target_gene_variants
    """
}
