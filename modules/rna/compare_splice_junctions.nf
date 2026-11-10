process COMPARE_SPLICE_JUNCTIONS {
    tag "${case_info.name}"

    publishDir "${params.outdir ?: 'results'}/rna/junctions/comparisons", mode: 'copy', overwrite: true

    input:
    tuple val(case_info), path(case_counts), path(control_counts)
    path compare_script

    output:
    tuple val(case_info.name),
        path("${case_info.name}.junction_comparison.tsv"),
        path("${case_info.name}.isopedia_splice_queries.bed"),
        path("${case_info.name}.introverse_lookup.tsv"),
        path("${case_info.name}.introverse_gene_list.csv"),
        emit: comparison

    script:
    """
    set -euo pipefail
    mkdir -p output

    # Compare junction read support between case and control samples.
    # This is a descriptive comparison, not a statistical test of differential splicing.

    python3 ${compare_script} --case ${case_counts} --control ${control_counts} --gene "${case_info.gene}" --min-support ${params.rna_junction_min_support ?: 1} --outdir output

    mv output/junction_comparison.tsv "${case_info.name}.junction_comparison.tsv"
    mv output/isopedia_splice_queries.bed "${case_info.name}.isopedia_splice_queries.bed"
    mv output/introverse_lookup.tsv "${case_info.name}.introverse_lookup.tsv"
    printf 'gene\\n%s\\n' "${case_info.gene}" > "${case_info.name}.introverse_gene_list.csv"
    """
}
