process RELATIVE_JUNCTION_USAGE {
    tag "${comparison_name}"
    label 'python'

    publishDir "${params.get('outdir') ?: params.get('results_dir') ?: 'results'}/rna/junctions/usage",
        mode: 'copy', overwrite: true

    input:
    tuple val(comparison_name), path(annotated_comparison_tsv)
    path usage_script
    val min_site_support

    output:
    tuple val(comparison_name),
          path("${comparison_name}.junction_relative_usage.tsv"),
          path("${comparison_name}.junction_usage_sites.tsv"),
          path("${comparison_name}.junction_usage_summary.tsv"),
          emit: usage

    script:
    """
    set -euo pipefail

    # Calculate each junction's relative share of observed read support
    # at its donor and acceptor sites in case and control samples.
    # These descriptive proportions are not percent-spliced-in.

    python3 ${usage_script} \
        --input ${annotated_comparison_tsv} \
        --outdir . \
        --prefix "${comparison_name}" \
        --min-site-support ${min_site_support}
    """
}
