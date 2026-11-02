process ANNOTATE_SPLICE_JUNCTIONS {

    tag "${comparison_name}:${gene}"
    label 'python'

    publishDir "${params.outdir ?: 'results'}/rna/junctions/annotations",
        mode: 'copy', overwrite: true

    input:
    tuple val(comparison_name), val(gene), path(comparison_tsv)
    path reference_gtf
    path annotation_script

    output:
    tuple val(comparison_name),
          path("${comparison_name}.annotated_junction_comparison.tsv"),
          path("${comparison_name}.annotation_summary.tsv"),
          emit: annotations

    script:
    """
    set -euo pipefail

    mkdir -p output

    python3 ${annotation_script} \
        --comparison ${comparison_tsv} \
        --reference-gtf ${reference_gtf} \
        --gene "${gene}" \
        --outdir output

    mv output/annotated_junction_comparison.tsv \
        "${comparison_name}.annotated_junction_comparison.tsv"

    mv output/annotation_summary.tsv \
        "${comparison_name}.annotation_summary.tsv"
    """
}
