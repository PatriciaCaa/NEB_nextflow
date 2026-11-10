process PLOT_JUNCTION_FIGURES {
    tag "${comparison_name}:${gene}"

    // Keep plots separate by comparison and respect the configured run root.
    publishDir { "${params.get('outdir') ?: params.get('results_dir') ?: 'results'}/rna/junctions/figures/${comparison_name}" },
        mode: 'copy', overwrite: true

    input:
    tuple val(comparison_name), val(gene), val(chrom),
          path(annotated_tsv), path(usage_tsv),
          path(case_chains), path(control_chains)
    path reference_gtf
    path exon_script
    path plot_script
    val rscript_executable

    output:
    tuple val(comparison_name),
          path("${comparison_name}.*.png"),
          path("${comparison_name}.*.pdf"),
          path("${comparison_name}.reference_union_exons.tsv"),
          emit: figures

    script:
    """
    set -euo pipefail
    test -s ${annotated_tsv}
    test -s ${usage_tsv}
    test -s ${case_chains}
    test -s ${control_chains}
    command -v "${rscript_executable}" >/dev/null

    # Merge overlapping reference exons across all annotated transcripts
    # of the target gene to provide a common genomic context for the plots.

    python3 ${exon_script} \
        --reference-gtf ${reference_gtf} \
        --gene "${gene}" \
        --chrom "${chrom}" \
        --output "${comparison_name}.reference_union_exons.tsv"

    "${rscript_executable}" ${plot_script} \
        ${usage_tsv} \
        ${annotated_tsv} \
        "${comparison_name}.reference_union_exons.tsv" \
        ${case_chains} \
        ${control_chains} \
        . \
        "${comparison_name}"

    # Arc and junction-chain figures are required.
    for suffix in junction_arcs junction_chains; do
        test -s "${comparison_name}.\${suffix}.png"
        test -s "${comparison_name}.\${suffix}.pdf"
    done

    # Relative usage is optional when no sites are comparable.
    # If generated, both PNG and PDF must be present.
    if [[ -e "${comparison_name}.relative_usage.png" ||
          -e "${comparison_name}.relative_usage.pdf" ]]; then
        test -s "${comparison_name}.relative_usage.png"
        test -s "${comparison_name}.relative_usage.pdf"
    fi
    """
}
