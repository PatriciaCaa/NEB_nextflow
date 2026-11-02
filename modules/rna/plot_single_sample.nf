process PLOT_SINGLE_SAMPLE {

    tag "${sample}"
    label 'rna_report'

    publishDir "${params.results_dir}/rna/characterisation/single_sample",
        mode: 'copy',
        saveAs: { name -> "${sample}/${name}" }


    input:
    tuple val(sample), path(summary_tsv)
    path plot_script
    val mane_id

    output:
    tuple val(sample), path("single_sample_figures"),
        emit: figures

    script:
    def mane_arg = mane_id ? "\"${mane_id}\"" : ""

    """
    set -euo pipefail

    mkdir -p single_sample_figures

    Rscript ${plot_script} \
        ${summary_tsv} \
        "${sample}" \
        single_sample_figures \
        ${mane_arg}
    """
}
