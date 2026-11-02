process PLOT_ALL_SAMPLES {

    tag "Transcript characterisation"
    label 'rna_report'

    publishDir "${params.results_dir}/rna/characterisation/all_samples",
        mode: 'copy'


    input:
    path summary_tsv
    path plot_script

    output:
    path "all_samples_figures", emit: figures

    script:
    """
    set -euo pipefail

    mkdir -p all_samples_figures

    Rscript ${plot_script} \
        ${summary_tsv} \
        all_samples_figures
    """
}
