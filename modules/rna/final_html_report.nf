process FINAL_HTML_REPORT {

    tag "${case_info.name}"

    label 'rna_report'

    cpus 2
    memory '8 GB'

    input:
    val case_info
    path combined_gtf
    path mane_transcripts
    path clinvar_variants
    path report_script

    output:
    tuple val(case_info.name),
          path("${case_info.name}.html"),
          path("${case_info.name}_transcripts_detected_plot.png"),
          path("${case_info.name}_proportion_plot.png"),
          path("${case_info.name}_transcript_structure_plot.png"),
          path("${case_info.name}_additional_info.csv"),
          emit: report_files

    script:
    def variant_names = (case_info.variant_name ?: []).join(',')
    def variant_positions = (case_info.variant_pos ?: []).join(',')
    def threshold = case_info.prop_threshold ?: 0.10

    """
    set -euo pipefail

    Rscript ${report_script} \
        --case_name "${case_info.name}" \
        --case_sample "${case_info.case_sample}" \
        --case_control "${case_info.case_control}" \
        --chrom "${case_info.chrom}" \
        --start_pos ${case_info.start_pos} \
        --end_pos ${case_info.end_pos} \
        --gene "${case_info.gene}" \
        --mane_id "${case_info.mane_id ?: ''}" \
        --variant_name "${variant_names}" \
        --variant_pos "${variant_positions}" \
        --prop_threshold ${threshold} \
        --merged_gtf ${combined_gtf} \
        --MANE_transcripts ${mane_transcripts} \
        --clinvar_variants ${clinvar_variants} \
        --output_html ${case_info.name}.html
    """
}
