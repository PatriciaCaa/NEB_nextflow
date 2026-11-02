process COMBINE_REQUANT_GTFS {

    tag "all_samples"

    cpus 1
    memory '2 GB'

    input:
    path gtfs

    output:
    path "final_project_output.gtf", emit: combined_gtf

    script:
    """
    set -euo pipefail

    : > final_project_output.gtf

    for gtf in *_mergedQuant.gtf; do
        if [ ! -f "\$gtf" ]; then
            echo "ERROR: No requantification GTF files found" >&2
            exit 1
        fi

        sample_base=\$(basename "\$gtf" _mergedQuant.gtf)

        awk -v SAMPLE="\$sample_base" '
            BEGIN { FS = OFS = sprintf("%c", 9) }

            /^#/ { print; next }

            NF >= 9 {
                if (substr(\$9, length(\$9), 1) != ";") {
                    \$9 = \$9 ";"
                }

                \$9 = \$9 " sample_id " sprintf("%c", 34) SAMPLE sprintf("%c", 34) ";"
                print
            }
        ' "\$gtf" >> final_project_output.gtf
    done
    """
}
