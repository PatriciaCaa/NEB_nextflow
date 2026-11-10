process STRINGTIE_MERGE {

    tag "all_samples"

    label 'rna_tools'

    cpus 4
    memory '8 GB'

    input:
    path gtfs
    path reference_gtf

    output:
    path "merged.gtf", emit: merged_gtf

    script:

    // Sort the input GTF files by filename for reproducible ordering.
    // Prepare a newline-separated list for StringTie2's merge operation.
    def gtf_list = gtfs.sort { gtf -> gtf.name }.collect { gtf -> gtf.name }.join('\n')

    """
    set -euo pipefail

    cat > gtf_list.txt <<'GTF_LIST'
${gtf_list}
GTF_LIST

    stringtie \
        --merge \
        -o merged.gtf \
        -G ${reference_gtf} \
        -p ${task.cpus} \
        gtf_list.txt
    """
}
