// Re-basecalling module for raw Oxford Nanopore adaptive-sampling data.
// Uses POD5 input and will perform higher-accuracy Dorado basecalling.

process BASECALL{


    label 'dorado'

    accelerator 1, type: 'nvidia'

    publishDir "${params.results_dir}/basecalling", mode: 'copy'

    input:
    path pod5_dir          // folder containing raw POD5 files
    val barcode_kit        // dorado barcode kit identifier
    val model              // dorado basecalling model

    output:
    path "basecalled.bam"

    // Assign barcode classifications during basecalling so that
    // demultiplexing can use the existing labels without reclassification.

    script:
    """
    dorado basecaller ${model} ${pod5_dir} \
    --kit-name ${barcode_kit} \
    > basecalled.bam

    """
}

// QC for basecalled Oxford Nanopore reads

process BASECALL_NANOSTAT_QC {

    label 'nanostat'

    cpus 2

    publishDir "${params.results_dir}/qc/basecalling", mode: 'copy'

    input:
    path basecalled_bam

    output:
    path "*.nanostat.txt", emit: nanostat

    script:
    """
    sample_id=\$(basename ${basecalled_bam} .bam)

    NanoStat \
        --bam ${basecalled_bam} \
        --threads ${task.cpus} \
        --name \${sample_id}.nanostat.txt
    """
}
