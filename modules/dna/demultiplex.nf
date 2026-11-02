// Demultiplexing module for barcoded Oxford Nanopore DNA reads.
// Assumes barcode classification has already been assigned during basecalling.

process DEMULTIPLEX {

    publishDir "${params.results_dir}/demultiplexing", mode: 'copy'

    label 'dorado'

    input:
    path basecalled_bam

    output:
    path "demultiplexed/*_barcode*.bam", emit: classified_bams
    path "demultiplexed/unclassified.bam", emit: unclassified_bam, optional: true

    script:
    """
    mkdir -p demultiplexed

    dorado demux \
        --no-classify \
        --output-dir demultiplexed \
        ${basecalled_bam}
    """
}

// QUality control step for the output of the demultiplexing stage
// Counts how many reads are in each barcode BAM
process DEMULTIPLEX_QC{

    label 'samtools'

    cpus 1

    publishDir "${params.results_dir}/qc/demultiplexing", mode: 'copy'
    
    input:
    path classified_bam

    output:
    path "*.read_count.txt", emit: read_count

    script:
    """
    sample_id=\$(basename ${classified_bam} .bam)

    samtools view -c ${classified_bam} \
        > \${sample_id}.read_count.txt

    """
}
