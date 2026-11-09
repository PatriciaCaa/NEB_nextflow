

process ALIGN {

    cpus 8
    publishDir "${params.results_dir}/alignment", mode: 'copy'

    label 'alignment'    

    input:
    path sample_bam
    path reference_fasta

    output:
    tuple path("*.aligned.sorted.bam"),
        path("*.aligned.sorted.bam.bai"),
        emit: aligned_files  

    // This script extracts a sample ID from the input BAM filename so each sample gets unique output filenames.
    script:
    """
    sample_id=\$(basename ${sample_bam} .bam)

    samtools fastq ${sample_bam} \
        | minimap2 -t ${task.cpus} -ax map-ont ${reference_fasta} - \
        | samtools sort -@ ${task.cpus} -o \${sample_id}.aligned.sorted.bam

    samtools index -@ ${task.cpus} \${sample_id}.aligned.sorted.bam
    """

}

// Assess alignment quality using samtools flagstat and stats,
// and calculate read depth across the target genomic regions.

process ALIGNMENT_QC{
    cpus 2
    publishDir "${params.results_dir}/qc/alignment", mode: 'copy'
    
    label 'samtools'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    path target_bed

    output:
    path "*.flagstat.txt", emit: flagstat
    path "*.stats.txt", emit: stats
    path "*.target.depth.txt", emit: target_depth

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    samtools flagstat \
        -@ ${task.cpus} \
        ${aligned_bam} \
        > \${sample_id}.flagstat.txt 

    samtools stats \
        ${aligned_bam} \
        > \${sample_id}.stats.txt

    samtools depth \
        -a \
        -b ${target_bed} \
        ${aligned_bam} \
        > \${sample_id}.target.depth.txt
    """

}
