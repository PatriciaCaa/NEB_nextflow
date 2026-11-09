

process SMALL_VARIANTS {

    cpus { threads }
    publishDir "${params.results_dir}/small_variants", mode: 'copy'

    label 'clair3'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    // Clair3 also needs the same reference genome the reads were aligned against
    path reference_fasta
    path reference_fai
    val clair3_model
    val threads
    path target_bed

    output:
    path "*.clair3/merge_output.vcf.gz", emit: vcf
    path "*.clair3/merge_output.vcf.gz.tbi", emit: vcf_index
    path "*.clair3", emit: clair3_dir

    /// Call small variants using Clair3 and the specified ONT-compatible model.
    // Restrict variant calling to the genomic regions defined in target_bed.
    // The Clair3 model directory is supplied through the clair3_model parameter.

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    run_clair3.sh \
        --bam_fn=${aligned_bam} \
        --ref_fn=${reference_fasta} \
        --model_path=${clair3_model} \
        --platform=ont \
        --threads=${threads} \
        --bed_fn=${target_bed} \
        --output=\${sample_id}.clair3
        
    """
}

process SMALL_VARIANTS_QC {

    publishDir "${params.results_dir}/qc/small_variants", mode: 'copy'
    label 'bcftools'

    input:
    path clair3_dir

    output:
    path "*.small_variants.stats.txt", emit: stats

    script:
    """
    sample_id=\$(basename ${clair3_dir} .clair3)

    bcftools stats ${clair3_dir}/merge_output.vcf.gz \
        > \${sample_id}.small_variants.stats.txt
    """
}

// Normalise variant representation against the reference genome,
// including left-alignment of indels and splitting multiallelic records.

process NORMALISE_SMALL_VARIANTS {

    publishDir "${params.results_dir}/small_variants", mode: 'copy'
    label 'bcftools'

    input:
    path clair3_dir
    path reference_fasta
    path reference_fai

    output:
    path "*.small_variants.normalised.vcf.gz", emit: normalised_vcf
    path "*.small_variants.normalised.vcf.gz.tbi", emit: normalised_vcf_index

    script:
    """
    sample_id=\$(basename ${clair3_dir} .clair3)

    bcftools norm \
        -f ${reference_fasta} \
        -m -any \
        -Oz \
        -o \${sample_id}.small_variants.normalised.vcf.gz \
        ${clair3_dir}/merge_output.vcf.gz

    bcftools index -t \
        \${sample_id}.small_variants.normalised.vcf.gz
    """
}
