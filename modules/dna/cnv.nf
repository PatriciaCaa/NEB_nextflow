// QDNAseq CNV calling uses Oxford Nanopore's wf-cnv implementation
// and its bundled run_qdnaseq.r script for adaptive-sampling / shallow-WGS data.
// QDNAseq already generates QC-related outputs such as bin-level data
// and diagnostic plots, so no separate QC process is required here.
// The corrected CNV calls VCF is compressed and indexed for downstream annotation.
// The segmented VCF is retained separately in the QDNAseq output directory.

process QDNASEQ_CNV {

    label 'wf_cnv'

    publishDir "${params.results_dir}/cnv/qdnaseq", mode: 'copy'

    cpus 1
    memory '16 GB'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    val genome_build
    val qdnaseq_bin_size

    output:
    path "*.qdnaseq.cnv.vcf.gz", emit: qdnaseq_vcf
    path "*.qdnaseq.cnv.vcf.gz.tbi", emit: qdnaseq_vcf_index
    path "qdnaseq_output/*", emit: qdnaseq_outputs

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    mkdir -p qdnaseq_output

    run_qdnaseq.r \
        --bam ${aligned_bam} \
        --out_prefix qdnaseq_output/\${sample_id} \
        --binsize ${qdnaseq_bin_size} \
        --reference ${genome_build}

    # Preserve the raw QDNAseq VCFs
    mv qdnaseq_output/\${sample_id}_calls.vcf \
        qdnaseq_output/\${sample_id}_calls.raw.vcf

    mv qdnaseq_output/\${sample_id}_segs.vcf \
        qdnaseq_output/\${sample_id}_segs.raw.vcf

    # Fix the main CNV calls VCF used downstream
    fix_qdnaseq_vcf.py \
        -i qdnaseq_output/\${sample_id}_calls.raw.vcf \
        -o \${sample_id}.qdnaseq.cnv.vcf \
        --sample_id \${sample_id}

    # Repair the segmented VCF and retain it
    fix_qdnaseq_vcf.py \
        -i qdnaseq_output/\${sample_id}_segs.raw.vcf \
        -o qdnaseq_output/\${sample_id}_segs.vcf \
        --sample_id \${sample_id}

    bgzip \${sample_id}.qdnaseq.cnv.vcf
    tabix -f -p vcf \${sample_id}.qdnaseq.cnv.vcf.gz
    """
}


// Calculate binned genome coverage for Spectre CNV calling.

process MOSDEPTH_CNV {

    label 'mosdepth'

    publishDir "${params.results_dir}/cnv/spectre/mosdepth", mode: 'copy'

    input:
    tuple path(aligned_bam), path(aligned_bai)
    val bin_size

    output:
    tuple path("*.regions.bed.gz"),
          path("*.regions.bed.gz.csi"),
          emit: coverage

    script:
    """
    sample_id=\$(basename ${aligned_bam} .aligned.sorted.bam)

    mosdepth \
        -n \
        --by ${bin_size} \
        \${sample_id} \
        ${aligned_bam}

    tabix -C -p bed \${sample_id}.regions.bed.gz
    """
}

process SPECTRE_CNV {

    publishDir "${params.results_dir}/cnv/spectre", mode: 'copy'

    label 'spectre'

    input:
    tuple path(coverage_bed), path(coverage_index)
    path reference_fasta

    output:
    path "*.spectre.cnv.vcf.gz", emit: spectre_vcf
    path "spectre_output/*", emit: spectre_outputs

    script:
    """
    sample_id=\$(basename ${coverage_bed} .regions.bed.gz)

    mkdir -p spectre_output

    spectre CNVCaller \
        --coverage ${coverage_bed} \
        --sample-id \${sample_id} \
        --output-dir spectre_output \
        --reference ${reference_fasta}

    raw_vcf=\$(find spectre_output -name "*.vcf.gz" | head -n 1) 

    cp "\${raw_vcf}" \${sample_id}.spectre.cnv.vcf.gz
    """
}


process SPECTRE_CNV_QC {

    publishDir "${params.results_dir}/qc/cnv/spectre", mode: 'copy'

    label 'bcftools'

    input:
    path spectre_vcf

    output:
    path "*.spectre.stats.txt", emit: stats

    script:
    """
    sample_id=\$(basename ${spectre_vcf} .spectre.cnv.vcf.gz)

    bcftools stats ${spectre_vcf} \
        > \${sample_id}.spectre.stats.txt
    """
}
