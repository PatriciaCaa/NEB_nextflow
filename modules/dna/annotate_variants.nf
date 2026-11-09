// Annotate small variants, structural variants and CNVs using Ensembl VEP.
// Compress and index the annotated VCF files for downstream analysis.

process ANNOTATE_SMALL_VARIANTS {

    tag "VEP small variant annotation"

    cpus 4

    label 'vep'

    publishDir "${params.results_dir}/annotation/small_variants", mode: 'copy'

    input:
    path vcf
    path reference_fasta
    val vep_cache
    val genome_build

    output:
    path "*.small_variants.vep.vcf.gz", emit: annotated_vcf
    path "*.small_variants.vep.vcf.gz.tbi", emit: annotated_vcf_index
    path "*.small_variants.vep_summary.html", emit: vep_stats

    script:
    """
    sample_id=\$(basename ${vcf} .small_variants.normalised.vcf.gz)

    vep \
        --input_file ${vcf} \
        --output_file STDOUT \
        --vcf \
        --compress_output bgzip \
        --cache \
        --offline \
        --dir_cache ${vep_cache} \
        --assembly ${genome_build} \
        --fasta ${reference_fasta} \
        --fork ${task.cpus} \
        --symbol \
        --hgvs \
        --canonical \
        --mane \
        --variant_class \
        --numbers \
        --domains \
        --protein \
        --biotype \
        --check_existing \
        --af_gnomade \
        --sift b \
        --polyphen b \
        --stats_file \${sample_id}.small_variants.vep_summary.html \
        --force_overwrite \
        > \${sample_id}.small_variants.vep.vcf.gz

    tabix -p vcf \${sample_id}.small_variants.vep.vcf.gz
    """
}

process ANNOTATE_STRUCTURAL_VARIANTS {

    tag "VEP structural variant annotation"

    cpus 4

    label 'vep'

    publishDir "${params.results_dir}/annotation/structural_variants", mode: 'copy'

    input:
    path vcf
    path reference_fasta
    val vep_cache
    val genome_build

    output:
    path "*.structural_variants.vep.vcf.gz", emit: annotated_sv_vcf
    path "*.structural_variants.vep.vcf.gz.tbi", emit: annotated_sv_vcf_index
    path "*.structural_variants.vep_summary.html", emit: sv_vep_stats

    script:
    """
    sample_id=\$(basename ${vcf} .sv.vcf.gz)

    vep \
        --input_file ${vcf} \
        --output_file STDOUT \
        --vcf \
        --compress_output bgzip \
        --cache \
        --offline \
        --dir_cache ${vep_cache} \
        --assembly ${genome_build} \
        --fasta ${reference_fasta} \
        --fork ${task.cpus} \
        --symbol \
        --canonical \
        --mane \
        --variant_class \
        --biotype \
        --check_existing \
        --stats_file \${sample_id}.structural_variants.vep_summary.html \
        --force_overwrite \
        > \${sample_id}.structural_variants.vep.vcf.gz

    tabix -p vcf \${sample_id}.structural_variants.vep.vcf.gz
    """
}

process ANNOTATE_SPECTRE_CNV {

    tag "VEP Spectre CNV annotation"

    cpus 4

    label 'vep'

    publishDir "${params.results_dir}/annotation/cnv/spectre", mode: 'copy'

    input:
    path vcf
    path reference_fasta
    val vep_cache
    val genome_build

    output:
    path "*.spectre.cnv.vep.vcf.gz", emit: annotated_cnv_vcf
    path "*.spectre.cnv.vep.vcf.gz.tbi", emit: annotated_cnv_vcf_index
    path "*.spectre.cnv.vep_summary.html", emit: cnv_vep_stats

    script:
    """
    sample_id=\$(basename ${vcf} .spectre.cnv.vcf.gz)

    vep \
        --input_file ${vcf} \
        --output_file STDOUT \
        --vcf \
        --compress_output bgzip \
        --cache \
        --offline \
        --dir_cache ${vep_cache} \
        --assembly ${genome_build} \
        --fasta ${reference_fasta} \
        --fork ${task.cpus} \
        --symbol \
        --canonical \
        --mane \
        --variant_class \
        --biotype \
        --check_existing \
        --stats_file \${sample_id}.spectre.cnv.vep_summary.html \
        --force_overwrite \
        > \${sample_id}.spectre.cnv.vep.vcf.gz

    tabix -p vcf \${sample_id}.spectre.cnv.vep.vcf.gz
    """
}

process ANNOTATE_QDNASEQ_CNV {

    tag "VEP QDNAseq CNV annotation"

    cpus 4

    label 'vep'

    publishDir "${params.results_dir}/annotation/cnv/qdnaseq", mode: 'copy'

    input:
    path vcf
    path reference_fasta
    val vep_cache
    val genome_build

    output:
    path "*.qdnaseq.cnv.vep.vcf.gz", emit: annotated_qdnaseq_vcf
    path "*.qdnaseq.cnv.vep.vcf.gz.tbi", emit: annotated_qdnaseq_vcf_index
    path "*.qdnaseq.cnv.vep_summary.html", emit: qdnaseq_vep_stats

    script:
    """
    sample_id=\$(basename ${vcf} .qdnaseq.cnv.vcf.gz)

    vep \
        --input_file ${vcf} \
        --output_file STDOUT \
        --vcf \
        --compress_output bgzip \
        --cache \
        --offline \
        --dir_cache ${vep_cache} \
        --assembly ${genome_build} \
        --fasta ${reference_fasta} \
        --fork ${task.cpus} \
        --symbol \
        --canonical \
        --mane \
        --variant_class \
        --biotype \
        --check_existing \
        --stats_file \${sample_id}.qdnaseq.cnv.vep_summary.html \
        --force_overwrite \
        > \${sample_id}.qdnaseq.cnv.vep.vcf.gz

    tabix -p vcf \${sample_id}.qdnaseq.cnv.vep.vcf.gz
    """
}
