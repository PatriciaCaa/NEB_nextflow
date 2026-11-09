// Orchestrates the long-read DNA analysis workflow,
// from basecalling and alignment to variant detection,
// annotation and gene-focused reporting.



// Import the DNA analysis processes defined in separate module files.
// These processes are invoked within DNA_WORKFLOW below.

include { BASECALL; BASECALL_NANOSTAT_QC } from '../modules/dna/basecall.nf'
include { ALIGN; ALIGNMENT_QC } from '../modules/dna/align.nf'
include { SMALL_VARIANTS; SMALL_VARIANTS_QC; NORMALISE_SMALL_VARIANTS } from '../modules/dna/small_variants.nf'
include { BASECALL; BASECALL_NANOSTAT_QC } from '../modules/dna/basecall.nf'
include { ALIGN; ALIGNMENT_QC } from '../modules/dna/align.nf'
include { SMALL_VARIANTS; SMALL_VARIANTS_QC; NORMALISE_SMALL_VARIANTS } from '../modules/dna/small_variants.nf'
include { STRUCTURAL_VARIANTS; STRUCTURAL_VARIANTS_QC } from '../modules/dna/structural_variants.nf'
include { DEMULTIPLEX; DEMULTIPLEX_QC } from '../modules/dna/demultiplex.nf'
include {
    QDNASEQ_CNV;
    MOSDEPTH_CNV;
    SPECTRE_CNV;
    SPECTRE_CNV_QC
} from '../modules/dna/cnv.nf'
include { 
    ANNOTATE_SMALL_VARIANTS; 
    ANNOTATE_STRUCTURAL_VARIANTS;
    ANNOTATE_SPECTRE_CNV;
    ANNOTATE_QDNASEQ_CNV 
} from '../modules/dna/annotate_variants.nf'
include { GENERATE_VARIANT_REPORT } from '../modules/dna/variant_report.nf'
include { EXTRACT_GENE_REGIONS } from '../modules/dna/gene_regions.nf'

workflow DNA_WORKFLOW {

    // Declare the inputs received from the parent workflow.
    take:
    pod5_dir
    barcode_kit
    model
    reference_fasta
    reference_fai
    annotation_gtf	
    clair3_model
    threads
    target_bed
    qdnaseq_bin_size
    spectre_bin_size
    cnv_caller
    genome_build
    vep_assembly
    vep_cache
    target_gene
    known_variant

    // Define the analysis steps and process connections.
    main:

    BASECALL(
        pod5_dir,
        barcode_kit,
        model
    )

    BASECALL_NANOSTAT_QC(
        BASECALL.out
    )

    DEMULTIPLEX(
        BASECALL.out
    )

    classified_bams_ch = DEMULTIPLEX.out.classified_bams.flatten()

    DEMULTIPLEX_QC(
        classified_bams_ch
    )

    ALIGN(
        classified_bams_ch,
        reference_fasta
    )

    ALIGNMENT_QC(
        ALIGN.out.aligned_files,
        target_bed
    )

    EXTRACT_GENE_REGIONS(
        target_gene,
        annotation_gtf,
        reference_fai,
        channel.value(file("${projectDir}/scripts/extract_gene_regions.py", checkIfExists: true))
    )

    SMALL_VARIANTS(
        ALIGN.out.aligned_files,
        reference_fasta,
        reference_fai,
        clair3_model,
        threads,
        target_bed
    )

    SMALL_VARIANTS_QC(
        SMALL_VARIANTS.out.clair3_dir
    )

    NORMALISE_SMALL_VARIANTS(
        SMALL_VARIANTS.out.clair3_dir,
        reference_fasta,
        reference_fai
    )

    ANNOTATE_SMALL_VARIANTS(
        NORMALISE_SMALL_VARIANTS.out.normalised_vcf,
        reference_fasta,
        vep_cache,
        vep_assembly
    )

    STRUCTURAL_VARIANTS(
        ALIGN.out.aligned_files,
        reference_fasta
    )

    ANNOTATE_STRUCTURAL_VARIANTS(
        STRUCTURAL_VARIANTS.out.sv_vcf,
        reference_fasta,
        vep_cache,
        vep_assembly
    )

    STRUCTURAL_VARIANTS_QC(
        STRUCTURAL_VARIANTS.out.sv_vcf
    )

    // Select the CNV analysis method and collect its annotated VCF
    // for use in the final variant report.

    def annotated_cnv_vcf_ch

    if (cnv_caller == 'qdnaseq') {

        QDNASEQ_CNV(
            ALIGN.out.aligned_files,
            genome_build,
            qdnaseq_bin_size  
        )

        ANNOTATE_QDNASEQ_CNV(
            QDNASEQ_CNV.out.qdnaseq_vcf,
            reference_fasta,
            vep_cache,
            vep_assembly
        )

        annotated_cnv_vcf_ch = ANNOTATE_QDNASEQ_CNV.out.annotated_qdnaseq_vcf

    } else if (cnv_caller == 'spectre') {

        MOSDEPTH_CNV(
            ALIGN.out.aligned_files,
            spectre_bin_size
        )

        SPECTRE_CNV(
            MOSDEPTH_CNV.out.coverage,
            reference_fasta
        )

        ANNOTATE_SPECTRE_CNV(
            SPECTRE_CNV.out.spectre_vcf,
            reference_fasta,
            vep_cache,
            vep_assembly
        )

        annotated_cnv_vcf_ch = ANNOTATE_SPECTRE_CNV.out.annotated_cnv_vcf

        SPECTRE_CNV_QC(
            SPECTRE_CNV.out.spectre_vcf
        )

    } else {

        error "Unsupported CNV caller: ${cnv_caller}. Use 'qdnaseq' or 'spectre'."
    }

    // Combine annotated small variants, structural variants and CNVs
    // with gene coordinates to generate the gene-focused report.

    GENERATE_VARIANT_REPORT(
        ANNOTATE_SMALL_VARIANTS.out.annotated_vcf,
        ANNOTATE_STRUCTURAL_VARIANTS.out.annotated_sv_vcf,
        annotated_cnv_vcf_ch,
        target_gene,
        EXTRACT_GENE_REGIONS.out.gene_bed,
        EXTRACT_GENE_REGIONS.out.exon_bed,
        known_variant,
        file("${projectDir}/scripts/build_variant_report.py", checkIfExists: true),
        file("${projectDir}/scripts/plot_gene_variants.py", checkIfExists: true)
    )
}
