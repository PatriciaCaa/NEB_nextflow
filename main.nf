
// This file will call the DNA and RNA subworkflows

params.pod5_dir = null // i need the file path
params.barcode_kit = 'SQK-NBD114-24' // Given by Kylie as possible kit (needs confirmation)
params.model = 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' // This is the re-basecalling model Kylie used.
params.reference_fasta = null
params.reference_fai = null
params.annotation_gtf = null // GTF annotation matching the reference genome assembly
params.clair3_model = 'r1041_e82_400bps_sup_v520' // Clair3 model matched to Dorado v5.2.0 SUP R10.4.1 E8.2 400bps data.
params.threads = 8 // This is just a defafult, I can change it depending on the HPC resources
params.target_bed = null
params.qdnaseq_bin_size = 500
params.cnv_caller = 'qdnaseq'
params.spectre_bin_size = 1000
params.genome_build = null // Must match the reference FASTA used for the specific study. For QDNAseq use hg19 or hg38.
params.results_dir = 'results'
params.vep_assembly = null
params.vep_cache = null

// DNA reporting
params.target_gene = null
params.known_variant = null



// RNA inputs and execution options.
params.analysis_mode = 'dna'
params.samples = null
params.cases = []
params.fwd_primer = null
params.rev_primer = null
params.reference_gtf = null
params.mane_transcripts = null
params.clinvar_variants = null
params.report_script = "${projectDir}/scripts/final_report_nextflow.R"
params.rna_preprocessing = null
params.rna_trim_mode = 'independent'
params.rna_minimap2_strand_mode = 'auto'
params.rna_report_mode = 'comparison'
params.mane_transcript_id = null

include { DNA_WORKFLOW } from './subworkflows/dna_workflow.nf'
include { RNA_WORKFLOW } from './subworkflows/rna_workflow.nf'
include { CAPTURE_RUN_METADATA } from './modules/dna/provenance/run_metadata.nf'


workflow {
    def metadataPreprocessing = null

    if (params.analysis_mode == 'dna') {
    DNA_WORKFLOW(
        params.pod5_dir, 
        params.barcode_kit, 
        params.model,
        params.reference_fasta,
        params.reference_fai,
	params.annotation_gtf,
        params.clair3_model,
        params.threads,
        params.target_bed,
        params.qdnaseq_bin_size,
        params.spectre_bin_size,
        params.cnv_caller,
        params.genome_build,
        params.vep_assembly,
        params.vep_cache,
        params.target_gene,
        params.known_variant

    )


    } else if (params.analysis_mode == 'rna') {

        def preprocessing = params.rna_preprocessing ?:
            ((params.rna_trim_mode ?: 'independent') == 'none'
                ? 'none' : 'cutadapt')

        def reportMode = params.rna_report_mode ?: 'comparison'

        if (!['cutadapt', 'pychopper', 'none'].contains(preprocessing)) {
            error "Unsupported RNA preprocessing: ${preprocessing}"
        }

        if (!['comparison', 'none'].contains(reportMode)) {
            error "Unsupported RNA report mode: ${reportMode}"
        }

        def required = [
            'samples', 'reference_fasta', 'reference_gtf'
        ]

        if (preprocessing == 'cutadapt') {
            required += ['fwd_primer', 'rev_primer']
        }

        if (reportMode == 'comparison') {
            required += [
                'cases', 'mane_transcripts',
                'clinvar_variants', 'report_script'
            ]
        }

        required.each { key ->
            if (!params[key]) {
                error "Missing required RNA parameter: ${key}"
            }
        }

        metadataPreprocessing = preprocessing

        samples_ch = channel.fromList(params.samples)
            .map { sample ->
                tuple(
                    sample.id.toString(),
                    file(sample.fastq, checkIfExists: true)
                )
            }

        cases_ch = channel.fromList(params.cases ?: [])

        RNA_WORKFLOW(
            samples_ch,
            params.fwd_primer
                ? channel.value(file(params.fwd_primer, checkIfExists: true))
                : channel.empty(),
            params.rev_primer
                ? channel.value(file(params.rev_primer, checkIfExists: true))
                : channel.empty(),
            channel.value(file(params.reference_fasta, checkIfExists: true)),
            channel.value(file(params.reference_gtf, checkIfExists: true)),
            cases_ch,
            params.mane_transcripts
                ? channel.value(file(params.mane_transcripts, checkIfExists: true))
                : channel.empty(),
            params.clinvar_variants
                ? channel.value(file(params.clinvar_variants, checkIfExists: true))
                : channel.empty(),
            params.report_script
                ? channel.value(file(params.report_script, checkIfExists: true))
                : channel.empty()
        )

    } else {
        error "Unsupported analysis mode: ${params.analysis_mode}"
    }


    // Shared execution metadata for DNA and RNA.
    def runMetadata = [
        analysis_id: workflow.sessionId.toString(),
        run_name: workflow.runName,
        start_time: workflow.start.toString(),
        nextflow_version: nextflow.version.toString(),
        pipeline_revision: workflow.commitId ?: 'unavailable',
        analysis_type: params.analysis_mode.toString().toUpperCase(),
        parameters: params
    ]

    if (params.analysis_mode == 'rna') {
        runMetadata.preprocessing_method = metadataPreprocessing
    }

    CAPTURE_RUN_METADATA(
        groovy.json.JsonOutput.prettyPrint(
            groovy.json.JsonOutput.toJson(runMetadata)
        )
    )

}
