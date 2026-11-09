// RNA subworkflow for long-read transcript analysis.
// Coordinates read preprocessing, alignment, splice-junction analysis,
// transcript assembly, comparative analysis and reporting.
// Incorporates RAPID-derived stages alongside newly developed modules.

include { NANOSTAT } from '../modules/rna/nanostat'
include { CUTADAPT } from '../modules/rna/cutadapt'
include { PYCHOPPER } from '../modules/rna/pychopper'
include { MINIMAP2 } from '../modules/rna/minimap2'
include { INDEX_BAM } from '../modules/rna/index_bam'
include { FLAIR_TRANSCRIPTOME } from '../modules/rna/flair_transcriptome'
include { EXTRACT_SPLICE_JUNCTIONS } from '../modules/rna/extract_splice_junctions'
include { COMPARE_SPLICE_JUNCTIONS } from '../modules/rna/compare_splice_junctions'
include { STRINGTIE_UNGUIDED } from '../modules/rna/stringtie_unguided'
include { STRINGTIE_MERGE } from '../modules/rna/stringtie_merge'
include { STRINGTIE_REQUANT } from '../modules/rna/stringtie_requant'
include { COMBINE_REQUANT_GTFS } from '../modules/rna/combine_requant_gtfs'
include { FINAL_HTML_REPORT } from '../modules/rna/final_html_report'
include { EXTRACT_GENE_TRANSCRIPT_SUMMARY } from '../modules/rna/extract_gene_transcript_summary'
include { PLOT_SINGLE_SAMPLE } from '../modules/rna/plot_single_sample'
include { PLOT_ALL_SAMPLES } from '../modules/rna/plot_all_samples'

include { ANNOTATE_SPLICE_JUNCTIONS } from '../modules/rna/annotate_splice_junctions'
include { RELATIVE_JUNCTION_USAGE } from '../modules/rna/relative_junction_usage'
include { PLOT_JUNCTION_FIGURES } from '../modules/rna/plot_junction_figures'

workflow RNA_WORKFLOW {

    take:
    samples_ch
    fwd_primer_ch
    rev_primer_ch
    reference_fasta_ch
    reference_gtf_ch
    cases_ch
    mane_transcripts_ch
    clinvar_variants_ch
    report_script_ch

    main:
    // QC can run independently of trimming.
    NANOSTAT(samples_ch)

    // Select the RNA read preprocessing method.
    // Use Cutadapt by default unless trimming is disabled.
    def preprocessing = params.rna_preprocessing ?:
        ((params.rna_trim_mode ?: 'independent') == 'none'
            ? 'none' : 'cutadapt')

    pychopper_stats_ch = channel.empty()
    pychopper_unclassified_ch = channel.empty()
    pychopper_rescued_ch = channel.empty()
    pychopper_report_ch = channel.empty()

    if (preprocessing == 'none') {
        reads_for_alignment_ch = samples_ch

    } else if (preprocessing == 'cutadapt') {
        CUTADAPT(samples_ch, fwd_primer_ch, rev_primer_ch)
        reads_for_alignment_ch = CUTADAPT.out.trimmed_reads

    } else if (preprocessing == 'pychopper') {
        PYCHOPPER(samples_ch)
        reads_for_alignment_ch = PYCHOPPER.out.trimmed_reads
        pychopper_stats_ch = PYCHOPPER.out.stats
        pychopper_unclassified_ch = PYCHOPPER.out.unclassified
        pychopper_rescued_ch = PYCHOPPER.out.rescued
        pychopper_report_ch = PYCHOPPER.out.report

    } else {
        error "Unsupported RNA preprocessing: ${preprocessing}"
    }

    MINIMAP2(
        reads_for_alignment_ch,
        reference_fasta_ch
    )

    // Select StringTie2, FLAIR3 or both transcript assemblers.
    // StringTie2 is the default when no assembler is specified.
    // Validate the selection before enabling the corresponding processes.

    def assembler = params.rna_assembler ?: 'stringtie'
    if (!(assembler in ['stringtie', 'flair', 'both'])) {
        error "Unsupported RNA assembler: ${assembler}"
    }
    def useStringtie = assembler in ['stringtie', 'both']
    def useFlair = assembler in ['flair', 'both']

    // Index once for junction analysis and/or FLAIR3.
    if (params.rna_junction_analysis != false || useFlair) {
        INDEX_BAM(MINIMAP2.out.sorted_bam)
    }

    // Primary internal-splicing analysis, independent of transcript assembly.
    // Runs on each aligned sample even if no comparison has been configured.
    junction_counts_ch = channel.empty()
    junction_chains_ch = channel.empty()
    junction_summaries_ch = channel.empty()
    junction_comparisons_ch = channel.empty()
    junction_annotations_ch = channel.empty()
    junction_relative_usage_ch = channel.empty()
    junction_figures_ch = channel.empty()

    if (params.rna_junction_analysis != false) {

        // Use the explicitly configured junction-analysis gene when available.
        // Otherwise, use the first case's gene, falling back to NEB.

        def junctionGene = params.rna_junction_gene ?:
            ((params.cases ?: []).size() > 0 ? params.cases[0].gene : 'NEB')

        EXTRACT_SPLICE_JUNCTIONS(
            INDEX_BAM.out.indexed_bam,
            reference_gtf_ch,
            channel.value(file("${projectDir}/scripts/rna/extract_splice_junctions.py", checkIfExists: true)),
            junctionGene.toString()
        )
        junction_counts_ch = EXTRACT_SPLICE_JUNCTIONS.out.counts
        junction_chains_ch = EXTRACT_SPLICE_JUNCTIONS.out.chains
        junction_summaries_ch = EXTRACT_SPLICE_JUNCTIONS.out.summaries

        // Match junction counts to the configured case and control sample IDs.
        // Pair the matching counts by comparison name to support multiple comparisons.
        def case_matches = cases_ch
            .combine(EXTRACT_SPLICE_JUNCTIONS.out.counts)
            .filter { info, sample, _counts -> info.case_sample.toString() == sample.toString() }
            .map { info, _sample, counts -> tuple(info.name.toString(), info, counts) }

        def control_matches = cases_ch
            .combine(EXTRACT_SPLICE_JUNCTIONS.out.counts)
            .filter { info, sample, _counts -> info.case_control.toString() == sample.toString() }
            .map { info, _sample, counts -> tuple(info.name.toString(), counts) }

        def paired_counts = case_matches.join(control_matches)
            .map { _name, info, case_counts, control_counts ->
                tuple(info, case_counts, control_counts)
            }

        COMPARE_SPLICE_JUNCTIONS(
            paired_counts,
            channel.value(file("${projectDir}/scripts/rna/compare_splice_junctions.py", checkIfExists: true))
        )
        junction_comparisons_ch = COMPARE_SPLICE_JUNCTIONS.out.comparison

        // Match each comparison to its configured target gene.
        annotation_inputs_ch = cases_ch
            .map { info -> tuple(info.name.toString(), info.gene.toString()) }
            .join(COMPARE_SPLICE_JUNCTIONS.out.comparison)
            .map { name, gene, comparison, _bed, _lookup, _gene_list ->
                tuple(name, gene, comparison)
            }

        ANNOTATE_SPLICE_JUNCTIONS(
            annotation_inputs_ch,
            reference_gtf_ch,
            channel.value(file("${projectDir}/scripts/rna/annotate_splice_junctions.py", checkIfExists: true))
        )

        junction_annotations_ch = ANNOTATE_SPLICE_JUNCTIONS.out.annotations

        // Annotated junctions carry gene strand information for donor/acceptor grouping.
        def usage_inputs_ch = ANNOTATE_SPLICE_JUNCTIONS.out.annotations
            .map { comparison_name, annotated_tsv, _annotation_summary ->
                tuple(comparison_name, annotated_tsv)
            }

        RELATIVE_JUNCTION_USAGE(
            usage_inputs_ch,
            channel.value(file(
                "${projectDir}/scripts/rna/calculate_relative_junction_usage.py",
                checkIfExists: true
            )),
            (params.get('rna_junction_min_site_support') ?: 10) as Integer
        )

        junction_relative_usage_ch = RELATIVE_JUNCTION_USAGE.out.usage

        // Generate optional read-level splice-junction figures when explicitly enabled.
        // These figures combine junction annotations, relative usage and observed
        // junction chains from the corresponding case and control samples.
        if (params.get('rna_junction_figures') == true) {
            def figure_meta_ch = cases_ch.map { info ->
                tuple(info.name.toString(), info.gene.toString(),
                      info.chrom.toString(), info.case_sample.toString(),
                      info.case_control.toString())
            }
            def figure_annotation_ch = ANNOTATE_SPLICE_JUNCTIONS.out.annotations
                .map { name, annotated, _summary -> tuple(name, annotated) }
            def figure_usage_ch = RELATIVE_JUNCTION_USAGE.out.usage
                .map { name, usage, _sites, _summary -> tuple(name, usage) }

            // Join each sample's own observed chains, not a Cartesian product.
            def figure_case_ch = cases_ch
                .combine(EXTRACT_SPLICE_JUNCTIONS.out.chains)
                .filter { info, sample, _chains ->
                    info.case_sample.toString() == sample.toString()
                }
                .map { info, _sample, chains -> tuple(info.name.toString(), chains) }
            def figure_control_ch = cases_ch
                .combine(EXTRACT_SPLICE_JUNCTIONS.out.chains)
                .filter { info, sample, _chains ->
                    info.case_control.toString() == sample.toString()
                }
                .map { info, _sample, chains -> tuple(info.name.toString(), chains) }

            def figure_inputs_ch = figure_meta_ch
                .join(figure_annotation_ch)
                .join(figure_usage_ch)
                .join(figure_case_ch)
                .join(figure_control_ch)
                .map { name, gene, chrom, _case, _control, anno, usage, case_chains, control_chains ->
                    tuple(name, gene, chrom, anno, usage, case_chains, control_chains)
                }

            PLOT_JUNCTION_FIGURES(
                figure_inputs_ch,
                reference_gtf_ch,
                channel.value(file(
                    "${projectDir}/scripts/rna/extract_reference_gene_exons.py",
                    checkIfExists: true
                )),
                channel.value(file(
                    "${projectDir}/scripts/rna/plot_junction_structure.R",
                    checkIfExists: true
                )),
                params.get('rna_plot_rscript') ?: 'Rscript'
            )
            junction_figures_ch = PLOT_JUNCTION_FIGURES.out.figures
        }
    }

    // Run FLAIR3 transcript assembly when 'flair' or 'both' is selected.
    // Initialise empty output channels so the workflow also supports
    // runs where FLAIR3 is disabled.

    flair_gtf_ch = channel.empty()
    flair_bed_ch = channel.empty()
    flair_fasta_ch = channel.empty()

    if (useFlair) {
        FLAIR_TRANSCRIPTOME(
            INDEX_BAM.out.indexed_bam,
            reference_fasta_ch,
            reference_gtf_ch
        )
        flair_gtf_ch = FLAIR_TRANSCRIPTOME.out.isoform_gtf
        flair_bed_ch = FLAIR_TRANSCRIPTOME.out.isoform_bed
        flair_fasta_ch = FLAIR_TRANSCRIPTOME.out.isoform_fasta
    }

    unguided_gtf_ch = channel.empty()
    merged_gtf_ch = channel.empty()
    requantified_ch = channel.empty()
    combined_gtf_ch = channel.empty()

    if (useStringtie) {
        // Per-sample unguided transcript assembly.
        STRINGTIE_UNGUIDED(
            MINIMAP2.out.sorted_bam
        )

        // Collect sample assemblies for a single merge.
        // The merge module sorts GTF filenames.
        unguided_gtfs_ch = STRINGTIE_UNGUIDED.out.unguided_gtf
            .map { _sample, gtf -> gtf }
            .collect()

        STRINGTIE_MERGE(
            unguided_gtfs_ch,
            reference_gtf_ch
        )

        // Requantify each sample using the merged annotation.
        STRINGTIE_REQUANT(
            MINIMAP2.out.sorted_bam,
            STRINGTIE_MERGE.out.merged_gtf.first()
        )

        requant_gtfs_ch = STRINGTIE_REQUANT.out.requantified
            .map { _sample, gtf, _abundance -> gtf }
            .collect()

        COMBINE_REQUANT_GTFS(requant_gtfs_ch)
        unguided_gtf_ch = STRINGTIE_UNGUIDED.out.unguided_gtf
        merged_gtf_ch = STRINGTIE_MERGE.out.merged_gtf
        requantified_ch = STRINGTIE_REQUANT.out.requantified
        combined_gtf_ch = COMBINE_REQUANT_GTFS.out.combined_gtf
    }

    // Optional gene-level transcript characterisation using StringTie2
    // and/or FLAIR3 assemblies, depending on the selected assembler.
    transcript_summaries_ch = channel.empty()
    single_sample_figures_ch = channel.empty()
    all_samples_figures_ch = channel.empty()

    if (params.get('rna_characterisation') == true) {

        def configuredCases = params.get('cases') ?: []

        def characterisationGene =
            params.get('rna_characterisation_gene') ?:
            params.get('target_gene') ?:
            (configuredCases.size() > 0 ? configuredCases[0].gene : null)

        if (!characterisationGene) {
            error 'RNA characterisation requires --rna_characterisation_gene or --target_gene'
        }

        // STRINGTIE_REQUANT emits (sample, GTF, abundance).
        // Each sample must be characterised using its own GTF.
        def summary_inputs_ch = channel.empty()

        // Generate single-sample and multi-sample figures from StringTie2 summaries.
        // FLAIR3 summaries are retained separately but are not used by these plots.
        if (useStringtie) {
            summary_inputs_ch = STRINGTIE_REQUANT.out.requantified
                .map { sample, gtf, _a ->
                    tuple(sample, 'stringtie', gtf)
                }
        }
        if (useFlair) {
            def f = FLAIR_TRANSCRIPTOME.out.isoform_gtf
                .map { sample, gtf ->
                    tuple(sample, 'flair', gtf)
                }
            summary_inputs_ch = useStringtie
                ? summary_inputs_ch.mix(f) : f
        }

        EXTRACT_GENE_TRANSCRIPT_SUMMARY(
            summary_inputs_ch,
            reference_gtf_ch,
            channel.value(file(
                "${projectDir}/scripts/rna/extract_gene_transcript_summary.py",
                checkIfExists: true
            )),
            characterisationGene.toString()
        )

        transcript_summaries_ch =
            EXTRACT_GENE_TRANSCRIPT_SUMMARY.out.summary

        if (useStringtie) {
            def stringtie_summaries_ch = transcript_summaries_ch
                .filter { sample, method, tsv ->
                    method == 'stringtie'
                }
            // Preserve the correct sample-to-summary pairing.
            def single_plot_inputs_ch = stringtie_summaries_ch
                .map { sample, _assembler, summary ->
                    tuple(sample, summary)
                }

            PLOT_SINGLE_SAMPLE(
                single_plot_inputs_ch,
                channel.value(file(
                    "${projectDir}/scripts/rna/plot_single_sample.R",
                    checkIfExists: true
                )),
                params.get('mane_transcript_id') ?: ''
            )

            single_sample_figures_ch = PLOT_SINGLE_SAMPLE.out.figures

            // One header followed by every sample's transcript records.
            def outputRoot =
                params.get('results_dir') ?: params.get('outdir') ?: 'results'

            // Combine sample summaries with consistent line endings.
            // Validate column structure and retain exactly one header.
            def combined_summary_ch = stringtie_summaries_ch
                .map { sample, _assembler, summary ->
                    tuple(sample.toString(), summary.text.readLines())
                }
                .collect(flat: false)
                .map { entries ->
                    if (!entries) {
                        error 'No transcript summaries available'
                    }

                    def ordered = entries.sort { a, b -> a[0] <=> b[0] }
                    def header = ordered[0][1][0]
                    def ncols = header.split('\\t', -1).size()
                    def records = []

                    ordered.each { item ->
                        def sample = item[0]
                        def lines = item[1]

                        if (!lines || lines[0] != header) {
                            error "Transcript summary header mismatch: ${sample}"
                        }

                        lines.drop(1).findAll { it.trim() }.each { line ->
                            if (line.split('\\t', -1).size() != ncols) {
                                error "Invalid transcript summary columns: ${sample}"
                            }
                            records << line
                        }
                    }

                    ([header] + records).join('\n') + '\n'
                }
                .collectFile(
                    name: 'combined_gene_transcript_summary.tsv',
                    sort: false,
                    storeDir: "${outputRoot}/rna/characterisation/summaries"
                )

            PLOT_ALL_SAMPLES(
                combined_summary_ch,
                channel.value(file(
                    "${projectDir}/scripts/rna/plot_all_samples.R",
                    checkIfExists: true
                ))
            )

            all_samples_figures_ch = PLOT_ALL_SAMPLES.out.figures
        }

    }

    // Generate the final HTML report when StringTie2 is enabled.
    // Reporting uses the combined StringTie2 GTF; FLAIR3-only runs
    // do not currently produce this report.
    report_files_ch = channel.empty()
    if (useStringtie) {
        FINAL_HTML_REPORT(
            cases_ch,
            COMBINE_REQUANT_GTFS.out.combined_gtf.first(),
            mane_transcripts_ch,
            clinvar_variants_ch,
            report_script_ch
        )
        report_files_ch = FINAL_HTML_REPORT.out.report_files
    }

    emit:
    nanostat_stats = NANOSTAT.out.stats
    trimmed_reads = reads_for_alignment_ch
    pychopper_stats = pychopper_stats_ch
    pychopper_unclassified = pychopper_unclassified_ch
    pychopper_rescued = pychopper_rescued_ch
    pychopper_reports = pychopper_report_ch
    sorted_bam = MINIMAP2.out.sorted_bam
    unguided_gtf = unguided_gtf_ch
    merged_gtf = merged_gtf_ch
    requantified = requantified_ch
    combined_gtf = combined_gtf_ch
    junction_counts = junction_counts_ch
    junction_chains = junction_chains_ch
    junction_summaries = junction_summaries_ch
    junction_comparisons = junction_comparisons_ch
    junction_annotations = junction_annotations_ch
    junction_relative_usage = junction_relative_usage_ch
    junction_figures = junction_figures_ch
    transcript_summaries = transcript_summaries_ch
    single_sample_figures = single_sample_figures_ch
    all_samples_figures = all_samples_figures_ch
    reports = report_files_ch
    flair_gtf = flair_gtf_ch
    flair_bed = flair_bed_ch
    flair_fasta = flair_fasta_ch
}
