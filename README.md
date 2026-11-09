# NEB Nextflow

**A modular Nextflow DSL2 workflow for Oxford Nanopore long-read DNA and RNA sequencing analysis.**

NEB Nextflow was developed during a research internship at the UK Dementia Research Institute, University of Cambridge, within the Ryten Lab.

The project was motivated by the challenges associated with investigating *NEB*, a large, structurally complex gene encoding nebulin. Its extensive exon structure, repetitive regions and alternative splicing present challenges for variant interpretation and transcript characterisation.

The workflow was designed to support gene-focused genomic analysis and the investigation of internal RNA splicing structures. Although motivated by *NEB*, the implementation is configurable for other genes of interest.

## Overview

The project comprises two principal workflows implemented using Nextflow DSL2.

### DNA analysis

The DNA workflow supports Oxford Nanopore sequencing analysis through:

- Basecalling and demultiplexing using Dorado.
- Reference alignment and coverage assessment using Minimap2 and supporting tools.
- Small-variant detection using Clair3.
- Structural-variant detection using Sniffles2.
- Copy-number analysis using QDNAseq and SPECTRE.
- Variant annotation using Ensembl VEP.
- Gene-focused reporting and variant visualisation.

The resulting outputs are intended to support the identification and characterisation of candidate genomic variants and inform subsequent experimental investigation.

### RNA analysis

The RNA workflow supports long-read transcriptomic analysis, with particular emphasis on internal splice-junction structures.

Its principal capabilities include:

- Optional read preprocessing using Cutadapt or Pychopper.
- Splice-aware alignment using Minimap2.
- Transcript reconstruction using StringTie2 and FLAIR3.
- Direct extraction of splice junctions and consecutive junction chains from aligned reads.
- Reference annotation and comparison of observed junctions between samples.
- Descriptive analysis of relative junction usage at alternative splice sites.
- Transcript characterisation, visualisation and reporting.

The junction-focused analysis provides read-level evidence of internal splicing patterns, complementing the interpretation of reconstructed transcript models.

The RNA workflow builds on analytical approaches established in [RAPID](https://github.com/KylieMontgomery/RAPID), developed by Kylie-Ann Montgomery. NEB Nextflow extends this framework through additional junction-focused analysis, alternative transcript reconstruction and integration with a separate DNA analysis workflow.

## Public Datasets and Evaluation

Development and computational testing were conducted using publicly available Oxford Nanopore sequencing datasets.

The original *NEB* research-case data were not analysed with this workflow. Instead, independent public datasets were selected to evaluate specific analytical components.

| Dataset | Target | Application |
|---|---|---|
| HG00276 long-read DNA data | *CYP2D6* | DNA variant detection, annotation and gene-focused reporting |
| [GSE287954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287954): SRR32119253 and SRR32119294 | *NEB* | RNA junction extraction, reference annotation, sample comparison and transcript characterisation |
| [ERR3577082](https://www.ebi.ac.uk/ena/browser/view/ERR3577082) and [ERR3577098](https://www.ebi.ac.uk/ena/browser/view/ERR3577098) | *CACNA1C* | Junction visualisation, FLAIR3 reconstruction and comparison with StringTie2 |

### DNA evaluation

The *CYP2D6* analysis demonstrated the execution of selected DNA processing, variant-detection, annotation and reporting components.

Gene-focused outputs included small variants, structural-variant information and genomic visualisations. The evaluation primarily assessed computational functionality rather than comprehensive variant-calling accuracy.

### RNA junction analysis

The *NEB*-focused analysis of the GSE287954 samples identified 188 distinct internal splice junctions, of which 173 matched the selected Ensembl reference annotation and 15 were unannotated relative to that reference.

The workflow also generated read-supported junction chains and descriptive relative-usage measurements. No alternative splice sites met the minimum support criteria for comparable relative-usage differences between the two samples.

These results demonstrate the functionality of the junction-analysis components without establishing statistically significant differential splicing.

### Transcript reconstruction

FLAIR3 was independently evaluated using the *CACNA1C* datasets ERR3577082 and ERR3577098.

Comparison with StringTie2 identified differences in reconstructed exon structures and shared internal splice-junction chains, illustrating the influence of assembler selection on transcript reconstruction.

FLAIR3 was subsequently integrated into the main RNA workflow, allowing users to select `stringtie`, `flair` or `both` through the `rna_assembler` parameter.

The updated integration passed Nextflow preview validation. A separate end-to-end execution of the final combined-assembler configuration remains outstanding.

## Usage

The workflow was developed primarily for Linux and high-performance computing environments, including the Cambridge Service for Data Driven Discovery (CSD3).

It requires Nextflow, appropriate sequencing inputs, reference resources and the software environments associated with the selected analytical processes.

A general invocation is:

```bash
nextflow run main.nf -profile local -params-file /path/to/analysis.json
```

Execution can also be configured for Slurm-based HPC environments using the `hpc` profile.

Configuration files must be adapted to the dataset and computing environment. Some existing settings are specific to the original CSD3 development infrastructure.

## Reproducibility

The workflow includes mechanisms to support reproducible computational analyses, including recording run parameters, software information, execution metadata and input file checksums.

Supporting scripts and Nextflow configuration files are provided to facilitate analysis tracking and subsequent development.


## Acknowledgements

Developed during a research internship at the UK Dementia Research Institute, University of Cambridge, within the Ryten Lab.

The RNA workflow builds upon analytical approaches established in [RAPID](https://github.com/KylieMontgomery/RAPID), developed by Kylie-Ann Montgomery.

The project benefited from scientific guidance provided by the Ryten Lab, publicly available sequencing datasets and the open-source bioinformatics tools used throughout development.
