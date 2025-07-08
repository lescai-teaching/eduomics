# eduomics: Output

## Introduction

This document describes the output produced by the pipeline. The eduomics pipeline generates educational genomic and transcriptomic datasets with known ground truth, organized into clearly structured directories for easy use in teaching bioinformatics workflows.

## Pipeline Output Directory Structure

The pipeline organizes outputs into separate directories based on simulation type:

```
results/
├── pipeline_info/                    # Pipeline execution information
├── dna_simulations/                  # DNA variant simulation results
│   └── [simulation_id]/
│       ├── [variant_name]/           # Simulated data with specific variant
│       └── references/               # Reference files for analysis
└── rna_simulations/                  # RNA differential expression results
    └── [simulation_id]/
        ├── [gene_set_name]/          # Simulated data with specific gene set
        └── references/               # Reference files for analysis
```

## DNA Simulation Outputs

### DNA Simulation Data Directory: `dna_simulations/[simulation_id]/[variant_name]/`

Each DNA simulation creates a directory named after the simulated pathogenic variant (e.g., `chr22-12345-A-T`). This directory contains:

#### Simulated Sequencing Data

- **`normal_1.fq.gz`** and **`normal_2.fq.gz`**: Paired-end FASTQ files for the normal/control sample
- **`disease_1.fq.gz`** and **`disease_2.fq.gz`**: Paired-end FASTQ files for the disease/case sample containing the pathogenic variant

#### Validation and Solution Files

- **`solution_[variant_name].txt`**: Contains the exact variant that was simulated (the "answer key")
- **`simulated_validated.vcf.gz`**: VCF file containing the validated simulated variant
- **`[simulation_id]_scenario.txt`**: AI-generated educational scenario providing biological context for the simulation

### DNA Reference Bundle: `dna_simulations/[simulation_id]/references/`

Contains all reference files needed to analyze the simulated DNA data:

#### Reference Genome Files

- **`[chromosome].fa`**: Subset reference genome FASTA file (e.g., chr22.fa)
- **`[chromosome].fa.fai`**: FASTA index file
- **`[chromosome].dict`**: Sequence dictionary for GATK tools

#### Capture/Target Regions

- **`[simulation_id].capture.bed.gz`**: Compressed BED file defining capture regions
- **`[simulation_id].capture.bed.gz.tbi`**: Tabix index for capture regions
- **`[simulation_id]_[chromosome]_target.bed`**: Target regions for the specific chromosome
- **`[simulation_id]_[chromosome]_target.pad50.bed`**: Target regions with 50bp padding
- **`[simulation_id]_[chromosome]_target.pad500.bed`**: Target regions with 500bp padding

#### Alignment Index

- **`bwa/`**: Directory containing BWA alignment index files
  - `[chromosome].amb`, `[chromosome].ann`, `[chromosome].bwt`, `[chromosome].pac`, `[chromosome].sa`

#### Known Variant Databases (Subset to Target Regions)

- **`[simulation_id]_gnomad.vcf.gz`** and **`[simulation_id]_gnomad.vcf.gz.tbi`**: gnomAD variants
- **`[simulation_id]_mills.vcf.gz`** and **`[simulation_id]_mills.vcf.gz.tbi`**: Mills indels
- **`[simulation_id]_1000g.vcf.gz`** and **`[simulation_id]_1000g.vcf.gz.tbi`**: 1000 Genomes variants
- **`[simulation_id]_dbsnp.vcf.gz`** and **`[simulation_id]_dbsnp.vcf.gz.tbi`**: dbSNP variants

## RNA Simulation Outputs

### RNA Simulation Data Directory: `rna_simulations/[simulation_id]/[gene_set_name]/`

Each RNA simulation creates a directory named after the first 5 differentially expressed genes (e.g., `GENE1_GENE2_GENE3_GENE4_GENE5`). This directory contains:

#### Simulated RNA-seq Data

- **`validated_reads/`**: Directory containing simulated RNA-seq FASTQ files
  - **`sample_01_1.fasta.gz`**, **`sample_01_2.fasta.gz`**: Paired-end reads for replicate 1, group 1
  - **`sample_02_1.fasta.gz`**, **`sample_02_2.fasta.gz`**: Paired-end reads for replicate 2, group 1
  - **`sample_0X_1.fasta.gz`**, **`sample_0X_2.fasta.gz`**: Additional replicates based on simulation parameters
  - Pattern continues for all replicates and groups as specified in the samplesheet

#### Differential Expression Analysis Results

- **`deseq2_results.tsv`**: Complete DESeq2 differential expression results table
- **`deseq2_de_genes.txt`**: List of significantly differentially expressed genes
- **`deseq2_tx2gene.tsv`**: Transcript-to-gene mapping file

#### Visualization and Quality Control

- **`deseq2_ma_plot.pdf`**: MA plot showing log fold change vs. mean expression
- **`deseq2_dispersion_plot.pdf`**: Dispersion estimates plot
- **`deseq2_count_plot.pdf`**: Normalized count plots for top genes
- **`deseq2_heatmap_plot.pdf`**: Heatmap of differentially expressed genes
- **`deseq2_pca_plot.pdf`**: Principal component analysis plot

#### Functional Enrichment Analysis

- **`enrichment_results.rds`**: R data file containing complete enrichment analysis results
- **`dotplot_BP.png`**: Dot plot for Biological Process GO terms
- **`dotplot_MF.png`**: Dot plot for Molecular Function GO terms
- **`dotplot_CC.png`**: Dot plot for Cellular Component GO terms
- **`cnetplot_BP.png`**: Concept network plot for Biological Process terms
- **`cnetplot_MF.png`**: Concept network plot for Molecular Function terms
- **`cnetplot_CC.png`**: Concept network plot for Cellular Component terms

#### Validation and Educational Materials

- **`validation_result.txt`**: Validation status of the simulation ("GOOD SIMULATION" if passed)
- **`[simulation_id]_scenario.txt`**: AI-generated educational scenario explaining the biological context

### RNA Reference Bundle: `rna_simulations/[simulation_id]/references/`

Contains all reference files needed to analyze the simulated RNA-seq data:

- **`gencode_transcripts_novers_[chromosome].fasta`**: Subset transcriptome FASTA file
- **`filtered_annotation.gff3`**: Subset gene annotation file for the target chromosome
- **`[chromosome].fa`**: Subset genome FASTA file
- **`salmon/`**: Directory containing Salmon index files for transcript quantification

## Pipeline Information: `pipeline_info/`

Contains technical information about the pipeline execution:

- **`execution_timeline_[timestamp].html`**: Timeline of process execution
- **`execution_report_[timestamp].html`**: Detailed execution report
- **`execution_trace_[timestamp].txt`**: Trace of all executed processes
- **`pipeline_dag_[timestamp].html`**: Directed acyclic graph of the pipeline
- **`nf_core_eduomics_software_versions.yml`**: Versions of all software tools used

## Understanding the Output Structure

### DNA Simulation Use Case

The DNA simulation outputs are designed for teaching variant calling workflows:

1. **Students receive**: Simulated FASTQ files and reference bundle
2. **Students perform**: Read alignment, variant calling, and annotation
3. **Validation**: Checks that the simulated pathogenic variant is indeed present in the reads of the "disease" sample and absent from the "normal" sample, ensuring that the simulation is suitable for teaching variant calling workflows.
4. **Educational context**: Provided through AI-generated scenarios

### RNA Simulation Use Case

The RNA simulation outputs are designed for teaching differential expression and enrichment analysis:

1. **Students receive**: Simulated RNA-seq FASTQ files and reference bundle
2. **Students perform**: Read quantification, differential expression analysis, and functional enrichment
3. **Validation**: Verifies that the simulated data produce a statistically significant and coherent differential expression pattern, with detectable and biologically plausible functional enrichment.
4. **Educational context**: Provided through AI-generated scenarios and visualization

### File Naming Conventions

- **Simulation ID**: User-defined identifier from the samplesheet
- **Variant Name**: Format `chr[N]-[position]-[ref]-[alt]` (e.g., `chr22-12345-A-T`)
- **Gene Set Name**: First 5 differentially expressed genes joined by underscores
- **Sample Names**: Follow pattern `sample_[group][replicate]_[read].fasta.gz`

### Quality Control and Validation

The pipeline includes built-in validation steps:

- **DNA simulations**: Verified that the injected variant is detectable in the simulated reads
- **RNA simulations**: Validated that the differential expression pattern meets quality thresholds
- Only simulations passing validation are included in the final output

This comprehensive output structure ensures that educators have everything needed to create effective bioinformatics learning experiences, while students receive realistic datasets with known ground truth for validation of their analytical approaches.
