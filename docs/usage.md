# eduomics: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/eduomics/usage](https://nf-co.re/eduomics/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The eduomics pipeline is designed to create realistic, educational genomic datasets that support problem-based learning through structured storylines. Rather than generating random data, eduomics produces biologically plausible datasets tailored to guide students through complete bioinformatics workflows, from raw sequencing data to interpretation of the results. This guide will walk you through setting up and running simulations for both DNA variant calling and RNA differential expression scenarios.

## Quick Start Tutorial

### Prerequisites

Before running the pipeline, ensure you have:

1. **Nextflow installed** (version ≥25.04.4)
2. **Container system** (Docker, Singularity, or Conda)
3. **Reference genome configured** (we recommend using `GATK.GRCh38`)

### Understanding the Input Samplesheet

The pipeline uses a CSV samplesheet to define simulation parameters. Here's the structure:

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
```

**Column Descriptions:**

- `id`: Unique identifier for your simulation (e.g. `simulation1`)
- `type`: Either `dna` for variant calling simulations or `rna` for differential expression simulations.
- `chromosome`: Target chromosome (e.g., `chr22`, `chr1`).
- `coverage`: Sequencing depth. Default is `30` for RNA-seq and `100` for DNA.
- `capture`: BED file URL for DNA capture regions (leave empty for RNA simulations).
- `reps`: Number of biological replicates per group.
- `groups`: Number of experimental groups (typically 2 for case/control).
- `simthreshold`: Simulation threshold for gene selection (leave empty for DNA simulations). Default is `0.3`.

## DNA Variant Simulation Tutorial

### Creating a DNA Simulation

1. **Prepare your samplesheet** (`dna_samplesheet.csv`):

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
my_dna_sim,dna,chr22,100,path/to/capture.bed,1,2,,
```

**Parameter Explanation:**

- `my_dna_sim`: Your simulation name.
- `dna`: DNA simulation mode.
- `chr22`: Focus on chromosome 22 (computationally efficient for teaching).
- `100`: 100x coverage (typical for exome sequencing).
- `capture`: path to exome capture BED file defining target regions.
- `1`: 1 replicate per group.
- `2`: 2 groups (e.g. normal vs. disease).
- Empty simthreshold field (not needed for DNA).

2. **Run the DNA simulation**:

```bash
nextflow run nf-core/eduomics \
    -profile docker \
    --input dna_samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir dna_results
```

### Overview of the DNA Simulation Workflow

1. **Reference Preparation**: Subsets the reference genome and annotation to your target chromosome and capture regions.
2. **Variant Selection**: Extracts pathogenic variants from ClinVar database within your target regions.
3. **Profile Generation**: Creates a sequencing profile from existing BAM files to ensure realistic read simulations.
4. **Read Simulation**: Uses SimuSCoP to generate paired-end FASTQ files containing the injected variants.
5. **Validation**: Performs variant calling to verify that injected variants can be detected.
6. **Scenario Generation**: Creates AI-powered educational scenarios explaining the biological context.

### Expected DNA Output Structure

```
dna_results/
├── dna_simulations/
│   └── my_dna_sim/
│       ├── chr22-12345-A-T/                 # Variant-specific folder
│       │   ├── normal_1.fq.gz               # Normal sample reads
│       │   ├── normal_2.fq.gz
│       │   ├── disease_1.fq.gz              # Disease sample reads
│       │   ├── disease_2.fq.gz
│       │   ├── simulated_validated.vcf      # Validated VCF
│       │   ├── solution_chr22-12345-A-T.txt # Variant-specific solution
│       │   └── my_dna_sim_scenario.txt      # Educational scenario
│       └── references/                      # Reference bundle
│           ├── reference.fa
│           ├── capture_regions.bed
│           ├── known_variants.vcf
│           └── bwa_index/
```

## RNAseq and Differential Expression Analysis Simulation Tutorial

### Creating an RNA Simulation

1. **Prepare your samplesheet** (`rna_samplesheet.csv`):

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
my_rna_sim,rna,chr22,30,,3,2,0.3
```

**Parameter Explanation:**

- `my_rna_sim`: Your simulation name.
- `rna`: RNA simulation mode.
- `chr22`: Focus on chromosome 22 (computationally efficient for teaching).
- `30`: 30x coverage (typical for RNA-seq).
- Empty capture field (not needed for RNA).
- `3`: 3 replicates per group.
- `2`: 2 groups (e.g. normal vs. disease).
- `0.3`: Jaccard index threshold used to construct a similarity network from Gene Ontology (GO) annotations.

2. **Run the RNA simulation**:

```bash
nextflow run nf-core/eduomics \
    -profile docker \
    --input rna_samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir rna_results
```

### Overview of the RNA Simulation Workflow

1. **Reference Preparation**: Subsets the reference genome and annotation to your target chromosome.
2. **Gene Selection**: Identifies genes suitable for differential expression based on Gene Ontology (GO) annotation.
3. **Count Matrix Generation**: Creates realistic count matrices with known differential expression patterns.
4. **Read Simulation**: Uses Polyester to generate FASTQ files matching the count matrices.
5. **Expression Quantification**: Runs Salmon to quantify transcript expression.
6. **Differential Analysis**: Performs DESeq2 analysis to identify differentially expressed (DE) genes.
7. **Functional Enrichment**: Conducts GO enrichment analysis on DE genes.
8. **Validation**: Ensures the simulation produces detectable differential expression and realistic enriched pathways.
9. **Scenario Generation**: Creates AI-powered educational scenarios explaining the biological context.

### Expected RNA Output Structure

```
rna_results/
├── rna_simulations/
│   └── my_rna_sim/
│       ├── GENE1_GENE2_GENE3_GENE4_GENE5/  # Top 5 DE genes folder
│       │   ├── validated_reads/
│       │   │   ├── sample_01_1.fasta.gz
│       │   │   ├── sample_01_2.fasta.gz
│       │   │   ├── sample_02_1.fasta.gz
│       │   │   ├── sample_02_2.fasta.gz
│       │   │   └── ...
│       │   ├── deseq2_results.tsv          # DE analysis results
│       │   ├── deseq2_de_genes.txt         # List of DE genes
│       │   ├── deseq2_ma_plot.pdf          # MA plot
│       │   ├── deseq2_dispersion_plot.pdf  # Dispersion plot
│       │   ├── deseq2_count_plot.pdf       # Count plot
│       │   ├── deseq2_heatmap_plot.pdf     # Heatmap
│       │   ├── deseq2_pca_plot.pdf         # PCA plot
│       │   ├── enrichment_results.rds      # GO enrichment
│       │   ├── dotplot_*.png               # Dotplot for BP, MF and CC
|       |   |__ cnetplot_*.png              # Cnetplot for BP, MF and CC
│       │   ├── validation_result.txt       # Validation status
│       │   └── my_rna_sim_scenario.txt     # Educational scenario
│       └── references/                     # Reference bundle
│           ├── transcripts.fa
│           ├── annotation.gff3
│           ├── genome.fa
|           |── tx2gene.tsv
│           └── salmon_index/
```

## Using the Simulated Data for Teaching

### DNA Simulations

The generated data can be used to teach:

1. **Quality Control**: FastQC analysis of the FASTQ files.
2. **Read Alignment**: BWA-MEM alignment to reference genome.
3. **Variant Calling**: GATK HaplotypeCaller workflow.
4. **Variant Annotation**: Using tools like SnpEff.
5. **Clinical Interpretation**: Analysis of the pathogenicity of detected variants.

### RNA Simulations

The generated data can be used to teach:

1. **Quality Control**: FastQC and MultiQC analysis.
2. **Quantification**: Salmon transcript quantification.
3. **Differential Expression and results visualisation**: DESeq2 analysis and results interpretation.
4. **Functional Analysis**: GO enrichment and pathway analysis.

## AI-powered Educational Scenarios

The Gemini Flash API is used to generate plausible clinical descriptions that offer biological context and realistic use cases for students to analyse the data. Each simulation provides:

- **Biological and clinical context** derived from variant and gene expression patterns.
- **Real-world clinical relevance** of the findings.
- **Engaging learning objectives** for students to analyse omics data within a structured problem-solving approach and connect molecular alterations with phenotypic outcomes.

These scenarios help in contextualising the simulations within a biologically and clinically relevant framework.

## Advanced Configuration

### Multiple Simulations

You can run multiple simulations in a single samplesheet:

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
dna_easy,dna,chr22,100,path/to/capture.bed,1,2,0.1,
dna_hard,dna,chr22,50,path/to/capture.bed,2,3,0.5,
rna_basic,rna,chr22,30,,3,2,0.2
rna_complex,rna,chr1,50,,5,3,0.4
```

> ⚠️ **Current limitation:**
>
> Defining multiple chromosomes per simulation type in the same CSV is **not currently supported** and will cause the workflow to fail.
> This feature will be implemented in a future release.

### Adjusting Simulation Complexity

**For Beginners:**

- Use `chr22` (smaller chromosome)
- Lower coverage (30-50x)
- Fewer replicates (2-3)

**For Advanced Users:**

- Use larger chromosomes (`chr1`, `chr2`)
- Higher coverage (100x+)
- More replicates (5+)

## Test Run & Debug

**Test Run:**

```bash
# Always test with the provided test dataset first:
nextflow run nf-core/eduomics -profile test,docker --outdir test_results
```

### Validation Failures

If simulations fail validation:

1. **Check the coverage**: Ensure sufficient coverage for variant detection (DNA) or expression quantification (RNA).
2. **Verify the capture region**: Ensure BED file format is correct and contains target regions.
3. **Adjust the similarity threshold**: Modify the `simthreshold` value to make similarity criteria more or less stringent.

## Best Practices

1. **Start Small**: Begin with small chromosomes (e.g. chr22) before moving to larger chromosomes.
2. **Test First**: Always run the test profile before your custom simulations.
3. **Plan Storage**: Simulations can generate several GB of data per sample.
4. **Document Parameters**: Keep track of simulation parameters for reproducibility.
5. **Validate Results**: Check that simulations meet your educational objectives.

## Getting Help

- **Pipeline Issues**: Open an issue on the [GitHub repository](https://github.com/nf-core/eduomics/issues).
- **Educational Applications**: Contact the development team for teaching-specific guidance.

## Next Steps

After running your simulations:

1. **Review the output structure** (see [output documentation](output.md)).
2. **Design your teaching workflow** using the generated data.
3. **Create assessment materials** based on the known ground truth.
4. **Share your educational scenarios** with the community.

The eduomics pipeline provides a foundation for problem-based learning where students can validate their analytical skills against known biological truth.
