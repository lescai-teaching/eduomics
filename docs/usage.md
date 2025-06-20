# nf-core/eduomics: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/eduomics/usage](https://nf-co.re/eduomics/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The eduomics pipeline is designed to create realistic, educational genomic datasets for teaching bioinformatics analysis. This guide will walk you through setting up and running simulations for both DNA variant calling and RNA differential expression scenarios.

## Quick Start Tutorial

### Prerequisites

Before running the pipeline, ensure you have:

1. **Nextflow installed** (version ≥24.04.2)
2. **Container system** (Docker, Singularity, or Conda)
3. **Reference genome configured** (we recommend using `GATK.GRCh38`)

### Step 1: Understanding the Input Samplesheet

The pipeline uses a CSV samplesheet to define simulation parameters. Here's the structure:

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
```

**Column Descriptions:**

- `id`: Unique identifier for your simulation
- `type`: Either `dna` for variant calling simulations or `rna` for differential expression simulations
- `chromosome`: Target chromosome (e.g., `chr22`, `chr1`)
- `coverage`: Sequencing depth (e.g., `30` for RNA-seq, `100` for DNA)
- `capture`: BED file URL for DNA capture regions (leave empty for RNA simulations)
- `reps`: Number of biological replicates per group
- `groups`: Number of experimental groups (typically 2 for case/control)
- `simthreshold`: Simulation threshold for gene selection (0.1-0.5 recommended)

### Step 2: DNA Variant Simulation Tutorial

#### Creating a DNA Simulation

1. **Prepare your samplesheet** (`dna_samplesheet.csv`):

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
my_dna_sim,dna,chr22,100,https://raw.githubusercontent.com/lescai-teaching/eduomics_testdata/refs/heads/main/dna/whole_chr22/Twist_exome_2.0_covered_chr22_500pad.bed,1,2,0.3
```

**Parameter Explanation:**

- `my_dna_sim`: Your simulation name
- `dna`: DNA simulation mode
- `chr22`: Focus on chromosome 22 (computationally efficient for teaching)
- `100`: 100x coverage (typical for exome sequencing)
- `capture`: URL to exome capture BED file defining target regions
- `1,2`: 1 replicate per group, 2 groups (normal vs. disease)
- `0.3`: 30% of genes will be considered for variant injection

2. **Run the DNA simulation**:

```bash
nextflow run nf-core/eduomics \
    -profile docker \
    --input dna_samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir dna_results
```

#### What the DNA Simulation Does

1. **Reference Preparation**: Subsets the reference genome and annotation to your target chromosome and capture regions
2. **Variant Selection**: Extracts pathogenic variants from ClinVar database within your target regions
3. **Profile Generation**: Creates a sequencing profile from existing BAM files to ensure realistic read characteristics
4. **Read Simulation**: Uses SimuSCoP to generate paired-end FASTQ files containing the injected variants
5. **Validation**: Performs variant calling to verify that injected variants can be detected
6. **Scenario Generation**: Creates AI-powered educational scenarios explaining the biological context

#### Expected DNA Output Structure

```
dna_results/
├── dna_simulations/
│   └── my_dna_sim/
│       ├── chr22-12345-A-T/          # Variant-specific folder
│       │   ├── normal_1.fq.gz        # Normal sample reads
│       │   ├── normal_2.fq.gz
│       │   ├── disease_1.fq.gz       # Disease sample reads
│       │   ├── disease_2.fq.gz
│       │   ├── simulated_validated.vcf # Validation VCF
│       │   ├── solution_chr22-12345-A-T.txt # Answer key
│       │   └── my_dna_sim_scenario.txt # Educational scenario
│       └── references/               # Reference bundle
│           ├── reference.fa
│           ├── capture_regions.bed
│           ├── known_variants.vcf
│           └── bwa_index/
```

### Step 3: RNA Differential Expression Simulation Tutorial

#### Creating an RNA Simulation

1. **Prepare your samplesheet** (`rna_samplesheet.csv`):

```csv
id,type,chromosome,coverage,reps,groups,simthreshold
my_rna_sim,rna,chr22,30,,3,2,0.3
```

**Parameter Explanation:**

- `my_rna_sim`: Your simulation name
- `rna`: RNA simulation mode
- `chr22`: Focus on chromosome 22
- `30`: 30x coverage (typical for RNA-seq)
- Empty capture field (not needed for RNA)
- `3,2`: 3 replicates per group, 2 groups
- `0.3`: 30% of genes will show differential expression

2. **Run the RNA simulation**:

```bash
nextflow run nf-core/eduomics \
    -profile docker \
    --input rna_samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir rna_results
```

#### What the RNA Simulation Does

1. **Transcriptome Preparation**: Subsets transcriptome references to your target chromosome
2. **Gene Selection**: Identifies genes suitable for differential expression based on functional annotations
3. **Count Matrix Generation**: Creates realistic count matrices with known differential expression patterns
4. **Read Simulation**: Uses Polyester to generate RNA-seq FASTQ files matching the count matrices
5. **Expression Quantification**: Runs Salmon to quantify transcript expression
6. **Differential Analysis**: Performs DESeq2 analysis to identify differentially expressed genes
7. **Functional Enrichment**: Conducts GO enrichment analysis on differentially expressed genes
8. **Validation**: Ensures the simulation produces detectable differential expression
9. **Scenario Generation**: Creates educational scenarios explaining the biological context

#### Expected RNA Output Structure

```
rna_results/
├── rna_simulations/
│   └── my_rna_sim/
│       ├── GENE1_GENE2_GENE3_GENE4_GENE5/  # Top 5 DE genes folder
│       │   ├── validated_reads/
│       │   │   ├── sample_01_1.fasta.gz    # Group 1 replicates
│       │   │   ├── sample_01_2.fasta.gz
│       │   │   ├── sample_02_1.fasta.gz
│       │   │   ├── sample_02_2.fasta.gz
│       │   │   └── ...
│       │   ├── deseq2_results.tsv          # DE analysis results
│       │   ├── deseq2_de_genes.txt         # List of DE genes
│       │   ├── deseq2_ma_plot.pdf          # MA plot
│       │   ├── deseq2_pca_plot.pdf         # PCA plot
│       │   ├── enrichment_results.rds      # GO enrichment
│       │   ├── dotplot_BP.png              # Enrichment plots
│       │   ├── validation_result.txt       # Validation status
│       │   └── my_rna_sim_scenario.txt     # Educational scenario
│       └── references/                     # Reference bundle
│           ├── transcripts.fa
│           ├── annotation.gff3
│           └── salmon_index/
```

### Step 4: Advanced Configuration

#### Multiple Simulations

You can run multiple simulations in a single samplesheet:

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
dna_easy,dna,chr22,100,https://example.com/capture.bed,1,2,0.1
dna_hard,dna,chr22,50,https://example.com/capture.bed,2,3,0.5
rna_basic,rna,chr22,30,,3,2,0.2
rna_complex,rna,chr1,50,,5,3,0.4
```

#### Adjusting Simulation Complexity

**For Beginners:**

- Use `chr22` (smaller chromosome)
- Lower coverage (30-50x)
- Fewer replicates (2-3)
- Lower simthreshold (0.1-0.2)

**For Advanced Users:**

- Use larger chromosomes (`chr1`, `chr2`)
- Higher coverage (100x+)
- More replicates (5+)
- Higher simthreshold (0.4-0.5)

#### Custom Capture Regions

For DNA simulations, you can provide your own capture BED file:

```csv
id,type,chromosome,coverage,capture,reps,groups,simthreshold
custom_panel,dna,chr17,150,/path/to/my_panel.bed,2,2,0.3
```

### Step 5: Using the Simulated Data for Teaching

#### For DNA Simulations

The generated data can be used to teach:

1. **Quality Control**: FastQC analysis of the FASTQ files
2. **Read Alignment**: BWA-MEM alignment to reference genome
3. **Variant Calling**: GATK HaplotypeCaller workflow
4. **Variant Annotation**: Using tools like VEP or ANNOVAR
5. **Clinical Interpretation**: Analyzing the pathogenicity of detected variants

**Teaching Workflow:**

```bash
# Students can practice this workflow:
# 1. Quality control
fastqc *.fq.gz

# 2. Alignment
bwa mem reference.fa sample_1.fq.gz sample_2.fq.gz | samtools sort > sample.bam

# 3. Variant calling
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.vcf

# 4. Compare with solution
diff sample.vcf solution_chr22-12345-A-T.txt
```

#### For RNA Simulations

The generated data can be used to teach:

1. **Quality Control**: FastQC and MultiQC analysis
2. **Quantification**: Salmon or Kallisto transcript quantification
3. **Differential Expression**: DESeq2 or edgeR analysis
4. **Functional Analysis**: GO enrichment and pathway analysis
5. **Visualization**: Creating plots and heatmaps

**Teaching Workflow:**

```bash
# Students can practice this workflow:
# 1. Quantification
salmon quant -i salmon_index -l A -1 sample_1.fasta.gz -2 sample_2.fasta.gz -o sample_quant

# 2. Import to R and run DESeq2
# (R code for differential expression analysis)

# 3. Compare results with provided solution
# Compare detected DE genes with deseq2_de_genes.txt
```

### Step 6: Troubleshooting

#### Common Issues

**Memory Requirements:**

```bash
# If you encounter memory issues, add:
export NXF_OPTS='-Xms1g -Xmx4g'
```

**Container Issues:**

```bash
# For Singularity users:
nextflow run nf-core/eduomics -profile singularity --singularity_pull_docker_container

# For Conda users:
nextflow run nf-core/eduomics -profile conda
```

**Test Run:**

```bash
# Always test with the provided test dataset first:
nextflow run nf-core/eduomics -profile test,docker --outdir test_results
```

#### Validation Failures

If simulations fail validation:

1. **Check coverage**: Ensure sufficient coverage for variant detection (DNA) or expression quantification (RNA)
2. **Verify capture regions**: Ensure BED file format is correct and contains target regions
3. **Adjust thresholds**: Lower simthreshold values for more conservative simulations

### Step 7: Educational Scenarios

Each simulation generates an AI-powered educational scenario that provides:

- **Biological context** for the simulated variants or expression changes
- **Clinical relevance** of the findings
- **Learning objectives** for the dataset
- **Expected outcomes** students should achieve

These scenarios help instructors frame the computational exercise within a meaningful biological context.

## Resource Requirements

### Minimum Requirements

- **CPU**: 4 cores
- **Memory**: 8 GB RAM
- **Storage**: 50 GB free space
- **Time**: 1-4 hours depending on simulation complexity

### Recommended Requirements

- **CPU**: 8+ cores
- **Memory**: 16+ GB RAM
- **Storage**: 100+ GB free space
- **Time**: 30 minutes - 2 hours

## Best Practices

1. **Start Small**: Begin with chr22 simulations before moving to larger chromosomes
2. **Test First**: Always run the test profile before your custom simulations
3. **Plan Storage**: Simulations can generate several GB of data per sample
4. **Document Parameters**: Keep track of simulation parameters for reproducibility
5. **Validate Results**: Check that simulations meet your educational objectives

## Getting Help

- **Pipeline Issues**: Open an issue on the [GitHub repository](https://github.com/nf-core/eduomics/issues)
- **Usage Questions**: Join the [nf-core Slack](https://nf-co.re/join/slack) #eduomics channel
- **Educational Applications**: Contact the development team for teaching-specific guidance

## Next Steps

After running your simulations:

1. **Review the output structure** (see [output documentation](output.md))
2. **Design your teaching workflow** using the generated data
3. **Create assessment materials** based on the known ground truth
4. **Share your educational scenarios** with the community

The eduomics pipeline provides a foundation for evidence-based bioinformatics education where students can validate their analytical skills against known biological truth.
