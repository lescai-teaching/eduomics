process POLYESTER_SIMULATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e2c5eb02f27d03de' :
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e85f578a0b8f70ca' }"

    input:
    tuple val(meta),  path(countmatrix)
    tuple val(meta2), path(foldchange)
    tuple val(meta3), path(filtered_txfasta)

    output:
    tuple val(meta), path('simulated_reads/*.fasta.gz'), emit: reads
    path 'versions.yml'                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reps = meta2.reps ?: 3
    def groups = meta2.groups ?: 2
    def repsList = (1..groups).collect { "${reps}" }.join(',')
    def numRepsString = "num_reps=c(${repsList})"
    def readData = countmatrix ? "countmat = readRDS('${countmatrix}')" : "fold_changes = readRDS('${foldchange}')"
    def readsPerTxFactor = meta.coverage ?: 30
    def simulation_function = countmatrix ?
        "simulate_experiment_countmat('${filtered_txfasta}', readmat=countmat, outdir='simulated_reads', cores = ${task.cpus}, error_model = 'illumina5', gzip = TRUE)" :
        "simulate_experiment('${filtered_txfasta}', reads_per_transcript=readspertx, ${numRepsString}, fold_changes=fold_changes, outdir='simulated_reads', cores = ${task.cpus}, error_model = 'illumina5', gzip = TRUE)"

    """
    cat <<EOF >simulate_reads.R
    #!/usr/bin/env R

    library(polyester)
    library(Biostrings)

    # Load the fasta
    fasta <- readDNAStringSet("${filtered_txfasta}")
    readspertx = round(${readsPerTxFactor} * width(fasta) / 100)

    # Load the count matrix or fold changes
    ${readData}

    ${simulation_function}
    EOF

    Rscript simulate_reads.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polyester: \$(Rscript -e 'cat(as.character(packageVersion("polyester")))')
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p simulated_reads
    touch simulated_reads/${prefix}_1.fasta.gz
    touch simulated_reads/${prefix}_2.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polyester: \$(Rscript -e 'cat(as.character(packageVersion("polyester")))')
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """
}
