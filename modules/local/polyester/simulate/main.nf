process POLYESTER_SIMULATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e2c5eb02f27d03de':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e85f578a0b8f70ca' }"

    input:
    tuple val(meta), path(countmatrix)
    tuple val(meta2), path(foldchange)
    path(txfasta)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reps = meta2.reps ?: meta2.reps : 3
    def groups = meta2.groups ?: meta2.groups : 2
    // Create a list that repeats 'reps' for the number of 'groups' times
    def repsList = (1..groups).collect { "${reps}" }.join(",")
    def numRepsString = "num_reps=c(${repsList})"
    def readData = countmatrix ? "countmat = readRDS('${countmatrix}')" : "fold_changes = readRDS('${foldchange}')"
    def simulation_function = (countmatrix
        ? "simulate_experiment_countmat(${txfasta}, readmat=countmat, outdir='simulated_reads')"
        : "simulate_experiment(${txfasta}, reads_per_transcript=readspertx, ${numRepsString}, fold_changes=fold_changes, outdir='simulated_reads')")

    """
    cat <<EOF >>simulate_reads.R
    #!/usr/bin/env R

    library(polyester)
    library(Biostrings)

    # Load the count matrix or fold changes
    ${readData}

    ${simulation_function}
    EOF

    Rscript simulate_reads.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polyester: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polyester: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
