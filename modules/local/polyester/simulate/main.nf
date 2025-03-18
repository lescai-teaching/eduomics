process POLYESTER_SIMULATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e2c5eb02f27d03de':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-polyester:e85f578a0b8f70ca' }"

    input:
    tuple val(meta), path(countmatrix)
    tuple val(meta), path(foldchange)
    path(txfasta)
    val(reps)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def simulation_args = countmatrix ? "readmat=${countmatrix}" : "fold_changes=${foldchange}"
    def simulation_function = countmatrix ? "simulate_experiment_countmat" : "simulate_experiment"

    """
    #!/usr/bin/env R

    library(polyester)
    library(Biostrings)

    ${simulation_function}(
        "${txfasta}",
        reps = c(${reps}, ${reps}),
        ${simulation_args},
        output_dir = "."
    )

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
