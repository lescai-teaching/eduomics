process COUNTMATRICES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:a28f97e8be230fb0':
        'community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:285d583c318ea12d' }"

    input:
    tuple val(meta),  path(filtered_txfasta)
    tuple val(meta2), path(transcriptData)
    tuple val(meta3), path(geneLists)

    output:
    tuple val(meta), path ("countMatrix_*.rds")    , emit: simcountMatrix
    tuple val(meta), path ("expected_*.rds")       , emit: simAnnords
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def simreps = meta.reps ?: 3
    def simgroups = meta.groups ?: 2

    """
    count_matrices.R \\
        '${meta.coverage}' \\
        '${simreps}' \\
        '${simgroups}' \\
        '${filtered_txfasta}' \\
        '${transcriptData}' \\
        '${geneLists}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch countMatrix_stub.rds
    touch expected_stub.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """
}
