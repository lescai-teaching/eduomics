process COUNTMATRICES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:59cdfcce7d992560':
        'community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:4645f56e7256e01f' }"

    input:
    tuple val(meta),  path(filtered_txfasta)
    tuple val(meta2), path(filtered_gff3)
    tuple val(meta3), path(geneList)

    output:
    tuple val(meta), path ("countMatrix_*.rds")    , emit: simcountMatrix
    tuple val(meta), path ("expected_*.rds")       , emit: simAnnords
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/count_matrices.R \\
        '${meta.coverage}' \\
        '${meta.reps}' \\
        '${meta.groups}' \\
        '${txfasta}' \\
        '${filtered_gff3}' \\
        '${geneList}'

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
