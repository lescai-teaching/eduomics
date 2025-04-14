process COUNTMATRICES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:59cdfcce7d992560':
        'community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:4645f56e7256e01f' }"

    input:
    val(meta)
    path(txfasta)
    path(gff3)
    path(geneList)

    output:
    path "${meta.id}_countMatrix_*.rds"    , emit: simcountMatrix
    path "${meta.id}_expected_*.rds"       , emit: simAnnords
    path "${meta.id}_expected_*.tsv"       , emit: simAnnotsv
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/count_matrices.R \\
        '${prefix}' \\
        '${meta.coverage}' \\
        '${meta.reps}' \\
        '${meta.groups}' \\
        '${txfasta}' \\
        '${gff3}' \\
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
    touch ${meta.id}_countMatrix_stub.rds
    touch ${meta.id}_expected_stub.rds
    touch ${meta.id}_expected_stub.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """
}
