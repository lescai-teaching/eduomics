process ENRICHMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_r-tidyverse:7874040bf6eb56b5':
        'community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_r-tidyverse:dad2b1c57b5305d6' }"

    input:
    tuple val(meta),  path(resdata)
    tuple val(meta2), path(tx2gene)

    output:
    tuple val(meta), path("enrichment_results.rds"), path("*.png")    , emit: enrichment_results
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    enrichment.R \\
        '${resdata}' \\
        '${tx2gene}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch enrichment_results.rds
    for ont in BP MF CC; do
        touch dotplot_\${ont}.png
        touch cnetplot_\${ont}.png
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
    END_VERSIONS
    """
}
