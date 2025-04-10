process ENRICHMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_r-tidyverse:d8812438cb8ebb59':
        'community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_r-tidyverse:2a42661ad4d31ae0' }"

    input:
    tuple val(meta), path(resdata)
    path tx2gene

    output:
    tuple val(meta), path("${prefix}_enrichment_results"), emit: results
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/enrichment.R \\
        '${resdata}' \\
        '${tx2gene}' \\
        '${prefix}_enrichment_results'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_enrichment_results
    touch ${prefix}_enrichment_results/enrichment_results.rds
    for ont in BP MF CC; do
        touch ${prefix}_enrichment_results/dotplot_\${ont}.png
        touch ${prefix}_enrichment_results/cnetplot_\${ont}.png
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
    END_VERSIONS
    """
}
