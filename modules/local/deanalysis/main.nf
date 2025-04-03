process DEANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:fccda8ad2d4d4c2d':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:3db5891c323469dc' }"

    input:
    tuple val(meta), path(quant)
    path (transcriptData)             // tibble resulting from subsetgff module

    output:
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_ma_plot.pdf")            , emit: ma_plot
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_dispersion_plot.pdf")    , emit: dispersion_plot
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_count_plot.pdf")         , emit: count_plot
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_results.tsv")            , emit: results
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_heatmap_plot.pdf")       , emit: heatmap_plot
    tuple val(meta), path("${meta.id}_deanalysis/deseq2_pca_plot.pdf")           , emit: pca_plot
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_deanalysis
    Rscript ${baseDir}/bin/de_analysis.R ${prefix} ${meta.reps} ${meta.groups} ${transcriptData} ${quant}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_deanalysis
    touch ${prefix}_deanalysis/deseq2_ma_plot.pdf
    touch ${prefix}_deanalysis/deseq2_dispersion_plot.pdf
    touch ${prefix}_deanalysis/deseq2_count_plot.pdf
    touch ${prefix}_deanalysis/deseq2_results.tsv
    touch ${prefix}_deanalysis/deseq2_heatmap_plot.pdf
    touch ${prefix}_deanalysis/deseq2_pca_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
