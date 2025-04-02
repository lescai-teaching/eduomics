process DEANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-tidyverse:9621ff2877385470':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:3db5891c323469dc' }"

    input:
    tuple val(meta), path(counts)
    path (transcripData)             // tibble resulting from subsetgff module

    output:
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_ma_plot.pdf")            , emit: ma_plot
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_dispersion_plot.pdf")    , emit: dispersion_plot
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_count_plot.pdf")         , emit: count_plot
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_results.tsv")            , emit: results
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_heatmap_plot.pdf")       , emit: heatmap_plot
    tuple val(meta), path("${prefix}_deanalysis/*.deseq2_pca_plot.pdf")           , emit: pca_plot
    path "versions.yml"                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/count_matrices.R ${meta.reps} ${meta.groups} ${gff3} ${gff3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('deseq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_deseq2_ma_plot.pdf
    touch ${prefix}_deseq2_dispersion_plot.pdf
    touch ${prefix}_deseq2_count_plot.pdf
    touch ${prefix}_deseq2_results.tsv
    touch ${prefix}_deseq2_heatmap_plot.pdf
    touch ${prefix}_deseq2_pca_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('deseq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
