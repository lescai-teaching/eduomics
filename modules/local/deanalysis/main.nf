process DESEQ2_SALMON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:fccda8ad2d4d4c2d':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:3db5891c323469dc' }"

    input:
    tuple val(meta), path(quant_dirs)
    path tx2gene

    output:
    tuple val(meta), path("${prefix}_deseq2_results"), emit: results
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript ${baseDir}/bin/de_analysis.R \\
        '${prefix}' \\
        '${meta.reps}' \\
        '${meta.groups}' \\
        '${tx2gene}' \\
        '${quant_dirs.join(',')}' \\
        '${prefix}_deseq2_results'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_deseq2_results
    touch ${prefix}_deseq2_results/deseq2_results.tsv
    touch ${prefix}_deseq2_results/ma_plot.pdf
    touch ${prefix}_deseq2_results/dispersion_plot.pdf
    touch ${prefix}_deseq2_results/count_plot.pdf
    touch ${prefix}_deseq2_results/heatmap_plot.pdf
    touch ${prefix}_deseq2_results/pca_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
