process DEANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:fccda8ad2d4d4c2d':
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-pheatmap_r-tidyverse:3db5891c323469dc' }"

    input:
    tuple val(meta), path(quant_dirs)
    tuple val(meta2), path(transcriptData)

    output:
    tuple val(meta), path("deseq2_results.tsv"), path("deseq2_de_genes.txt"), path("*.pdf")    , emit: deseq2_results
    tuple val(meta), path("deseq2_tx2gene.tsv")                                                , emit: deseq2_tx2gene
    path "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/de_analysis.R \\
        '${meta.reps}' \\
        '${meta.groups}' \\
        '${transcriptData}' \\
        '${quant_dirs.join(',')}'

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
    touch deseq2_results.tsv
    touch deseq2_tx2gene.tsv
    touch deseq2_de_genes.txt
    touch deseq2_ma_plot.pdf
    touch deseq2_dispersion_plot.pdf
    touch deseq2_count_plot.pdf
    touch deseq2_heatmap_plot.pdf
    touch deseq2_pca_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        bioconductor-tximport: \$(Rscript -e "cat(as.character(packageVersion('tximport')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
