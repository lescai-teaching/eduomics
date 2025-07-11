process SUBSETGFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/lescailab/gff-parsing:r-base-4.3.2_bioconductor-annotationdbi-1.62.2_bioconductor-org.hs.eg.db-3.17.0--d87518872d6a6346':
        'ghcr.io/lescailab/gff-parsing:r-base-4.3.2_bioconductor-annotationdbi-1.62.2_bioconductor-org.hs.eg.db-3.17.0--2784b1cf02b4cfd3' }"

    input:
    val(meta)
    path(gff3)

    output:
    tuple val(meta), path("filtered_annotation.gff3")        , emit: filtered_annotation
    tuple val(meta), path("valid_gene_lists.rds")            , emit: geneLists
    tuple val(meta), path("list_gene_association.tsv")       , emit: genes_list_association
    tuple val(meta), path("transcript_data.rds")             , emit: transcript_data
    tuple val(meta), path("subsetgff_parsing_log.txt")       , emit: subsetgff_parsing_log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    subset_gff.R ${meta.chromosome} ${meta.simthreshold} ${gff3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-annotationdbi: \$(Rscript -e "cat(as.character(packageVersion('AnnotationDbi')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-igraph: \$(Rscript -e "cat(as.character(packageVersion('igraph')))")
        r-httr: \$(Rscript -e "cat(as.character(packageVersion('httr')))")
        r-jsonlite: \$(Rscript -e "cat(as.character(packageVersion('jsonlite')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
        r-matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch filtered_annotation.gff3
    touch valid_gene_lists.rds
    touch list_gene_association.tsv
    touch filtered_transcript_data.rds
    touch subsetgff_parsing_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-annotationdbi: \$(Rscript -e "cat(as.character(packageVersion('AnnotationDbi')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-igraph: \$(Rscript -e "cat(as.character(packageVersion('igraph')))")
        r-httr: \$(Rscript -e "cat(as.character(packageVersion('httr')))")
        r-jsonlite: \$(Rscript -e "cat(as.character(packageVersion('jsonlite')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
        r-matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
    END_VERSIONS
    """
}
