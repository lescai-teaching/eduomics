process SUBSETGFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-biostrings_bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_pruned:35831cd366df77bf':
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-biostrings_bioconductor-clusterprofiler_bioconductor-org.hs.eg.db_pruned:ced2353e22e6ba1f' }"

    input:
    val(meta)
    path(gff3)

    output:
    tuple val(meta), path("${meta.id}_filtered_annotation_novers_${meta.chromosome}.gff3"),    emit: filtered_annotation
    tuple val(meta), path("valid_gene_lists.rds"),                                             emit: geneLists
    tuple val(meta), path("list_gene_association.tsv"),                                        emit: genes_list_association
    tuple val(meta), path("transcript_data.rds"),                                              emit: transcript_data
    tuple val(meta), path("${meta.id}_tx2gene_${meta.chromosome}.tsv"),                        emit: tx2gene
    tuple val(meta), path("subsetgff_parsing_log.txt"),                                        emit: subsetgff_parsing_log
    path "versions.yml",                                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    subset_gff.R ${meta.id} ${meta.chromosome} ${meta.simthreshold} ${gff3}

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
    touch ${meta.id}_filtered_annotation_novers_${meta.chromosome}.gff3
    touch valid_gene_lists.rds
    touch list_gene_association.tsv
    touch transcript_data.rds
    touch ${meta.id}_tx2gene_${meta.chromosome}.tsv
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
