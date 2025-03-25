process SUBSETGFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    val(meta)
    path(gff3)

    output:
    path "valid_gene_lists.rds"     , emit: geneLists
    path "transcript_data.rds"      , emit: transcriptData
    path "parsing_log.txt"          , emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    subsetgff.R ${meta.chromosome} ${gff3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-annotationdbi: \$(Rscript -e "cat(as.character(packageVersion('AnnotationDbi')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-igraph: \$(Rscript -e "cat(as.character(packageVersion('igraph')))")
        r-s3: \$(Rscript -e "cat(as.character(packageVersion('s3')))")
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
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
        bioconductor-annotationdbi: \$(Rscript -e "cat(as.character(packageVersion('AnnotationDbi')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-igraph: \$(Rscript -e "cat(as.character(packageVersion('igraph')))")
        r-s3: \$(Rscript -e "cat(as.character(packageVersion('s3')))")
        r-httr: \$(Rscript -e "cat(as.character(packageVersion('httr')))")
        r-jsonlite: \$(Rscript -e "cat(as.character(packageVersion('jsonlite')))")
        bioconductor-clusterprofiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
        r-matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
    END_VERSIONS
    """
}
