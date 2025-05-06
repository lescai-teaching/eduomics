process SUBSETFASTATX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-rtracklayer_r-tidyverse:f420c00b549a4380':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-rtracklayer_r-tidyverse:e36d3b6eec6fd274' }"

    input:
    val(meta)
    path(txfasta)
    path(gff3)

    output:
    tuple val(meta), path ("gencode_transcripts_noversion.fasta")    , emit: filtered_txfasta
    tuple val(meta), path ("subsetfastatx_parsing_log.txt")          , emit: log
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/subset_fastatx.R ${meta.chromosome} ${txfasta} ${gff3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch gencode_transcripts_noversion.fasta
    touch subsetfastatx_parsing_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
