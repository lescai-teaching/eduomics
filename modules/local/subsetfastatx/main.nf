process SUBSETFASTATX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:a28f97e8be230fb0':
        'community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:285d583c318ea12d' }"

    input:
    path(txfasta)
    tuple val(meta), path(transcript_data)

    output:
    tuple val(meta), path ("gencode_transcripts_novers_*.fasta"), emit: filtered_txfasta
    tuple val(meta), path ("subsetfastatx_parsing_log.txt")     , emit: subsetfastatx_parsing_log
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    subset_fastatx.R ${meta.chromosome} ${txfasta} ${transcript_data}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
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
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
