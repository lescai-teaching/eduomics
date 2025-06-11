process SUBVAR {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_tabix:ef572dd1df9a29c1' :
        'community.wave.seqera.io/library/htslib_tabix:9596dbe66ce87412'}"

    input:
    tuple val(meta), path(vcf)
    path bed

    output:
    tuple val(meta), path('*_benign.vcf')        , optional: true, emit: benign_vcf
    tuple val(meta), path('*_pathogenic.vcf')    , optional: true, emit: pathogenic_vcf
    tuple val(meta), path('*_ontarget.vcf.gz')   , optional: true, emit: selected_vcf
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tabix -p vcf ${vcf}
    tabix -h -R ${bed} ${vcf} | bgzip -c > ${prefix}_ontarget.vcf.gz
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Benign" -e "_benign" | grep -e "#" -e "multiple_submitters" >${prefix}_benign.vcf
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Pathogenic" -e "_pathogenic" | grep -e "#" -e "multiple_submitters" | grep -e "#" -e "nonsense" >${prefix}_pathogenic.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_benign.vcf
    touch ${prefix}_pathogenic.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """
}
