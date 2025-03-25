process SUBVAR {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2' :
        'biocontainers/htslib:1.20--h5efdd21_2'}"

    input:
    tuple val(meta), path(vcf)
    path bed

    output:
    tuple val(meta), path('*_benign.vcf.gz')        , optional: true, emit: benign_vcf
    tuple val(meta), path('*_pathogenic.vcf.gz')    , optional: true, emit: pathogenic_vcf
    tuple val(meta), path('*_selected.vcf')         , optional: true, emit: selected_vcf
    path 'versions.yml'                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tabix -p vcf ${vcf}
    tabix -h -R ${bed} ${vcf} | bgzip -c > ${prefix}_ontarget.vcf.gz
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Benign" -e "_benign" | bgzip -c > ${prefix}_benign.vcf.gz
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Pathogenic" -e "_pathogenic" | bgzip -c > ${prefix}_pathogenic.vcf.gz

    # Additional filtering steps
    zcat ${prefix}_ontarget.vcf.gz | grep -v "no_assertion_criteria_provided" | grep -e "^#" -e "reviewed_by_expert_panel" > ${prefix}_selected.vcf
    zcat ${prefix}_ontarget.vcf.gz | grep -v "^#" | grep -e "multiple_submitters" | grep "nonsense" >> ${prefix}_selected.vcf



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version 2>&1 | sed 's/^.*tabix //')
        bgzip: \$(bgzip --version 2>&1 | sed 's/^.*bgzip //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_benign.vcf.gz
    touch ${prefix}_pathogenic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo "1.20")
        bgzip: \$(echo "1.20")
    END_VERSIONS
    """
}
