
process SIMUSCOP_SEQTOPROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/simuscop:1.1.2--a03a6e546b386805':
        'community.wave.seqera.io/library/simuscop:1.1.2--251f0a0bfca0fa85' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(vcf)
    tuple path(capture)
    tuple path(fasta)


    output:
    tuple val(meta), path("*.profile"), emit: profile
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // seqToProfile requires VCF files to be uncompressed
    def vcf_command = vcf.name.endsWith('.gz') ? "gunzip -c ${vcf} > ${vcf.baseName} && vcf_file=${vcf.baseName}" : "vcf_file=${vcf}"

    """
    ${vcf_command}
    seqToProfile \\
    -b $bam \\
    -v \$vcf_file \\
    -t $capture \\
    -r $fasta \\
    $args >${prefix}.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simuscop: \$(seqToProfile -h 2>&1 | grep "^Version" | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simuscop: \$(seqToProfile -h 2>&1 | grep "^Version" | sed 's/Version: //')
    END_VERSIONS
    """
}
