
process SIMUSCOP_SIMUREADS {
    tag "$meta2.simulatedvar"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/simuscop:2.0.0--46cea91693c644bb':
        'community.wave.seqera.io/library/simuscop:2.0.0--44718e83828ef52d' }"

    input:
    tuple val(meta), path(profile)
    tuple path(fasta), path(fai)
    tuple val(meta2), path(variantstoinject)
    path(capture)

    output:
    tuple val(meta2), path("simulated_reads/*.fq.gz"), emit: reads
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coverage = meta.coverage ? "${meta.coverage}" : "50"
    """
    cat <<EOF >abundance_casecontrol.txt
    0	1
    1	0
    EOF

    cat <<EOF >simulation.config
    ref=${fasta}
    profile=${profile}
    variation=${variantstoinject}
    abundance=abundance_casecontrol.txt
    target=${capture}
    name=normal,disease
    output=simulated_reads
    threads=${task.cpus}
    coverage=${coverage}
    layout=PE
    EOF

    simuReads simulation.config

    for read in \$(ls simulated_reads/*.fq)
    do
    gzip \$read
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simuscop: \$(seqToProfile -h 2>&1 | grep "^Version" | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p simulated_reads
    touch simulated_reads/normal_1.fq.gz
    touch simulated_reads/normal_2.fq.gz
    touch simulated_reads/case_1.fq.gz
    touch simulated_reads/case_2.fq.gz

    touch simulation.config
    touch abundance_casecontrol.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simuscop: \$(seqToProfile -h 2>&1 | grep "^Version" | sed 's/Version: //')
    END_VERSIONS
    """
}
