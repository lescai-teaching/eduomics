process AISCENARIOS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_google-genai:842241e3cb24ce31':
        'community.wave.seqera.io/library/pip_google-genai:657c5e158292b26a' }"

    input:
    tuple val(meta), val(variant)
    tuple val(meta2), val(genes)

    output:
    tuple val(meta), path("*_scenario.txt"), emit: scenario
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def script_command = variant ? "generate_variant_scenarios.py --variant ${variant} --output ${prefix}_scenario.txt" : "generate_differentialexpression_scenarios.py --genes ${genes} --output ${prefix}_scenario.txt"

    """
    $script_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aiscenarios: 1.0.0)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def script_command = variant ? "./generate_variant_scenarios --variant ${variant} --output ${prefix}_scenario.txt" : "generate_differentialexpression_scenario --genes ${genes} --output ${prefix}_scenario.txt"

    """
    touch ${prefix}_scenario.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aiscenarios: 1.0.0)
    END_VERSIONS
    """
}
