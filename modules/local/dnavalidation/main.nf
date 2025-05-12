process DNAVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/bash:5.2.21--67f53ad0451dfdce' :
    'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(reads)

    output:
    tuple val(meta), path("dna_${meta.simulatedvar}_validation")      , optional: true, emit: dna_validated_results
    tuple val(meta), path("dna_${meta.simulatedvar}_validation/*.vcf"), optional: true, emit: vcf_file
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant = meta.simulatedvar ?: ''
    """
    variantpos=\$(echo "${variant}" | cut -d"-" -f2)
    if grep -q "\${variantpos}" ${vcf} && [ "${meta.simulatedvar}" = "${meta2.simulatedvar}" ]; then
        mkdir -p dna_${variant}_validation
        echo "${variant}" > dna_${variant}_validation/solution_${variant}.txt
        cp ${vcf} dna_${variant}_validation/
        cp ${reads} dna_${variant}_validation/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def variant = meta.simulatedvar ?: 'chr22-1234-A-T'
    """
    mkdir -p dna_${variant}_validation
    touch dna_${variant}_validation/simulated_validated.vcf
    touch dna_${variant}_validation/simulated_validated_read_1.fastq.gz
    touch dna_${variant}_validation/simulated_validated_read_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.0.17
    END_VERSIONS
    """
}
