process DNAVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/bash:5.2.21--67f53ad0451dfdce' :
    'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"

    input:
    tuple val(meta), path(vcf), path(reads)

    output:
    tuple val(meta), path("dna_${meta.simulatedvar}_validation"), optional: true, emit: dna_validated_results
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant = meta.simulatedvar ?: ''
    """
    variantpos=\$(echo "${variant}" | cut -d"-" -f2)
    if zcat ${vcf} | grep -v \"#\" | grep -q "\${variantpos}"; then
        mkdir -p dna_${variant}_validation
        echo "${variant}" > dna_${variant}_validation/solution_${variant}.txt
        cp ${vcf} dna_${variant}_validation/
        for read in ${reads}; do
            cp \$read dna_${variant}_validation/
        done
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
    touch dna_${variant}_validation/solution_${variant}.txt
    touch dna_${variant}_validation/simulated_validated.vcf
    for read in ${reads}; do
        touch dna_${variant}_validation/\$(basename \$read)
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.0.17
    END_VERSIONS
    """
}
