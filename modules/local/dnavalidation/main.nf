process DNAVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/bash_grep_gzip:68663f683f82a964' :
    'community.wave.seqera.io/library/bash_grep_gzip:9781f12aeda63f1e' }"

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
    variantpos=\$(cut -d"-" -f2 <<<"${variant}")
    result_dir="dna_${variant}_validation"

    # check whether the variant position is present in the VCF
    if gzip -cd ${vcf} | grep -v '^#' | grep -q "\${variantpos}"; then
        mkdir -p "\${result_dir}"
        echo "${variant}" > "\${result_dir}/solution_${variant}.txt"
        cp ${vcf} "\${result_dir}/"
        for read in ${reads}; do
            cp "\$read" "\${result_dir}/"
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
    touch dna_${variant}_validation/simulated_validated.vcf.gz
    for read in ${reads}; do
        touch dna_${variant}_validation/\$(basename \$read)
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.0.17
    END_VERSIONS
    """
}
