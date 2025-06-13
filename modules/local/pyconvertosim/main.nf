process PYCONVERTOSIM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.9.21--7652b035af009a82':
        'community.wave.seqera.io/library/python:3.9.21--8c83a010dbf906d0' }"

    input:
    tuple val(meta), path(vcf_benign)
    tuple val(meta2), path(vcf_pathogenic)

    output:
    tuple val(meta), path("*_base_variation.txt") , emit: base_variation
    tuple val(meta), path("*_patho_variation.txt"), emit: patho_variation
    tuple val(meta), path("*_simvar_*.txt")       , emit: combined_variations
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    convert_vcf_to_variation.py \
        -i ${vcf_benign} \
        -o ${prefix}_base_variation_unsorted.txt \
        -n True

    convert_vcf_to_variation.py \
        -i ${vcf_pathogenic} \
        -o ${prefix}_patho_variation_unsorted.txt

    sort -k4 -n ${prefix}_base_variation_unsorted.txt > ${prefix}_base_variation.txt
    sort -k4 -n ${prefix}_patho_variation_unsorted.txt > ${prefix}_patho_variation.txt

    # Process pathogenic variants one by one
    counter=0
    while IFS= read -r line
    do
        counter=\$((counter+1))
        variant=\$(echo -e \"\$line\" | awk -F'\\t' '{print \$3 "-" \$4 "-" \$5 "-" \$6}')
        cp ${prefix}_base_variation.txt ${prefix}_simvar_\${variant}_\${counter}_unsorted.txt
        echo "\$line" >> ${prefix}_simvar_\${variant}_\${counter}_unsorted.txt
        sort -k4 -n ${prefix}_simvar_\${variant}_\${counter}_unsorted.txt > ${prefix}_simvar_\${variant}_\${counter}.txt
    done < ${prefix}_patho_variation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_base_variation.txt
    touch ${prefix}_patho_variation.txt
    touch ${prefix}_simvar_chr22-12345-A-T_1.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
