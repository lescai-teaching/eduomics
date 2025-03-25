process PYCONVERTOSIM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py27_0':
        'biocontainers/biopython:1.70--np112py27_0' }"

input:
    tuple val(meta), path(vcf_benign)
    tuple val(meta2), path(vcf_pathogenic)

    output:
    tuple val(meta), path("*_base_variation.txt") , emit: base_variation
    tuple val(meta), path("*_patho_variation.txt"), emit: patho_variation
    tuple val(meta), path("*_variation_*.txt")    , emit: combined_variations
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python ${baseDir}/bin/convert_vcf_to_variation.py \
        -i ${vcf_benign} \
        -o ${prefix}_base_variation.txt \
        -n True

    python ${baseDir}/bin/convert_vcf_to_variation.py \
        -i ${vcf_pathogenic} \
        -o ${prefix}_patho_variation.txt

    sort -k4 -n ${prefix}_base_variation.txt > tmp && mv tmp ${prefix}_base_variation.txt
    sort -k4 -n ${prefix}_patho_variation.txt > tmp && mv tmp ${prefix}_patho_variation.txt

    # Process pathogenic variants one by one
    counter=0
    while IFS= read -r line
    do
        counter=\$((counter+1))
        cp ${prefix}_base_variation.txt ${prefix}_variation_\${counter}.txt
        echo "\$line" >> ${prefix}_variation_\${counter}.txt
        sort -k4 -n ${prefix}_variation_\${counter}.txt > tmp && mv tmp ${prefix}_variation_\${counter}.txt
    done < ${prefix}_patho_variation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_base_variation.txt
    touch ${prefix}_patho_variation.txt
    touch ${prefix}_variation_1.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
