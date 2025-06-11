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
    # Check if VCF uses 'chr' prefix in chromosome names (more robust approach)Add commentMore actions
    if zcat ${vcf} | grep -v "^#" | head -1 | cut -f1 | grep -q "^chr"; then
        vcf_has_chr="yes"
    else
        vcf_has_chr="no"
    fi

    # Check if BED uses 'chr' prefix in chromosome names
    if head -1 ${bed} | cut -f1 | grep -q "^chr"; then
        bed_has_chr="yes"
    else
        bed_has_chr="no"
    fi

    echo "DEBUG: VCF has chr prefix: \$vcf_has_chr"
    echo "DEBUG: BED has chr prefix: \$bed_has_chr"

    # Create adjusted BED file if needed
    if [[ "\$vcf_has_chr" == "yes" && "\$bed_has_chr" == "no" ]]; then
        echo "DEBUG: Adding chr prefix to BED file"
        awk 'BEGIN{OFS="\\t"} {if(\$1 !~ /^chr/) \$1="chr"\$1; print}' ${bed} > adjusted_${bed}
        bed_file="adjusted_${bed}"
    elif [[ "\$vcf_has_chr" == "no" && "\$bed_has_chr" == "yes" ]]; then
        echo "DEBUG: Removing chr prefix from BED file"
        awk 'BEGIN{OFS="\\t"} {gsub(/^chr/, "", \$1); print}' ${bed} > adjusted_${bed}
        bed_file="adjusted_${bed}"
    else
        echo "DEBUG: No BED file adjustment needed"
        bed_file="${bed}"
    fi

    echo "DEBUG: Using bed file: \$bed_file"

    # Index the VCF file
    tabix -p vcf ${vcf}

    # Extract variants using the adjusted BED fileAdd commentMore actions
    tabix -h -R \$bed_file ${vcf} | bgzip -c > ${prefix}_ontarget.vcf.gz

    # Extract benign variants
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Benign" -e "_benign" | grep -e "#" -e "multiple_submitters" > ${prefix}_benign.vcf

    # Extract pathogenic variants
    zcat ${prefix}_ontarget.vcf.gz | grep -e "^#" -e "Pathogenic" -e "_pathogenic" | grep -e "#" -e "multiple_submitters" | grep -e "#" -e "nonsense" > ${prefix}_pathogenic.vcf


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
