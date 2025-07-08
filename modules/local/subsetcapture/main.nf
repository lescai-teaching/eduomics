process SUBSETCAPTURE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedtools_htslib:62c15a7abca4d6b2':
        'community.wave.seqera.io/library/bedtools_htslib:04fa82f202d44d67' }"

    input:
    val(meta)
    path target_chrom_size
    path capture_bed

    output:
    tuple val(meta), path("*.capture.bed.gz"), emit: capture_bed_gz
    tuple val(meta), path("*.capture.bed.gz.tbi"), emit: capture_bed_index
    tuple val(meta), path("*_target.bed"), emit: target_bed
    tuple val(meta), path("*_target.pad50.bed"), emit: target_bed_pad50
    tuple val(meta), path("*_target.pad500.bed"), emit: target_bed_pad500
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Compress and index origin capture BED
    bgzip -c ${capture_bed} > ${prefix}.capture.bed.gz
    tabix -p bed ${prefix}.capture.bed.gz

    # Subset by chrom
    tabix -p bed ${prefix}.capture.bed.gz ${meta.chromosome} > ${prefix}_${meta.chromosome}_target.bed


    #padding
    bedtools slop -i ${prefix}_${meta.chromosome}_target.bed -g ${target_chrom_size} -b 500 > ${prefix}_${meta.chromosome}_target.pad500.bed
    bedtools slop -i ${prefix}_${meta.chromosome}_target.bed -g ${target_chrom_size} -b 50 > ${prefix}_${meta.chromosome}_target.pad50.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.capture.bed.gz
    touch ${prefix}.capture.bed.gz.tbi
    touch ${prefix}_${meta.chromosome}_target.bed
    touch ${prefix}_${meta.chromosome}_target.pad50.bed
    touch ${prefix}_${meta.chromosome}_target.pad500.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
