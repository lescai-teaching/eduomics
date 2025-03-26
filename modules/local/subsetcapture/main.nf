process SUBSETCAPTURE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedtools_htslib:46c2f90c5e5c04ad':
        'community.wave.seqera.io/library/bedtools_htslib:d780123314824c66' }"

    input:
    tuple val(meta), val(chromosome)
    path target_chrom_size
    path capture_bed

    output:
    //tuple val(meta), path("*_target.chrom.sizes"), emit: target_chrom_size
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
    tabix -p bed ${prefix}.capture.bed.gz ${chromosome} > ${prefix}_${chromosome}_target.bed


    #padding
    bedtools slop -i ${prefix}_${chromosome}_target.bed -g ${target_chrom_size} -b 500 > ${prefix}_${chromosome}_target.pad500.bed
    bedtools slop -i ${prefix}_${chromosome}_target.bed -g ${target_chrom_size} -b 50 > ${prefix}_${chromosome}_target.pad50.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subsetcapture: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subsetcapture: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
