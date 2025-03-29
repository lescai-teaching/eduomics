include { PYCONVERTOSIM      } from '../../../modules/local/pyconvertosim/main'
include { SIMUSCOP_SIMUREADS } from '../../../modules/local/simuscop/simureads/main'

workflow PROFILE_SIMULATE_VARS_FASTQ {

    take:
    vcf_benign // channel: [mandatory] [[meta], vcf_benign]
    vcf_patho  // channel: [mandatory] [[meta], vcf_patho]
    simprofile // channel: [mandatory] [[meta], profile]
    fasta_fai  // channel: [mandatory] [[meta], fasta, fai]
    capture    // channel: [mandatory] [capture_500pad]

    main:

    ch_versions = Channel.empty()

    PYCONVERTOSIM( vcf_benign, vcf_patho )
    ch_versions = ch_versions.mix(PYCONVERTOSIM.out.versions)

    variants_to_inject = PYCONVERTOSIM.out.combined_variations.flatten()


    SIMUSCOP_SIMUREADS(
        simprofile,
        fasta_fai,
        variants_to_inject,
        capture
    )



    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

