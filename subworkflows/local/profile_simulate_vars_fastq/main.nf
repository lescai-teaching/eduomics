include { PYCONVERTOSIM      } from '../../../modules/local/pyconvertosim/main'
include { SIMUSCOP_SIMUREADS } from '../../../modules/local/simuscop/simureads/main'

workflow PROFILE_SIMULATE_VARS_FASTQ {

    take:
    vcf_benign // channel: [mandatory] [[meta], vcf_benign]
    vcf_patho  // channel: [mandatory] [[meta], vcf_patho]
    simprofile // channel: [mandatory] [[meta], profile]
    fasta_fai  // channel: [mandatory] [fasta, fai]
    capture    // channel: [mandatory] [capture_500pad]

    main:

    ch_versions = Channel.empty()

    PYCONVERTOSIM( vcf_benign, vcf_patho )
    ch_versions = ch_versions.mix(PYCONVERTOSIM.out.versions.ifEmpty([]))

    variants_to_inject = PYCONVERTOSIM.out.combined_variations.flatMap { m, files ->
                                                files.collect { file ->
                                                def var = file.getName().split('_')[2]
                                                def newmeta = m + [simulatedvar: "${var}"]
                                                return [newmeta, file]
                                                } }
    ch_variants_to_inject = params.istest
        ? variants_to_inject.take(params.test_limit)
        : variants_to_inject

    variants_to_inject.dump(tag: 'variants to inject')
    ch_variants_to_inject.dump(tag: 'variants selected for injection')

    SIMUSCOP_SIMUREADS(
        simprofile,
        fasta_fai,
        ch_variants_to_inject,
        capture
    )
    ch_versions = ch_versions.mix(SIMUSCOP_SIMUREADS.out.versions.ifEmpty([]))


    emit:
    simreads = SIMUSCOP_SIMUREADS.out.reads    // channel: [ val(meta), [ "control_1.fq.gz", "control_2.fq.gz", "case_1.fq.gz", "case_2.fq.gz" ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

