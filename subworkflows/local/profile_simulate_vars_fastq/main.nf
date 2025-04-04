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
    ch_versions = ch_versions.mix(PYCONVERTOSIM.out.versions)

    variants_to_inject = PYCONVERTOSIM.out.combined_variations.flatMap { m, files ->
                                                files.collect { file ->
                                                def var = file.getName().split('_')[2]
                                                def newmeta = m + [simulatedvar: "${var}"]
                                                return [newmeta, file]
                                                } }


    SIMUSCOP_SIMUREADS(
        simprofile,
        fasta_fai,
        variants_to_inject,
        capture
    )
    ch_versions = ch_versions.mix(SIMUSCOP_SIMUREADS.out.versions)

    simulated_reads_ch = SIMUSCOP_SIMUREADS.out.reads
            .map{ m, files ->
                    def grouped = files.groupBy { file ->
                        file.name.replaceFirst(/_[12]\.fq\.gz$/, '')
                        }
                    grouped.collect { sampleName, group ->
                        def updatedMeta = m + [sample: sampleName]
                        [updatedMeta, group.sort()]
                        }
            }
            .flatMap{
                m, reads ->
                [m, reads]
            }


    emit:
    simreads = simulated_reads_ch    // channel: [ val(meta), [ "reads_1.fq.gz", "reads_2.fq.gz" ] ]
    versions = ch_versions           // channel: [ versions.yml ]
}

