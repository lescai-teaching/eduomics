include { WGSIM                          } from '../../../modules/nf-core/wgsim/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MD     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_RECAL  } from '../../../modules/nf-core/samtools/index/main'
//include { SAMTOOLS_SORT                  } from '../../../modules/nf-core/samtools/sort/main'
include { BWA_MEM                        } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR         } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { SIMUSCOP_SEQTOPROFILE          } from '../../../modules/local/simuscop/seqtoprofile/main'

workflow FASTA_WGSIM_TO_PROFILE {

    take:
    fasta       // channel: [mandatory] [ val(meta), [ fasta ] ]
    fai         // channel: [mandatory] [ val(meta), [fai] ]
    dict        // channel: [mandatory] [ val(meta), [dict] ]
    bwa_index   // channel: [mandatory] [ val(meta), [bwa_index] ]
    dbsnp       // channel: [mandatory] [ [dbsnp] ]
    dbsnp_tbi   // channel: [mandatory] [ [dbsnp_tbi] ]
    mills       // channel: [mandatory] [ [mills] ]
    mills_tbi   // channel: [mandatory] [ [mills_tbi] ]
    capture     // channel: [mandatory] [ [capture] ]

    main:

    ch_versions = Channel.empty()

    // simulate reads on the chosen fasta
    WGSIM(fasta)

    // align simulated reads to their reference
    BWA_MEM( WGSIM.out.fastq, bwa_index, fasta, true )

    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // ----> GATK4 MINI WORKFLOW to call normal variants <-----

    GATK4_MARKDUPLICATES(
        BWA_MEM.out.bam,
        fasta.map{ meta, it -> [ it ] },
        fai.map{ meta, it -> [ it ] },
        )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    empty_intervals_ch = Channel.value([[]])

    INDEX_MD(GATK4_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(INDEX_MD.out.versions)

    bam_for_recal = GATK4_MARKDUPLICATES.out.bam
                        .combine(INDEX_MD.out.bai, by: 0)
                        .combine(empty_intervals_ch)

    known_sites_all = dbsnp.mix(mills).collect()
    known_sites_all_tbi = dbsnp_tbi.mix(mills_tbi).collect()

    GATK4_BASERECALIBRATOR(
        bam_for_recal,
        fasta,
        fai,
        dict,
        known_sites_all.map{ it -> [[id:'sites'], it] },
        known_sites_all_tbi.map{ it -> [[id:'sites'], it] }
        )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    bam_for_applybqsr = BWA_MEM.out.bam
        .combine(SAMTOOLS_INDEX.out.bai, by: 0)
        .combine(GATK4_BASERECALIBRATOR.out.table, by: 0)
        .combine(empty_intervals_ch)

    GATK4_APPLYBQSR(
        bam_for_applybqsr,
        fasta.map{ meta, it -> [ it ] },
        fai.map{ meta, it -> [ it ] },
        dict.map{ meta, it -> [ it ] }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    INDEX_RECAL(GATK4_APPLYBQSR.out.bam)

    empty_models_ch = Channel.value([[]])

    bam_for_calling = GATK4_APPLYBQSR.out.bam
        .combine(INDEX_RECAL.out.bai, by: 0)
        .combine(empty_intervals_ch)
        .combine(empty_models_ch)


    GATK4_HAPLOTYPECALLER(
        bam_for_calling,
        fasta,
        fai,
        dict,
        dbsnp.map{ it -> [[id:'test'], it] },
        dbsnp_tbi.map{ it -> [[id:'test'], it] }
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)


    SIMUSCOP_SEQTOPROFILE(
        GATK4_APPLYBQSR.out.bam,
        GATK4_HAPLOTYPECALLER.out.vcf,
        capture,
        fasta.map{ meta, it -> [ it ] }
    )
    ch_versions = ch_versions.mix(SIMUSCOP_SEQTOPROFILE.out.versions)



    emit:
    bam      = GATK4_APPLYBQSR.out.bam              // channel: [ val(meta), [ bam ] ]
    vcf      = GATK4_HAPLOTYPECALLER.out.vcf        // channel: [ val(meta), [ vcf ] ]
    profile  = SIMUSCOP_SEQTOPROFILE.out.profile    // channel: [ val(meta), [ profile ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}

