include { BWA_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MD     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_RECAL  } from '../../../modules/nf-core/samtools/index/main'
/// to add the module to compare gVCFs and .solutions files

    take:
    reads_ch       // channel: [ val(meta), [ reads ] ]
    reference_ch   // channel: [ val(meta), [ fasta, fai, dict ] ]
    bwa_index      // channel: [ val(meta), [ bwa_index ] ] // from??? check previous modules/subworkflows!!!!
    known_sites_ch // channel: [ val(meta), [ known_sites_vcf, known_sites_tbi ] ]
    dbsnp          // channel: [ val(meta), [ dbsnp_vcf ] ]
    dbsnp_tbi      // channel: [ val(meta), [ dbsnp_tbi ] ]


workflow VALIDATE_FASTQ_TO_VCF {

    take:
    reads_ch       // channel: [ val(meta), [ reads ] ]
    reference_ch   // channel: [ val(meta), reference_fasta ]
    known_sites_ch // channel: [ val(meta), known_sites_vcf ]

    main:
    //Alignment
    BWA_MEM(reads_ch, bwa_index, reference_ch.map{ meta, fasta, fai, dict -> [ fasta ] }, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

// marking duplicates
    GATK4_MARKDUPLICATES(
        BWA_MEM.out.bam,
        reference_ch.map{ meta, fasta, fai, dict -> [ fasta ] },
        reference_ch.map{ meta, fasta, fai, dict -> [ fai ] }
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    INDEX_MD(GATK4_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(INDEX_MD.out.versions)

    // BQSR
    //take bam/bai
    bam_for_recal = GATK4_MARKDUPLICATES.out.bam
        .combine(INDEX_MD.out.bai, by: 0)
        .combine(Channel.value([[]]))

    GATK4_BASERECALIBRATOR(
        bam_for_recal,
        reference_ch,
        known_sites_ch.map{ meta, vcf, tbi -> [ [ id: 'known_sites' ], vcf ] },
        known_sites_ch.map{ meta, vcf, tbi -> [ [ id: 'known_sites' ], tbi ] }
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    // Apply BQSR
    bam_for_applybqsr = GATK4_MARKDUPLICATES.out.bam
        .combine(INDEX_MD.out.bai, by: 0)
        .combine(GATK4_BASERECALIBRATOR.out.table, by: 0)
        .combine(Channel.value([[]]))

    GATK4_APPLYBQSR(
        bam_for_applybqsr,
        reference_ch.map{ meta, fasta, fai, dict -> [ fasta ] },
        reference_ch.map{ meta, fasta, fai, dict -> [ fai ] },
        reference_ch.map{ meta, fasta, fai, dict -> [ dict ] }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    INDEX_RECAL(GATK4_APPLYBQSR.out.bam)
    ch_versions = ch_versions.mix(INDEX_RECAL.out.versions)

    // HaplotypeCaller, joint variant calling!!!
    GATK4_HAPLOTYPECALLER(GATK4_APPLYBQSR.out.bam, reference_ch)

    // Import Genomics DB
    GATK4_GENOMICSDBIMPORT(GATK4_HAPLOTYPECALLER.out.vcf.map{ meta, vcf -> [ [ id: 'genomicsdb' ], vcf ] })
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)


    GATK4_GENOTYPEGVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb,
        reference_ch
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    emit:
    bam             = GATK4_APPLYBQSR.out.bam        // channel: [ val(meta), [ bam ] ]
    vcf             = GATK4_GENOTYPEGVCFS.out.vcf    // channel: [ val(meta), [ vcf ] ]
    vcf_index       = GATK4_GENOTYPEGVCFS.out.tbi    // channel: [ val(meta), [ tbi ] ]
    versions        = ch_versions                    // channel: [ versions.yml ]


    emit:
    final_vcf = GATK4_GENOTYPEGVCFS.out.vcf
    final_vcf_index = GATK4_GENOTYPEGVCFS.out.vcf_index
    // compare .solution and obtained gvcf --> variants and also GT, e.g. 0/0, 0|0
}
