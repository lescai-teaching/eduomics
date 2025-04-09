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
// to add the module to compare gVCFs and .solutions files


workflow VALIDATE_FASTQ_TO_VCF {

    take:
    fastq       // channel: [ val(meta), [ fastq ] ] // from wgsim
    fasta       // channel: [mandatory] [ val(meta), [ fasta ] ]
    fai         // channel: [mandatory] [ val(meta), [fai] ]
    dict        // channel: [mandatory] [ val(meta), [dict] ]
    bwa_index   // channel: [mandatory] [ val(meta), [bwa_index] ]
    target_bed  // channel: [mandatory] [ [target_bed] ]
    dbsnp       // channel: [mandatory] [ [dbsnp] ]
    dbsnp_tbi   // channel: [mandatory] [ [dbsnp_tbi] ]
    mills       // channel: [mandatory] [ [mills] ]
    mills_tbi   // channel: [mandatory] [ [mills_tbi] ]
    //known_sites_ch // channel: [ val(meta), known_sites_vcf ], use different channels for mills and dbsnp

    main:
    //Alignment
    BWA_MEM( reads_ch, bwa_index, fasta, true )
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

    known_sites_vcf_ch = dbsnp.mix(mills).collect()
    known_sites_tbi_ch = dbsnp_tbi.mix(mills_tbi).collect()

    GATK4_BASERECALIBRATOR(
        bam_for_recal,
        //reference_ch,
        fasta,
        fai,
        dict,
        known_sites_vcf_ch.map{ meta, vcf, tbi -> [ [ id: 'known_sites' ], vcf ] },
        known_sites_tbi_ch.map{ meta, vcf, tbi -> [ [ id: 'known_sites' ], tbi ] }
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    // Apply BQSR
    bam_for_applybqsr = GATK4_MARKDUPLICATES.out.bam
        .combine(INDEX_MD.out.bai, by: 0)
        .combine(GATK4_BASERECALIBRATOR.out.table, by: 0)
        .combine(Channel.value([[]]))

    GATK4_APPLYBQSR(
        bam_for_applybqsr,
        //reference_ch.map{ meta, fasta, fai, dict -> [ fasta ] },
        //reference_ch.map{ meta, fasta, fai, dict -> [ fai ] },
        //reference_ch.map{ meta, fasta, fai, dict -> [ dict ] }
        fasta.map{ meta, it -> [ it ] },
        fai.map{ meta, it -> [ it ] },
        dict.map{ meta, it -> [ it ] }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    INDEX_RECAL(GATK4_APPLYBQSR.out.bam)
    ch_versions = ch_versions.mix(INDEX_RECAL.out.versions)

    // HaplotypeCaller, joint variant calling!!!
    bam_for_calling = GATK4_APPLYBQSR.out.bam
        .combine(INDEX_RECAL.out.bai, by: 0)
        .combine(empty_intervals_ch) // check for intervals or target capture?
        .combine(empty_models_ch)

    GATK4_HAPLOTYPECALLER(GATK4_APPLYBQSR.out.bam, fasta, fai, dict,  dbsnp.map{ it -> [[id:'test'], it] },dbsnp_tbi.map{ it -> [[id:'test'], it] })



    ch_vcf          = GATK4_HAPLOTYPECALLER.out.vcf.collect{it[1]}.toList()
    ch_index        = GATK4_HAPLOTYPECALLER.out.tbi.collect{it[1]}.toList()
    ch_dict_gendb   = ch_dict.map{meta, dict -> return dict}.toList()

    ch_gendb_input  = Channel.of([id:genomicsdb])
        .combine(ch_vcf)
        .combine(ch_index)
        .combine(ch_gendb_intervals)
        .combine(ch_dict_gendb)
        .map{meta, vcf, tbi, interval, dict -> [meta, vcf, tbi, interval, [], dict]}

    GATK4_GENOMICSDBIMPORT ( ch_gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())


    GATK4_GENOTYPEGVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb,
        ///reference_ch, GATK4_GENOTYPEGVCFS has 6 inputs!
        tuple val(meta), path(input), path(gvcf_index), path(intervals), path(intervals_index)
        fasta,
        fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    emit:
    vcf             = GATK4_GENOTYPEGVCFS.out.vcf    // channel: [ val(meta), [ vcf ] ]
    vcf_index       = GATK4_GENOTYPEGVCFS.out.tbi    // channel: [ val(meta), [ tbi ] ]
    versions        = ch_versions                    // channel: [ versions.yml ]
    // compare .solution and obtained gvcf --> variants and also GT, e.g. 0/0, 0|0
}
