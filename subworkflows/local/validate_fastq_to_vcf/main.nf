include { BWA_MEM                        } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR         } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT         } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS            } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MD     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_RECAL  } from '../../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_SORT                  } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS           as MERGE_GENOTYPEGVCFS       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS           as MERGE_VQSR                } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// to add the module to compare gVCFs and .solutions files


workflow VALIDATE_FASTQ_TO_VCF {

    take:
    simulated_reads_ch       // channel: [ val(meta), [ "reads_1.fq.gz", "reads_2.fq.gz" ] ]
    fasta       // channel: [mandatory] [ val(meta), [ fasta ] ]
    fai         // channel: [mandatory] [ val(meta), [fai] ]
    dict        // channel: [mandatory] [ val(meta), [dict] ]
    bwa_index   // channel: [mandatory] [ val(meta), [bwa_index] ]
    target_bed  // channel: [mandatory] [ [target_bed] ]
    dbsnp       // channel: [mandatory] [ [dbsnp] ]
    dbsnp_tbi   // channel: [mandatory] [ [dbsnp_tbi] ]
    dbsnp_vqsr
    known_indels //mills
    known_indels_tbi
    known_indels_vqsr
    known_snps
    known_snps_tbi //
    known_snps_vqsr
    resource_indels_vcf //?
    resource_indels_tbi //?
    resource_snps_vcf //?
    resource_snps_tbi //?

    main:
    ch_versions = Channel.empty()

    //Alignment
    BWA_MEM( simulated_reads_ch, bwa_index, fasta, true )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // marking duplicates
    GATK4_MARKDUPLICATES(BWA_MEM.out.bam,
        fasta.map{ meta, it -> [ it ] },
        fai.map{ meta, it -> [ it ] })
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    empty_intervals_ch = Channel.value([[]])

    INDEX_MD(GATK4_MARKDUPLICATES.out.bam)

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


    ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> vcf }.collect()
    ch_index = GATK4_HAPLOTYPECALLER.out.tbi.map { meta, tbi -> tbi }.collect()

    ch_gendb_input  = Channel.value()
        .combine(ch_vcf)
        .combine(ch_index)
        .combine(target_bed)
        .combine(Channel.value([]))  // empty interval_value
        .combine(Channel.value([]))  // empty workspace
        //.map{meta, gvcf, tbi, interval -> [meta, vcf, tbi, interval, [], dict]}
        .map{ meta, gvcf, tbi, intervals -> [ [ id:'joint_variant_calling', intervals_name:intervals.baseName, num_intervals:meta.num_intervals ], gvcf, tbi, intervals ] }
        .groupTuple(by:3) //join on interval file
        .map{ meta_list, gvcf, tbi, intervals ->
            // meta is now a list of [meta1, meta2] but they are all the same. So take the first element.
            [ meta_list[0], gvcf, tbi, intervals, [], [] ]


    GATK4_GENOMICSDBIMPORT ( ch_gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    // Joint genotyping performed using GenotypeGVCFs
    // Sort vcfs called by interval within each VCF

    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fai, dict, dbsnp.map{ it -> [ [:], it ] }, dbsnp_tbi.map{ it -> [ [:], it ] })
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)
    gvcf_to_merge = BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [ meta.subMap('num_intervals') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'haplotypecaller' ], vcf ]}.groupTuple()

    MERGE_GENOTYPEGVCFS(gvcf_to_merge, dict)

    vqsr_input = MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true)
    indels_resource_label = known_indels_vqsr.mix(dbsnp_vqsr).collect()
    snps_resource_label = known_snps_vqsr.mix(dbsnp_vqsr).collect()

    // Recalibrate INDELs and SNPs separately
    VARIANTRECALIBRATOR_INDEL(
        vqsr_input,
        resource_indels_vcf,
        resource_indels_tbi,
        indels_resource_label,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    VARIANTRECALIBRATOR_SNP(
        vqsr_input,
        resource_snps_vcf,
        resource_snps_tbi,
        snps_resource_label,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    //Prepare SNPs and INDELs for ApplyVQSR
    // Step 1. : ApplyVQSR to SNPs
    // Step 2. : Use ApplyVQSR_SNP output and run ApplyVQSR_INDEL. This avoids duplicate entries in the vcf as described here: https://hpc.nih.gov/training/gatk_tutorial/vqsr.html

    // Join results of variant recalibration into a single channel tuple
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_snp = vqsr_input.join(VARIANTRECALIBRATOR_SNP.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_SNP.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    GATK4_APPLYVQSR_SNP(
        vqsr_input_snp,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })

    // Join results of ApplyVQSR_SNP and use as input for Indels to avoid duplicate entries in the result
    // Rework meta for variantscalled.csv and annotation tools
    vqsr_input_indel = GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi).map{ meta, vcf, tbi -> [ meta + [ id:'joint_variant_calling' ], vcf, tbi ]}
        .join(VARIANTRECALIBRATOR_INDEL.out.recal, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.idx, failOnDuplicate: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.tranches, failOnDuplicate: true)
        .map{ meta, vcf, tbi, recal, index, tranche -> [ meta + [ id:'recalibrated_joint_variant_calling' ], vcf, tbi, recal, index, tranche ] }

    GATK4_APPLYVQSR_INDEL(
        vqsr_input_indel,
        fasta.map{ meta, fasta -> [ fasta ] },
        fai.map{ meta, fai -> [ fai ] },
        dict.map{ meta, dict -> [ dict ] })


    // The following is an ugly monster to achieve the following:
    // When MERGE_GENOTYPEGVCFS and GATK4_APPLYVQSR are run, then use output from APPLYVQSR
    // When MERGE_GENOTYPEGVCFS and NOT GATK4_APPLYVQSR , then use the output from MERGE_GENOTYPEGVCFS

    merge_vcf_for_join = MERGE_GENOTYPEGVCFS.out.vcf.map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
    merge_tbi_for_join = MERGE_GENOTYPEGVCFS.out.tbi.map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

    // Remap for both to have the same key, if ApplyBQSR is not run, the channel is empty --> populate with empty elements
    vqsr_vcf_for_join = GATK4_APPLYVQSR_INDEL.out.vcf.ifEmpty([[:], []]).map{meta, vcf -> [[id: 'joint_variant_calling'] , vcf]}
    vqsr_tbi_for_join = GATK4_APPLYVQSR_INDEL.out.tbi.ifEmpty([[:], []]).map{meta, tbi -> [[id: 'joint_variant_calling'] , tbi]}

    // Join on metamap
    // If both --> meta, vcf_merged, vcf_bqsr
    // If not VQSR --> meta, vcf_merged, []
    // if the second is empty, use the first
    genotype_vcf = merge_vcf_for_join.join(vqsr_vcf_for_join, remainder: true).map{
        meta, joint_vcf, recal_vcf ->

        vcf_out = recal_vcf ?: joint_vcf

        [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"haplotypecaller"], vcf_out]
    }

    genotype_index = merge_tbi_for_join.join(vqsr_tbi_for_join, remainder: true).map{
        meta, joint_tbi, recal_tbi ->

        tbi_out = recal_tbi ?: joint_tbi

        [[id:"joint_variant_calling", patient:"all_samples", variantcaller:"haplotypecaller"], tbi_out]
    }

    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    versions = versions.mix(VARIANTRECALIBRATOR_SNP.out.versions)
    versions = versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    emit:
    genotype_index  // channel: [ val(meta), [ tbi ] ]
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]

    versions        // channel: [ versions.yml ]



//    emit:
//    vcf             = GATK4_GENOTYPEGVCFS.out.vcf    // channel: [ val(meta), [ vcf ] ]
//    vcf_index       = GATK4_GENOTYPEGVCFS.out.tbi    // channel: [ val(meta), [ tbi ] ]
//    versions        = ch_versions                    // channel: [ versions.yml ]
    // create .solution with deleterious variant
}

