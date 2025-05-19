include { BWA_MEM                                                } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES                                   } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR                                 } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                                        } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER                                  } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT                                 } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { SAMTOOLS_INDEX                                         } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MD                             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_RECAL                          } from '../../../modules/nf-core/samtools/index/main'


// to add the module to compare gVCFs and .solutions files


workflow VALIDATE_FASTQ_TO_VCF {

    take:
    simulated_reads_ch       // channel: [ val(meta), [ "reads_1.fq.gz", "reads_2.fq.gz" ] ]
    fasta                    // channel: [mandatory] [ val(meta), [ fasta ] ]
    fai                      // channel: [mandatory] [ val(meta), [fai] ]
    dict                     // channel: [mandatory] [ val(meta), [dict] ]
    bwa_index                // channel: [mandatory] [ val(meta), [bwa_index] ]
    intervals               // channel: [mandatory] [ [intervals] ]
    dbsnp                    // channel: [mandatory] [ [dbsnp] ]
    dbsnp_tbi                // channel: [mandatory] [ [dbsnp_tbi] ]
    known_snps               // channel: [mandatory] [ [known_snps] ]
    known_snps_tbi           // channel: [mandatory] [ [known_snps_tbi] ]


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

    known_sites_all = dbsnp.mix(known_snps).collect()
    known_sites_all_tbi = dbsnp_tbi.mix(known_snps_tbi).collect()

GATK4_BASERECALIBRATOR(
    bam_for_recal,
    fasta,
    fai,
    dict,
    known_sites_all.map{ it -> [[id:'sites'], it] },
    known_sites_all_tbi.map{ it -> [[id:'sites'], it] }
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    bam_for_applybqsr = GATK4_MARKDUPLICATES.out.bam
        .combine(INDEX_MD.out.bai, by: 0)
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
    ch_dict_gendb   = dict.map{meta, dict -> return dict}.toList()
    //see https://github.com/nf-core/modules/blob/e053975682e471dfed382fe6cdb29a71406f6630/subworkflows/nf-core/bam_create_som_pon_gatk/main.nf#L6
    ch_dict_gendb.view { "ch_dict_gendb: $it" }

//ch_gendb_input = Channel.value([id: 'test'])
//    .combine(ch_vcf)
//    .combine(ch_index)
//    .combine(intervals)
//    .combine(ch_dict_gendb)
//    .map { meta, vcf, tbi, interval, dict ->
//        [
//            meta,           // val(meta)
//            vcf.flatten(),  // path(vcf)
//            tbi.flatten(),  // path(tbi)
//            interval,       // path(interval_file)
//            [],             // val(interval_value)
//            []              // path(wspace)
//        ]
//    }

    ch_gendb_input = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)
        .combine(intervals)
        .map{ meta, vcf, tbi, intervals ->
            [
                [ id:'joint_variant_calling', intervals_name:intervals[0].baseName, num_intervals:1 ],
                [vcf],
                [tbi],
                intervals[1],
                [],
                []
            ]
        }

    ch_gendb_input.view { "ch_gendb_input: $it" }


    GATK4_GENOMICSDBIMPORT ( ch_gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    // Joint genotyping performed using GenotypeGVCFs
    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fai, dict, dbsnp.map{ it -> [ [:], it ] }, dbsnp_tbi.map{ it -> [ [:], it ] })
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)


    emit:
    genotype_vcf    = GATK4_GENOTYPEGVCFS.out.vcf    // channel: [ val(meta), [ vcf ] ]
    genotype_index  = GATK4_GENOTYPEGVCFS.out.tbi    // channel: [ val(meta), [ tbi ] ]
    versions        = ch_versions                    // channel: [ versions.yml ]

}
