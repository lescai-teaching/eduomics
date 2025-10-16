include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MD     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_RECAL  } from '../../../modules/nf-core/samtools/index/main'
include { BWA_MEM                        } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR         } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT         } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS            } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { DNAVALIDATION                  } from '../../../modules/local/dnavalidation/main'


workflow FASTQ_VARIANT_TO_VALIDATION {

    take:
    reads       // channel: [mandatory] [ val(meta), [ reads ]   ] <- NB: this channel contains 4 reads 2 * case + 2 * control / modified in sub-workflow
    fasta       // channel: [mandatory] [ val(meta), [ fasta ]   ]
    fai         // channel: [mandatory] [ val(meta), [fai]       ]
    dict        // channel: [mandatory] [ val(meta), [dict]      ]
    bwa_index   // channel: [mandatory] [ val(meta), [bwa_index] ]
    dbsnp       // channel: [mandatory] [ [dbsnp]                ]
    dbsnp_tbi   // channel: [mandatory] [ [dbsnp_tbi]            ]
    mills       // channel: [mandatory] [ [mills]                ]
    mills_tbi   // channel: [mandatory] [ [mills_tbi]            ]
    capture     // channel: [mandatory] [ [capture]              ]

    main:

    ch_versions = Channel.empty()

    simulated_reads_ch = reads
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
                m, samplereads ->
                [m, samplereads]
            }
    simulated_reads_ch.dump(tag: 'sim reads by sample')

    // align simulated reads to their reference
    BWA_MEM( simulated_reads_ch, bwa_index, fasta, true )

    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

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
        // pair each BAM with its index using join to avoid cartesian products
        .join(INDEX_MD.out.bai)
        .map { meta, bam, bai -> [meta, bam, bai] }
        .combine(empty_intervals_ch)
    bam_for_recal.dump(tag: 'bam for recal')

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

    bam_for_applybqsr = GATK4_MARKDUPLICATES.out.bam
        // match BAM, index and recalibration table by sample
        .join(INDEX_MD.out.bai)
        .join(GATK4_BASERECALIBRATOR.out.table)
        .map { meta, bam, bai, table -> [meta, bam, bai, table] }
        .combine(capture)

    bam_for_applybqsr.dump(tag: 'bam for applybqsr')

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
        .join(INDEX_RECAL.out.bai)
        .map { meta,bam,bai -> [meta, bam, bai] }
        .combine(capture)
        .combine(empty_models_ch)

    bam_for_calling.dump(tag: 'bam for calling')

    GATK4_HAPLOTYPECALLER(
        bam_for_calling,
        fasta,
        fai,
        dict,
        dbsnp.map{ it -> [[id:'test'], it] },
        dbsnp_tbi.map{ it -> [[id:'test'], it] }
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    // even if the samples in this execution will all be about the same variant
    // we need to collect the VCFs and we cannot just group by meta because
    // the sample ID will be different

    GATK4_HAPLOTYPECALLER.out.vcf.dump(tag: 'single sample VCF')

    gendb_input = GATK4_HAPLOTYPECALLER.out.vcf
    .join(GATK4_HAPLOTYPECALLER.out.tbi)
    .map { meta, vcf, tbi ->
        def key = meta.simulatedvar
        return [key, meta, vcf, tbi]
    }
    .groupTuple(by: 0)
    .combine(capture)
    .map { key, meta_list, vcf_list, tbi_list, interval ->
        return [meta_list[0], vcf_list, tbi_list, interval, [], []]
    }

    gendb_input.dump(tag: 'genDB channel')

    // Import GVCFs into GenomicsDB
    GATK4_GENOMICSDBIMPORT(
        gendb_input,
        false, // run_intlist
        false, // run_updatewspace
        false  // input_map
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    GATK4_GENOMICSDBIMPORT.out.genomicsdb.dump(tag: 'genomicsDB out')

    // Joint genotyping
    GATK4_GENOTYPEGVCFS(
        GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, db -> [ meta, db, [], [], [] ] },
        fasta,
        fai,
        dict,
        dbsnp.map{ it -> [[id:'test'], it] },
        dbsnp_tbi.map{ it -> [[id:'test'], it] }
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    reads.dump(tag: 'original reads')

    // the sample key in a multi-sample VCF file doesn't make sense
    // and more importantly is not present in the original case/control reads
    // therefore must be removed

    validation_ch = GATK4_GENOTYPEGVCFS.out.vcf
        .map{ meta, vcf ->
            def newmeta = meta.findAll { it.key != 'sample' }
            return [newmeta, vcf]
        }
        .join(reads)

    validation_ch.dump(tag: 'validation CH')

    DNAVALIDATION( validation_ch )
    ch_versions = ch_versions.mix(DNAVALIDATION.out.versions)

    DNAVALIDATION.out.dna_validated_results.dump(tag: 'dnavalidation output')

    scenarios_ch = DNAVALIDATION.out.dna_validated_results
        .map { meta, results ->
            return [meta, meta.simulatedvar, meta.simulatedgene]
        }

    scenarios_ch.dump(tag: 'scenario CH')

    emit:
    // the following channels could be empty as a result of optional emission
    // must be handled in main before they are passed to the aiscenarios module
    simulation = DNAVALIDATION.out.dna_validated_results  // channel: [ val(meta), path(validated_results_folder/*) ]
    scenario   = scenarios_ch                             // channel: [ val(meta), val(variant), val(gene)          ]
    versions   = ch_versions                              // channel: [ versions.yml                                ]

}
