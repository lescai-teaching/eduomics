include { WGSIM                          } from '../../../modules/nf-core/wgsim/main'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                  } from '../../../modules/nf-core/samtools/sort/main'
include { BWA_INDEX                      } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                        } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR         } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_APPLYBQSR                } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { SIMUSCOP_SEQTOPROFILE          } from '../../../modules/local/simuscop/seqtoprofile/main'

workflow FASTA_WGSIM_TO_PROFILE {

    take:
    fasta       // channel: [mandatory] [ val(meta), [ fasta ] ]
    fai         // channel: [optional]  [ val(meta), [fai] ]
    dict        // channel: [optional]  [ val(meta), [dict] ]
    bwa_index   // channel: [optional]  [ val(meta), [bwaindex] ]
    dbsnp       // channel: [mandatory] [ [dbsnp], [dbsnp_tbi] ]
    mills       // channel: [mandatory] [ [mills], [mills_tbi] ]

    main:

    ch_versions = Channel.empty()

    // simulate reads on the chosen fasta
    WGSIM(ch_fasta)

    if(!fai){
        SAMTOOLS_FAIDX ( fasta )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    // sets fasta index to correct input
    fastaindex = fai ?: SAMTOOLS_FAIDX.out.fai

    if(!dict){
        GATK4_CREATESEQUENCEDICTIONARY ( fasta )
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }
    // sets dictionary to correct input
    fastadict = dict ?: GATK4_CREATESEQUENCEDICTIONARY.out.dict

    if(!bwa_index){
        BWAMEM_INDEX ( fasta )
        ch_versions = ch_versions.mix(BWAMEM_INDEX.out.versions)
    }
    // sets bwaindex to correct input
    bwaindex    = bwa_index ?: BWAMEM_INDEX.out.index


    // align simulated reads to their reference
    BWA_MEM( WGSIM.out.fastq, bwaindex, fasta, true )

    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    GATK4_MARKDUPLICATES()


    GATK4_BASERECALIBRATOR()


    GATK4_APPLYBQSR()

    GATK4_HAPLOTYPECALLER()





    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

