include { SUBSETFASTATX      } from '../../../modules/local/subsetfastatx/main'
include { SUBSETGFF          } from '../../../modules/local/subsetgff/main'
include { SALMON_INDEX       } from '../../../modules/nf-core/salmon/index/main'

workflow PREPARE_RNA_GENOME {

    take:

    ch_fasta   // channel: [ path(fasta) ]
    ch_txfasta // channel: [ path(txfasta) ]
    ch_gff3    // channel: [ path(gff3) ]
    ch_meta    // channel: [ val(meta) ]

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    SUBSETFASTATX ( ch_meta, ch_txfasta, ch_gff3 )
    ch_versions = ch_versions.mix(SUBSETFASTATX.out.versions.first())
    ch_log      = ch_log.mix(SUBSETFASTATX.out.log)

    SUBSETGFF ( ch_meta, ch_gff3 )
    ch_versions = ch_versions.mix(SUBSETGFF.out.versions.first())
    ch_log      = ch_log.mix(SUBSETGFF.out.log)

    SALMON_INDEX ( ch_fasta, SUBSETFASTATX.out.fasta )
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())

    emit:
    txfasta_index      = SALMON_INDEX.out.index          // channel: [ path(txfasta_index) ]
    txfasta            = SUBSETFASTATX.out.fasta         // channel: [ path(txfasta) ]
    gene_lists         = SUBSETGFF.out.geneLists         // channel: [ path(geneLists) ]
    transcript_data    = SUBSETGFF.out.transcriptData    // channel: [ path(transcriptData) ]
    log_files          = ch_log                          // channel: [ path(log_files) ]
    versions           = ch_versions                     // channel: [ versions.yml ]
}

