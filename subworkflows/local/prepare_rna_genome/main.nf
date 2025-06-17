include { SUBSETFASTATX                              } from '../../../modules/local/subsetfastatx/main'
include { SUBSETGFF                                  } from '../../../modules/local/subsetgff/main'
include { SALMON_INDEX                               } from '../../../modules/nf-core/salmon/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SUBSET    } from '../../../modules/nf-core/samtools/faidx/main'


workflow PREPARE_RNA_GENOME {

    take:
    ch_meta          // channel: [ val(meta)         ]
    ch_gff3          // channel: [ path(gff3)        ]
    ch_txfasta       // channel: [ path(txfasta)     ]
    ch_genomefasta   // channel: [ path(genomefasta) ]

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    SUBSETGFF ( ch_meta, ch_gff3 )
    ch_versions = ch_versions.mix(SUBSETGFF.out.versions)
    ch_log      = ch_log.mix(SUBSETGFF.out.subsetgff_parsing_log)

    SUBSETFASTATX ( ch_txfasta, SUBSETGFF.out.filtered_transcript_data )
    ch_versions = ch_versions.mix(SUBSETFASTATX.out.versions)
    ch_log      = ch_log.mix(SUBSETFASTATX.out.subsetfastatx_parsing_log)

    // Add the meta to the genome fasta
    // Samtools faidx requires an input with the meta
    ch_genomefasta_with_meta = ch_meta.combine(ch_genomefasta)
        .map { meta, fasta -> [meta, fasta] }

    SAMTOOLS_FAIDX_SUBSET ( ch_genomefasta_with_meta, [[],[]], [] )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SUBSET.out.versions)

    // Remove the meta from the filtered txfasta and from the filtered genome fasta
    // Salmon index requires inputs without a meta
    ch_filtered_txfasta_nometa = SUBSETFASTATX.out.filtered_txfasta
        .map { meta, fasta -> fasta }

    ch_filtered_genomefasta_nometa = SAMTOOLS_FAIDX_SUBSET.out.fa
        .map { meta, fasta -> fasta }

    SALMON_INDEX ( ch_filtered_genomefasta_nometa, ch_filtered_txfasta_nometa )
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    ch_rna_bundle = SUBSETFASTATX.out.filtered_txfasta
        .join(SUBSETGFF.out.filtered_annotation)
        .join(SAMTOOLS_FAIDX_SUBSET.out.fa)
        .combine(SALMON_INDEX.out.index)
        .map { it ->
            [it[0], it[1..-1]]
        }

    emit:
    filtered_annotation      = SUBSETGFF.out.filtered_annotation.collect()      // channel: [ val(meta), path(filtered_annotation)                                                                     ]
    gene_lists               = SUBSETGFF.out.geneLists                          // channel: [ val(meta), path(geneLists)                                                                               ]
    gene_list_association    = SUBSETGFF.out.genes_list_association             // channel: [ val(meta), path(gene_list_association)                                                                   ]
    filtered_transcript_data = SUBSETGFF.out.filtered_transcript_data.collect() // channel: [ val(meta), path(filtered_transcript_data)                                                                ]
    filtered_txfasta         = SUBSETFASTATX.out.filtered_txfasta.collect()     // channel: [ val(meta), path(filtered_txfasta)                                                                        ]
    filtered_genomefasta     = SAMTOOLS_FAIDX_SUBSET.out.fa.collect()           // channel: [ val(meta), path(filtered_genomefasta)                                                                    ]
    txfasta_index            = SALMON_INDEX.out.index.collect()                 // channel: [ path(salmon_index)                                                                                       ]
    log_files                = ch_log                                           // channel: [ val(meta), path(log_files)                                                                               ]
    rna_bundle               = ch_rna_bundle.collect()                          // channel: [ val(meta), [path(filtered_txfasta), path(filtered_gff3), path(gfiltered_genomefasta), path(salmonindex)] ]
    versions                 = ch_versions                                      // channel: [ versions.yml                                                                                             ]
}

