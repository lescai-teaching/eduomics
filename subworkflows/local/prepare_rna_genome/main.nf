include { SUBSETFASTATX      } from '../../../modules/local/subsetfastatx/main'
include { SUBSETGFF          } from '../../../modules/local/subsetgff/main'
include { SALMON_INDEX       } from '../../../modules/nf-core/salmon/index/main'

workflow PREPARE_RNA_GENOME {

    take:
    ch_fasta   // channel: [ path(fasta)   ]
    ch_txfasta // channel: [ path(txfasta) ]
    ch_gff3    // channel: [ path(gff3)    ]
    ch_meta    // channel: [ val(meta)     ]

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    SUBSETFASTATX ( ch_meta, ch_txfasta, ch_gff3 )
    ch_versions = ch_versions.mix(SUBSETFASTATX.out.versions.first())
    ch_log      = ch_log.mix(SUBSETFASTATX.out.log)

    SUBSETGFF ( ch_meta, ch_gff3 )
    ch_versions = ch_versions.mix(SUBSETGFF.out.versions.first())
    ch_log      = ch_log.mix(SUBSETGFF.out.log)

    // Remove the meta from the filtered txfasta
    // Salmon index requires inputs without a meta
    ch_filtered_txfasta_nometa = SUBSETFASTATX.out.filtered_txfasta
        .map { meta, fasta -> fasta }

    SALMON_INDEX ( ch_fasta, ch_filtered_txfasta_nometa )
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())

    emit:
    txfasta_index           = SALMON_INDEX.out.index                 // channel: [ path(salmon_index)                       ]
    filtered_txfasta        = SUBSETFASTATX.out.filtered_txfasta     // channel: [ val(meta), path(filtered_txfasta)        ]
    gene_lists              = SUBSETGFF.out.geneLists                // channel: [ val(meta), path(geneLists)               ]
    genes_list_associations = SUBSETGFF.out.genes_list_association   // channel: [ val(meta), path(genes_list_associations) ]. This is a TSV file containing the association between a specific list and its genes. Necessary to add genes in meta by coupling the genes with a specific countMatrix
    filtered_gff3           = SUBSETGFF.out.filtered_gff3            // channel: [ val(meta), path(gff3_filtered)           ]
    log_files               = ch_log                                 // channel: [ val(meta), path(log_files)               ]
    versions                = ch_versions                            // channel: [ versions.yml                             ]
}

