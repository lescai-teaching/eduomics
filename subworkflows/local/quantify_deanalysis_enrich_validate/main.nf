
include { SALMON_QUANT         } from '../../../modules/nf-core/salmon/quant/main'
include { DEANALYSIS           } from '../../../modules/local/deanalysis/main'
include { ENRICHMENT           } from '../../../modules/local/enrichment/main'
include { RNASEQVALIDATION     } from '../../../modules/local/rnaseqvalidation/main'

workflow QUANTIFY_DEANALYSIS_ENRICH_VALIDATE {

    take:
    ch_simreads           // channel: [ val(meta), path(simreads)        ]
    ch_index              // channel: [ path(salmon_index)               ]
    ch_filteredgff3       // channel: [ val(meta), path(filteredgff3)    ]
    ch_filteredtxfasta    // channel: [ val(meta), path(filteredtxfasta) ]
    ch_alignment_mode     // channel: [ val(alignment_mode)              ]
    ch_libtype            // channel: [ val(libtype)                     ]
    ch_transcriptData     // channel: [ val(meta), path(transcriptData)  ]. Utilised by the analysis module to produce the tx2gene

    main:

    ch_versions = Channel.empty()

    // Remove the meta from filteredgff3 and filteredtxfasta
    // Salmon quant module requires input with no meta
    ch_filteredgff3_nometa = ch_filteredgff3
                .map { meta, gff3 -> gff3 }

    ch_filteredtxfasta_nometa = ch_filteredtxfasta
                .map { meta, txfasta -> txfasta }

    SALMON_QUANT ( ch_simreads, ch_index, ch_filteredgff3_nometa, ch_filteredtxfasta_nometa, ch_alignment_mode, ch_libtype )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    DEANALYSIS ( SALMON_QUANT.out.results, ch_transcriptData )
    ch_versions = ch_versions.mix(DEANALYSIS.out.versions.first())

    ENRICHMENT ( DEANALYSIS.out.deseq2_results, DEANALYSIS.out.deseq2_tx2gene )
    ch_versions = ch_versions.mix(ENRICHMENT.out.versions.first())

    RNASEQVALIDATION ( DEANALYSIS.out.deseq2_results, ch_simreads, ENRICHMENT.out.enrichment_results, DEANALYSIS.out.deseq2_tx2gene)
    ch_versions = ch_versions.mix(RNASEQVALIDATION.out.versions.first())

    emit:
    // SALMON_QUANT outputs
    salmon_results      = SALMON_QUANT.out.results                               // channel: [ val(meta), path(quant_dir)         ]
    salmon_json_info    = SALMON_QUANT.out.json_info                             // channel: [ val(meta), path(json_info)         ], optional
    salmon_lib_format   = SALMON_QUANT.out.lib_format_counts                     // channel: [ val(meta), path(lib_format_counts) ], optional

    // DEANALYSIS outputs
    deseq2_results      = DEANALYSIS.out.deseq2_results                           // channel: [ val(meta), path(deseq2_results.tsv), path(deseq2_de_genes.txt), path(*.pdf) ]
    deseq2_tx2gene      = DEANALYSIS.out.deseq2_tx2gene                           // channel: [ val(meta), path(deseq2_tx2gene.tsv)                                         ]

    // ENRICHMENT outputs
    enrichment_results  = ENRICHMENT.out.enrichment_results                      // channel: [ val(meta), path(enrichment_results.rds), path(*.png) ]

    // RNASEQVALIDATION outputs
    rnaseq_validated_results  = RNASEQVALIDATION.out.rnaseq_validated_results    // channel: [ val(meta), path(rnaseq_validation) ], optional

    versions            = ch_versions                                            // channel: [ versions.yml ]
}

