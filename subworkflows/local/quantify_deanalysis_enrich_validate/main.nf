
include { SALMON_QUANT         } from '../../../modules/nf-core/salmon/quant/main'
include { DEANALYSIS           } from '../../../modules/local/deanalysis/main'
include { ENRICHMENT           } from '../../../modules/local/enrichment/main'
include { RNASEQVALIDATION     } from '../../../modules/local/rnaseqvalidation/main'

workflow QUANTIFY_DEANALYSIS_ENRICH_VALIDATE {

    take:
    ch_simreads                   // channel: [ val(meta), path(simreads)                          ]
    ch_index                      // channel: [ path(salmon_index)                                 ]
    ch_filteredgff3               // channel: [ val(meta), path(filteredgff3)                      ]
    ch_filteredtxfasta            // channel: [ val(meta), path(filteredtxfasta)                   ]
    ch_alignment_mode             // value  : [ boolean                                            ]
    ch_libtype                    // value  : [ val(libtype)                                       ]
    ch_filtered_transcriptData    // channel: [ val(meta), path(filtered_transcriptData)           ]

    main:

    ch_versions = Channel.empty()

    // Modify `meta.id` to make it unique for each sample based on filenames
    ch_simreads_modified = ch_simreads.map { meta, reads ->
        def new_meta = meta.clone()                                  // clone the original meta to mantain it. It is required in the following steps
        def base_name = reads[0].getBaseName()                       // e.g. "sample_01_1"
        def sample = base_name.split('_')[0..1].join('_')            // e.g. "sample_01"
        new_meta.original_id = meta.id
        new_meta.id = sample                                         // update meta.id with sample name
        tuple(new_meta, reads)
    }

    ch_simreads_modified.dump(tag: "rna simreads modified")

    // Strip away meta for filtered gff3 and transcript fasta
    ch_filteredgff3_nometa = ch_filteredgff3.map { meta, gff3 -> gff3 }
    ch_filteredtxfasta_nometa = ch_filteredtxfasta.map { meta, txfasta -> txfasta }

    // Run quantification
    SALMON_QUANT (
        ch_simreads_modified,
        ch_index,
        ch_filteredgff3_nometa,
        ch_filteredtxfasta_nometa,
        ch_alignment_mode,
        ch_libtype
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.ifEmpty([]).first())

    SALMON_QUANT.out.results.dump(tag: 'salmon quant results')

    // Restore original `meta.id` for downstream DE analysis (grouping by meta)
    ch_deanalysis_input = SALMON_QUANT.out.results.map { meta, quant_dir ->
        def cleaned_meta = meta.clone()
        cleaned_meta.id = meta.original_id
        cleaned_meta.remove('original_id')
        cleaned_meta.remove('sample')
        tuple(cleaned_meta, quant_dir)
    }
    .groupTuple(by: 0)

    ch_deanalysis_input.dump(tag: 'quant dirs')

    // Run differential expression analysis
    DEANALYSIS ( ch_deanalysis_input, ch_filtered_transcriptData )
    ch_versions = ch_versions.mix(DEANALYSIS.out.versions.ifEmpty([]))

    // Extract only deseq2_results.tsv for enrichment module
    ch_deseq2_results_tsv_only = DEANALYSIS.out.deseq2_results.map { meta, tsv, file1, list1 ->
        tuple(meta, tsv)
    }

    // Run enrichment analysis
    ENRICHMENT ( ch_deseq2_results_tsv_only, DEANALYSIS.out.deseq2_tx2gene )
    ch_versions = ch_versions.mix(ENRICHMENT.out.versions.ifEmpty([]))

    // Perform final validation of results
    RNASEQVALIDATION (
        ch_simreads,
        DEANALYSIS.out.deseq2_results,
        ENRICHMENT.out.enrichment_results,
        DEANALYSIS.out.deseq2_tx2gene
    )
    ch_versions = ch_versions.mix(RNASEQVALIDATION.out.versions.ifEmpty([]))

    emit:
    // SALMON_QUANT outputs
    salmon_results      = SALMON_QUANT.out.results                               // channel: [ val(meta), path(quant_dir)         ]
    salmon_json_info    = SALMON_QUANT.out.json_info                             // channel: [ val(meta), path(json_info)         ], optional
    salmon_lib_format   = SALMON_QUANT.out.lib_format_counts                     // channel: [ val(meta), path(lib_format_counts) ], optional

    // DEANALYSIS outputs
    deseq2_results      = DEANALYSIS.out.deseq2_results                           // channel: [ val(meta), path(deseq2_results), path(deseq2_de_genes), path(*.pdf) ]
    deseq2_tx2gene      = DEANALYSIS.out.deseq2_tx2gene                           // channel: [ val(meta), path(deseq2_tx2gene)                                     ]

    // ENRICHMENT outputs
    enrichment_results  = ENRICHMENT.out.enrichment_results                      // channel: [ val(meta), path(enrichment_results), path(*.png) ]

    // RNASEQVALIDATION outputs
    rnaseq_validated_results  = RNASEQVALIDATION.out.rnaseq_validated_results    // channel: [ val(meta), path(rnaseq_validation) ], optional

    versions            = ch_versions                                            // channel: [ versions.yml ]
}

