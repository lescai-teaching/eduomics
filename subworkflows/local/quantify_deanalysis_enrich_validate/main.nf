
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
    ch_tx2gene                    // channel: [ val(meta), path(tx2gene)                           ]

    main:

    ch_versions = Channel.empty()

    ch_simreads.dump(tag: "rna simreads")

    // Modify `meta.id` to make it unique for each sample based on filenames
    ch_simreads_modified = ch_simreads.map{ m, files ->
                def grouped = files.groupBy { file ->
                    file.name.replaceFirst(/_[12]\.fasta\.gz$/, '')
                    }
                grouped.collect { sampleName, group ->
                    def updatedMeta = m.clone()
                    updatedMeta.original_id = m.id
                    updatedMeta.id = sampleName
                    [updatedMeta, group.sort()]
                    }
            }
            .flatMap()

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
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

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
    DEANALYSIS ( ch_deanalysis_input, ch_tx2gene )
    ch_versions = ch_versions.mix(DEANALYSIS.out.versions)

    DEANALYSIS.out.deseq2_results.dump(tag: 'deseq2 results')

    // Extract only deseq2_results.tsv for enrichment module
    ch_deseq2_results_tsv_only = DEANALYSIS.out.deseq2_results.map { meta, tsv, txt, pdfs ->
        tuple(meta, tsv)
    }

    ch_deseq2_results_tsv_only.dump(tag: 'enrichment input')

    // Run enrichment analysis
    ENRICHMENT ( ch_deseq2_results_tsv_only, ch_tx2gene )
    ch_versions = ch_versions.mix(ENRICHMENT.out.versions)

    ENRICHMENT.out.enrichment_results.dump(tag: 'enrichment results')

    // Join channels
    ch_validation_input = ch_simreads
    .join(DEANALYSIS.out.deseq2_results, by: 0)
    .join(ENRICHMENT.out.enrichment_results, by: 0)

    ch_validation_input.dump(tag: 'validation CH')

    // Perform final validation of results
    RNASEQVALIDATION ( ch_validation_input )
    ch_versions = ch_versions.mix(RNASEQVALIDATION.out.versions)

    RNASEQVALIDATION.out.rnaseq_validated_results.dump(tag: 'rnaseq validation output')

    ch_scenario = RNASEQVALIDATION.out.rnaseq_validated_results
        .map { m, folder ->
            return [m, false, m.genes]
        }

    ch_scenario.dump(tag: 'scenario CH')

    emit:
    // SALMON_QUANT outputs
    salmon_results            = SALMON_QUANT.out.results                          // channel: [ val(meta), path(quant_dir)         ]
    salmon_json_info          = SALMON_QUANT.out.json_info                        // channel: [ val(meta), path(json_info)         ], optional
    salmon_lib_format         = SALMON_QUANT.out.lib_format_counts                // channel: [ val(meta), path(lib_format_counts) ], optional

    // DEANALYSIS outputs
    deseq2_results            = DEANALYSIS.out.deseq2_results                     // channel: [ val(meta), path(deseq2_results), path(deseq2_de_genes), path(*.pdf) ]

    // ENRICHMENT outputs
    enrichment_results        = ENRICHMENT.out.enrichment_results                 // channel: [ val(meta), path(enrichment_results), path(*.png) ]

    // RNASEQVALIDATION outputs
    rnaseq_validated_results  = RNASEQVALIDATION.out.rnaseq_validated_results    // channel: [ val(meta), path(rnaseq_validation) ], optional

    scenario                  = ch_scenario                                      // channel: [ val(meta), false, val(genes) ]

    versions                  = ch_versions                                      // channel: [ versions.yml ]
}

