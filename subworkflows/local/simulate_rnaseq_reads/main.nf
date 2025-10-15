include { COUNTMATRICES      } from '../../../modules/local/countmatrices/main'
include { POLYESTER_SIMULATE } from '../../../modules/local/polyester/simulate/main'

workflow SIMULATE_RNASEQ_READS {

    take:
    ch_filtered_txfasta            // channel: [ val(meta), path(filtered_txfasta)        ]
    ch_transcriptData              // channel: [ val(meta), path(transcriptData)          ]
    ch_genelists                   // channel: [ val(meta), path(genelists)               ]
    ch_gene_list_association       // channel: [ val(meta), path(gene_list_association)   ]

    main:

    ch_versions = Channel.empty()

    ch_gene_list_association.dump(tag: 'gene list association')

    // Generate count matrices
    COUNTMATRICES(ch_filtered_txfasta, ch_transcriptData, ch_genelists)
    ch_versions = ch_versions.mix(COUNTMATRICES.out.versions)

    // Build gene map
    genesMap = ch_gene_list_association
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def list_name = row[1].gene_list
            def genes = row[1].genes.split(',')
            [list_name, genes]
        }
    genesMap.dump(tag: 'genes map processed')

    // Associate the genes to countmatrices and define a new meta
    ch_matrices_with_genes = COUNTMATRICES.out.simcountMatrix
        .combine(genesMap)
        .flatMap { meta, path_list, list_key, genes ->
            path_list
                .findAll { path -> list_key == "list${path.name.replaceAll(/\D+/, '')}" }
                .collect { filterd_path ->
                def newmeta = meta + [genes: genes.join(',')]
                [newmeta, filterd_path] }
        }
    ch_matrices_with_genes.dump(tag: 'count matrices with genes')

    ch_fold_change = Channel.value([[id: 'null'], []]) // this simulation is currently not implemented

    ch_matrices_with_genes_limited = params.istest
        ? ch_matrices_with_genes.take(params.test_limit)
        : ch_matrices_with_genes

    ch_matrices_with_genes.view { "genes: ${it[0]}" }

    // Simulate the reads
    POLYESTER_SIMULATE(ch_matrices_with_genes_limited, ch_fold_change, ch_filtered_txfasta)
    ch_versions = ch_versions.mix(POLYESTER_SIMULATE.out.versions)

    emit:
    countMatrix    = COUNTMATRICES.out.simcountMatrix  // channel: [ val(meta),    path(countMatrix)    ]
    expectedLog2FC = COUNTMATRICES.out.simAnnords      // channel: [ val(meta),    path(expectedLog2FC) ]
    simreads       = POLYESTER_SIMULATE.out.reads      // channel: [ val(newmeta), path(simreads)       ]. The newmeta will look like this [ id: 'test', chromosome: 'chr22', coverage: 30, reps: 3, groups: 2, simthreshold: 0.3, genes: 'A,B,C' ]
    versions       = ch_versions                       // channel: [ versions.yml                       ]
}

