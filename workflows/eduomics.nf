/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                    } from 'plugin/nf-schema'
include { softwareVersionsToYAML              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText              } from '../subworkflows/local/utils_nfcore_eduomics_pipeline'

include { PREPARE_RNA_GENOME                  } from '../subworkflows/local/prepare_rna_genome'
include { SIMULATE_RNASEQ_READS               } from '../subworkflows/local/simulate_rnaseq_reads'
include { QUANTIFY_DEANALYSIS_ENRICH_VALIDATE } from '../subworkflows/local/quantify_deanalysis_enrich_validate'

include { AISCENARIOS                         } from '../modules/local/aiscenarios'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EDUOMICS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input


    main:

    ch_versions  = Channel.empty()
    ch_scenarios = Channel.empty()

    ch_samplesheet
    .branch { m, c ->
        dna: m.type == "dna"
        rna: m.type == "rna"
    }
    .set{ input_bytype_ch }

    // CREATING CHANNELS FROM REFERENCE FILES
    ch_fasta   = Channel.fromPath(params.fasta)
    ch_fai     = Channel.fromPath(params.fai)
    ch_txfasta = Channel.fromPath(params.txfasta)
    ch_gff3    = Channel.fromPath(params.gff3)


    PREPARE_RNA_GENOME(
        input_bytype_ch.rna.map { meta, capture -> meta},
        ch_gff3,
        ch_txfasta,
        ch_fasta
    )

    ch_foldchange = Channel.empty()

    SIMULATE_RNASEQ_READS(
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        PREPARE_RNA_GENOME.out.filtered_transcript_data,
        PREPARE_RNA_GENOME.out.gene_lists,
        ch_foldchange,
        PREPARE_RNA_GENOME.out.gene_list_association
    )

    ch_libtype = Channel.value(params.salmon_libtype)

    QUANTIFY_DEANALYSIS_ENRICH_VALIDATE(
        SIMULATE_RNASEQ_READS.out.simreads,
        PREPARE_RNA_GENOME.out.txfasta_index,
        PREPARE_RNA_GENOME.out.filtered_annotation,
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        params.salmon_alignmode,
        ch_libtype,
        PREPARE_RNA_GENOME.out.filtered_transcript_data
    )

    ch_rna_scenario = QUANTIFY_DEANALYSIS_ENRICH_VALIDATE.out.deseq2_tx2gene
        .map { m, tx ->
            return [m, false, m.genes]
        }
    ch_scenarios = ch_scenarios.mix(ch_rna_scenario)

    AISCENARIOS(ch_scenarios)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'eduomics_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
