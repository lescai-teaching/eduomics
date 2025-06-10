#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/eduomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/eduomics
    Website: https://nf-co.re/eduomics
    Slack  : https://nfcore.slack.com/channels/eduomics
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.preview.output = true

include { EDUOMICS                     } from './workflows/eduomics'
include { SUBSET_REFERENCES_TO_TARGETS } from './subworkflows/local/subset_references_to_targets'
include { PREPARE_RNA_GENOME           } from './subworkflows/local/prepare_rna_genome'
include { PIPELINE_INITIALISATION      } from './subworkflows/local/utils_nfcore_eduomics_pipeline'
include { PIPELINE_COMPLETION          } from './subworkflows/local/utils_nfcore_eduomics_pipeline'
include { getGenomeAttribute           } from './subworkflows/local/utils_nfcore_eduomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// getting attributes for genome files
params.fasta       = getGenomeAttribute('fasta')
params.gnomad      = getGenomeAttribute('germline_resource')
params.mills       = getGenomeAttribute('known_indels')
params.vcf1000g    = getGenomeAttribute('known_snps')
params.dbsnp       = getGenomeAttribute('dbsnp')
params.clinvar     = getGenomeAttribute('clinvar_vcf')
params.gff3        = getGenomeAttribute('gff3')
params.txfasta     = getGenomeAttribute('txfasta')
params.bwa_index   = getGenomeAttribute('bwa_index')



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_EDUOMICS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    samplesheet
    .branch { m, c ->
        dna: m.type == "dna"
        rna: m.type == "rna"
    }
    .set{ input_bytype_ch }

    // CREATING CHANNELS FROM REFERENCE FILES
    ch_fasta       = Channel.fromPath(params.fasta)
    ch_txfasta     = Channel.fromPath(params.txfasta)
    ch_gff3        = Channel.fromPath(params.gff3)
    ch_capture_bed = input_bytype_ch.dna.map { meta, capture -> capture }
    ch_gnomad      = Channel.fromPath(params.gnomad)
    ch_mills       = Channel.fromPath(params.mills)
    ch_1000g       = Channel.fromPath(params.vcf1000g)
    ch_dbsnp       = Channel.fromPath(params.dbsnp)
    ch_clinvar     = Channel.fromPath(params.clinvar)

    //
    // WORKFLOW: Run pipeline
    //
    EDUOMICS (
        input_bytype_ch.dna,
        input_bytype_ch.rna,
        ch_fasta,
        ch_txfasta,
        ch_gff3,
        ch_capture_bed,
        ch_gnomad,
        ch_mills,
        ch_1000g,
        ch_dbsnp,
        ch_clinvar
    )

    emit:
    versions                 = EDUOMICS.out.versions                   // channel: [ path(versions.yml)                          ]
    fastq_validated_variants = EDUOMICS.out.fastq_validated_variants   // channel: [ val(meta), path(validated_results_folder/*) ]
    rnaseq_validated_reads   = EDUOMICS.out.rnaseq_validated_reads     // channel: [ val(meta), path(rnaseq_validation)          ]
    dnabundle                = EDUOMICS.out.dnabundle                  // channel: [ val(meta), [all references bundle] ]
    rnabundle                = EDUOMICS.out.rnabundle                  // channel: [ val(meta), [path(txfasta), path(gff3), path(salmonindex)] ]
    scenario_description     = EDUOMICS.out.scenario_description       // channel: [ val(meta), path(scenario.txt)               ]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_EDUOMICS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )

    publish:
    // software versions
    softwareversions = NFCORE_EDUOMICS.out.versions
    //simulation results
    dnasimulation    = NFCORE_EDUOMICS.out.fastq_validated_variants
    rnasimulation    = NFCORE_EDUOMICS.out.rnaseq_validated_reads
    scenario         = NFCORE_EDUOMICS.out.scenario_description
    // dna reference and bundle needed for the analysis of simulated DNA data
    dnabundle        = NFCORE_EDUOMICS.out.dnabundle
    // rna reference and bundle needed for the analysis of simulated RNA data
    rnabundle        = NFCORE_EDUOMICS.out.rnabundle

}

output {

    softwareversions {
        path "${params.outdir}/pipeline_info"
    }

    dnasimulation {
        path { meta, files ->
            "dna_simulations/${meta.id}/${meta.simulatedvar}"
        }
    }

    dnabundle {
        path { meta, files ->
            "dna_simulations/${meta.id}/references"
        }
    }

    rnasimulation {
        path { meta, files ->
            def simfolder = meta.genes.take(16)
            "rna_simulations/${meta.id}/${simfolder}"
        }
    }

    rnabundle {
        path { meta, files ->
            "rna_simulations/${meta.id}/references"
        }
    }

    scenario {
        path { meta, text ->
            if (meta.type == "dna"){
                "dna_simulations/${meta.id}/${meta.simulatedvar}"
            }
            else {
                def simfolder = meta.genes.take(16)
                "rna_simulations/${meta.id}/${simfolder}"
            }
        }
    }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
