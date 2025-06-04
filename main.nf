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

    //
    // WORKFLOW: Run pipeline
    //
    EDUOMICS (
        samplesheet
    )

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
    softwareversions = EDUOMICS.out.versions
    //simulation results
    dnasimulation    = EDUOMICS.out.fastq_validated_variants
    rnasimulation    = EDUOMICS.out.rnaseq_validated_reads
    scenario         = EDUOMICS.out.scenario_description
    // dna reference needed for the analysis of simulated DNA data
    dnareference     = SUBSET_REFERENCES_TO_TARGETS.out.target_fa
    dnareferenceindex = SUBSET_REFERENCES_TO_TARGETS.out.target_fai
    // gatk bundke needed for the analysis of simulated DNA data



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
