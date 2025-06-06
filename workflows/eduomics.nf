/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                    } from 'plugin/nf-schema'
include { softwareVersionsToYAML              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText              } from '../subworkflows/local/utils_nfcore_eduomics_pipeline'

// RNA subworkflows
include { PREPARE_RNA_GENOME                  } from '../subworkflows/local/prepare_rna_genome'
include { SIMULATE_RNASEQ_READS               } from '../subworkflows/local/simulate_rnaseq_reads'
include { QUANTIFY_DEANALYSIS_ENRICH_VALIDATE } from '../subworkflows/local/quantify_deanalysis_enrich_validate'

// DNA subworkflows
include { SUBSET_REFERENCES_TO_TARGETS        } from '../subworkflows/local/subset_references_to_targets'
include { FASTA_WGSIM_TO_PROFILE              } from '../subworkflows/local/fasta_wgsim_to_profile'
include { PROFILE_SIMULATE_VARS_FASTQ         } from '../subworkflows/local/profile_simulate_vars_fastq'
include { FASTQ_VARIANT_TO_VALIDATION         } from '../subworkflows/local/fastq_variant_to_validation'

// common module
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
    ch_fasta       = Channel.fromPath(params.fasta)
    ch_txfasta     = Channel.fromPath(params.txfasta)
    ch_gff3        = Channel.fromPath(params.gff3)
    ch_capture_bed = Channel.fromPath(params.capture)
    ch_gnomad      = Channel.fromPath(params.gnomad)
    ch_mills       = Channel.fromPath(params.mills)
    ch_1000g       = Channel.fromPath(params.vcf1000g)
    ch_dbsnp       = Channel.fromPath(params.dbsnp)
    ch_clinvar     = Channel.fromPath(params.clinvar)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNA ANALYSIS BRANCH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    PREPARE_RNA_GENOME(
        input_bytype_ch.rna.map { meta, capture -> meta},
        ch_gff3,
        ch_txfasta,
        ch_fasta
    )

    ch_versions = ch_versions.mix(PREPARE_RNA_GENOME.out.versions)

    ch_foldchange = Channel.empty()

    SIMULATE_RNASEQ_READS(
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        PREPARE_RNA_GENOME.out.filtered_transcript_data,
        PREPARE_RNA_GENOME.out.gene_lists,
        ch_foldchange,
        PREPARE_RNA_GENOME.out.gene_list_association
    )

    ch_versions = ch_versions.mix(SIMULATE_RNASEQ_READS.out.versions)

    // Set to true in test.config to run the test with a smaller dataset
    if (params.istest) {
        SIMULATE_RNASEQ_READS.out.simreads = SIMULATE_RNASEQ_READS.out.simreads.take(5)
    }

    QUANTIFY_DEANALYSIS_ENRICH_VALIDATE(
        SIMULATE_RNASEQ_READS.out.simreads,
        PREPARE_RNA_GENOME.out.txfasta_index,
        PREPARE_RNA_GENOME.out.filtered_annotation,
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        params.salmon_alignmode,
        params.salmon_libtype,
        PREPARE_RNA_GENOME.out.filtered_transcript_data
    )

    ch_versions = ch_versions.mix(QUANTIFY_DEANALYSIS_ENRICH_VALIDATE.out.versions)

    ch_rna_scenario = QUANTIFY_DEANALYSIS_ENRICH_VALIDATE.out.deseq2_tx2gene
        .map { m, tx ->
            return [m, false, m.genes]
        }
    ch_scenarios = ch_scenarios.mix(ch_rna_scenario)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DNA ANALYSIS BRANCH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    SUBSET_REFERENCES_TO_TARGETS(
        input_bytype_ch.dna.map { meta, capture -> meta},
        ch_fasta,
        params.get_sizes_bool,
        ch_capture_bed,
        ch_gnomad,
        ch_mills,
        ch_1000g,
        ch_dbsnp,
        ch_clinvar
        )

    ch_versions = ch_versions.mix(SUBSET_REFERENCES_TO_TARGETS.out.versions)


    FASTA_WGSIM_TO_PROFILE(
        SUBSET_REFERENCES_TO_TARGETS.out.target_fa,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fai,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dict,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bwa_index,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bed
    )

    ch_versions = ch_versions.mix(FASTA_WGSIM_TO_PROFILE.out.versions)

    ch_fasta_fai = SUBSET_REFERENCES_TO_TARGETS.out.target_fa
            .join(SUBSET_REFERENCES_TO_TARGETS.out.target_fai)
            .map { meta, fasta, fai -> [fasta, fai] }

    PROFILE_SIMULATE_VARS_FASTQ(
        SUBSET_REFERENCES_TO_TARGETS.out.clinvar_benign_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.clinvar_pathogenic_vcf,
        FASTA_WGSIM_TO_PROFILE.out.profile,
        ch_fasta_fai,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bed_pad500
    )

    ch_versions = ch_versions.mix(PROFILE_SIMULATE_VARS_FASTQ.out.versions)

    // Set to true in test.config to run the test with a smaller dataset
    if (params.istest) {
        PROFILE_SIMULATE_VARS_FASTQ.out.simreads = PROFILE_SIMULATE_VARS_FASTQ.out.simreads.take(5)
    }

    FASTQ_VARIANT_TO_VALIDATION(
        PROFILE_SIMULATE_VARS_FASTQ.out.simreads,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fa,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fai,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dict,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bwa_index,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bed
    )

    ch_versions = ch_versions.mix(FASTQ_VARIANT_TO_VALIDATION.out.versions)

    ch_dna_scenario = FASTQ_VARIANT_TO_VALIDATION.out.scenario
        .map { m, var ->
            return [m, m.var, false]
        }
    ch_scenarios = ch_scenarios.mix(ch_dna_scenario)


    AISCENARIOS(ch_scenarios)

    ch_versions = ch_versions.mix(AISCENARIOS.out.versions)


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
    versions                 = ch_collated_versions                                               // channel: [ path(versions.yml)                          ]
    fastq_validated_variants = FASTQ_VARIANT_TO_VALIDATION.out.simulation                         // channel: [ val(meta), path(validated_results_folder/*) ]
    rnaseq_validated_reads   = QUANTIFY_DEANALYSIS_ENRICH_VALIDATE.out.rnaseq_validated_results   // channel: [ val(meta), path(rnaseq_validation)          ]
    scenario_description     = AISCENARIOS.out.scenario                                           // channel: [ val(meta), path(scenario.txt)               ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
