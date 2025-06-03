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
    ch_fai         = Channel.fromPath(params.fai)
    ch_txfasta     = Channel.fromPath(params.txfasta)
    ch_gff3        = Channel.fromPath(params.gff3)
    ch_capture_bed = Channel.fromPath(params.capture_bed)
    ch_gnomad      = Channel.fromPath(params.gnomad)
    ch_mills       = Channel.fromPath(params.mills)
    ch_1000g       = Channel.fromPath(params.vcf1000g)
    ch_dbsnp       = Channel.fromPath(params.dbsnp)
    ch_clinvar     = Channel.fromPath(params.clinvar)
    ch_bwa_index   = Channel.fromPath(params.bwa_index)

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

    ch_foldchange = Channel.empty()

    SIMULATE_RNASEQ_READS(
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        PREPARE_RNA_GENOME.out.filtered_transcript_data,
        PREPARE_RNA_GENOME.out.gene_lists,
        ch_foldchange,
        PREPARE_RNA_GENOME.out.gene_list_association
    )

    QUANTIFY_DEANALYSIS_ENRICH_VALIDATE(
        SIMULATE_RNASEQ_READS.out.simreads,
        PREPARE_RNA_GENOME.out.txfasta_index,
        PREPARE_RNA_GENOME.out.filtered_annotation,
        PREPARE_RNA_GENOME.out.filtered_txfasta,
        params.salmon_alignmode,
        params.salmon_libtype,
        PREPARE_RNA_GENOME.out.filtered_transcript_data
    )

    ch_rna_scenario = QUANTIFY_DEANALYSIS_ENRICH_VALIDATE.out.deseq2_tx2gene
        .map { m, tx ->
            return [m, false, m.genes]
        }
    ch_scenarios = ch_scenarios.mix(ch_rna_scenario)

    AISCENARIOS(ch_scenarios)

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

    FASTA_WGSIM_TO_PROFILE(
        SUBSET_REFERENCES_TO_TARGETS.out.target_fa,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fai,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dict,
        ch_bwa_index,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bed
    )

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

    FASTQ_VARIANT_TO_VALIDATION(
        PROFILE_SIMULATE_VARS_FASTQ.out.simreads,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fa,
        SUBSET_REFERENCES_TO_TARGETS.out.target_fai,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dict,
        ch_bwa_index,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_dbsnp_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_vcf,
        SUBSET_REFERENCES_TO_TARGETS.out.target_mills_tbi,
        SUBSET_REFERENCES_TO_TARGETS.out.target_bed
    )

    ch_dna_scenario = FASTQ_VARIANT_TO_VALIDATION.out.scenario
        .map { m, var ->
            return [m, m.var, false]
        }
    ch_scenarios = ch_scenarios.mix(ch_dna_scenario)

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
