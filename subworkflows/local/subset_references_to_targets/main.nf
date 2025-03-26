include { GATK4_SELECTVARIANTS     } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { SAMTOOLS_FAIDX           } from '../../../modules/nf-core/samtools/faidx/main'
include { SUBVAR                   } from '../../../modules/local/subvar/main'
include { SUBSETCAPTURE            } from '../../../modules/local/subsetcapture/main'

workflow SUBSET_REFERENCES_TO_TARGETS {

    take:
    ch_fasta         // channel: [ val(meta), [ fasta ] ]
    ch_fai           // channel: [ val(meta), [ fai ] ]
    ch_get_sizes     // channel: [ get_sizes ]
    ch_capture_bed   // channel: [ capture_bed ]
    ch_chromosome    // channel: [ val(meta), val(chromosome) ]
    ch_vcf_idx_intervals // channel: [ val(meta), [ vcf, vcf_idx, intervals ] ]
    ch_clinvar_vcf   // channel: [ val(meta), path(vcf) ]

    main:

    ch_versions = Channel.empty()

    // Index the reference genome if not already indexed
    SAMTOOLS_FAIDX ( ch_fasta.combine(ch_fai), [], false )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Subset capture regions
    SUBSETCAPTURE ( ch_chromosome, SAMTOOLS_FAIDX.out.sizes, ch_capture_bed )
    ch_versions = ch_versions.mix(SUBSETCAPTURE.out.versions)

    // Select variants based on intervals
    GATK4_SELECTVARIANTS ( ch_vcf_idx_intervals )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)

    // Subset clinvar variants
    SUBVAR ( ch_clinvar_vcf, SUBSETCAPTURE.out.target_bed )
    ch_versions = ch_versions.mix(SUBVAR.out.versions)

    emit:
    fa                = SAMTOOLS_FAIDX.out.fa              // channel:  [ val(meta), [ fa, fasta ] ]
    fai               = SAMTOOLS_FAIDX.out.fai              // channel: [ val(meta), [ fai ] ]
    sizes             = SAMTOOLS_FAIDX.out.sizes            // channel: [ val(meta), [ sizes ] ]
    gzi               = SAMTOOLS_FAIDX.out.gzi              // channel: [ val(meta), [ gzi ] ]
    capture_bed_gz    = SUBSETCAPTURE.out.capture_bed_gz    // channel: [ val(meta), [ capture_bed_gz ] ]
    capture_bed_index = SUBSETCAPTURE.out.capture_bed_index // channel: [ val(meta), [ capture_bed_index ] ]
    target_bed        = SUBSETCAPTURE.out.target_bed        // channel: [ val(meta), [ target_bed ] ]
    target_bed_pad50  = SUBSETCAPTURE.out.target_bed_pad50  // channel: [ val(meta), [ target_bed_pad50 ] ]
    target_bed_pad500 = SUBSETCAPTURE.out.target_bed_pad500 // channel: [ val(meta), [ target_bed_pad500 ] ]
    selected_vcf      = GATK4_SELECTVARIANTS.out.vcf        // channel: [ val(meta), [ vcf ] ]
    selected_tbi      = GATK4_SELECTVARIANTS.out.tbi        // channel: [ val(meta), [ tbi ] ]
    benign_vcf        = SUBVAR.out.benign_vcf               // channel: [ val(meta), [ benign_vcf ] ]
    pathogenic_vcf    = SUBVAR.out.pathogenic_vcf           // channel: [ val(meta), [ pathogenic_vcf ] ]
    selected_vcf      = SUBVAR.out.selected_vcf             // channel: [ val(meta), [ selected_vcf ] ]

    versions = ch_versions                                  // channel: [ versions.yml ]
}
