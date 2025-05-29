include { SUBVAR                                              } from '../../../modules/local/subvar/main'
include { SUBSETCAPTURE                                       } from '../../../modules/local/subsetcapture/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SUBSET             } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_INDEX              } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SIZES              } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_GNOMAD } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_MILLS  } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_1000G  } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_DBSNP  } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_CREATESEQUENCEDICTIONARY                      } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'



workflow SUBSET_REFERENCES_TO_TARGETS {

    take:
    ch_meta          // channel: [ val(meta)   ]
    ch_fasta         // channel: [ path(fasta) ]
    ch_fai           // channel: [ path(fai)   ]
    ch_get_sizes     // channel: [ get_sizes   ] // boolean
    ch_capture_bed   // channel: [ path(capture_bed) ]
    ch_gnomad        // channel: [ path(vcf), path(vcf_idx) ]
    ch_mills         // channel: [ path(vcf), path(vcf_idx) ]
    ch_1000g         // channel: [ path(vcf), path(vcf_idx) ]
    ch_dbsnp         // channel: [ path(vcf), path(vcf_idx) ]
    ch_clinvar_vcf   // channel: [ path(vcf) ]

    main:

    ch_versions = Channel.empty()

    ch_fasta_with_meta = ch_meta.combine(ch_fasta).map { meta, fasta -> [meta, fasta] }
    ch_fai_with_meta = ch_meta.combine(ch_fai).map { meta, fai -> [meta, fai] }

    SAMTOOLS_FAIDX_SUBSET ( ch_fasta_with_meta, [[],[]], ch_get_sizes )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SUBSET.out.versions)

    SAMTOOLS_FAIDX_INDEX (SAMTOOLS_FAIDX_SUBSET.out.fa, [[],[]], ch_get_sizes )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_INDEX.out.versions)

    SAMTOOLS_FAIDX_SIZES ( SAMTOOLS_FAIDX_SUBSET.out.fa, SAMTOOLS_FAIDX_INDEX.out.fai, true )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SIZES.out.versions)

    GATK4_CREATESEQUENCEDICTIONARY ( SAMTOOLS_FAIDX_SUBSET.out.fa )
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    //Subset capture regions by chrom
    SUBSETCAPTURE (
        ch_meta,
        SAMTOOLS_FAIDX_SIZES.out.sizes.map { meta, sizes -> sizes },
        ch_capture_bed
        )

    ch_versions = ch_versions.mix(SUBSETCAPTURE.out.versions)

    // Subset gatk bundle vcf
    GATK4_SELECTVARIANTS_GNOMAD ( SUBSETCAPTURE.out.target_bed.combine(ch_gnomad).map { meta, target_bed, vcf, vcf_idx ->
        [meta, vcf, vcf_idx, target_bed]
    })
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_GNOMAD.out.versions)

    GATK4_SELECTVARIANTS_MILLS ( SUBSETCAPTURE.out.target_bed.combine(ch_mills).map { meta, target_bed, vcf, vcf_idx ->
        [meta, vcf, vcf_idx, target_bed]
    })
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_MILLS.out.versions)

    GATK4_SELECTVARIANTS_1000G ( SUBSETCAPTURE.out.target_bed.combine(ch_1000g).map { meta, target_bed, vcf, vcf_idx ->
        [meta, vcf, vcf_idx, target_bed]
    })
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_1000G.out.versions)

    GATK4_SELECTVARIANTS_DBSNP (SUBSETCAPTURE.out.target_bed.combine(ch_dbsnp).map { meta, target_bed, vcf, vcf_idx ->
        [meta, vcf, vcf_idx, target_bed]
    })
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_DBSNP.out.versions)

    // Subset clinvar variants
    SUBVAR (
    ch_meta.combine(ch_clinvar_vcf).map { meta, vcf -> [meta, vcf] },
    SUBSETCAPTURE.out.target_bed.map { meta, target_bed -> target_bed }
    )

    ch_versions = ch_versions.mix(SUBVAR.out.versions)


    emit:
    target_fa         = SAMTOOLS_FAIDX_SUBSET.out.fa                                 // channel:  [ val(meta), [ fa, fasta ] ]
    target_fai        = SAMTOOLS_FAIDX_INDEX.out.fai                                 // channel: [ val(meta), [ fai ] ]
    target_sizes      = SAMTOOLS_FAIDX_SIZES.out.sizes                               // channel: [ val(meta), [ sizes ] ]
    capture_bed_gz    = SUBSETCAPTURE.out.capture_bed_gz.map { meta, bed -> bed }    // channel: [ [ capture_bed_gz ] ]
    capture_bed_index = SUBSETCAPTURE.out.capture_bed_index.map { meta, bed -> bed } // channel: [ [ capture_bed_index ] ]
    target_bed        = SUBSETCAPTURE.out.target_bed.map { meta, bed -> bed }        // channel: [ [ target_bed ] ]
    target_bed_pad50  = SUBSETCAPTURE.out.target_bed_pad50.map { meta, bed -> bed }  // channel: [ [ target_bed_pad50 ] ]
    target_bed_pad500 = SUBSETCAPTURE.out.target_bed_pad500.map { meta, bed -> bed } // channel: [ [ target_bed_pad500 ] ]
    target_gnomad_vcf = GATK4_SELECTVARIANTS_GNOMAD.out.vcf.map { meta, vcf -> vcf } // channel: [ [ vcf ] ]
    target_gnomad_tbi = GATK4_SELECTVARIANTS_GNOMAD.out.tbi.map { meta, tbi -> tbi } // channel: [ [ tbi ] ]
    target_mills_vcf  = GATK4_SELECTVARIANTS_MILLS.out.vcf.map { meta, vcf -> vcf }  // channel: [ [ vcf ] ]
    target_mills_tbi  = GATK4_SELECTVARIANTS_MILLS.out.tbi.map { meta, tbi -> tbi }  // channel: [ [ tbi ] ]
    target_1000g_vcf  = GATK4_SELECTVARIANTS_1000G.out.vcf.map { meta, vcf -> vcf }  // channel: [ [ vcf ] ]
    target_1000g_tbi  = GATK4_SELECTVARIANTS_1000G.out.tbi.map { meta, tbi -> tbi }  // channel: [ [ tbi ] ]
    target_dbsnp_vcf  = GATK4_SELECTVARIANTS_DBSNP.out.vcf.map { meta, vcf -> vcf }  // channel: [ [ vcf ] ]
    target_dbsnp_tbi  = GATK4_SELECTVARIANTS_DBSNP.out.tbi.map { meta, tbi -> tbi }  // channel: [ [ tbi ] ]
    clinvar_benign_vcf      = SUBVAR.out.benign_vcf                                  // channel: [ val(meta), [ benign_vcf ] ]
    clinvar_pathogenic_vcf  = SUBVAR.out.pathogenic_vcf                              // channel: [ val(meta), [ pathogenic_vcf ] ]
    clinvar_selected_vcf    = SUBVAR.out.selected_vcf                                // channel: [ val(meta), [ selected_vcf ] ]
    target_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    versions = ch_versions                                                           // channel: [ versions.yml ]
}



