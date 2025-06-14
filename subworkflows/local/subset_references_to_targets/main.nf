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
include { BWA_INDEX                                           } from '../../../modules/nf-core/bwa/index'




workflow SUBSET_REFERENCES_TO_TARGETS {

    take:
    ch_meta          // channel: [ val(meta)   ]
    ch_fasta         // channel: [ path(fasta) ]
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

    SAMTOOLS_FAIDX_SUBSET ( ch_fasta_with_meta, [[],[]], ch_get_sizes )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SUBSET.out.versions)

    SAMTOOLS_FAIDX_INDEX (SAMTOOLS_FAIDX_SUBSET.out.fa, [[],[]], ch_get_sizes )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_INDEX.out.versions)

    SAMTOOLS_FAIDX_SIZES ( SAMTOOLS_FAIDX_SUBSET.out.fa, SAMTOOLS_FAIDX_INDEX.out.fai, true )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SIZES.out.versions)

    BWA_INDEX( SAMTOOLS_FAIDX_SUBSET.out.fa )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

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

    ch_dna_bundle = SAMTOOLS_FAIDX_SUBSET.out.fa
        .join(SAMTOOLS_FAIDX_INDEX.out.fai)
        .join(SUBSETCAPTURE.out.target_bed)
        .join(SUBSETCAPTURE.out.target_bed_pad50)
        .join(GATK4_SELECTVARIANTS_GNOMAD.out.vcf)
        .join(GATK4_SELECTVARIANTS_GNOMAD.out.tbi)
        .join(GATK4_SELECTVARIANTS_MILLS.out.vcf)
        .join(GATK4_SELECTVARIANTS_MILLS.out.tbi)
        .join(GATK4_SELECTVARIANTS_1000G.out.vcf)
        .join(GATK4_SELECTVARIANTS_1000G.out.tbi)
        .join(GATK4_SELECTVARIANTS_DBSNP.out.vcf)
        .join(GATK4_SELECTVARIANTS_DBSNP.out.tbi)
        .join(GATK4_CREATESEQUENCEDICTIONARY.out.dict)
        .join(BWA_INDEX.out.index)
        .map { it ->
            [it[0], it[1..-1]]
        }


    emit:
    target_fa               = SAMTOOLS_FAIDX_SUBSET.out.fa.collect()                                 // channel: [ val(meta), [ fasta ]               ]
    target_fai              = SAMTOOLS_FAIDX_INDEX.out.fai.collect()                                 // channel: [ val(meta), [ fai ]                 ]
    target_sizes            = SAMTOOLS_FAIDX_SIZES.out.sizes.collect()                               // channel: [ val(meta), [ sizes ]               ]
    target_bwa_index        = BWA_INDEX.out.index.collect()                                          // channel: [ val(meta), [ bwa ]                 ]
    capture_bed_gz          = SUBSETCAPTURE.out.capture_bed_gz.map    { meta, bed -> bed }.collect() // channel: [ [ capture_bed_gz ]                 ]
    capture_bed_index       = SUBSETCAPTURE.out.capture_bed_index.map { meta, bed -> bed }.collect() // channel: [ [ capture_bed_index ]              ]
    target_bed              = SUBSETCAPTURE.out.target_bed.map        { meta, bed -> bed }.collect() // channel: [ [ target_bed ]                     ]
    target_bed_pad50        = SUBSETCAPTURE.out.target_bed_pad50.map  { meta, bed -> bed }.collect() // channel: [ [ target_bed_pad50 ]               ]
    target_bed_pad500       = SUBSETCAPTURE.out.target_bed_pad500.map { meta, bed -> bed }.collect() // channel: [ [ target_bed_pad500 ]              ]
    target_gnomad_vcf       = GATK4_SELECTVARIANTS_GNOMAD.out.vcf.map { meta, vcf -> vcf }.collect() // channel: [ [ vcf ]                            ]
    target_gnomad_tbi       = GATK4_SELECTVARIANTS_GNOMAD.out.tbi.map { meta, tbi -> tbi }.collect() // channel: [ [ tbi ]                            ]
    target_mills_vcf        = GATK4_SELECTVARIANTS_MILLS.out.vcf.map  { meta, vcf -> vcf }.collect() // channel: [ [ vcf ]                            ]
    target_mills_tbi        = GATK4_SELECTVARIANTS_MILLS.out.tbi.map  { meta, tbi -> tbi }.collect() // channel: [ [ tbi ]                            ]
    target_1000g_vcf        = GATK4_SELECTVARIANTS_1000G.out.vcf.map  { meta, vcf -> vcf }.collect() // channel: [ [ vcf ]                            ]
    target_1000g_tbi        = GATK4_SELECTVARIANTS_1000G.out.tbi.map  { meta, tbi -> tbi }.collect() // channel: [ [ tbi ]                            ]
    target_dbsnp_vcf        = GATK4_SELECTVARIANTS_DBSNP.out.vcf.map  { meta, vcf -> vcf }.collect() // channel: [ [ vcf ]                            ]
    target_dbsnp_tbi        = GATK4_SELECTVARIANTS_DBSNP.out.tbi.map  { meta, tbi -> tbi }.collect() // channel: [ [ tbi ]                            ]
    clinvar_benign_vcf      = SUBVAR.out.benign_vcf.collect()                                        // channel: [ val(meta), [ benign_vcf ]          ]
    clinvar_pathogenic_vcf  = SUBVAR.out.pathogenic_vcf.collect()                                    // channel: [ val(meta), [ pathogenic_vcf ]      ]
    clinvar_selected_vcf    = SUBVAR.out.selected_vcf.collect()                                      // channel: [ val(meta), [ selected_vcf ]        ]
    target_dict             = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()                      // channel: [ val(meta), [ dict ]                ]
    versions                = ch_versions                                                            // channel: [ versions.yml                       ]
    dna_bundle              = ch_dna_bundle.collect()                                                // channel: [ val(meta), [all references bundle] ]
}
