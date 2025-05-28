include { SUBVAR                                   } from '../../../modules/local/subvar/main'
include { SUBSETCAPTURE                            } from '../../../modules/local/subsetcapture/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SUBSET  } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_INDEX   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SIZES   } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_GNOMAD } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_MILLS  } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_1000G  } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_DBSNP  } from '../../../modules/nf-core/gatk4/selectvariants/main'


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

    // Subset the reference genome, index and get chrom sizes
    //SAMTOOLS_FAIDX ( ch_fasta.combine(ch_fai), [], false )
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions

    //SAMTOOLS_FAIDX_SUBSET ( ch_fasta, [[id:'null'], []], false )
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SUBSET.out.versions)
    // false: do not generate chrom.sizes file

    ch_fasta_with_meta = ch_meta.combine(ch_fasta).map { meta, fasta -> [meta, fasta] }
    ch_fai_with_meta = ch_meta.combine(ch_fai).map { meta, fai -> [meta, fai] }

    SAMTOOLS_FAIDX_SUBSET ( ch_fasta_with_meta, [[],[]], false )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SUBSET.out.versions)
    // [[id:'null'], []]

    //SAMTOOLS_FAIDX_INDEX ( SAMTOOLS_FAIDX_SUBSET.out.fa, [], false )
    // index fasta after subsetting

    SAMTOOLS_FAIDX_INDEX (SAMTOOLS_FAIDX_SUBSET.out.fa, [[],[]], false )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_INDEX.out.versions)

    SAMTOOLS_FAIDX_SIZES ( SAMTOOLS_FAIDX_SUBSET.out.fa, SAMTOOLS_FAIDX_INDEX.out.fai, true )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SIZES.out.versions)
    // generating the .sizes from .fai index

    //Subset capture regions by chrom
    SUBSETCAPTURE (
        ch_meta,
        SAMTOOLS_FAIDX_SIZES.out.sizes.map{ map, sizes -> [sizes]},
        ch_capture_bed
        )

    ch_versions = ch_versions.mix(SUBSETCAPTURE.out.versions)

    // Select variants based on intervals (target capture BED) - GATK BUNDLE
    // TO FIX: to do for several VCF files!
    //GATK4_SELECTVARIANTS ( ch_vcf_idx_intervals )
    GATK4_SELECTVARIANTS_GNOMAD ( ch_gnomad.combine(SUBSETCAPTURE.out.target_bed.map{ map, target_bed -> [target_bed]}) )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_GNOMAD.out.versions)


    out_gnomad_vcf = GATK4_SELECTVARIANTS_GNOMAD.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz$/, "_${meta.chromosome}_gnomad.vcf.gz") ]
    }
    out_gnomad_tbi = GATK4_SELECTVARIANTS_GNOMAD.out.vcf.map { meta, tbi ->
        [ meta, tbi.getName().replaceAll(/.vcf.gz.tbi$/, "_${meta.chromosome}_gnomad.vcf.gz.tbi") ]
    }

    GATK4_SELECTVARIANTS_MILLS  ( ch_mills.combine(SUBSETCAPTURE.out.target_bed.map{ map, target_bed -> [target_bed]}) )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_MILLS.out.versions)

    out_mills_vcf = GATK4_SELECTVARIANTS_MILLS.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz$/, "_${meta.chromosome}_mills.vcf.gz") ]
    }
    out_mills_tbi = GATK4_SELECTVARIANTS_MILLS.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz.tbi$/, "_${meta.chromosome}_mills.vcf.gz.tbi") ]
    }

    GATK4_SELECTVARIANTS_1000G  ( ch_1000g.combine(SUBSETCAPTURE.out.target_bed.map{ map, target_bed -> [target_bed]}) )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_1000G.out.versions)


    out_1000G_vcf = GATK4_SELECTVARIANTS_1000G.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz$/, "_${meta.chromosome}_1000G.vcf.gz") ]
    }
    out_1000G_tbi = GATK4_SELECTVARIANTS_1000G.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz.tbi$/, "_${meta.chromosome}_1000G.vcf.gz.tbi") ]
    }

    GATK4_SELECTVARIANTS_DBSNP  ( ch_dbsnp.combine(SUBSETCAPTURE.out.target_bed.map{ map, target_bed -> [target_bed]}) )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_DBSNP.out.versions)

    out_dbsnp_vcf = GATK4_SELECTVARIANTS_DBSNP.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz$/, "_${meta.chromosome}_dbsnp.vcf.gz") ]
    }
    out_dbsnp_tbi = GATK4_SELECTVARIANTS_DBSNP.out.vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/.vcf.gz.tbi$/, "_${meta.chromosome}_dbsnp.vcf.gz.tbi") ]
    }

    // Subset clinvar variants
    SUBVAR ( ch_clinvar_vcf,
    SUBSETCAPTURE.out.target_bed.map{ map, target_bed -> [target_bed]})
    ch_versions = ch_versions.mix(SUBVAR.out.versions)
    out_clinvar_benign_vcf = SUBVAR.out.benign_vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/_benign.vcf$/, "_${meta.chromosome}_clinvar_benign.vcf") ]
    }
    out_pathogenic_benign_vcf = SUBVAR.out.pathogenic_vcf.map { meta, vcf ->
        [ meta, vcf.getName().replaceAll(/_pathogenic.vcf$/, "_${meta.chromosome}_clinvar_pathogenic.vcf") ]
    }
    emit:
    target_fa         = SAMTOOLS_FAIDX_SUBSET.out.fa        // channel:  [ val(meta), [ fa, fasta ] ]
    target_fai        = SAMTOOLS_FAIDX_INDEX.out.fai        // channel: [ val(meta), [ fai ] ]
    target_sizes      = SAMTOOLS_FAIDX_SIZES.out.sizes      // channel: [ val(meta), [ sizes ] ]
    //gzi             = SAMTOOLS_FAIDX.out.gzi              // channel: [ val(meta), [ gzi ] ]
    capture_bed_gz    = SUBSETCAPTURE.out.capture_bed_gz    // channel: [ val(meta), [ capture_bed_gz ] ]
    capture_bed_index = SUBSETCAPTURE.out.capture_bed_index // channel: [ val(meta), [ capture_bed_index ] ]
    target_bed        = SUBSETCAPTURE.out.target_bed        // channel: [ val(meta), [ target_bed ] ]
    target_bed_pad50  = SUBSETCAPTURE.out.target_bed_pad50  // channel: [ val(meta), [ target_bed_pad50 ] ]
    target_bed_pad500 = SUBSETCAPTURE.out.target_bed_pad500 // channel: [ val(meta), [ target_bed_pad500 ] ]
    //target_gatk_vcf      = GATK4_SELECTVARIANTS.out.vcf   // channel: [ val(meta), [ vcf ] ]
    //target_gatk_tbi      = GATK4_SELECTVARIANTS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    target_gnomad_vcf = out_gnomad_vcf // channel: [ val(meta), [ vcf ] ]
    target_gnomad_tbi = out_gnomad_tbi // channel: [ val(meta), [ tbi ] ]
    target_mills_vcf  = out_mills_vcf  // channel: [ val(meta), [ vcf ] ]
    target_mills_tbi  = out_mills_tbi  // channel: [ val(meta), [ tbi ] ]
    target_1000g_vcf  = out_1000G_vcf  // channel: [ val(meta), [ vcf ] ]
    target_1000g_tbi  = out_1000G_tbi  // channel: [ val(meta), [ tbi ] ]
    target_dbsnp_vcf  = out_dbsnp_vcf  // channel: [ val(meta), [ vcf ] ]
    target_dbsnp_tbi  = out_dbsnp_tbi  // channel: [ val(meta), [ tbi ] ]
    clinvar_benign_vcf        = out_clinvar_benign_vcf       // channel: [ val(meta), [ benign_vcf ] ]
    clinvar_pathogenic_vcf    = out_pathogenic_benign_vcf   // channel: [ val(meta), [ pathogenic_vcf ] ]
    clinvar_selected_vcf      = SUBVAR.out.selected_vcf     // channel: [ val(meta), [ selected_vcf ] ]
    versions = ch_versions                                  // channel: [ versions.yml ]
}
