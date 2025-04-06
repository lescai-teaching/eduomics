include { BWA_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'


workflow VALIDATE_FASTQ_TO_VCF {

    take:
    reads_ch       // channel: [ val(meta), [ reads ] ]
    reference_ch   // channel: [ val(meta), reference_fasta ]
    known_sites_ch // channel: [ val(meta), known_sites_vcf ]

    main:
    // BWA Align
    BWA_MEM(reads_ch, reference_ch)

    // Mark Duplicates
    GATK4_MARKDUPLICATES(BWA_MEM.out.bam)

    // Base Quality Score Recalibration (BQSR)
    GATK4_BASERECALIBRATOR(GATK4_MARKDUPLICATES.out.bam, reference_ch, known_sites_ch)

    // Apply BQSR
    GATK4_APPLYBQSR(MARKDUPLICATES.out.bam, BQSR.out.recalibration_table, reference_ch)

    // HaplotypeCaller
    GATK4_HAPLOTYPECALLER(GATK4_APPLYBQSR.out.bam, reference_ch)

    // Import Genomics DB
    GATK4_GENOMICSDBIMPORT(GATK4_HAPLOTYPECALLER.out.gvcf)

    // Genotype GVCF
    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.genomics_db, reference_ch)

    emit:
    final_vcf = GATK4_GENOTYPEGVCFS.out.vcf
    final_vcf_index = GATK4_GENOTYPEGVCFS.out.vcf_index
}
