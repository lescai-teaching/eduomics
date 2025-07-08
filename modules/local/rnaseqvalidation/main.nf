process RNASEQVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-dose:4.0.0--5dccd2d8d274537d':
        'community.wave.seqera.io/library/bioconductor-dose:4.0.0--57136ca31aaf6317' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(deseq2_results_tsv), path(deseq2_de_genes_txt), path(deseq2_pdf)
    tuple val(meta3), path(enrichment_rds), path(enrichment_png)
    tuple val(meta4), path(deseq2_tx2gene)

    output:
    tuple val(meta), path("rnaseq_validation")    , optional: true, emit: rnaseq_validated_results
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rnaseqvalidation.R ${enrichment_rds}

    if [ "\$(cat validation_result.txt)" == "GOOD SIMULATION" ]; then
        mkdir -p rnaseq_validation/validated_reads

        for read_file in ${reads}; do
            cp "\${read_file}" rnaseq_validation/validated_reads
        done

        cp "${deseq2_results_tsv}" rnaseq_validation/
        cp "${deseq2_de_genes_txt}" rnaseq_validation/

        for pdf_file in ${deseq2_pdf}; do
            cp "\${pdf_file}" rnaseq_validation/
        done

        cp "${enrichment_rds}" rnaseq_validation/

        for png_file in ${enrichment_png}; do
            cp "\${png_file}" rnaseq_validation/
        done

        cp "${deseq2_tx2gene}" rnaseq_validation/

        cp validation_result.txt rnaseq_validation/

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-DOSE: \$(Rscript -e "cat(as.character(packageVersion('DOSE')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p rnaseq_validation/validated_reads
    touch rnaseq_validation/validated_reads/${prefix}_1.fasta.gz
    touch rnaseq_validation/validated_reads/${prefix}_2.fasta.gz
    touch rnaseq_validation/deseq2_results.tsv
    touch rnaseq_validation/deseq2_de_genes.txt
    touch rnaseq_validation/deseq2_ma_plot.pdf
    touch rnaseq_validation/deseq2_dispersion_plot.pdf
    touch rnaseq_validation/deseq2_count_plot.pdf
    touch rnaseq_validation/deseq2_heatmap_plot.pdf
    touch rnaseq_validation/deseq2_pca_plot.pdf
    touch rnaseq_validation/enrichment_results.rds
    for ont in BP MF CC; do
        touch rnaseq_validation/dotplot_\${ont}.png
        touch rnaseq_validation/cnetplot_\${ont}.png
    done
    touch rnaseq_validation/validation_result.txt
    touch rnaseq_validation/deseq2_tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-DOSE: \$(Rscript -e "cat(as.character(packageVersion('DOSE')))")
    END_VERSIONS
    """
}
