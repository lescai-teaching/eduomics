process RNASEQVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-dose:4.0.0--5dccd2d8d274537d':
        'community.wave.seqera.io/library/bioconductor-dose:4.0.0--57136ca31aaf6317' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(deseq2_results_tsv), path(deseq2_tx2gene_tsv), path(deseq2_de_genes_txt), path(deseq2_pdf)
    tuple val(meta3), path(enrichment_rds), path(enrichment_png)

    output:
    tuple val(meta), path("rnaseq_validation/validated_reads/*.fasta.gz"), path("rnaseq_validation/deseq2_results.tsv"), path("rnaseq_validation/deseq2_tx2gene.tsv"), path("rnaseq_validation/deseq2_de_genes.txt"), path("rnaseq_validation/*.pdf"), path("rnaseq_validation/enrichment_results.rds"), path("rnaseq_validation/*.png"), path("rnaseq_validation/validation_result.txt")    , optional: true, emit: rnaseq_validated_results
    path "versions.yml"                                                                                                                                                                                                                                                                                                                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${baseDir}/bin/rnaseqvalidation.R ${enrichment_rds}

    if [ "\$(cat validation_result.txt)" == "GOOD SIMULATION" ]; then
        mkdir -p rnaseq_validation/validated_reads

        for read_file in ${reads}; do
            cp "\${read_file}" rnaseq_validation/validated_reads
        done

        cp "${deseq2_results_tsv}" rnaseq_validation/
        cp "${deseq2_tx2gene_tsv}" rnaseq_validation/
        cp "${deseq2_de_genes_txt}" rnaseq_validation/
        cp "${deseq2_pdf}" rnaseq_validation/

        cp "${enrichment_rds}" rnaseq_validation/
        if [ -n "${enrichment_png}" ]; then
            cp ${enrichment_png} rnaseq_validation/
        fi

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
    touch rnaseq_validation/deseq2_tx2gene.tsv
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-DOSE: \$(Rscript -e "cat(as.character(packageVersion('DOSE')))")
    END_VERSIONS
    """
}
