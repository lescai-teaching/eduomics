process RNASEQVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-dose:4.0.0--5dccd2d8d274537d':
        'community.wave.seqera.io/library/bioconductor-dose:4.0.0--57136ca31aaf6317' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(deseq2_results)
    tuple val(meta3), path(enrichment_results)

    output:
    tuple val(meta), path("rnaseq_validation/validated_reads/*.fasta.gz"), path("rnaseq_validation/deseq2_results.tsv"), path("rnaseq_validation/deseq2_tx2gene.tsv"), path("rnaseq_validation/deseq2_de_genes.txt"), path("rnaseq_validation/*.pdf"), path("rnaseq_validation/enrichment_results.rds"), path("rnaseq_validation/*.png"), path("rnaseq_validation/validation_result.txt")    , optional: true, emit: rnaseq_validated_results
    path "versions.yml"                                                                                                                                                                                                                                                                                                                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat <<EOF > validate_enrichment.R
    #!/usr/bin/env Rscript

    library(DOSE)

    enrichment_results <- readRDS("${enrichment_results}")
    validate_enrichment <- function(enrichment_results) {
        categories <- intersect(names(enrichment_results), c("BP", "MF", "CC"))

        for (cat in categories) {
            if (dim(enrichment_results[[cat]])[1] > 3) {
                write("GOOD SIMULATION", file = "validation_result.txt")
                return(invisible(NULL))
            }
        }

        write("SIMULATION NOT GOOD", file = "validation_result.txt")
        return(invisible(NULL))
    }

    validate_enrichment(enrichment_results)
    EOF

    Rscript validate_enrichment.R

    if [ "\$(cat validation_result.txt)" == "GOOD SIMULATION" ]; then
        mkdir -p rnaseq_validation/validated_reads

        for read_file in ${reads}; do
            cp $\{read_file} rnaseq_validation/validated_reads
        done

        for f in ${deseq2_results}; do
            cp "\$f" rnaseq_validation/;
        done

        for f in ${enrichment_results}; do
            cp "\$f" rnaseq_validation/;
        done

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
    touch validated_reads/${prefix}_1.fasta.gz
    touch validated_reads/${prefix}_2.fasta.gz
    touch deseq2_results.tsv
    touch deseq2_tx2gene.tsv
    touch deseq2_de_genes.txt
    touch deseq2_ma_plot.pdf
    touch deseq2_dispersion_plot.pdf
    touch deseq2_count_plot.pdf
    touch deseq2_heatmap_plot.pdf
    touch deseq2_pca_plot.pdf
    touch enrichment_results.rds
    for ont in BP MF CC; do
        touch dotplot_\${ont}.png
        touch cnetplot_\${ont}.png
    done
    touch validated_de_genes.txt
    touch validation_result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-DOSE: \$(Rscript -e "cat(as.character(packageVersion('DOSE')))")
    END_VERSIONS
    """
}
