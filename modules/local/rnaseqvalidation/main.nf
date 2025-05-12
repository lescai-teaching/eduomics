process RNASEQVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-dose:4.0.0--5dccd2d8d274537d':
        'community.wave.seqera.io/library/bioconductor-dose:4.0.0--57136ca31aaf6317' }"

    input:
    tuple val(meta),  path(enrichment_results)
    tuple val(meta2), path(reads)
    tuple val(meta3), path(de_genes)

    output:
    tuple val(meta), path("validated_reads/*.fasta.gz")    , optional: true, emit: valid_reads
    tuple val(meta), path("validated_de_genes.txt")        , optional: true, emit: valid_de_genes
    tuple val(meta), path("validation_result.txt")         , emit: validation_result
    path "versions.yml"                                    , emit: versions

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
        mkdir -p validated_reads
        for read_file in ${reads}; do
            ln -s "\$(readlink -f "\$read_file")" validated_reads/
        done
        cp ${de_genes} validated_de_genes.txt
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
    mkdir -p validated_reads
    touch validated_reads/${prefix}_1.fasta.gz
    touch validated_reads/${prefix}_2.fasta.gz
    touch validated_de_genes.txt
    touch validation_result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-DOSE: \$(Rscript -e "cat(as.character(packageVersion('DOSE')))")
    END_VERSIONS
    """
}
