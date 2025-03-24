process SUBSETFASTATX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-rtracklayer_r-tidyverse:f420c00b549a4380':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-rtracklayer_r-tidyverse:e36d3b6eec6fd274' }"

    input:
    val(meta)
    path(txfasta)
    path(gff3)

    output:
    path "*.fasta"                  , emit: fasta
    path "parsing_log.txt"          , emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    #### Load the libraries ####
    library(rtracklayer)
    library(tidyverse)
    library(Biostrings)


    #### Input selection ####
    chromosome_of_interest <- "${meta.chromosome}"
    log_file <- "parsing_log.txt"


    #### Parse GFF file ####
    gff_data <- import.gff3("${gff3}")


    #### Filter transcripts of interest ####

    # Filter for transcripts on the specified chromosome
    selected_transcripts <- gff_data[gff_data\$type == "transcript" & seqnames(gff_data) == chromosome_of_interest]
    length(unique(selected_transcripts\$transcript_id))

    # Convert to tibble and filter for protein-coding transcripts; select relevant columns.
    transcript_data <- as_tibble(as.data.frame(selected_transcripts)) %>%
    filter(gene_type == "protein_coding") %>%
    dplyr::select(transcript_id, gene_id, gene_name, gene_type)

    # Clean transcript and gene IDs by removing version numbers
    transcript_data <- transcript_data %>%
    mutate(
        transcript_id = gsub("\\\\.\\\\d+\$", "", transcript_id, perl = TRUE),
        gene_id = gsub("\\\\.\\\\d+\$", "", gene_id, perl = TRUE)
    )

    cat("1) Extracted", nrow(transcript_data),
        "protein-coding transcripts from", paste0(chromosome_of_interest),
        "\n", file = log_file, append = TRUE)


    #### MATCH FASTA FILE NEEDED FOR SIMULATIONS ####

    # Decompress and read the DNA sequences from the FASTA file
    fasta_all = readDNAStringSet("${txfasta}")

    cat("\n2) Number of transcripts in FASTA before filtering", length(fasta_all),
        "\n", file = log_file, append = TRUE)

    # Function to check if a given FASTA entry corresponds to an annotated transcript
    match_fasta_name <- function(fasta_name, parsed_gff_tibble){
    # Split the FASTA header by "|"
    data = unlist(stringr::str_split(fasta_name, "\\|"))
    # Extract the transcript ID, removing the version suffix (e.g., ".1")
    transcript_id = gsub("\\\\.\\\\d+\$","", data[[1]], perl=T)
    # Check if the transcript ID is present in the parsed GFF data
    if (transcript_id %in% parsed_gff_tibble\$transcript_id){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
    }

    # Filter FASTA sequences that match the transcript IDs in `transcript_data`
    fasta_annotated <- fasta_all[unlist(lapply(names(fasta_all), match_fasta_name, transcript_data)),]

    # Modify sequence names by removing version numbers and additional details
    names(fasta_annotated) <- gsub("\\\\.\\\\d+\\\\|.*\$","", names(fasta_annotated), perl=T)

    cat("\n3) Number of transcripts in FASTA after filtering: ", length(fasta_annotated),
        "\n", file = log_file, append = TRUE)


    #### Save the filtered FASTA sequences using the specified chromosome of interest ####
    writeXStringSet(fasta_annotated, paste0("gencode_transcripts_noversion_", chromosome_of_interest, ".fasta"), format = "fasta")

    # Generate versions.yml file
    writeLines(
    c(
        '"${task.process}":',
        paste0("    bioconductor-rtracklayer: ", packageVersion("rtracklayer")),
        paste0("    bioconductor-biostrings: ", packageVersion("Biostrings")),
        paste0("    r-tidyverse: ", packageVersion("tidyverse"))
    ),
    "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gencode_transcripts_noversion_${meta.chromosome}.fasta
    touch parsing_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
