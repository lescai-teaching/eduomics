#!/usr/bin/env R

#### Load the libraries ####
library(rtracklayer)
library(tidyverse)
library(Biostrings)

argv <- commandArgs(trailingOnly = TRUE)

#### Input selection ####
chromosome_of_interest <- argv[1]
fasta <- argv[2]
gff3 <- argv[3]
log_file <- "parsing_log.txt"


#### Parse GFF file ####
gff_data <- import.gff3(gff3)


#### Filter transcripts of interest ####

# Filter for transcripts on the specified chromosome
selected_transcripts <- gff_data[gff_data$type == "transcript" & seqnames(gff_data) == chromosome_of_interest]
length(unique(selected_transcripts$transcript_id))

# Convert to tibble and filter for protein-coding transcripts; select relevant columns.
transcript_data <- as_tibble(as.data.frame(selected_transcripts)) %>%
filter(gene_type == "protein_coding") %>%
dplyr::select(transcript_id, gene_id, gene_name, gene_type)

# Clean transcript and gene IDs by removing version numbers
transcript_data <- transcript_data %>%
mutate(
    transcript_id = gsub("\\.\\d+$", "", transcript_id, perl = TRUE),
    gene_id = gsub("\\.\\d+$", "", gene_id, perl = TRUE)
)

cat("1) Extracted", nrow(transcript_data),
    "protein-coding transcripts from", paste0(chromosome_of_interest),
    "\n", file = log_file, append = TRUE)


#### MATCH FASTA FILE NEEDED FOR SIMULATIONS ####

# Decompress and read the DNA sequences from the FASTA file
fasta_all = readDNAStringSet(fasta)

cat("\n2) Number of transcripts in FASTA before filtering", length(fasta_all),
    "\n", file = log_file, append = TRUE)

# Function to check if a given FASTA entry corresponds to an annotated transcript
match_fasta_name <- function(fasta_name, parsed_gff_tibble){
# Split the FASTA header by "|"
data = unlist(stringr::str_split(fasta_name, "\|"))
# Extract the transcript ID, removing the version suffix (e.g., ".1")
transcript_id = gsub("\\.\\d+$","", data[[1]], perl=T)
# Check if the transcript ID is present in the parsed GFF data
if (transcript_id %in% parsed_gff_tibble$transcript_id){
    return(TRUE)
}
else{
    return(FALSE)
}
}

# Filter FASTA sequences that match the transcript IDs in `transcript_data`
fasta_annotated <- fasta_all[unlist(lapply(names(fasta_all), match_fasta_name, transcript_data)),]

# Modify sequence names by removing version numbers and additional details
names(fasta_annotated) <- gsub("\\.\\d+\\|.*$","", names(fasta_annotated), perl=T)

cat("\n3) Number of transcripts in FASTA after filtering: ", length(fasta_annotated),
    "\n", file = log_file, append = TRUE)


#### Save the filtered FASTA sequences using the specified chromosome of interest ####
writeXStringSet(fasta_annotated, paste0("gencode_transcripts_noversion_", chromosome_of_interest, ".fasta"), format = "fasta")
