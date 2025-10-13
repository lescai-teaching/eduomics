#!/usr/bin/env Rscript

#### Load the libraries ####
library(tidyverse)
library(Biostrings)


#### Input selection ####

argv <- commandArgs(trailingOnly = TRUE)

id <- argv[1]
chromosome_of_interest <- argv[2]
fasta <- argv[3]
transcript_data <- argv[4]
log_file <- "subsetfastatx_parsing_log.txt"


#### Import the filtered transcript data ####
# This file containing protein-coding transcripts from the simulated chromosome
transcript_data <- readRDS(transcript_data)

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
    data = unlist(stringr::str_split(fasta_name, fixed("|")))
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
writeXStringSet(fasta_annotated, paste0(id, "_filtered_transcripts_novers_", chromosome_of_interest, ".fasta"), format = "fasta")
