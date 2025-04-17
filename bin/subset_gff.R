#!/usr/bin/env R

#### Load the libraries ####
library(rtracklayer)
library(igraph)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)
library(clusterProfiler)
library(Matrix)


#### Input selection ####

argv <- commandArgs(trailingOnly = TRUE)

chromosome_of_interest <- argv[1]
simthreshold <- as.numeric(argv[2])
gff3 <- argv[3]
log_file <- "subsetgff_parsing_log.txt"


#### Import the GFF file ####
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

cat("1) Annotation parsing\nExtracted", nrow(transcript_data),
    "protein-coding transcripts from", paste0(chromosome_of_interest),
    "\n", file = log_file, append = TRUE)


#### RETRIEVE GO ANNOTATIONS USING org.Hs.eg.db ####

# Use gene names as keys (SYMBOL) to retrieve GO annotations.
unique_genes <- unique(transcript_data$gene_name)
go_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = unique_genes,
                                        columns = c("GOALL", "ONTOLOGYALL"),
                                        keytype = "SYMBOL")
# Filter to exclude lacking ontologies
go_annotations <- subset(go_annotations, !is.na(ONTOLOGYALL))
go_annotations <- unique(go_annotations[, c("SYMBOL", "GOALL", "ONTOLOGYALL")])

cat("\n2) Number of Gene Ontologies\nRetrieved", nrow(go_annotations),
    "unique GO annotations (ALL) for", length(unique_genes), "genes",
    "\n", file = log_file, append = TRUE)


#### CREATE TRANSCRIPT-TO-GO MAPPING ####

# Merge transcript data with GO annotations using gene_name
transcript_go <- transcript_data %>%
inner_join(go_annotations, by = c("gene_name" = "SYMBOL"))

# Build a list mapping each transcript to its set of GO terms
transcript_go_list <- split(transcript_go$GOALL, transcript_go$transcript_id)
transcript_go_list <- transcript_go_list[sapply(transcript_go_list, length) > 0]


#### BUILD THE NETWORK ####

# Get unique transcript IDs and GO terms
transcript_ids <- names(transcript_go_list)
all_go_terms <- unique(unlist(transcript_go_list))

# Create a binary incidence matrix: rows=transcripts, columns=GO terms
binary_mat <- matrix(0, nrow = length(transcript_ids), ncol = length(all_go_terms),
                    dimnames = list(transcript_ids, all_go_terms))

# Fill the matrix: set entry = 1 if transcript is annotated with the GO term
for (i in seq_along(transcript_ids)) {
binary_mat[transcript_ids[i], transcript_go_list[[transcript_ids[i]]]] <- 1
}
# convert to a sparse matrix for memory efficiency (there will be lots of places without matches):
binary_mat <- Matrix(binary_mat, sparse = TRUE)

# Compute the intersection matrix: each (i,j) entry is the number of shared GO terms
# tcrossprod -> is the matrix cross product %*% of the transposed: this is used to quickly
# get the intersection between genes and GO terms they have
intersection <- tcrossprod(binary_mat)

# Compute the number of GO terms per transcript (row sums)
go_counts <- rowSums(binary_mat)

# Compute the union matrix using: union = count[i] + count[j] - intersection[i,j]
# 'outer' computes the pairwise sums, then subtract the intersection.
union_matrix <- outer(go_counts, go_counts, "+") - intersection

# Compute the Jaccard similarity matrix
jaccard_matrix <- intersection / union_matrix

# Set diagonal to NA (or 0) to ignore self-comparison
diag(jaccard_matrix) <- NA

# Now, we can extract transcript pairs that meet the similarity threshold
# We'll extract the upper triangle only (to avoid duplicates)
# this is the parameter we set already at the beginning
similarity_threshold <-  simthreshold
idx <- which(jaccard_matrix >= similarity_threshold & upper.tri(jaccard_matrix), arr.ind = TRUE)

# Create an edge list from these indices
edges <- data.frame(
transcript1 = rownames(jaccard_matrix)[idx[, 1]],
transcript2 = colnames(jaccard_matrix)[idx[, 2]],
weight = jaccard_matrix[idx],
stringsAsFactors = FALSE
)

# Build an undirected, weighted graph of transcripts.
g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = transcript_ids)

cat("\n3) Number of edges\nConstructed an edge list with",
    nrow(edges), "edges based on vectorised GO similarity",
    "\n", file = log_file, append = TRUE)


#### GRAPH-BASED CLUSTERING (LOUVAIN ALGORITHM) ####

communities <- cluster_louvain(g, weights = E(g)$weight)

# now we can add the memberships calculated to the original tibble
# where we have one transcript per row

# remember 'communities' is the result from cluster_louvain and membership_vector is:
membership_vector <- membership(communities)

# Create a tibble from the membership vector
membership_tb <- tibble(
transcript_id = names(membership_vector),
membership = as.integer(membership_vector)
)

# Define the minimum and maximum number of DE genes
chr_genes <- length(unique(transcript_data$gene_name))
chr_genes
percent_de_genes <- 0.1

# now left join with the original transcript_data tibble and filter out genes with NA in the membership column
transcript_data <- transcript_data %>%
left_join(membership_tb, by = "transcript_id")


#### PHENOTYPE DATA FROM MONARCH ####

# This step has been silenced because it is very time-consuming and for the first release is not necessary

# get_phenotypes_for_gene <- function(gene_id) {
#   # Ensure gene_id is in CURIE format. If missing ":", assume it's an ENSG id and prepend "ENSEMBL:"
#   parse_url <- "&limit=20&offset=0"
#   # Construct the query URL using the new API endpoint.
#   base_url <- "https://api-v3.monarchinitiative.org/v3/api/search?q=ENSEMBL%3A"
#   query_url <- paste0(base_url, gene_id, parse_url)
#
#   # Query the API.
#   res <- GET(query_url)
#   if (res$status_code != 200) {
#     warning("Query failed for gene: ", gene_id)
#     return(NA_character_)
#   }
#
#   # Parse the JSON response
#   res_text <- content(res, as = "text", encoding = "UTF-8")
#   json_data <- fromJSON(res_text, flatten = TRUE)
#
#   # Extract the 'has_phenotype_label' from all returned docs
#   phenos_description <- paste0(json_data$items$has_phenotype_label[[1]], collapse = ";")
#   phenos_hpo <- paste0(json_data$items$has_phenotype[[1]], collapse = ";")
#
#   if (length(phenos_description) == 0) {
#     return(NA_character_)
#   }
#
#   # Return unique phenotype labels
#   return(c(phenos_hpo,phenos_description))
# }
#
# transcript_data_pheno <- transcript_data %>%
#   rowwise() %>%
#   mutate(phenotypes = list(get_phenotypes_for_gene(gene_id))) %>%
#   ungroup() %>%
#   mutate(
#     phenotype_hpo = map_chr(phenotypes, ~ .x[1]),
#     phenotype_description = map_chr(phenotypes, ~ .x[2])
#   ) %>%
#   dplyr::select(-phenotypes)

# Tidy the data before enrichment
transcript_data <- transcript_data %>%
filter(!is.na(membership)) %>%
group_by(membership) %>%
mutate(
    gene_number = n_distinct(gene_name),
) %>%
ungroup() %>%
arrange(gene_number) %>%
relocate(gene_number, .before = membership) %>%
filter(gene_number > 3 & gene_number < (chr_genes * percent_de_genes))


#### Enrichment with EnrichGO ####

# Construct the gene list
gene_lists <- transcript_data %>%
group_by(membership) %>%
summarise(
    unique_genes = list(unique(gene_name)),
    .groups = "drop"
) %>%
mutate(
    gene_list = paste0("list", row_number())
)

cat("\n4) Number of Gene Lists\nRetrieved", nrow(gene_lists),
    "gene lists from", paste0(chromosome_of_interest),
    "\n", file = log_file, append = TRUE)


# EnrichGO function
executeGO <- function(sig_genes, universe, ontology) {
# Perform GO enrichment analysis
res <- enrichGO(
    gene = sig_genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)
return(res)
}

# Initialise an empty list to store significant results (both positive and negative)
significant_results <- list()

write("\n5) GO enrichment", file = log_file, append = TRUE)
for (i in 1:nrow(gene_lists)) {
sig_genes <- gene_lists$unique_genes[[i]]  # Extract the unique genes for the group
universe <- unique_genes                   # Define the universe

# Perform enrichment for BP, MF, and CC
ego_bp <- executeGO(sig_genes, universe, "BP")
ego_mf <- executeGO(sig_genes, universe, "MF")
ego_cc <- executeGO(sig_genes, universe, "CC")

# Check if the results are NULL, and if so, assign "No enrichment"
enrichment_results <- list(
    BP = if (!is.null(ego_bp) && dim(ego_bp)[1] >= 3) ego_bp else "No enrichment",
    MF = if (!is.null(ego_mf) && dim(ego_mf)[1] >= 3) ego_mf else "No enrichment",
    CC = if (!is.null(ego_cc) && dim(ego_cc)[1] >= 3) ego_cc else "No enrichment"
)

# Add the enrichment results to significant_results with the gene list ID
significant_results[[paste("list", i)]] <- enrichment_results

# Determine and print final message
enriched_categories <- names(enrichment_results)[enrichment_results != "No enrichment"]

if (length(enriched_categories) > 0) {
    msg <- paste("Category", paste(enriched_categories, collapse = ", "), "enriched for gene list", i)
} else {
    msg <- paste("No enrichment for gene list", i)
}

# Save the message in the log and print it
write(msg, file = log_file, append = TRUE)
}


##### Extract final gene lists ####

# Initialise an empty list to store valid gene lists
valid_gene_lists <- list()

write("\n6) Valid gene sets", file = log_file, append = TRUE)
for (i in 1:length(significant_results)) {

# Create a variable to store enriched categories for the current gene list
enriched_categories <- c()

# Check for significant results in 'BP', 'MF', and 'CC' categories
# If the result is not a character and has at least one row, it is considered significant
if (!is.character(significant_results[[i]][["BP"]]) && nrow(significant_results[[i]][["BP"]]) > 0) {
    enriched_categories <- c(enriched_categories, "BP")
}

if (!is.character(significant_results[[i]][["MF"]]) && nrow(significant_results[[i]][["MF"]]) > 0) {
    enriched_categories <- c(enriched_categories, "MF")
}

if (!is.character(significant_results[[i]][["CC"]]) && nrow(significant_results[[i]][["CC"]]) > 0) {
    enriched_categories <- c(enriched_categories, "CC")
}

# If there are any enriched categories, store the gene list in valid_gene_lists
if (length(enriched_categories) > 0) {
    valid_gene_lists[[paste0("list", i)]] <- gene_lists$unique_genes[[i]]

    # Create a message for enrichment and write it to the log file
    msg <- paste("Gene set", i, "has enriched categories:", paste(enriched_categories, collapse = ", "))
    write(msg, file = log_file, append = TRUE)
} else {
    # Create a message for lists with no significant enrichment and write it to the log file
    msg <- paste("No enrichment for gene set", i)
    write(msg, file = log_file, append = TRUE)
}
}

# Tibble for valid gene lists
valid_gene_lists_df <- tibble(
gene_list = names(valid_gene_lists),
genes = sapply(valid_gene_lists, function(x) paste(x, collapse = ","))
)


#### Save the resulting key files ####
write.table(valid_gene_lists_df, file = "list_gene_association.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(valid_gene_lists, "valid_gene_lists.rds")
saveRDS(transcript_data, "filtered_gff3.rds")
