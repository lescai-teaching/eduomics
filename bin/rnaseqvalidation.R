#!/usr/bin/env Rscript

#### Load the library ####
library(DOSE)


#### Input selection ####
args <- commandArgs(trailingOnly = TRUE)

enrichment_rds <- args[1]


#### Load the input file ####
enrichment_results <- readRDS(enrichment_rds)


#### Loop to validate the simulation based on the GO enrichment ####
# At least one GO category (BP or MF or CC) with at least >= 3 enriched pathways
validate_enrichment <- function(enrichment_results) {
    categories <- intersect(names(enrichment_results), c("BP", "MF", "CC"))

    for (cat in categories) {
        if (dim(enrichment_results[[cat]])[1] >= 3) {
            write("GOOD SIMULATION", file = "validation_result.txt")
            return(invisible(NULL))
        }
    }

    write("SIMULATION NOT GOOD", file = "validation_result.txt")
    return(invisible(NULL))
}

validate_enrichment(enrichment_results)
