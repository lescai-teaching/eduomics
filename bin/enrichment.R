#!/usr/bin/env Rscript

#### Load the library ####
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

#### Input selection ####
argv <- commandArgs(trailingOnly = TRUE)

deseq2_resdata <- argv[1]
tx2gene <- argv[2]

# Read the input files
resdata <- read_tsv(deseq2_resdata)
tx2gene <- read_tsv(tx2gene)

# Extract significant genes
sig_genes <- resdata$gene[which(resdata$padj < 0.05)]

#### Enrich GO analysis #####
ontologies <- c("BP", "MF", "CC")
enrichment_results <- list()

# Enrich GO loop
for (ont in ontologies) {
    ego <- tryCatch(
        enrichGO(
            gene = sig_genes,
            universe = unique(tx2gene$gene_id),
            OrgDb = org.Hs.eg.db,
            keyType = 'ENSEMBL',
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05
        ),
        error = function(e) NULL
    )

    enrichment_results[[ont]] <- ego

    if (!is.null(ego)) {
        # Dotplot for each ontology
        tryCatch(
            {
                file_name <- paste0("dotplot_", ont, ".png")
                png(file_name)
                print(dotplot(ego, showCategory = 10))
                dev.off()
            },
            error = function(e) NULL
        )

        # Cnetplot for each ontology
        tryCatch(
            {
                file_name <- paste0("cnetplot_", ont, ".png")
                png(file_name)
                print(cnetplot(ego, foldChange = resdata$log2FoldChange[which(resdata$padj < 0.05)]))
                dev.off()
            },
            error = function(e) NULL
        )
    }
}

# Save the enrichment results
saveRDS(enrichment_results, file = "enrichment_results.rds")
