#### Load the library ####
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)


#### Input selection ####

argv <- commandArgs(trailingOnly = TRUE)

deseq2_resdata <- argv[1]
deseq2_tx2gene <- argv[2]
outdir <- argv[3]

dir.create(outdir, showWarnings = FALSE)


# Read the input files
resdata = readRDS(deseq2_resdata)
tx2gene = read_tsv(deseq2_tx2gene) %>%
  dplyr::select(transcript_id, gene_id)

# Extract significant genes
sig_genes <- resdata$gene[which(resdata$padj < 0.05)]


#### Enrich GO analysis #####

ontologies <- c("BP", "MF", "CC")
enrichment_results <- list()

#### Enrich GO loop ####
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

  if (!is.null(ego) && any(ego@result$p.adjust < 0.05, na.rm = TRUE)) {
    ego@result <- ego@result[ego@result$p.adjust < 0.05, , drop = FALSE]
    enrichment_results[[ont]] <- ego

    # Dotplot for each ontology
    tryCatch(
      {
        file_name <- file.path(outdir, paste0("dotplot_", ont, ".png"))
        png(file_name)
        print(dotplot(ego, showCategory = 10))
        dev.off()
      },
      error = function(e) NULL
    )

    # Cnetplot for each ontology
    tryCatch(
      {
        file_name <- file.path(outdir, paste0("cnetplot_", ont, ".png"))
        png(file_name)
        print(cnetplot(ego, foldChange = resdata$log2FoldChange[which(resdata$padj < 0.05)]))
        dev.off()
      },
      error = function(e) NULL
    )
  }
}

# Save the enrichment results
saveRDS(enrichment_results, file = file.path(outdir, "enrichment_results.rds"))
