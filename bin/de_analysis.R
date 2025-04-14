#!/usr/bin/env Rscript

#### Load the libraries ####
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)


#### Input selection ####

argv <- commandArgs(trailingOnly = TRUE)

replica <- as.numeric(argv[1])
group <- as.numeric(argv[2])
tx2gene <- argv[3]
quant_dirs <- strsplit(argv[4], ",")[[1]]


#### Load the input files ####
group_labels <- c("control", "case")
dataset <- as_tibble(expand.grid(replica = 1:replica, group = 1:group) %>%
                       mutate(sample = paste0("sample_0", row_number()),
                              condition = factor(group, levels = 1:2, labels = group_labels)) %>%
                       dplyr::select(sample, condition))


tx2gene <- readRDS(tx2gene) %>%
  dplyr::select(transcript_id, gene_id)

write.table(tx2gene, file = "deseq2_tx2gene.tsv", sep = "\t", row.names = FALSE)


#### Load .quant files from Salmon  ####
files <- sapply(dataset$sample, function(sample) {
    quant_dir <- quant_dirs[which(grepl(sample, quant_dirs))]
    file.path(quant_dir, "quant.sf")
})
names(files) <- dataset$sample

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)


##### Differential expression analysis #####

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)

# prefilter min counts >10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# make sure base level is control
dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds)


#### Extract the results #####
res <- results(dds)
resOrdered <- res[order(res$pvalue),]

pdf("deseq2_ma_plot.pdf")
plotMA(res, ylim=c(-3,3))
dev.off()

pdf("deseq2_dispersion_plot.pdf")
plotDispEsts(dds)
dev.off()

pdf("deseq2_count_plot.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()


#### Save the results of the analysis ####
resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)

write_tsv(resdata, "deseq2_results.tsv")


#### Clustering ####
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("condition")])

pdf("deseq2_heatmap_plot.pdf")
pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)
dev.off()

pdf("deseq2_pca_plot.pdf")
plotPCA(ntd, intgroup=c("condition"))
dev.off()
