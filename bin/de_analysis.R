#!/usr/bin/env Rscript

#### Load the library ####
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)

argv <- commandArgs(trailingOnly = TRUE)


#### Input selection ####
replica <- as.numeric(argv[1])
group <- as.numeric(argv[2])
gff3 <- argv[3]
counts <- argv[4]


#### Creation of input files ####

# Dataset with samples and associated conditions
group_labels <- c("control", "case")
dataset <- as_tibble(expand.grid(replica = 1:replica, group = 1:group) %>%
  mutate(sample = paste0("sample_", row_number()),
         condition = group_labels[group]) %>%
  dplyr::select(sample, condition))

# tx2gene
tx2gene <- readRDS(gff3) %>%
  dplyr::select(transcript_id, gene_id)


#### Load .quant files from Salmon  ####

files <- list.files(counts, pattern = "quant.sf", recursive = T, full.names = T)
names <- basename(dirname(files))
names(files) <- names

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)


#### Prefilter min counts >10 ####

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### make sure base level is control
dds$condition <- relevel(dds$condition, ref = "control")


##### Differential expression analysis #####

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
