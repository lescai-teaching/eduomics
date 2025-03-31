#!/usr/bin/env Rscript

#### Load the library ####
library(tidyverse)
library(Biostrings)

argv <- commandArgs(trailingOnly = TRUE)


#### Input selection ####
meta_id <- argv[1]
coverage <- as.numeric(argv[2])
replica <- as.numeric(argv[3])
group <- as.numeric(argv[4])
fasta <- argv[5]
gff3 <- argv[6]
geneList <- argv[7]


#### Load the input file ####
fasta_annotated = readDNAStringSet(fasta)
annotation_data = readRDS(gff3)
geneList = readRDS(geneList)


# coverage ----> reads per transcript = transcriptlength/readlength * coverage
# here all transcripts will have ~equal FPKM
readspertx = round(coverage * width(fasta_annotated) / 100)

if(any(readspertx <= 0)){ stop("Reads_per_transcript contains zero or negative values") }


### BASE COUNT MATRIX
num_timepoints = replica * group
countmat = matrix(readspertx, nrow=length(fasta_annotated), ncol=num_timepoints)


### Create a list of possible fold changes
### cannot put fold changes too low in negative ==> risk to create negative read counts
### for transcripts with low reads per tx -> this of course gives an error
up_changes <- c(4, 6, 8, 10, 14)
down_changes <- c(0.20, 0.40, 0.60, 0.80, 0.90)
varchange <- c(1, 1.2, 1.5, 2, 0.9, 0.8, 0.6)


introduce_var <- function(countmatrix, varget, replica, group){
  base = countmatrix
  num_cols = replica * group
  for (i in 1:num_cols) {
    countmatrix[,i] <- base[,i] * sample(varget, length(base[,i]), replace = T)
  }
  return(countmatrix)
}

introduce_fold_change <- function(countmatrix, index_tx, fold_changes, replica, group){
  base = countmatrix
  num_cols = replica * group
  case_cols = ((replica + 1):num_cols)
  for (i in case_cols) {
    countmatrix[index_tx,i] <- sample(fold_changes, length(index_tx), replace = T) * base[index_tx,i]
  }
  return(countmatrix)
}

createCountMatrices <- function(meta_id, basematrix, variancechange, up_changes,
                                down_changes, genesSelectionList, annotationData,
                                replica, group) {
  countMatrices = list()
  total_samples <- replica * group
  half_samples <- total_samples / 2
  control_indices <- seq(1, half_samples)
  case_indices <- seq(half_samples + 1, total_samples)

  for (selectionListName in names(genesSelectionList)){
    writeLines(paste0("creating count matrix for ", selectionListName))
    selectionList = genesSelectionList[[selectionListName]]
    index <- which(annotationData$gene_name %in% selectionList)
    index_genes_upregulated <- index[1:round(length(index)/2, 0)]
    index_genes_downregulated <- index[(round(length(index)/2, 0) + 1):length(index)]
    simulationMatrix <- basematrix
    simulationMatrix <- introduce_var(simulationMatrix, variancechange, replica, group)
    simulationMatrix <- introduce_fold_change(simulationMatrix, index_genes_upregulated, up_changes, replica, group)
    simulationMatrix <- introduce_fold_change(simulationMatrix, index_genes_downregulated, down_changes, replica, group)
    rownames(simulationMatrix) <- rownames(annotationData)
    simulationAnnotation <- annotationData[match(selectionList, annotationData$gene_name), ]
    index_in_simulationMatrix <- match(simulationAnnotation$gene_name, annotationData$gene_name)
    simulationAnnotation$expMeanCASE <- rowMeans(simulationMatrix[index_in_simulationMatrix, case_indices])
    simulationAnnotation$expMeanCTRL <- rowMeans(simulationMatrix[index_in_simulationMatrix, control_indices])
    simulationAnnotation$expLog2FC <- log2(simulationAnnotation$expMeanCASE) - log2(simulationAnnotation$expMeanCTRL)
    cleanSelectionListName <- gsub(" ", "", selectionListName)
    saveRDS(simulationMatrix, paste0(meta_id, "_countMatrix_", cleanSelectionListName, ".rds"))
    saveRDS(simulationAnnotation, paste0(meta_id, "_expected_", cleanSelectionListName, ".rds"))
    write_tsv(simulationAnnotation, paste0(meta_id, "_expected_", cleanSelectionListName, ".tsv"))
    countMatrices[[selectionListName]] <- simulationMatrix
    write(selectionListName, file="datasets_simulated.txt", append=TRUE)
  }

  return(countMatrices)
}


all_simulated_counts <- createCountMatrices(meta_id, countmat, varchange, up_changes,
                                            down_changes, geneList, annotation_data, replica, group)
