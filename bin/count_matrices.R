#!/usr/bin/env Rscript

#### Load the library ####
library(tidyverse)

argv <- commandArgs(trailingOnly = TRUE)

#### Input selection ####
fasta <- argv[1]
gene_list <- argv[2]

fasta_annotated <- readDNAStringSet(fasta)

# ~30x coverage ----> reads per transcript = transcriptlength/readlength * 30
# here all transcripts will have ~equal FPKM
readspertx = round(30 * width(fasta_annotated) / 100)

if(any(readspertx <= 0)){error('Reads_per_transcript contains zero or negative values')}

# base count matrix
num_timepoints = 6
countmat = matrix(readspertx, nrow=length(fasta_annotated), ncol=num_timepoints)

# create a list of possible fold changes to sample from
change <- c(4, 6, 8, 10, 14)
varchange <- c(1, 1.2, 1.5, 2, 0.9, 0.8, 0.6)

#### Functions ####
introduce_var <- function(countmatrix, varget){
  base = countmatrix
  for (i in 1:6) {
    countmatrix[,i] <- base[,i] * sample(varget, length(base[,i]), replace = T)
  }
  return(countmatrix)
}

introduce_fold_change <- function(countmatrix, index_tx, fold_changes){
  base = countmatrix
  for (i in 4:6) {
    countmatrix[index_tx,i] <- sample(fold_changes, length(index_tx), replace = T) * base[index_tx,i]
  }
  return(countmatrix)
}

createCountMatrices <- function(basematrix, variancechange, foldchanges, genesSelectionList, annotationData){
  countMatrices = list()
  for (selectionListName in names(genesSelectionList)){
    writeLines(paste0("creating count matrix for ", selectionListName))
    selectionList = genesSelectionList[[selectionListName]]
    simulationMatrix <- basematrix
    simulationMatrix <- introduce_var(simulationMatrix, variancechange)
    simulationMatrix <- introduce_fold_change(simulationMatrix, selectionList, foldchanges)
    simulationAnnotation <- annotationData[selectionList,]
    simulationAnnotation$expMeanCASE <- rowMeans(simulationMatrix[selectionList,4:6])
    simulationAnnotation$expMeanCTRL <- rowMeans(simulationMatrix[selectionList,1:3])
    simulationAnnotation$expLog2FC <- log2(simulationAnnotation$expMeanCASE) - log2(simulationAnnotation$expMeanCTRL)
    saveRDS(simulationMatrix, paste0("simulatedCountMatrix_",selectionListName,".RData"))
    saveRDS(simulationAnnotation, paste0("expected_",selectionListName,".RData"))
    write_tsv(simulationAnnotation, paste0("expected_",selectionListName,".tsv"))
    countMatrices[[selectionListName]] <- simulationMatrix
    write(selectionListName,file="datasets_simulated.txt",append=TRUE)
  }
  return(countMatrices)
}

all_simulated_counts <- createCountMatrices(countmat, varchange, change, gene_list, fasta_annotated)
