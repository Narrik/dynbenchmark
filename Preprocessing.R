# Load libraries
library(dyno)
library(dynbenchmark)
library(dyneval)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(TENxPBMCData)

#df <- read_csv("pbmcd8k.csv")
sce <- TENxPBMCData(dataset = "pbmc8k")
colnames(sce) <- colnames(df)[-1]

# Filter cell outliers
qcstats <- perCellQCMetrics(sce)
filtered <- quickPerCellQC(qcstats)
sce <- sce[, !filtered$discard]

# Normalize and transform into log2 the filtered counts
sce <- logNormCounts(sce)
# Keep prop=0.1 Highly Variable Genes
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)
sce_counts <- counts(sce)[hvg,]
sce_logcounts <- logcounts(sce)[hvg,]

dim(sce)
dim(sce_counts)
dim(sce_logcounts)

# Transform from SCE to matrix
countsMatrix <- sce_counts %>% as.matrix() %>% t
logCountsMatrix <- sce_logcounts %>% as.matrix() %>% t

sce_counts <- counts(sce)[hvg,]
sce_logcounts <- logcounts(sce)[hvg,]

# Transform from SCE to matrix
countsMatrix <- sce_counts %>% as.matrix() %>% t
logCountsMatrix <- sce_logcounts %>% as.matrix() %>% t

dataset<- wrap_expression(
  counts = countsMatrix,
  expression = logCountsMatrix
)
