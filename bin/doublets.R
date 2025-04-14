#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
path_sce <- args[3]

cat(
    "\nSample: ", sample_id,
    "\nInput SCE: ", path_sce
)

seed = 42
set.seed(seed)

# Load the libraries
suppressPackageStartupMessages({
  library(scDblFinder)
})

# Load SCEs
cat("\n\nLoading cells SCE...\n")
sce <- readRDS(path_sce)
print(sce)

clusters_igraph <- scran::quickCluster(
  sce,
  method = "igraph"
  )

doublet_scores <- scDblFinder(
  sce,
  clusters = factor(clusters_igraph),
  returnType = "table"
)

doublet_calls <-
  doubletThresholding(
    doublet_scores,
    method = "griffiths",
    returnType = "call"
  )

# Assign doublet counts
sce$scDblFinder_calls <- doublet_calls

cat("\n\nWriting doublet-labbeled SCE:\n")
print(sce)

# Write the output
saveRDS(sce, paste0(sample_id, "_singlets.sce"))

cat("\nDone!")