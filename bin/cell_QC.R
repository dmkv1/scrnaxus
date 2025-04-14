#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
path_sce <- args[3]

# Load the libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
})

seed = 42
set.seed(seed)

# Load SCEs
sce <- readRDS(path_sce)

# Work


# Extract metrics
metrics <- list(
  sample_id = sample_id  
)
writeLines(jsonlite::toJSON(metrics, pretty=TRUE), paste0(sample_id, "_QC_metrics.json"))

# Write the SCE
cat("\n\nWriting clean SCE:\n\n")
print(sce)
saveRDS(sce, paste0(sample_id, "_clean.sce"))
cat("\nDone!\n")

# Print session info and collect garbage before exiting
print(sessionInfo())
gc()