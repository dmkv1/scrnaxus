#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
path_sce <- args[3]

# Load the libraries
suppressPackageStartupMessages({
  library(scDblFinder)
})

seed = 42
set.seed(seed)

# Load SCEs
sce <- readRDS(path_sce)

# Annotate the doublets
sce <- scDblFinder(
  sce,
  clusters = FALSE
)

# Extract metrics
n_singlets <- sum(sce$scDblFinder.class == "singlet")
n_doublets <- sum(sce$scDblFinder.class == "doublet")
doublet_rate <- n_doublets / (n_singlets + n_doublets)

metrics <- list(
  sample_id = sample_id,
  n_singlets = n_singlets,
  n_doublets = n_doublets,
  doublet_rate = doublet_rate
)
writeLines(jsonlite::toJSON(metrics, pretty=TRUE), paste0(sample_id, "_doublet_metrics.json"))

# Write the SCE
cat("\n\nWriting doublet-annotated SCE:\n\n")
print(sce)
saveRDS(sce, paste0(sample_id, "_singlets.sce"))
cat("\nDone!\n")

# Print session info and collect garbage before exiting
print(sessionInfo())
gc()