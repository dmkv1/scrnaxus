#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
path_sce_cells <- args[3]
path_sce_droplets <- args[4]

cat(
    "\nSample:", sample_id,
    "\nSCE cells path:", path_sce_cells
)

# Load the libraries
suppressPackageStartupMessages({
    library(dplyr)
})

cat("\n\nLoading cells SCE...\n")

# Load the counts into an SCE
sce_cells <- readRDS(path_sce_cells)
print(sce_cells)

# Mock script
sce_cells <- sce_cells[, 1:10]

cat("\n\nWriting cells SCE:\n\n")
print(sce_cells)

# Write the output
saveRDS(sce_cells, paste0(sample_id, "_decont.sce"))

cat("\nDone!")