#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
path_sce_cells <- args[3]
path_sce_droplets <- args[4]

cat(
    "\nSample: ", sample_id,
    "\nCells SCE: ", path_sce_cells,
    "\nDroplets SCE: ", path_sce_cells
)

seed = 42
set.seed(seed)

# Load the libraries
suppressPackageStartupMessages({
  library(decontX)
})

# Load SCEs
cat("\n\nLoading cells SCE...\n")
sce_cells <- readRDS(path_sce_cells)
print(sce_cells)

cat("\n\nLoading droplets SCE...\n")
sce_droplets <- readRDS(path_sce_droplets)
print(sce_droplets)

# Conduct RNA decontamination
decont_sce <- decontX(sce_cells, background = sce_droplets)

cat("\n\nWriting decontaminated SCE:\n")
print(decont_sce)

# Write the output
saveRDS(decont_sce, paste0(sample_id, "_decont.sce"))

cat("\nDone!")