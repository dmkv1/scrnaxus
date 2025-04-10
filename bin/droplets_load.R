#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Load the libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(DropletUtils)
})

# Get the sample data from NextFlow
sample_id <- args[1]
expected_cells <- as.numeric(args[2])
patient_id <- args[3]
timepoint <- args[4]
compartment <- args[5]
counts_dir <- args[6]


cat(
    "\nSample:", sample_id,
    "\nExpected cells:", expected_cells,
    "\nPatient:", patient_id,
    "\nTimepoint:", timepoint,
    "\nCompartment:", compartment,
    "\nPath:", counts_dir
)

sample_path <- file.path(counts_dir, "/Gene/raw")

cat("\n\nLoading counts...\n")

# Load the counts into an SCE
sce <- DropletUtils::read10xCounts(
    sample_path,
    sample.names = sample_id,
    col.names = TRUE
)

cat("\nInitial SCE:\n\n")
print(sce)

# Calculate library sizes for each droplet
libSizes <- colSums(assays(sce)$counts)

# Remove droplets that are completely empty
sce <- sce[, libSizes > 0]
cat("\n\nRemoved droplets with zero library sizes.")

# Assing sample column values
sce[["Patient"]] <- patient_id
sce[["Timepoint"]] <- timepoint
sce[["Compartment"]] <- compartment

colnames(sce) <- paste(colData(sce)[["Barcode"]],
    colData(sce)[["Sample"]],
    sep = "-"
)

cat("\n\nWriting SCE:\n\n")
print(sce)

# Write the output
saveRDS(sce, paste0(sample_id, "_unfiltered.sce"))
cat("\nDone!")