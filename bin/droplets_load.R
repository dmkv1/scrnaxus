#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Load the libraries
library(DropletUtils)
library(SingleCellExperiment)

# Get the sample data from NextFlow
sample_id <- args[1]
expected_cells <- as.numeric(args[2])
patient_id <- args[3]
timepoint <- args[4]
compartment <- args[5]
counts_dir <- args[6]


cat(
    "\nSample:", sample_id,
    "\nExpected cell:", expected_cells,
    "\nPatient:", patient_id,
    "\nTimepoints:", timepoint,
    "\nCompartment:", compartment,
    "\nPath:", counts_dir
)


sample_path <- file.path(counts_dir, "/Gene/raw")

# Load the counts into an SCE
sce <- DropletUtils::read10xCounts(
    sample_path,
    sample.names = sample_id,
    col.names = TRUE
)

# Calculate library sizes for each droplet
libSizes <- colSums(assays(sce)$counts)

# Remove droplets that are completely empty
sce <- sce[, libSizes > 0]

# Assing sample column values
sce[["Patient"]] <- patient_id
sce[["Timepoint"]] <- timepoint
sce[["Compartment"]] <- compartment

colnames(sce) <- paste(colData(sce)[["Barcode"]],
    colData(sce)[["Sample"]],
    sep = "-"
)

# Write the output
saveRDS(sce, "${sample_id}_unfiltered.sce")
