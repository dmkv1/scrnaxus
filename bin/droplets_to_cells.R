#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]

# Get the sample data from NextFlow
sample_id <- args[2]
expected_cells <- as.numeric(args[3])
patient_id <- args[4]
timepoint <- args[5]
compartment <- args[6]
counts_dir <- args[7]

cat(
  "\nSample:", sample_id,
  "\nExpected cells:", expected_cells,
  "\nPatient:", patient_id,
  "\nTimepoint:", timepoint,
  "\nCompartment:", compartment,
  "\nPath:", counts_dir
)

# Load the libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
})

###
# Counts to droplets
###

cat("\n\nLoading counts...\n")

# Load the counts into an SCE
sample_path <- file.path(counts_dir, "/Gene/raw")
sce_unfiltered <- DropletUtils::read10xCounts(
  sample_path,
  sample.names = sample_id,
  col.names = TRUE
)

# Assing sample column values
sce_unfiltered[["Patient"]] <- patient_id
sce_unfiltered[["Timepoint"]] <- timepoint
sce_unfiltered[["Compartment"]] <- compartment

# Avoid barcode collision
colnames(sce_unfiltered) <- paste(colData(sce_unfiltered)[["Barcode"]],
                                  colData(sce_unfiltered)[["Sample"]],
                                  sep = "-"
)

cat("\nInitial unfiltered SCE:\n\n")
print(sce_unfiltered)

# Calculate library sizes for each droplet
libSizes <- colSums(assays(sce_unfiltered)$counts)

# Remove droplets that are completely empty
sce_droplets <- sce_unfiltered[, libSizes > 0]
cat("\n\nRemoved droplets with zero library sizes.")

cat("\n\nWriting droplets SCE:\n\n")
print(sce_droplets)

# Write the output
saveRDS(sce_droplets, paste0(sample_id, "_droplets.sce"))

### 
# Droplets to cells
###
total_droplets <- dim(sce_droplets)[2]

# Calculate empty droplet statistics
set.seed(seed)
droplet_stats_raw <- DropletUtils::emptyDropsCellRanger(
  counts(sce_droplets),
  n.expected.cells = expected_cells,
  niters = 50000)

# Filter the droplets using set thresholds
UMI_thresh = 500
FDR_thresh = 0.1

droplet_stats <- as.data.frame(droplet_stats_raw)
droplet_stats <- droplet_stats[!is.na(droplet_stats$FDR), ]

droplet_stats$is_cell <- FALSE
droplet_stats$is_cell[droplet_stats$FDR < FDR_thresh &
                        droplet_stats$Total > UMI_thresh] <- TRUE

# Gather and write metrics
n_cells <- sum(droplet_stats$is_cell == TRUE)
n_cell_dropouts <- sum(droplet_stats$is_cell == FALSE)
n_cell_candidates <- n_cells + n_cell_dropouts

metrics <- list(
  sample_id = sample_id,
  n_droplets = total_droplets,
  n_cell_candidates = n_cell_candidates,
  UMI_threshold = UMI_thresh,
  FDR_threshold = FDR_thresh,
  n_cells = n_cells,
  cells_ratio = n_cells / n_cell_candidates
)
writeLines(jsonlite::toJSON(metrics, pretty=TRUE), paste0(sample_id, "_droplet_metrics.json"))

# Filter the SCE by the cell barcodes
sce_cells <- sce_droplets[, rownames(droplet_stats[droplet_stats$is_cell, ])]

# Write the SCE
cat("\n\nWriting cells SCE:\n\n")
print(sce_cells)
saveRDS(sce_cells, paste0(sample_id, "_cells.sce"))
cat("\nDone!\n")

# Print session info and collect garbage before exiting
print(sessionInfo())
gc()