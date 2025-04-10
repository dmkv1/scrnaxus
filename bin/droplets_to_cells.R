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

###
# Counts to droplets
###

# Load the libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(DropletUtils)
    library(dplyr)
    library(ggplot2)
})

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

droplet_stats <- as.data.frame(droplet_stats_raw) %>% 
  filter(!is.na(FDR)) %>% 
  mutate(
    is_cell = case_when(
      FDR < FDR_thresh & Total > UMI_thresh ~ TRUE,
      TRUE ~ FALSE
    )
  )

droplets_table <- table(droplet_stats$is_cell)

subtitle_text <- paste0(
  "Initial droplets: ", total_droplets,
  "\nCell candidates: ", droplets_table["TRUE"] + droplets_table["FALSE"],
  "\nUMI threshold: ", UMI_thresh,
  "\nFDR threshold: ", FDR_thresh,
  "\nCells: ", droplets_table["TRUE"],
  "\nNot cells: ", droplets_table["FALSE"]
)

FDR_plot <- droplet_stats %>% 
  mutate(
    FDR = case_when(is.na(FDR) ~ 1, TRUE ~ FDR)
  ) %>% 
  ggplot(., aes(x = Total, y = -log10(FDR))) +
  geom_point(aes(color = is_cell)) +
  geom_hline(yintercept = -log10(FDR_thresh), color = "firebrick", linetype = "dashed") +
  geom_vline(xintercept = UMI_thresh, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(trans = "pseudo_log") +
  labs(
    title = sample_id,
    caption = subtitle_text,
    x = "Total UMI count for each barcode, log2 scale") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  paste0("FDRplot_", sample_id, ".png"),
  FDR_plot,
  device = "png", bg = "white", dpi = 300,
  width = 8.1, height = 5.6, units = "in")

# Filter the SCE by the cell barcodes
sce_cells <- sce_droplets[, rownames(filter(droplet_stats, is_cell))]

cat("\n\nWriting cells SCE:\n\n")
print(sce_cells)

# Write the output
saveRDS(sce_cells, paste0(sample_id, "_cells.sce"))

cat("\nDone!")