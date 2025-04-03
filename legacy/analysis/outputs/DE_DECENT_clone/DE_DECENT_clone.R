#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(scater))
suppressMessages(library(DECENT))

# Path to project folder
wd <- "/media/data/NGS/Projects/scrnaseq-clonal2"
source(file.path(wd, "scRNAseq/analysis/scripts/colors.R"))

output.path <- file.path(wd, "scRNAseq/analysis/outputs/DE_DECENT_clone")

iteration_input = "6_MCL_CNV"

sce.dir <- file.path(wd, "scRNAseq/analysis/outputs", iteration_input)
sce.MCL <- readRDS(file = file.path(sce.dir, "sce.rds"))

patients <- c("P069")
for (current_patient in patients) {
  sce <- sce.MCL[, sce.MCL$patient == current_patient]
  sce <- sce[, sce$timepoint == "DG"]
  sce <- sce[, !is.na(sce$CNV.clone)]
  
  # down-sampling for P069
  if (current_patient == "P069") {
    sce <- sce[, which((
      colnames(sce) %in% sample(colnames(sce), size = dim(sce)[2] * 0.5) &
        sce$CNV.clone == "DGex"
    ) | sce$CNV.clone == "DGs")]
  }
  
  sce$CNV.clone <- factor(sce$CNV.clone, levels = c("DGex", "DGs"))
  
  cat(paste0("\nProcessing ", current_patient, "\n"))
  
  mat <- counts(sce)
  
  cat(paste0("Initial counts matrix dims:\n"))
  cat(paste0("Genes: "), dim(mat)[1], "\n")
  cat(paste0("Cells: "), dim(mat)[2], "\n")
  
  # at least > 3% non-zero counts
  mat <- mat[(rowSums(mat == 0) / ncol(mat)) * 100 > 3, ]
  
  # > 5 non-zero counts
  mat <- mat[rowSums(mat != 0) > 5,]
  
  cat(paste0("Filtered counts matrix dims:\n"))
  cat(paste0("Genes: "), dim(mat)[1], "\n")
  cat(paste0("Cells: "), dim(mat)[2], "\n")
  
  output.path.sample <- file.path(output.path, current_patient)
  dir.create(output.path.sample, recursive = T, showWarnings = F)
  
  de.table <- decent(
    data.obs = as.matrix(mat),
    X = ~ sce$CNV.clone,
    use.spikes = F,
    normalize = "ML",
    CE.range = c(0.02, 0.1),
    # specify the range of the ranked random capture efficiency
    parallel = TRUE,
    n.cores = 64,
    dir = output.path.sample
  )
  
  output.filename <- paste0("DE_clone_", current_patient, ".csv")
  output.filepath <- file.path(output.path.sample, output.filename)
  
  write_csv(de.table, output.filepath)
  cat(paste0("\nWritten ", output.filename, "\n"))
}
