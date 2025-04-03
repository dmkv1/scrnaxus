#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(openxlsx)))
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(scater)))
suppressMessages(suppressWarnings(library(infercnv)))
options(scipen = 100)

# path to output folder
inferCNV.output.path <- args[[1]]
print(paste("Saving results to", inferCNV.output.path))

# Path to project folder
wd <- "/media/data/NGS/Projects/Clonal2"
iteration = "4_all_cells_labbeled"
sce.dir <- file.path(wd, "scRNAseq/analysis/outputs", iteration)

patients <- c("P009", "P022", "P027", "P087", "P069")
timepoints <- c("DG", "REL")
compartments <- c("PBMC", "BM", "LN", "TON", "GUT", "ASC")

samples.df <-
  read.xlsx(file.path(sce.dir, "cell_metadata.xlsx")) %>%
  mutate(
    sample_norep = str_remove(sample, "_rep[0-9]"),
    patient = factor(patient, levels = patients),
    timepoint = factor(timepoint, levels = timepoints),
    compartment = factor(compartment, levels = compartments)
  )

sce <- readRDS(file = file.path(sce.dir, "sce.rds"))

sce$cell_type_manual <- colData(sce) %>%
  as.data.frame() %>%
  mutate(cell_type_manual = case_when(
    cell_type_manual == "MCL" ~ paste("MCL", timepoint, compartment),
    TRUE ~ cell_type_manual
  )) %>% pull(cell_type_manual)

k_obs_groups.list <- list(
  "P009" = 5,
  "P022" = 2,
  "P027" = 3,
  "P069" = 3,
  "P087" = 6
)

for (current_sample in unique(samples.df$patient)) {
  sample.sce <- sce[, sce$patient == current_sample]
  
  sample_annotations <- data.frame(Barcode = colnames(sample.sce),
                                   Cell_type = sample.sce$cell_type_manual)
  
  inferCNV.sample.path <-
    file.path(inferCNV.output.path, current_sample)
  dir.create(inferCNV.sample.path, showWarnings = F)
  annotations_filename <-
    paste0(current_sample, ".sampleAnnotation.tsv")
  write_tsv(
    sample_annotations,
    file = file.path(inferCNV.sample.path, annotations_filename),
    col_names = FALSE,
    quote = "none"
  )
  
  cell_groups <- unique(sample.sce$cell_type_manual)
  ref_groups <- cell_groups[!str_detect(cell_groups, "MCL")]
  
  inferCNVobj <- CreateInfercnvObject(
    as.matrix(counts(sample.sce)),
    file.path(inferCNV.sample.path, annotations_filename),
    gene_order_file = "/media/data/NGS/refs/inferCNV/hg38_gencode_v27.txt",
    ref_group_names = ref_groups
  )
  
  inferCNVobj <- infercnv::run(
    inferCNVobj,
    cutoff = 0.1,
    out_dir = inferCNV.sample.path,
    resume_mode = TRUE,
    
    cluster_by_groups = FALSE,
    cluster_references=FALSE,
    denoise = TRUE,
    
    k_obs_groups = k_obs_groups.list[[current_sample]],
    
    HMM = TRUE,
    HMM_type = 'i6',
    analysis_mode = "subclusters",
    hclust_method='ward.D2',
    
    tumor_subcluster_partition_method = "random_trees",
    tumor_subcluster_pval = 0.01,
    
    write_phylo = TRUE,
    
    no_prelim_plot = TRUE,
    output_format = "png",
    png_res = 600,
    
    num_threads = 16
  )
  
  print(paste("Sample", current_sample, "done!\n"))
}

for (current_sample in unique(samples.df$patient)) {
  inferCNV.sample.path <-
    file.path(inferCNV.output.path, current_sample)
  file.copy(
    file.path(inferCNV.sample.path, "infercnv.png"),
    file.path(inferCNV.output.path, paste0(current_sample, ".png"))
  )
}
