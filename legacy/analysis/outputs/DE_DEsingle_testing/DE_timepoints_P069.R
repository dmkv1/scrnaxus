suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(scater))
suppressMessages(library(DEsingle))

# Parallelization
suppressMessages(library(BiocParallel))
bpp <- MulticoreParam(workers = 48)

# Path to project folder
wd <- "/media/data/NGS/Projects/Clonal2/"
source(file.path(wd, "scRNAseq/analysis/scripts/colors.R"))

output.path <- file.path(wd, "scRNAseq/analysis/outputs/DE_testing")

iteration_input = "4_all_cells_labbeled"
sce.dir <- file.path(wd, "scRNAseq/analysis/outputs", iteration_input)
sce.labelled <- readRDS(file = file.path(sce.dir, "sce.rds"))

sce.labelled$timepoint <- factor(sce.labelled$timepoint, levels = c("REL", "DG"))

patients <- sce.labelled$patient %>% levels()

DE.timepoints.results <- list()
for(current_patient in c("P069")){
  sce <- sce.labelled[, sce.labelled$patient == current_patient]
  
  # sampling
  #sce <- sce[which(colnames(sce) %in% sample(colnames(sce), size = dim(sce)[2] / 20)), ]
  
  DE.timepoints.results[[current_patient]] <- DEsingle(counts = sce,
                                                       group = sce$timepoint, 
                                                       parallel = FALSE)
  DE.timepoints.results[[current_patient]] %>% 
    rownames_to_column("gene") %>% 
    write.xlsx(.,
               file = file.path(output.path,
                                paste0("DE_timepoint_", current_patient, ".xlsx"))
    )
  
  gc()
}


