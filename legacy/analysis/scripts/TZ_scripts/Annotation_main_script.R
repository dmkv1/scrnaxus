
library(tidyverse)
library(cowplot)
library(openxlsx)

# SC methods
library(scran)
library(scater)
library(bluster)
library(Seurat)

# Annotations
library(SingleR)
library(tricycle)
library(harmony)

#to load Single cell experiment
sce <- readRDS(file = "/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/outputs/4_all_cells_labbeled/sce.rds")

#to make seurat object from Single cell experiment 
sce.to.seu <- sce
rownames(sce.to.seu) <- make.unique(rowData(sce.to.seu)[["Symbol"]])
seu <- SeuratObject::as.Seurat(sce.to.seu, counts = "counts", data = "logcounts")

# to add nFeature, nCount and percent.mt into metadata
seu$orig.ident <- seu$Sample
seu@assays[["RNA"]] <- seu@assays[["originalexp"]] 
seu[["percent_mt"]] <- PercentageFeatureSet(object =seu, pattern = "^MT-")
seu$nCount_RNA <- colSums(seu@assays$RNA@counts)
seu$nFeature_RNA <- colSums(seu@assays$RNA@counts > 0)

# to check the filtering
VlnPlot(seu , features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, group.by = "orig.ident")

# subset of cells with acceptable percent.mt i.e. less than 15%
seu_split <- SplitObject(seu, split.by= "orig.ident")

for (i in names(seu_split)) {
seu_split[[i]] <-
subset(x = seu_split[[i]], subset = nCount_RNA > 300 &nFeature_RNA > 200 & nFeature_RNA & percent_mt < 15)
}

seu_f  <- merge(x=seu_split [[1]], y=seu_split [2:length(seu_split)])

# to check the filtering
VlnPlot(seu_f , features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, group.by = "orig.ident")

## Table of metric summaries
library(dplyr)
library(tibble)
dplyr::bind_cols(
  tibble::tibble(metric = c("nCount_RNA", "nFeature_RNA", "percent_mt")),
  dplyr::bind_rows(summary(seu_f_sct$nCount_RNA), summary(seu_f_sct$nFeature_RNA), summary(seu_f_sct$percent_mt))
)


# data integration
seu_f_sct <- SCTransform(seu_f , verbose = TRUE, conserve.memory =T)
seu_f_sct <- RunPCA(object = seu_f_sct , npcs = 50)

# to check samples on UMAP before bath effect correction
seu_f_sct <- RunUMAP(seu_f_sct , dims = 1:30)
seu_f_sct  <- FindNeighbors(seu_f_sct, dims = 1:14)
seu_f_sct  <- FindClusters(seu_f_sct, resolution =0.2)
UMAPPlot(seu_f_sct ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(seu_f_sct ,group.by="seurat_clusters", label =T)+NoLegend()
UMAPPlot(seu_f_sct ,group.by="cell_type_manual", label =T)+NoLegend()+UMAPPlot(seu_f_sct ,group.by="seurat_clusters", label =T)
UMAPPlot(seu_f_sct ,group.by="seurat_clusters", label =T)+ FeaturePlot(seu_f_sct, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(seu_f_sct ,group.by="seurat_clusters", cols = c("25"="red", "23"="blue"),  label =T)+NoLegend()


# calculate marker genes for each cluster
DefaultAssay(seu_f_sct) <- "SCT"
DefaultAssay(seu_f_sct) <- "RNA"
seu_f_sct <- NormalizeData(seu_f_sct, scale.factor = 10000, verbose = TRUE)
seu_f_sct_markers <- FindAllMarkers(object = seu_f_sct, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)


# to run harmony (it doesn't help for identification of main cells types - B cells co-cluster with MCL cells, better to not use batch effect correction for main UMAP)
# seu_f_sct_harm  <- RunHarmony(seu_f_sct , group.by.vars = c("orig.ident"), assay.use =  "SCT")
# ElbowPlot(seu_f_sct_harm, reduction = "harmony", ndims = 50)
# seu_f_sct_harm  <- RunUMAP(seu_f_sct_harm , reduction = "harmony", dims = 1:14)
# seu_f_sct_harm  <- FindNeighbors(seu_f_sct_harm , reduction = "harmony", dims = 1:14)
# seu_f_sct_harm  <- FindClusters(seu_f_sct_harm, resolution =0.2)
# 
# UMAPPlot(seu_f_sct_harm ,group.by="seurat_clusters", label =T)#+NoLegend()
# UMAPPlot(seu_f_sct_harm ,group.by="seurat_clusters", label =T)+NoLegend()+ FeaturePlot(seu_f_sct_harm, features =c("JCHAIN","MS4A1","CD3D","FOXP3","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
# UMAPPlot(seu_f_sct_harm ,group.by="cell_type_manual", label =T)#+NoLegend()
# UMAPPlot(seu_f_sct_harm ,group.by="seurat_clusters", label =T)+NoLegend()+ FeaturePlot(seu_f_sct_harm, features =c("IGLL1","CD38","CD34","IL7R","CD22","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
# 
# # calculate marker genes for each cluster
# DefaultAssay(seu_f_sct_harm) <- "SCT"
# DefaultAssay(seu_f_sct_harm) <- "RNA"
# seu_f_sct_harm <- NormalizeData(seu_f_sct_harm, scale.factor = 10000, verbose = TRUE)
# seu_f_sct_harm_markers <- FindAllMarkers(object = seu_f_sct_harm, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)

# cells filtered out based on high percent.mt will get annotation lowQ
lowQ <- setdiff(Cells(seu),Cells(seu_f_sct))

# identification of proliferating cells
gene_list_prolif <- list(union(cc.genes$s.genes, cc.genes$g2m.genes))
seu_f_sct <- AddModuleScore(seu_f_sct, features = gene_list_prolif, name = "prolif_signature")
FeaturePlot(seu_f_sct, reduction = "umap", features = "prolif_signature1", raster = FALSE)

VlnPlot(seu_f_sct, features = "prolif_signature1") +
  NoLegend() &
  geom_hline(yintercept=0.06)

# first subset all normal non-lyphoma cells which form separate cluster
Idents(seu_f_sct) <- seu_f_sct$seurat_clusters
non_tumor_cells <- subset(seu_f_sct, idents = c("4","8","1","10","13","16"))

# second, identify non-lignant cells in other clusters (B and IGLL1)
###identification of non tumor B-cells ################
Bcell <- subset(seu_f_sct, idents = c("0","14"))
Bcell <- SCTransform(Bcell , verbose = TRUE, conserve.memory =T)
#Bcell [["SCT"]]@var.features <- Bcell [["SCT"]]@var.features[!(Bcell [["SCT"]]@var.features %in% VDJTgenes_filter_out)]
#sobj_filtered_sct [["SCT"]]@var.features <- sobj_filtered_sct [["SCT"]]@var.features[!(sobj_filtered_sct [["SCT"]]@var.features %in% filtered_genes$VDJ_HS)]
Bcell <- RunPCA(object = Bcell , npcs = 50)

Bcell  <- RunHarmony(Bcell , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(Bcell, reduction = "harmony", ndims = 50)
Bcell  <- RunUMAP(Bcell , reduction = "harmony", dims = 1:10)
Bcell  <- FindNeighbors(Bcell , reduction = "harmony", dims = 1:10)
Bcell  <- FindClusters(Bcell, resolution =0.2)

UMAPPlot(Bcell ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(Bcell ,group.by="seurat_clusters", label =T)+NoLegend()+ FeaturePlot(Bcell, features =c("JCHAIN","MS4A1","CD3D","FOXP3","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(Bcell ,group.by="seurat_clusters", cols = c("2"="red", "6"="blue"),  label =T)+NoLegend()

non_tumor_B_cells <- subset(Bcell, idents = c("2","6","9"))


###identificaion of (potentially) non tumor IGLL1 positive cells
IGLL1_cells <- subset(seu_f_sct, idents = c("3"))
IGLL1_cells <- SCTransform(IGLL1_cells , verbose = TRUE, conserve.memory =T)
#IGLL1_cells [["SCT"]]@var.features <- IGLL1_cells [["SCT"]]@var.features[!(IGLL1_cells [["SCT"]]@var.features %in% VDJTgenes_filter_out)]
#sobj_filtered_sct [["SCT"]]@var.features <- sobj_filtered_sct [["SCT"]]@var.features[!(sobj_filtered_sct [["SCT"]]@var.features %in% filtered_genes$VDJ_HS)]
IGLL1_cells <- RunPCA(object = IGLL1_cells , npcs = 50)

#IGLL1_cells  <- RunHarmony(IGLL1_cells , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(IGLL1_cells, ndims = 50)
IGLL1_cells  <- RunUMAP(IGLL1_cells , dims = 1:10)
IGLL1_cells  <- FindNeighbors(IGLL1_cells , dims = 1:10)
IGLL1_cells  <- FindClusters(IGLL1_cells, resolution =0.2)

UMAPPlot(IGLL1_cells ,group.by="seurat_clusters", label =T)#+NoLegend()
UMAPPlot(IGLL1_cells ,group.by="seurat_clusters", label =T)+NoLegend()+ FeaturePlot(IGLL1_cells, features =c("JCHAIN","MS4A1","CD3D","FOXP3","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "IGLL1","CCND1","LAG3"))

non_tumor_IGLL1_cells <- subset(IGLL1_cells, idents = c("5"))

# merge all non-tumor cells into one sobj
non_tumor_sobj  <- merge(x=non_tumor_cells, y= non_tumor_B_cells)
non_tumor_sobj <- merge(x=non_tumor_sobj, y= non_tumor_IGLL1_cells)

# pull barcodes of all MCL cell
MCL_cells <- setdiff(Cells(seu_f_sct),Cells(non_tumor_sobj))

# since in PBMC and tissue are different cell types, these two compartments need to be separated.
PBMC_samples <- grep("PBMC", unique(non_tumor_sobj$Sample), value = TRUE)
Idents(non_tumor_sobj) <- non_tumor_sobj$Sample
PBMC_samples_sobj <- subset(non_tumor_sobj, idents = PBMC_samples)
tissue_samples_sobj <- subset(non_tumor_sobj, idents = PBMC_samples, invert = TRUE)

## from now use scripts in R Script files: PBMC_samples and tissue_samples

# after the annotation of non-tumor cells is finished, make final annotation file of all cells in the dataset

final_annotation <- seu@meta.data
final_annotation$sample_barcode <- rownames(final_annotation)
final_annotation <- dplyr::select(final_annotation, "sample" = "orig.ident", sample_barcode, compartment)
final_annotation$cluster_name <- "cluster"

# to add annotations for lowQ cells = "low_quality_cell"
final_annotation[lowQ, "cluster_name"] <- "low_quality_cell"

# to add annotations for MCL cells = malignant
final_annotation[MCL_cells, "cluster_name"] <- "malignant"

# to add annotations from PBMC_seu sobj
PBMC_data_frame <- as.data.frame(PBMC_seu@meta.data)
rows_to_change_PBMC <- match(rownames(PBMC_data_frame), final_annotation$sample_barcode)
final_annotation$cluster_name[rows_to_change_PBMC] <- PBMC_data_frame$cluster

# to add annotations from tissue_seu sobj
tissue_data_frame <- as.data.frame(tissue_seu@meta.data)
rows_to_change_tissue <- match(rownames(tissue_data_frame), final_annotation$sample_barcode)
final_annotation$cluster_name[rows_to_change_tissue] <- tissue_data_frame$cluster

final_annotation$cluster_name[final_annotation$cluster_name == "MCL"] <- "malignant"

as.data.frame(table(final_annotation$cluster_name)) %>% dplyr::arrange(-Freq)
write.csv(final_annotation, "/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/final_annotation.csv")

final_annotation_cell_states <- final_annotation
final_annotation_cell_states$cell_type <- NA
final_annotation_cell_states$cell_state_1 <- NA
final_annotation_cell_states$cell_state_2 <- NA

cell_state_names <- read.csv("/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/cell_states_names.csv")
library(plyr)
final_annotation_cell_states$cell_type <- plyr::mapvalues(final_annotation_cell_states$cluster_name, from = cell_state_names$cluster_name, to = cell_state_names$cell_type)
final_annotation_cell_states$cell_state_1 <- plyr::mapvalues(final_annotation_cell_states$cluster_name, from = cell_state_names$cluster_name, to = cell_state_names$cell_state_1)
final_annotation_cell_states$cell_state_2 <- plyr::mapvalues(final_annotation_cell_states$cluster_name, from = cell_state_names$cluster_name, to = cell_state_names$cell_state_c2)

final_annotation_cell_states$cluster_name[final_annotation_cell_states$cluster_name == "unswitched_memory_B-cell "] <- "unswitched_memory_B-cell"
final_annotation_cell_states$cluster_name[final_annotation_cell_states$cluster_name == "intermediate_monocyte "] <- "intermediate_monocyte"
final_annotation_cell_states$cluster_name[final_annotation_cell_states$cluster_name == "atypical_memory_B-cell "] <- "atypical_memory_B-cell"
final_annotation_cell_states$cluster_name[final_annotation_cell_states$cluster_name == "transitional_stage_B-cell "] <- "transitional_stage_B-cell"

as.data.frame(table(final_annotation_cell_states$cluster_name)) %>% dplyr::arrange(-Freq)
as.data.frame(table(final_annotation_cell_states$cell_type)) %>% dplyr::arrange(-Freq)

final_annotation <- read.csv("/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/final_annotation.csv")
meta <- dplyr::select(final_annotation, sample_barcode, cluster_name, cell_type, cell_state_1, cell_state_2)
rownames(meta) <- meta$sample_barcode
seu_f_sct <- AddMetaData(seu_f_sct, metadata = meta)
seu_f_sct_wo_doublets <- subset(seu_f_sct, cell_type != "doublet")
seu_f_sct_wo_doublets_malignant <- subset(seu_f_sct_wo_doublets, cell_type != "malignant")

# Percentual distribution of predicted cell types accross individual samples (orig.ident)
ggplot(seu_f_sct_wo_doublets@meta.data, aes(x=orig.ident, fill=cell_type)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot(seu_f_sct_wo_doublets_malignant@meta.data, aes(x=orig.ident, fill=cell_type)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


as.data.frame(table(final_annotation$cell_type)) %>% dplyr::arrange(-Freq)


# Extract metadata
seu_f_sct_wo_doublets_malignant_df <- seu_f_sct_wo_doublets_malignant@meta.data

# Calculate frequencies of cell types for each sample
cell_type_frequencies <- seu_f_sct_wo_doublets_malignant_df %>%
  group_by(orig.ident, cell_type) %>%
  summarise(count = n()) %>%
  mutate(frequency = count / sum(count))

write.csv(cell_type_frequencies, "/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/cell_type_frequencies.csv")

# Plot the frequencies
ggplot(cell_type_frequencies, aes(x = orig.ident, y = frequency, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Frequency", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# to rename variables in specific rows
rows_to_change <- match(sobj$cell_barcode, df$cell_barcode)
df$cell_type[rows_to_change] <- sobj$cell_type



