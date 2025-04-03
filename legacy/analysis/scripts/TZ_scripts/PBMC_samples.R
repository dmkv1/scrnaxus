
#MAIN UMAP
# data integration
PBMC_seu <- SCTransform(PBMC_samples_sobj , verbose = TRUE, conserve.memory =T)
PBMC_seu <- RunPCA(object = PBMC_seu , npcs = 50)

PBMC_seu  <- RunHarmony(PBMC_seu , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_seu, reduction = "harmony", ndims = 50)
PBMC_seu  <- RunUMAP(PBMC_seu , reduction = "harmony", dims = 1:10)
PBMC_seu  <- FindNeighbors(PBMC_seu , reduction = "harmony", dims = 1:10)
PBMC_seu  <- FindClusters(PBMC_seu, resolution =0.2)
UMAPPlot(PBMC_seu ,group.by="orig.ident", label =T)
UMAPPlot(PBMC_seu ,group.by="seurat_clusters", label =T)#+NoLegend()
UMAPPlot(PBMC_seu ,group.by="cell_type_manual", label =T)+NoLegend()
UMAPPlot(PBMC_seu ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_seu, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_seu ,group.by="seurat_clusters", cols = c("15"="red", "17"="blue"),  label =T)+NoLegend()


# calculate marker genes for each cluster
DefaultAssay(PBMC_seu) <- "SCT"
DefaultAssay(PBMC_seu) <- "RNA"
PBMC_seu <- NormalizeData(PBMC_seu, scale.factor = 10000, verbose = TRUE)
PBMC_seu_markers <- FindAllMarkers(object = PBMC_seu, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)







#PBMC_B_cells_c147
PBMC_B_cells_c147 <- subset(PBMC_seu, idents = c("1","4","7"))
PBMC_B_cells_c147 <- SCTransform(PBMC_B_cells_c147 , verbose = TRUE, conserve.memory =T)
PBMC_B_cells_c147 <- RunPCA(object = PBMC_B_cells_c147 , npcs = 50)

PBMC_B_cells_c147  <- RunHarmony(PBMC_B_cells_c147 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_B_cells_c147, reduction = "harmony", ndims = 50)
PBMC_B_cells_c147  <- RunUMAP(PBMC_B_cells_c147 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147  <- FindNeighbors(PBMC_B_cells_c147 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147  <- FindClusters(PBMC_B_cells_c147, resolution = 1.5)
UMAPPlot(PBMC_B_cells_c147 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#B-cells (PBMC: transitional stage: MME(CD10), LEF1; atypical memory B cells: FCLR3,CXCR3,ITGAX (ZBTB32, ITGAX, FCLR5,FCGR2B,CD72), memory: CD27, naive: IL4R,TCL1A, activated: TFRC )
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("TFRC","ZBTB32","ITGAX","FCRL5","FCRL3","CXCR3","STMN1","IL4R","TCL1A","CD27","PSAP","RHOB","TNFRSF1B","FCGR2B","CD72","ISG15"))

#dendritic cells
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_B_cells_c147 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_B_cells_c147) <- "SCT"
DefaultAssay(PBMC_B_cells_c147) <- "RNA"
PBMC_B_cells_c147 <- NormalizeData(PBMC_B_cells_c147, scale.factor = 10000, verbose = TRUE)
PBMC_B_cells_c147_markers <- FindAllMarkers(object = PBMC_B_cells_c147, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#PBMC_B_cells_c147_c23
PBMC_B_cells_c147_c23 <- subset(PBMC_B_cells_c147, idents = c("2","3"))
PBMC_B_cells_c147_c23 <- SCTransform(PBMC_B_cells_c147_c23 , verbose = TRUE, conserve.memory =T)
PBMC_B_cells_c147_c23 <- RunPCA(object = PBMC_B_cells_c147_c23 , npcs = 50)

PBMC_B_cells_c147_c23  <- RunHarmony(PBMC_B_cells_c147_c23 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_B_cells_c147_c23, reduction = "harmony", ndims = 50)
PBMC_B_cells_c147_c23  <- RunUMAP(PBMC_B_cells_c147_c23 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147_c23  <- FindNeighbors(PBMC_B_cells_c147_c23 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147_c23  <- FindClusters(PBMC_B_cells_c147_c23, resolution =1.0)
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
VlnPlot(PBMC_B_cells_c147_c23, group.by = "seurat_clusters", features = "TCL1A")
#dendritic cells
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_B_cells_c147_c23 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c23, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_B_cells_c147_c23) <- "SCT"
DefaultAssay(PBMC_B_cells_c147_c23) <- "RNA"
PBMC_B_cells_c147_c23 <- NormalizeData(PBMC_B_cells_c147_c23, scale.factor = 10000, verbose = TRUE)
PBMC_B_cells_c147_c23_markers <- FindAllMarkers(object = PBMC_B_cells_c147_c23, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#PBMC_B_cells_c147_c4813
PBMC_B_cells_c147_c4813 <- subset(PBMC_B_cells_c147, idents = c("4","8","13"))
PBMC_B_cells_c147_c4813 <- SCTransform(PBMC_B_cells_c147_c4813 , verbose = TRUE, conserve.memory =T)
PBMC_B_cells_c147_c4813 <- RunPCA(object = PBMC_B_cells_c147_c4813 , npcs = 50)

PBMC_B_cells_c147_c4813  <- RunHarmony(PBMC_B_cells_c147_c4813 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_B_cells_c147_c4813, reduction = "harmony", ndims = 50)
PBMC_B_cells_c147_c4813  <- RunUMAP(PBMC_B_cells_c147_c4813 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147_c4813  <- FindNeighbors(PBMC_B_cells_c147_c4813 , reduction = "harmony", dims = 1:10)
PBMC_B_cells_c147_c4813  <- FindClusters(PBMC_B_cells_c147_c4813, resolution =1.5)
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_B_cells_c147_c4813 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_B_cells_c147_c4813, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_B_cells_c147_c4813) <- "SCT"
DefaultAssay(PBMC_B_cells_c147_c4813) <- "RNA"
PBMC_B_cells_c147_c4813 <- NormalizeData(PBMC_B_cells_c147_c4813, scale.factor = 10000, verbose = TRUE)
PBMC_B_cells_c147_c4813_markers <- FindAllMarkers(object = PBMC_B_cells_c147_c4813, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)




#PBMC_myeloid_cells_c56
PBMC_myeloid_cells_c56 <- subset(PBMC_seu, idents = c("5","6"))
PBMC_myeloid_cells_c56 <- SCTransform(PBMC_myeloid_cells_c56 , verbose = TRUE, conserve.memory =T)
PBMC_myeloid_cells_c56 <- RunPCA(object = PBMC_myeloid_cells_c56 , npcs = 50)

PBMC_myeloid_cells_c56  <- RunHarmony(PBMC_myeloid_cells_c56 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_myeloid_cells_c56, reduction = "harmony", ndims = 50)
PBMC_myeloid_cells_c56  <- RunUMAP(PBMC_myeloid_cells_c56 , reduction = "harmony", dims = 1:10)
PBMC_myeloid_cells_c56  <- FindNeighbors(PBMC_myeloid_cells_c56 , reduction = "harmony", dims = 1:10)
PBMC_myeloid_cells_c56  <- FindClusters(PBMC_myeloid_cells_c56, resolution =1.5)
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","STMN1"))
#T-cells
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_myeloid_cells_c56 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_myeloid_cells_c56, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_myeloid_cells_c56) <- "SCT"
DefaultAssay(PBMC_myeloid_cells_c56) <- "RNA"
PBMC_myeloid_cells_c56 <- NormalizeData(PBMC_myeloid_cells_c56, scale.factor = 10000, verbose = TRUE)
PBMC_myeloid_cells_c56_markers <- FindAllMarkers(object = PBMC_myeloid_cells_c56, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#PBMC_TNK_cells_c0238
PBMC_TNK_cells_c0238 <- subset(PBMC_seu, idents = c("0","2","3","8"))
PBMC_TNK_cells_c0238 <- SCTransform(PBMC_TNK_cells_c0238 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238 <- RunPCA(object = PBMC_TNK_cells_c0238 , npcs = 50)

PBMC_TNK_cells_c0238  <- RunHarmony(PBMC_TNK_cells_c0238 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238  <- RunUMAP(PBMC_TNK_cells_c0238 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238  <- FindNeighbors(PBMC_TNK_cells_c0238 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238  <- FindClusters(PBMC_TNK_cells_c0238, resolution =1.0)
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238) <- "RNA"
PBMC_TNK_cells_c0238 <- NormalizeData(PBMC_TNK_cells_c0238, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#PBMC_TNK_cells_c0238_c11
PBMC_TNK_cells_c0238_c11 <- subset(PBMC_TNK_cells_c0238, idents = c("11"))
PBMC_TNK_cells_c0238_c11 <- SCTransform(PBMC_TNK_cells_c0238_c11 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238_c11 <- RunPCA(object = PBMC_TNK_cells_c0238_c11 , npcs = 50)

PBMC_TNK_cells_c0238_c11  <- RunHarmony(PBMC_TNK_cells_c0238_c11 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238_c11, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238_c11  <- RunUMAP(PBMC_TNK_cells_c0238_c11 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c11  <- FindNeighbors(PBMC_TNK_cells_c0238_c11 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c11  <- FindClusters(PBMC_TNK_cells_c0238_c11, resolution =1.5)
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238_c11 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c11, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238_c11) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238_c11) <- "RNA"
PBMC_TNK_cells_c0238_c11 <- NormalizeData(PBMC_TNK_cells_c0238_c11, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_c11_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238_c11, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#PBMC_TNK_cells_c0238_c010
PBMC_TNK_cells_c0238_c010 <- subset(PBMC_TNK_cells_c0238, idents = c("0","10"))
PBMC_TNK_cells_c0238_c010 <- SCTransform(PBMC_TNK_cells_c0238_c010 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238_c010 <- RunPCA(object = PBMC_TNK_cells_c0238_c010 , npcs = 50)

PBMC_TNK_cells_c0238_c010  <- RunHarmony(PBMC_TNK_cells_c0238_c010 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238_c010, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238_c010  <- RunUMAP(PBMC_TNK_cells_c0238_c010 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c010  <- FindNeighbors(PBMC_TNK_cells_c0238_c010 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c010  <- FindClusters(PBMC_TNK_cells_c0238_c010, resolution =2.0)
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238_c010 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c010, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238_c010) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238_c010) <- "RNA"
PBMC_TNK_cells_c0238_c010 <- NormalizeData(PBMC_TNK_cells_c0238_c010, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_c010_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238_c010, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#PBMC_TNK_cells_c0238_c1389
PBMC_TNK_cells_c0238_c1389 <- subset(PBMC_TNK_cells_c0238, idents = c("1","3","8","9"))
PBMC_TNK_cells_c0238_c1389 <- SCTransform(PBMC_TNK_cells_c0238_c1389 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238_c1389 <- RunPCA(object = PBMC_TNK_cells_c0238_c1389 , npcs = 50)

PBMC_TNK_cells_c0238_c1389  <- RunHarmony(PBMC_TNK_cells_c0238_c1389 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238_c1389, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238_c1389  <- RunUMAP(PBMC_TNK_cells_c0238_c1389 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c1389  <- FindNeighbors(PBMC_TNK_cells_c0238_c1389 , reduction = "harmony", dims = 1:10)
PBMC_TNK_cells_c0238_c1389  <- FindClusters(PBMC_TNK_cells_c0238_c1389, resolution =1.0)
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238_c1389 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c1389, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238_c1389) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238_c1389) <- "RNA"
PBMC_TNK_cells_c0238_c1389 <- NormalizeData(PBMC_TNK_cells_c0238_c1389, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_c1389_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238_c1389, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#PBMC_TNK_cells_c0238_c24567
PBMC_TNK_cells_c0238_c24567 <- subset(PBMC_TNK_cells_c0238, idents = c("2","4","5","6","7"))
PBMC_TNK_cells_c0238_c24567 <- SCTransform(PBMC_TNK_cells_c0238_c24567 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238_c24567 <- RunPCA(object = PBMC_TNK_cells_c0238_c24567 , npcs = 50)

PBMC_TNK_cells_c0238_c24567  <- RunHarmony(PBMC_TNK_cells_c0238_c24567 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238_c24567, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238_c24567  <- RunUMAP(PBMC_TNK_cells_c0238_c24567 , reduction = "harmony", dims = 1:14)
PBMC_TNK_cells_c0238_c24567  <- FindNeighbors(PBMC_TNK_cells_c0238_c24567 , reduction = "harmony", dims = 1:14)
PBMC_TNK_cells_c0238_c24567  <- FindClusters(PBMC_TNK_cells_c0238_c24567, resolution =1.5)
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238_c24567 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238_c24567) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238_c24567) <- "RNA"
PBMC_TNK_cells_c0238_c24567 <- NormalizeData(PBMC_TNK_cells_c0238_c24567, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_c24567_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238_c24567, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#PBMC_TNK_cells_c0238_c24567_c9
PBMC_TNK_cells_c0238_c24567_c9 <- subset(PBMC_TNK_cells_c0238_c24567, idents = c("9"))
PBMC_TNK_cells_c0238_c24567_c9 <- SCTransform(PBMC_TNK_cells_c0238_c24567_c9 , verbose = TRUE, conserve.memory =T)
PBMC_TNK_cells_c0238_c24567_c9 <- RunPCA(object = PBMC_TNK_cells_c0238_c24567_c9 , npcs = 50)

PBMC_TNK_cells_c0238_c24567_c9  <- RunHarmony(PBMC_TNK_cells_c0238_c24567_c9 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(PBMC_TNK_cells_c0238_c24567_c9, reduction = "harmony", ndims = 50)
PBMC_TNK_cells_c0238_c24567_c9  <- RunUMAP(PBMC_TNK_cells_c0238_c24567_c9 , reduction = "harmony", dims = 1:14)
PBMC_TNK_cells_c0238_c24567_c9  <- FindNeighbors(PBMC_TNK_cells_c0238_c24567_c9 , reduction = "harmony", dims = 1:14)
PBMC_TNK_cells_c0238_c24567_c9  <- FindClusters(PBMC_TNK_cells_c0238_c24567_c9, resolution =1.0)
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMH","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#dnT
UMAPPlot(PBMC_TNK_cells_c0238_c24567_c9 ,group.by="seurat_clusters", label =T)+ FeaturePlot(PBMC_TNK_cells_c0238_c24567_c9, features =c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1", "CD8A","CD3D"))


# calculate marker genes for each cluster
DefaultAssay(PBMC_TNK_cells_c0238_c24567_c9) <- "SCT"
DefaultAssay(PBMC_TNK_cells_c0238_c24567_c9) <- "RNA"
PBMC_TNK_cells_c0238_c24567_c9 <- NormalizeData(PBMC_TNK_cells_c0238_c24567_c9, scale.factor = 10000, verbose = TRUE)
PBMC_TNK_cells_c0238_c24567_c9_markers <- FindAllMarkers(object = PBMC_TNK_cells_c0238_c24567_c9, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)




cluster_export <- read.csv(file = "/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/cluster_export_file_PBMC.csv", sep = ",", header = TRUE)
is_valid = TRUE
collected_cells <- character(0)

merged = data.frame()
for (i in 1:nrow(cluster_export)) {
  cluster_number = cluster_export[i,"cluster_number"]
  cluster_name = cluster_export[i,"cluster_name"]
  subset = cluster_export[i, "subset"]
  
  cell<-WhichCells(get(subset), idents = cluster_number)
  cell<- data.frame(cell)
  cell$cluster<-cluster_name
  cell$number<-cluster_number
  cell$object<-subset
  
  merged<-rbind(merged,cell)
}

colnames(merged) <- c("cell","cluster","cluster_number","original_object")
rownames(merged) <- merged$cell
merged_PBMC <- merged[,c(2,3,4)]
PBMC_seu <- AddMetaData(PBMC_seu, metadata = merged_PBMC)

UMAPPlot(PBMC_seu, group.by = "cluster", raster = FALSE, label = T, repel = T)
DotPlot(PBMC_seu, group.by = "cluster", features = c("SLC4A10", "CD3E","CD8A","CD4","FOXP3","SELL","CCR7", "CXCL13","GNLY","KLRF1", "MS4A1","TCL1A","CD27","VPREB1","CD38","PLD4", "IGHM","COCH","FCRL5","SIGLEC6", "COTL1", "MZB1","LYZ","C1QA","FCN1","LAMP3","CLEC9A","CDC1","CLEC10A", "STMN1","GZMK", "GZMH"), cluster.idents=T, cols = c("red", "green")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(x = "gene", y = "cell type")

setdiff(Cells(PBMC_seu), merged$cell)




# calculation of effector memory score for CD4 cells
# subsetting of CD4 memory cells
CD4_memory_T <- grep("CD4_effector_memory|CD4_central_memory", unique(PBMC_seu$cluster), value = TRUE)
Idents(PBMC_seu) <- PBMC_seu$cluster
CD4_memory_T_sobj <- subset(PBMC_seu, idents = CD4_memory_T)

# CD4 effector memory gene set score calculation
CD4_EM_genes <- list(c("CCL5","LYAR","KLRB1","GZMK","IL7R","CST7","CXCR4","IL32"))
CD4_memory_T_sobj <- AddModuleScore(object=CD4_memory_T_sobj, features=CD4_EM_genes, ctrl=100, name='gene_set_score_CD4_EM')
VlnPlot(CD4_memory_T_sobj, features = "gene_set_score_CD4_EM1")+ geom_hline(yintercept = 0.52)

# to rename variables in specific rows
expression_matrix <- as.data.frame(CD4_memory_T_sobj@meta.data)
cells_above_threshold <- rownames(expression_matrix)[expression_matrix[,"gene_set_score_CD4_EM1"] > 0.52]
CD4_memory_T_sobj$cluster[cells_above_threshold] <- "CD4_effector_memory_T-cell"
Idents(CD4_memory_T_sobj) <- CD4_memory_T_sobj$cluster
VlnPlot(CD4_memory_T_sobj, features = "gene_set_score_CD4_EM1")+ geom_hline(yintercept = 0.52)

# to add the annotation "CD4_effector_memory_T-cell" into PBMC_seu$cluster
as.data.frame(table(PBMC_seu$cluster)) %>% dplyr::arrange(-Freq)
PBMC_seu$cluster[cells_above_threshold] <- "CD4_effector_memory_T-cell"
as.data.frame(table(PBMC_seu$cluster)) %>% dplyr::arrange(-Freq)


# calculating cytotoxic gene score on all T-cells
Tcell <- grep("T-cell|Treg", unique(PBMC_seu$cluster), value = TRUE)
Idents(PBMC_seu) <- PBMC_seu$cluster
T_sobj <- subset(PBMC_seu, idents = Tcell)
T_sobj <- AddModuleScore(T_sobj, features = list(c("GNLY","NKG7","GZMH","PRF1","GZMB","GZMA","GZMM","KLRG1","KLRC1","KLRK1","CRTAM","IFNG")), name = "CTL_signature")
VlnPlot(T_sobj, features = "CTL_signature1")+ geom_hline(yintercept = 0)

# setting the threshold for identification of CD4 cytotoxic T-cells
CTL_signature_threshold <- 0

# subsetting all CD4 T-cells  from T-sobj 
CD4_T_cells <- grep("CD4", unique(T_sobj$cluster), value=T)
T_sobj_CD4 <- subset(T_sobj, idents = CD4_T_cells)

# to annotate CD4 T-cells above thereshold as CD4 cytotoxic T-cells
CD4_expression_matrix <- as.data.frame(T_sobj_CD4@meta.data)
CD4_cells_above_threshold <- rownames(CD4_expression_matrix)[CD4_expression_matrix[,"CTL_signature1"] > 0]
T_sobj_CD4$cluster[CD4_cells_above_threshold] <- "CD4_cytotoxic_T-cell"
Idents(T_sobj_CD4) <- T_sobj_CD4$cluster
VlnPlot(T_sobj_CD4, features = "CTL_signature1")+ geom_hline(yintercept = 0)

# to add the annotation "CD4_cytotoxic_T-cell" into PBMC_seu$cluster
as.data.frame(table(PBMC_seu$cluster)) %>% dplyr::arrange(-Freq)
PBMC_seu$cluster[CD4_cells_above_threshold] <- "CD4_cytotoxic_T-cell"
as.data.frame(table(PBMC_seu$cluster)) %>% dplyr::arrange(-Freq)

# to remove both temporary sobjects
rm(T_sobj)
rm(T_sobj_CD4)
rm(PBMC_seu_copy)

