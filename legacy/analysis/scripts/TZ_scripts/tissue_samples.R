
#MAIN UMAP
# data integration
tissue_seu <- SCTransform(tissue_samples_sobj , verbose = TRUE, conserve.memory =T)
tissue_seu <- RunPCA(object = tissue_seu , npcs = 50)

tissue_seu  <- RunHarmony(tissue_seu , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_seu, reduction = "harmony", ndims = 50)
tissue_seu  <- RunUMAP(tissue_seu , reduction = "harmony", dims = 1:10)
tissue_seu  <- FindNeighbors(tissue_seu , reduction = "harmony", dims = 1:10)
tissue_seu  <- FindClusters(tissue_seu, resolution =0.2)
UMAPPlot(tissue_seu ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_seu ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_seu ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_seu, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_seu ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_seu, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_seu, features =c("JCHAIN","MS4A1","CD3D","CD4","FCER1G","CD40LG","CD8A","CD8B","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_seu ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_seu, features =c("SOX4","MPO","RNASE2","MS4A3","HBB","S100A4","CCL5","GSTP1","PLEK", "APOC1","CCND1","MS4A2"))


# calculate marker genes for each cluster
DefaultAssay(tissue_seu) <- "SCT"
DefaultAssay(tissue_seu) <- "RNA"
tissue_seu <- NormalizeData(tissue_seu, scale.factor = 10000, verbose = TRUE)
tissue_seu_markers <- FindAllMarkers(object = tissue_seu, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)



#tissue_B_cells_c3
tissue_B_cells_c3 <- subset(tissue_seu, idents = c("3"))
tissue_B_cells_c3 <- SCTransform(tissue_B_cells_c3 , verbose = TRUE, conserve.memory =T)
tissue_B_cells_c3 <- RunPCA(object = tissue_B_cells_c3 , npcs = 50)

tissue_B_cells_c3  <- RunHarmony(tissue_B_cells_c3 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_B_cells_c3, reduction = "harmony", ndims = 50)
tissue_B_cells_c3  <- RunUMAP(tissue_B_cells_c3 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3  <- FindNeighbors(tissue_B_cells_c3 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3  <- FindClusters(tissue_B_cells_c3, resolution =1.5)
UMAPPlot(tissue_B_cells_c3 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#pre_B-cells
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("RAG1","RAG2","UHRF1" ,"BLNK","AGPS","CYGB","UMODL1","EBF1","MME","VPREB1","MS4A1","STMN1"))

#dendritic cells
UMAPPlot(tissue_B_cells_c3 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))


# calculate marker genes for each cluster
DefaultAssay(tissue_B_cells_c3) <- "SCT"
DefaultAssay(tissue_B_cells_c3) <- "RNA"
tissue_B_cells_c3 <- NormalizeData(tissue_B_cells_c3, scale.factor = 10000, verbose = TRUE)
tissue_B_cells_c3_markers <- FindAllMarkers(object = tissue_B_cells_c3, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#tissue_B_cells_c3_c1357
tissue_B_cells_c3_c1357 <- subset(tissue_B_cells_c3, idents = c("1", "3","5","7"))
tissue_B_cells_c3_c1357 <- SCTransform(tissue_B_cells_c3_c1357 , verbose = TRUE, conserve.memory =T)
tissue_B_cells_c3_c1357 <- RunPCA(object = tissue_B_cells_c3_c1357 , npcs = 50)

tissue_B_cells_c3_c1357  <- RunHarmony(tissue_B_cells_c3_c1357 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_B_cells_c3_c1357, reduction = "harmony", ndims = 50)
tissue_B_cells_c3_c1357  <- RunUMAP(tissue_B_cells_c3_c1357 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3_c1357  <- FindNeighbors(tissue_B_cells_c3_c1357 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3_c1357  <- FindClusters(tissue_B_cells_c3_c1357, resolution =1.5)
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#pre_B-cells
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("RAG1","RAG2","VPREB1" ,"IGLL1","SOX4","IGHM","PAX5","EBF1","CD24","VPREB1","MS4A1","STMN1"))
#dendritic cells
UMAPPlot(tissue_B_cells_c3_c1357 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c1357, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))


# calculate marker genes for each cluster
DefaultAssay(tissue_B_cells_c3_c1357) <- "SCT"
DefaultAssay(tissue_B_cells_c3_c1357) <- "RNA"
tissue_B_cells_c3_c1357 <- NormalizeData(tissue_B_cells_c3_c1357, scale.factor = 10000, verbose = TRUE)
tissue_B_cells_c3_c1357_markers <- FindAllMarkers(object = tissue_B_cells_c3_c1357, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#tissue_B_cells_c3_c468
tissue_B_cells_c3_c468 <- subset(tissue_B_cells_c3, idents = c("4","6","8"))
tissue_B_cells_c3_c468 <- SCTransform(tissue_B_cells_c3_c468 , verbose = TRUE, conserve.memory =T)
tissue_B_cells_c3_c468 <- RunPCA(object = tissue_B_cells_c3_c468 , npcs = 50)

tissue_B_cells_c3_c468  <- RunHarmony(tissue_B_cells_c3_c468 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_B_cells_c3_c468, reduction = "harmony", ndims = 50)
tissue_B_cells_c3_c468  <- RunUMAP(tissue_B_cells_c3_c468 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3_c468  <- FindNeighbors(tissue_B_cells_c3_c468 , reduction = "harmony", dims = 1:10)
tissue_B_cells_c3_c468  <- FindClusters(tissue_B_cells_c3_c468, resolution =1.5)
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3_c468 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_B_cells_c3_c468 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c468, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c468, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c468, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c468, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_B_cells_c3_c468 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_B_cells_c3_c468, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))


# calculate marker genes for each cluster
DefaultAssay(tissue_B_cells_c3_c468) <- "SCT"
DefaultAssay(tissue_B_cells_c3_c468) <- "RNA"
tissue_B_cells_c3_c468 <- NormalizeData(tissue_B_cells_c3_c468, scale.factor = 10000, verbose = TRUE)
tissue_B_cells_c3_c468_markers <- FindAllMarkers(object = tissue_B_cells_c3_c468, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)







#tissue_myeloid_cells_c7
tissue_myeloid_cells_c7 <- subset(tissue_seu, idents = c("7"))
tissue_myeloid_cells_c7 <- SCTransform(tissue_myeloid_cells_c7 , verbose = TRUE, conserve.memory =T)
tissue_myeloid_cells_c7 <- RunPCA(object = tissue_myeloid_cells_c7 , npcs = 50)

tissue_myeloid_cells_c7  <- RunHarmony(tissue_myeloid_cells_c7 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_myeloid_cells_c7, reduction = "harmony", ndims = 50)
tissue_myeloid_cells_c7  <- RunUMAP(tissue_myeloid_cells_c7 , reduction = "harmony", dims = 1:10)
tissue_myeloid_cells_c7  <- FindNeighbors(tissue_myeloid_cells_c7 , reduction = "harmony", dims = 1:10)
tissue_myeloid_cells_c7  <- FindClusters(tissue_myeloid_cells_c7, resolution =1.5)
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_myeloid_cells_c7 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_myeloid_cells_c7 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_myeloid_cells_c7 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))


# calculate marker genes for each cluster
DefaultAssay(tissue_myeloid_cells_c7) <- "SCT"
DefaultAssay(tissue_myeloid_cells_c7) <- "RNA"
tissue_myeloid_cells_c7 <- NormalizeData(tissue_myeloid_cells_c7, scale.factor = 10000, verbose = TRUE)
tissue_myeloid_cells_c7_markers <- FindAllMarkers(object = tissue_myeloid_cells_c7, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#tissue_myeloid_cells_c7_c14
tissue_myeloid_cells_c7_c14 <- subset(tissue_myeloid_cells_c7, idents = c("1","4"))
tissue_myeloid_cells_c7_c14 <- SCTransform(tissue_myeloid_cells_c7_c14 , verbose = TRUE, conserve.memory =T)
tissue_myeloid_cells_c7_c14 <- RunPCA(object = tissue_myeloid_cells_c7_c14 , npcs = 50)

tissue_myeloid_cells_c7_c14  <- RunHarmony(tissue_myeloid_cells_c7_c14 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_myeloid_cells_c7_c14, reduction = "harmony", ndims = 50)
tissue_myeloid_cells_c7_c14  <- RunUMAP(tissue_myeloid_cells_c7_c14 , reduction = "harmony", dims = 1:10)
tissue_myeloid_cells_c7_c14  <- FindNeighbors(tissue_myeloid_cells_c7_c14 , reduction = "harmony", dims = 1:10)
tissue_myeloid_cells_c7_c14  <- FindClusters(tissue_myeloid_cells_c7_c14, resolution =1.5)
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_myeloid_cells_c7_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_myeloid_cells_c7_c14, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))


# calculate marker genes for each cluster
DefaultAssay(tissue_myeloid_cells_c7_c14) <- "SCT"
DefaultAssay(tissue_myeloid_cells_c7_c14) <- "RNA"
tissue_myeloid_cells_c7_c14 <- NormalizeData(tissue_myeloid_cells_c7_c14, scale.factor = 10000, verbose = TRUE)
tissue_myeloid_cells_c7_c14_markers <- FindAllMarkers(object = tissue_myeloid_cells_c7_c14, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






##################################
#tissue_TNK_cells_c012456
tissue_TNK_cells_c012456 <- subset(tissue_seu, idents = c("0","1","2","4","5","6"))
tissue_TNK_cells_c012456 <- SCTransform(tissue_TNK_cells_c012456 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456 <- RunPCA(object = tissue_TNK_cells_c012456 , npcs = 50)

tissue_TNK_cells_c012456  <- RunHarmony(tissue_TNK_cells_c012456 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456  <- RunUMAP(tissue_TNK_cells_c012456 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456  <- FindNeighbors(tissue_TNK_cells_c012456 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456  <- FindClusters(tissue_TNK_cells_c012456, resolution =1.5)
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("LYZ","MS4A1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","TPSAB1","KLRB1"))
UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("CD4","CD8A","CD3D","FGFBP2","CD40LG","KLRC3","SELL","XCL1","KLRD1","TRDC","TPSAB1","KLRB1"))

UMAPPlot(tissue_TNK_cells_c012456 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456, features =c("SOX4","MPO","RNASE2","MS4A3","HBB","S100A4","CCL5","GSTP1","PLEK", "APOC1","CCND1","MS4A2"))


# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456) <- "RNA"
tissue_TNK_cells_c012456 <- NormalizeData(tissue_TNK_cells_c012456, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#tissue_TNK_cells_c012456_c14
tissue_TNK_cells_c012456_c14 <- subset(tissue_TNK_cells_c012456, idents = c("14"))
tissue_TNK_cells_c012456_c14 <- SCTransform(tissue_TNK_cells_c012456_c14 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c14 <- RunPCA(object = tissue_TNK_cells_c012456_c14 , npcs = 50)

tissue_TNK_cells_c012456_c14  <- RunHarmony(tissue_TNK_cells_c012456_c14 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c14, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c14  <- RunUMAP(tissue_TNK_cells_c012456_c14 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c14  <- FindNeighbors(tissue_TNK_cells_c012456_c14 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c14  <- FindClusters(tissue_TNK_cells_c012456_c14, resolution =2.0)
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))

UMAPPlot(tissue_TNK_cells_c012456_c14 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c14, features =c("SOX4","VPREB1","IGLL1","MS4A3","HBB","S100A4","CCL5","GSTP1","PLEK", "APOC1","CCND1","MS4A2"))

# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c14) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c14) <- "RNA"
tissue_TNK_cells_c012456_c14 <- NormalizeData(tissue_TNK_cells_c012456_c14, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c14_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c14, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)








#tissue_TNK_cells_c012456_c491316
tissue_TNK_cells_c012456_c491316 <- subset(tissue_TNK_cells_c012456, idents = c("4","9","13","16"))
tissue_TNK_cells_c012456_c491316 <- SCTransform(tissue_TNK_cells_c012456_c491316 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c491316 <- RunPCA(object = tissue_TNK_cells_c012456_c491316 , npcs = 50)

tissue_TNK_cells_c012456_c491316  <- RunHarmony(tissue_TNK_cells_c012456_c491316 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c491316, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c491316  <- RunUMAP(tissue_TNK_cells_c012456_c491316 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c491316  <- FindNeighbors(tissue_TNK_cells_c012456_c491316 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c491316  <- FindClusters(tissue_TNK_cells_c012456_c491316, resolution =2.0)
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T) +UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#CD4 T-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("CXCR3","NKG7","CD3D","CD4","CXCL13","CCL5","CD8A","GZMK","CCR6","TBX21","S100A4","KLRB1"))
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("CXCR3","NKG7","CD3D","CD4","SELL","TCF7","CCR7","LEF1","CCR6","KLRF1","S100A4","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))
#CD4 TCM
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"))

#CD4 TEM
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3", "COTL1", "NKG7"))
UMAPPlot(tissue_TNK_cells_c012456_c491316 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316, features =c("IL7R", "CCL5", "GNLY", "GZMK", "GZMH", "SELL", "KLRB1", "TCF7", "CCR7", "LEF1", "COTL1", "NKG7"))


# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c491316) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c491316) <- "RNA"
tissue_TNK_cells_c012456_c491316 <- NormalizeData(tissue_TNK_cells_c012456_c491316, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c491316_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c491316, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)






#tissue_TNK_cells_c012456_c491316_c1231112
tissue_TNK_cells_c012456_c491316_c1231112 <- subset(tissue_TNK_cells_c012456_c491316, idents = c("1","2","3","11","12"))
tissue_TNK_cells_c012456_c491316_c1231112 <- SCTransform(tissue_TNK_cells_c012456_c491316_c1231112 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c491316_c1231112 <- RunPCA(object = tissue_TNK_cells_c012456_c491316_c1231112 , npcs = 50)

tissue_TNK_cells_c012456_c491316_c1231112  <- RunHarmony(tissue_TNK_cells_c012456_c491316_c1231112 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c491316_c1231112, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c491316_c1231112  <- RunUMAP(tissue_TNK_cells_c012456_c491316_c1231112 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c491316_c1231112  <- FindNeighbors(tissue_TNK_cells_c012456_c491316_c1231112 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c491316_c1231112  <- FindClusters(tissue_TNK_cells_c012456_c491316_c1231112, resolution =2.0)
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T) +UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#CD4 T-cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("CXCR3","NKG7","CD3D","CD4","CXCL13","CCL5","CD8A","GZMK","CCR6","TBX21","S100A4","KLRB1"))
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("CXCR3","NKG7","CD3D","CD4","SELL","TCF7","CCR7","LEF1","CCR6","KLRF1","S100A4","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))

#CD4 TEM
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3", "COTL1", "NKG7"))
UMAPPlot(tissue_TNK_cells_c012456_c491316_c1231112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c491316_c1231112, features =c("IL7R", "CCL5", "FYB1", "GZMK", "CCR7", "GZMH", "CD8A", "CD3E", "CD4", "SELL", "COTL1", "NKG7"))

# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c491316_c1231112) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c491316_c1231112) <- "RNA"
tissue_TNK_cells_c012456_c491316_c1231112 <- NormalizeData(tissue_TNK_cells_c012456_c491316_c1231112, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c491316_c1231112_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c491316_c1231112, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)








#tissue_TNK_cells_c012456_c1112
tissue_TNK_cells_c012456_c1112 <- subset(tissue_TNK_cells_c012456, idents = c("11","12"))
tissue_TNK_cells_c012456_c1112 <- SCTransform(tissue_TNK_cells_c012456_c1112 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c1112 <- RunPCA(object = tissue_TNK_cells_c012456_c1112 , npcs = 50)

tissue_TNK_cells_c012456_c1112  <- RunHarmony(tissue_TNK_cells_c012456_c1112 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c1112, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c1112  <- RunUMAP(tissue_TNK_cells_c012456_c1112 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c1112  <- FindNeighbors(tissue_TNK_cells_c012456_c1112 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c1112  <- FindClusters(tissue_TNK_cells_c012456_c1112, resolution =1.0)
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="orig.ident", label =T)+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c1112 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c1112, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c1112) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c1112) <- "RNA"
tissue_TNK_cells_c012456_c1112 <- NormalizeData(tissue_TNK_cells_c012456_c1112, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c1112_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c1112, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)







#tissue_TNK_cells_c012456_c717
tissue_TNK_cells_c012456_c717 <- subset(tissue_TNK_cells_c012456, idents = c("7","17"))
tissue_TNK_cells_c012456_c717 <- SCTransform(tissue_TNK_cells_c012456_c717 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c717 <- RunPCA(object = tissue_TNK_cells_c012456_c717 , npcs = 50)

tissue_TNK_cells_c012456_c717  <- RunHarmony(tissue_TNK_cells_c012456_c717 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c717, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c717  <- RunUMAP(tissue_TNK_cells_c012456_c717 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c717  <- FindNeighbors(tissue_TNK_cells_c012456_c717 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c717  <- FindClusters(tissue_TNK_cells_c012456_c717, resolution =1.0)
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="cell_type_fine", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("GZMH","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c717 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c717, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c717) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c717) <- "RNA"
tissue_TNK_cells_c012456_c717 <- NormalizeData(tissue_TNK_cells_c012456_c717, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c717_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c717, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





#tissue_TNK_cells_c012456_c23810
tissue_TNK_cells_c012456_c23810 <- subset(tissue_TNK_cells_c012456, idents = c("2","3","8","10"))
tissue_TNK_cells_c012456_c23810 <- SCTransform(tissue_TNK_cells_c012456_c23810 , verbose = TRUE, conserve.memory =T)
tissue_TNK_cells_c012456_c23810 <- RunPCA(object = tissue_TNK_cells_c012456_c23810 , npcs = 50)

tissue_TNK_cells_c012456_c23810  <- RunHarmony(tissue_TNK_cells_c012456_c23810 , group.by.vars = c("orig.ident"), assay.use =  "SCT")
ElbowPlot(tissue_TNK_cells_c012456_c23810, reduction = "harmony", ndims = 50)
tissue_TNK_cells_c012456_c23810  <- RunUMAP(tissue_TNK_cells_c012456_c23810 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c23810  <- FindNeighbors(tissue_TNK_cells_c012456_c23810 , reduction = "harmony", dims = 1:10)
tissue_TNK_cells_c012456_c23810  <- FindClusters(tissue_TNK_cells_c012456_c23810, resolution =3.0)
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="orig.ident", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="cell_type_manual", label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="cell_type_fine", cols = c("Central memory CD8 T cells" = "red"), label =T)#+NoLegend()
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("JCHAIN","MS4A1","CD3D","IGLL1","FCER1G","FGFBP2","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", cols = c("3"="red", "8"="blue"),  label =T)+NoLegend()

#plasma cells
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("JCHAIN","MS4A1","MZB1","IGLL1","IGHA1","IGHM","FCN1","C1QA","STMN1", "LAMP3","CCND1","LAG3"))
#B-cells
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("JCHAIN","MS4A1","CD3D","TCL1A","CD27","RGS2","IGHM","COCH","STMN1","S100A10","SELL","CCND1"))
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("RGS13","MS4A1","CD3D","LMO2","SUGCT","PRDX1","IGHM","COCH","STMN1","S100A10","EGR1","LAG3"))
#dendritic cells
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("LYZ","MS4A1","CD3D","C1QA","LILRA4","IRF8","LAMP3","CLEC9A","STMN1","CLEC10A","SLC38A1","CD1C"))
#mono-mac
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("LYZ","MS4A1","CD3D","C1QA","FCN1","S100A8","IL1B","FCGR3A","S100A9","CLEC10A","SLC38A1","CD1C"))
#T-cells
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("LYZ","MS4A1","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","SLC4A10","KLRB1"))
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("CCL5","NKG7","CD3D","CD4","CXCL13","FOXP3","CD8A","GZMK","GNLY","KLRF1","LAG3","KLRB1"))
#NK cells
UMAPPlot(tissue_TNK_cells_c012456_c23810 ,group.by="seurat_clusters", label =T)+ FeaturePlot(tissue_TNK_cells_c012456_c23810, features =c("KLRF1","NCAM1","CD3D","FGFBP2","ITGAE","KLRC3","SELL","XCL1","KLRD1","TRDC","PRF1","KLRB1"))


# calculate marker genes for each cluster
DefaultAssay(tissue_TNK_cells_c012456_c23810) <- "SCT"
DefaultAssay(tissue_TNK_cells_c012456_c23810) <- "RNA"
tissue_TNK_cells_c012456_c23810 <- NormalizeData(tissue_TNK_cells_c012456_c23810, scale.factor = 10000, verbose = TRUE)
tissue_TNK_cells_c012456_c23810_markers <- FindAllMarkers(object = tissue_TNK_cells_c012456_c23810, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)





cluster_export <- read.csv(file = "/media/data/NGS/Projects/scrnaseq-clonal2/scRNAseq/analysis/scripts/TZ_scripts/cluster_export_file_tissue.csv", sep = ",", header = TRUE)
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
merged_tissue <- merged[,c(2,3,4)]
tissue_seu <- AddMetaData(tissue_seu, metadata = merged_tissue)

UMAPPlot(tissue_seu, group.by = "cluster", raster = FALSE, label = T, repel = T)
DotPlot(tissue_seu, group.by = "cluster", features = c("SLC4A10", "CD3E","CD8A","CD4","FOXP3","SELL","CCR7", "CXCL13","GNLY","KLRF1", "MS4A1","TCL1A","CD27","VPREB1","CD38","IGHM","COCH","FCRL5", "MZB1","LYZ","C1QA","FCN1","LAMP3","CLEC9A","CDC1","CLEC10A", "STMN1","GZMK", "GZMH"), cluster.idents=T, cols = c("red", "green")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(x = "gene", y = "cell type")

setdiff(Cells(tissue_seu), merged$cell)








# calculation of effector memory score for CD4 cells
# subsetting of CD4 memory cells
CD4_memory_T <- grep("CD4_effector_memory|CD4_central_memory", unique(tissue_seu$cluster), value = TRUE)
Idents(tissue_seu) <- tissue_seu$cluster
CD4_memory_T_sobj <- subset(tissue_seu, idents = CD4_memory_T)

# CD4 effector memory gene set score calculation
CD4_EM_genes <- list(c("CCL5","LYAR","KLRB1","GZMK","IL7R","CST7","CXCR4","IL32"))
CD4_memory_T_sobj <- AddModuleScore(object=CD4_memory_T_sobj, features=CD4_EM_genes, ctrl=100, name='gene_set_score_CD4_EM')
VlnPlot(CD4_memory_T_sobj, features = "gene_set_score_CD4_EM1")+ geom_hline(yintercept = 0.52)

# to rename variables in specific rows
expression_matrix <- as.data.frame(CD4_memory_T_sobj@meta.data)
cells_above_threshold <- rownames(expression_matrix)[expression_matrix[,"gene_set_score_CD4_EM1"] > 0.52]
cells_below_threshold <- rownames(expression_matrix)[expression_matrix[,"gene_set_score_CD4_EM1"] < 0.52]
CD4_memory_T_sobj$cluster[cells_above_threshold] <- "CD4_effector_memory_T-cell"
CD4_memory_T_sobj$cluster[cells_below_threshold] <- "CD4_central_memory_T-cell"
Idents(CD4_memory_T_sobj) <- CD4_memory_T_sobj$cluster
VlnPlot(CD4_memory_T_sobj, features = "gene_set_score_CD4_EM1")+ geom_hline(yintercept = 0.52)

# to add the annotation "CD4_effector_memory_T-cell" into tissue_seu$cluster
as.data.frame(table(tissue_seu$cluster)) %>% dplyr::arrange(-Freq)
tissue_seu$cluster[cells_above_threshold] <- "CD4_effector_memory_T-cell"
tissue_seu$cluster[cells_below_threshold] <- "CD4_central_memory_T-cell"
as.data.frame(table(tissue_seu$cluster)) %>% dplyr::arrange(-Freq)


# calculating cytotoxic gene score on all T-cells
Tcell <- grep("T-cell|Treg", unique(tissue_seu$cluster), value = TRUE)
Idents(tissue_seu) <- tissue_seu$cluster
T_sobj <- subset(tissue_seu, idents = Tcell)
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

# to add the annotation "CD4_cytotoxic_T-cell" into tissue_seu$cluster
as.data.frame(table(tissue_seu$cluster)) %>% dplyr::arrange(-Freq)
tissue_seu$cluster[CD4_cells_above_threshold] <- "CD4_cytotoxic_T-cell"
as.data.frame(table(tissue_seu$cluster)) %>% dplyr::arrange(-Freq)

# to remove both temporary sobjects
rm(T_sobj)
rm(T_sobj_CD4)
rm(tissue_seu_copy)


