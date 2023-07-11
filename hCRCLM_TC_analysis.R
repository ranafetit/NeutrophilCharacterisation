#isolating Tcells from h_CRCLM dataset based on markers in paper
#CD4 for CD4 T cells
##FOXP3 for Tregs
###CD8A for CD8 T cells
####SLC4A10 for MAIT cells, 

DimPlot(hCRC_s, reduction = "umap", label = T)

FeaturePlot(hCRC_s, features = c("CD4","CD8A","FOXP3","SLC4A10"))

#subsetting Tcells based on marker expression reported in paper
hCRC_CD4_TC_s <-subset(x = hCRC_s, subset = CD4 > 3)
hCRC_CD8_TC_s <-subset(x = hCRC_s, subset = CD8A > 3)
hCRC_FOXP3_TC_s <-subset(x = hCRC_s, subset = FOXP3 > 3)
hCRC_mait_TC_s <-subset(x = hCRC_s, subset = SLC4A10 > 3)

DimPlot(hCRCLM_TC, reduction = "umap", label = T)

#merge objects
hCRCLM_TC <- merge(x = hCRC_CD4_TC_s, y = list(hCRC_CD8_TC_s, hCRC_FOXP3_TC_s,hCRC_mait_TC_s))

###Repeat isolating cells with threshhold >1,5 not >3 to include all cells

#isolating Tcells from h_CRCLM dataset based on markers in paper
#CD4 for CD4 T cells
##FOXP3 for Tregs
###CD8A for CD8 T cells
####SLC4A10 for MAIT cells, 

DimPlot(hCRC_s, reduction = "umap", label = T)

FeaturePlot(hCRC_s, features = c("CD4","CD8A","FOXP3","SLC4A10"))

#subsetting Tcells based on marker expression reported in paper
hCRC_CD4_TC_s <-subset(x = hCRC_s, subset = CD4 > 1.5)
hCRC_CD8_TC_s <-subset(x = hCRC_s, subset = CD8A > 1.5)
hCRC_FOXP3_TC_s <-subset(x = hCRC_s, subset = FOXP3 > 1.5)
hCRC_mait_TC_s <-subset(x = hCRC_s, subset = SLC4A10 > 1.5)


#merge objects
hCRCLM_TC <- merge(x = hCRC_CD4_TC_s, y = list(hCRC_CD8_TC_s, hCRC_FOXP3_TC_s,hCRC_mait_TC_s))

#*****#
#recluster to analyse

hCRCLM_TC <- FindVariableFeatures(hCRCLM_TC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRCLM_TC), 10)
plot1 <- VariableFeaturePlot(hCRCLM_TC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hCRCLM_TC <- ScaleData(hCRCLM_TC)
hCRCLM_TC <- RunPCA(hCRCLM_TC, features = VariableFeatures(object = hCRCLM_TC))
VizDimLoadings(hCRCLM_TC, dims = 1:2, reduction = "pca")
DimPlot(hCRCLM_TC, reduction = "pca")
ElbowPlot(hCRCLM_TC)

hCRCLM_TC <- FindNeighbors(hCRCLM_TC, dims = 1:20)#SELECT FIRST 20 PCS
hCRCLM_TC <- FindClusters(hCRCLM_TC, resolution = 0.5)
hCRCLM_TC <- RunUMAP(hCRCLM_TC, dims = 1:20)
DimPlot(hCRCLM_TC, reduction = "umap", label = T)

FeaturePlot(hCRCLM_TC, features = c("CD4","CD8A","FOXP3","SLC4A10"))

saveRDS(hCRCLM_TC, file="hCRCLM_TC_th1.5.rds")

# Rename identity classes
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `13` = "CD4 T cell")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `17` = "Treg")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `7` = "CD8 T cell")

Idents(hCRCLM_TC)
hCRCLM_TC$Cell.type <- Idents(hCRCLM_TC)

DimPlot(hCRCLM_TC, reduction = "umap", group.by = "Cell.type", label = T, label.size = 5)

saveRDS(hCRCLM_TC, file="hCRCLM_TC.rds")
