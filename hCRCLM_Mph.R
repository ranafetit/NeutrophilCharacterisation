#markers for cell populations as in paper

#LYZ for myeloid cells

##MKI67+ proliferation macrophages:MKI67
##S100A8+ monocytes:S100A8, S100A9, and S100A12 
##HSP+ monocytes:HSPA1B, HSPH1, HSPA1A, and HSPB1
##cDC2 cells: CLEC10A
##CLEC9A+ cDC1 cells: CLEC9A
##MRC1+ CCL18+ M2-like macrophages: MRC1 and CCL18
##SPP1+ macrophages: SPP1
## CXCL10+ M1-like macrophages: CXCL10

library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(tidyverse)

DimPlot(hCRC_s, reduction = "umap", label = T, raster = FALSE)
FeaturePlot(hCRC_s, features = c("CD68", "SPP1"), raster = FALSE)
FeaturePlot(hCRC_s, features = c("CD68", "LYZ"), raster = FALSE)
FeaturePlot(hCRC_s, features = c("CD68", "SPP1","CXCL10"), raster = FALSE)
FeaturePlot(hCRC_s, features = c("HSPH1", "S100A8"), raster = FALSE)

FeaturePlot(hCRC_s, features = c("LYZ", "CD68","SPP1"), raster = FALSE)

#SELECT CLUSTERS 7,19,17 AS MACROPHAGES
hCRC_Mph_s <-subset(x = hCRC_s, idents = c("7","17","19"))

#recluster to analyse

hCRC_Mph_s <- FindVariableFeatures(hCRC_Mph_s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRC_Mph_s), 10)
plot1 <- VariableFeaturePlot(hCRC_Mph_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hCRC_Mph_s <- ScaleData(hCRC_Mph_s)
hCRC_Mph_s <- RunPCA(hCRC_Mph_s, features = VariableFeatures(object = hCRC_Mph_s))
VizDimLoadings(hCRC_Mph_s, dims = 1:2, reduction = "pca")
DimPlot(hCRC_Mph_s, reduction = "pca")
ElbowPlot(hCRC_Mph_s)

hCRC_Mph_s <- FindNeighbors(hCRC_Mph_s, dims = 1:20)#SELECT FIRST 20 PCS
hCRC_Mph_s <- FindClusters(hCRC_Mph_s, resolution = 0.5)
hCRC_Mph_s <- RunUMAP(hCRC_Mph_s, dims = 1:20)
DimPlot(hCRC_Mph_s, reduction = "umap", label = T)

FeaturePlot(hCRC_Mph_s, features = c("SPP1","MKI67","CXCL10","MRC1","CCL18"))

#M1_Mph_Markers
FeaturePlot(hCRC_Mph_s, features = c("CD80","CD86","NOS2","IL1A","IL6","IL1B","TLR2","TLR4"))
FeaturePlot(hCRC_Mph_s, features = c("CD80","CD40","CCR7", "CCL19","CD163"))
FeaturePlot(hCRC_Mph_s, features = c("IL1B","TNF","SOCS3","TLR2"))


#M2_Mph_Markers
FeaturePlot(hCRC_Mph_s, features = c("CSF1R","PPARG","ARG1","CD163","PDL2","CD301","MRC1"))

saveRDS(hCRC_Mph_s, file="hCRC_Mph_s.rds")
FeaturePlot(hCRC_Mph_s, features = c("MKI67"),max.cutoff = 2)

hCRC_Mph.cluster.markers <- FindAllMarkers(hCRC_Mph_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(hCRC_Mph.cluster.markers, file = "hCRC_Mph.cluster.markers.csv")

library(tidyverse)
hCRC_Mph.cluster.markers.top10<- hCRC_Mph.cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(hCRC_Mph.cluster.markers.top10, file = "hCRC_Mph.cluster.markersTop10.csv")


# Rename identity classes
##M1_Mph:3,4,6,10,14,
##M2_Mph:0,9,15,16,7,1
##Spp1_Mph:5,2,8, 
##MKi67_Mph: 11

### ??12,13,17
hCRC_Mph_s <- RenameIdents(object = hCRC_Mph_s, `14` = "M1-like-Mph")
hCRC_Mph_s <- RenameIdents(object = hCRC_Mph_s, `1` = "M2-like-Mph")
hCRC_Mph_s <- RenameIdents(object = hCRC_Mph_s, `8` = "SPP1+Mph")
hCRC_Mph_s <- RenameIdents(object = hCRC_Mph_s, `11` = "MKI67+Mph")

Idents(hCRC_Mph_s)
hCRC_Mph_s$Cell.type <- Idents(hCRC_Mph_s)

DimPlot(hCRC_Mph_s, reduction = "umap", group.by = "Cell.type", label = T)

hCRC_Mph_s1<-subset(x = hCRC_Mph_s, idents = c("12", "13", "17"), invert = TRUE)
DimPlot(hCRC_Mph_s1, reduction = "umap", group.by = "Cell.type", label = T, label.size = 5)


saveRDS(hCRC_Mph_s1, file="hCRC_Mph_s1.rds")

DimPlot(hCRC_neut_s, reduction = "umap", group.by = "seurat_clusters", label = T)
