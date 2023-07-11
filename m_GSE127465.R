#GSE127465_mouse_counts_normalized_matrix


genes <- read_tsv("GSE127465_gene_names_mouse_28205.tsv.gz", col_names = FALSE)
cells <- read_tsv("GSE127465_mouse_cell_metadata_15939x12.tsv.gz", col_names = T)
cells<-as.data.frame(cells)

library(tibble)
cells <- tibble::rownames_to_column(cells, "ID")

m_GSE127465 <- ReadMtx(
  mtx = "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz", features = "GSE127465_gene_names_mouse_28205.tsv.gz",
  cells = "GSE127465_mouse_cell_metadata_15939x12.tsv.gz", feature.column = 1, skip.cell = 1, mtx.transpose = T, cell.column = 1
)

m_GSE127465_s <- CreateSeuratObject(counts = m_GSE127465)

cell_IDs <- CellsByIdentities(m_GSE127465_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract tissue type from meta data
Tissue <- as.data.frame(cells[, c(1,2)])

library(tidyverse)
#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="ID")

#assign metadata to cells in seurat object
m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = Tissue,
  col.name = 'Tissue.idents'
)
head(x = m_GSE127465_s[[]])

#repeat for cell_type
cell_type <-as.data.frame(cells[, c(1,10)])
cell_type<-cell_type %>% remove_rownames %>% column_to_rownames(var="ID")

m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = cell_type,
  col.name = 'cell.types'
)
head(x = m_GSE127465_s[[]])

#repeat for cell_subtype
cell_subtype <-as.data.frame(cells[, c(1,11)])
cell_subtype<-cell_subtype %>% remove_rownames %>% column_to_rownames(var="ID")

m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = cell_subtype,
  col.name = 'cell.subtypes'
)
head(x = m_GSE127465_s[[]])

#Downstream analysis for clustering
m_GSE127465_s <- FindVariableFeatures(m_GSE127465_s, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(m_GSE127465_s), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m_GSE127465_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(m_GSE127465_s)
m_GSE127465_s <- ScaleData(m_GSE127465_s, features = all.genes)

#dim reduction
m_GSE127465_s <- RunPCA(m_GSE127465_s, features = VariableFeatures(object = m_GSE127465_s))
DimPlot(m_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")

#selecting PCs
m_GSE127465_s <- JackStraw(m_GSE127465_s, num.replicate = 100)
m_GSE127465_s <- ScoreJackStraw(m_GSE127465_s, dims = 1:20)
JackStrawPlot(m_GSE127465_s, dims = 1:20)

#clustering

m_GSE127465_s <- FindNeighbors(m_GSE127465_s, dims = 1:10)
m_GSE127465_s <- FindClusters(m_GSE127465_s, resolution = 0.5)
m_GSE127465_s <- RunUMAP(m_GSE127465_s, dims = 1:20)
DimPlot(m_GSE127465_s, reduction = "umap")
DimPlot(m_GSE127465_s, reduction = "umap", group.by = "cell.types")
saveRDS(m_GSE127465_s, file = "m_GSE127465_s.rds")


DimPlot(m_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")

#set active identity to cell type
m_GSE127465_s <- SetIdent(m_GSE127465_s, value = m_GSE127465_s@meta.data$cell.types)

#subset clusters that resemble neutrophils
mNeut_GSE127465_s <- subset(x = m_GSE127465_s, idents = "Neutrophils")

#recluster neutrophils
mNeut_GSE127465_s <- FindVariableFeatures(mNeut_GSE127465_s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(mNeut_GSE127465_s), 10)
Neut_plot1 <- VariableFeaturePlot(mNeut_GSE127465_s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(mNeut_GSE127465_s)
mNeut_GSE127465_s <- ScaleData(mNeut_GSE127465_s, features = neut.all.genes)
mNeut_GSE127465_s <- RunPCA(mNeut_GSE127465_s, features = VariableFeatures(object = mNeut_GSE127465_s))

DimPlot(mNeut_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")
ElbowPlot(mNeut_GSE127465_s)

#choose firts 15 PCAs

mNeut_GSE127465_s <- FindNeighbors(mNeut_GSE127465_s, dims = 1:15)
mNeut_GSE127465_s <- FindClusters(mNeut_GSE127465_s, resolution = 0.5)
mNeut_GSE127465_s <- RunUMAP(mNeut_GSE127465_s, dims = 1:15)
DimPlot(mNeut_GSE127465_s, reduction = "umap", label = T)

DimPlot(mNeut_GSE127465_s, reduction = "umap", label = T, group.by = "cell.subtypes")

mNeut.markers <- FindAllMarkers(mNeut_GSE127465_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mNeut.markers_Top10<- mNeut.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

mNeut.markers_Top5<- mNeut.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(mNeut.markers_Top10, "mNeut.markers_Top10.csv")

#plot canonical markers reported in paper

FeaturePlot(mNeut_GSE127465_s, features = c("Mmp8"), min.cutoff = "q10", max.cutoff = "q90")

#NETosis genes,
FeaturePlot(mNeut_GSE127465_s, features = c("G0s2","Mme","Bst1","Mpo","Creb5","Entpd4","Cyp4f3","Fpr2","Hpse","Selp","Ceacam3","Il17a","Vnn3"), min.cutoff = "q10", max.cutoff = "q90")

library(ggplot2)
N_subset_markers <-mNeut.markers_Top5$gene
DoHeatmap(mNeut_GSE127465_s, features = N_subset_markers, size = 3)+ scale_fill_gradientn(colors = c("white", "grey", "red"))


FeaturePlot(mNeut_GSE127465_s, features = c("Fgl2","Ptgs2","Csf3r","Gm2a"), min.cutoff = "q10", max.cutoff = "q90")


# store the current identities in a new column of meta.data called Cluster.idents
mNeut_GSE127465_s$Cluster.idents <- Idents(mNeut_GSE127465_s)
head(x = mNeut_GSE127465_s[[]])
Idents(mNeut_GSE127465_s) <- "cell.subtypes"
#extract top 10 markers for the subtypes in the paper
mNeut.subtypes.markers <- FindAllMarkers(mNeut_GSE127465_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mNeut.subtypes.markers<- mNeut.subtypes.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(mNeut.subtypes.markers, "mNeut.subtypes.markers.csv")

Idents(mNeut_GSE127465_s) <- "Cluster.idents"

