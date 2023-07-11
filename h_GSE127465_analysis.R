#non-small cell lung tumor and blood 7 patients
#GSE127465_human_counts_normalized_matrix
#adding meta data to h_seurat object

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

genes <- read_tsv("GSE127465_gene_names_human_41861.tsv.gz", col_names = FALSE)
cells <- read_tsv("GSE127465_human_cell_metadata_54773x25.tsv.gz", col_names = T)
cells<-as.data.frame(cells)
library(tibble)
cells <- tibble::rownames_to_column(cells, "ID")

h_GSE127465 <- ReadMtx(
  mtx = "GSE127465_human_counts_normalized_54773x41861.mtx.gz", features = "GSE127465_gene_names_human_41861.tsv.gz",
  cells = "GSE127465_human_cell_metadata_54773x25.tsv.gz", feature.column = 1, skip.cell = 1, mtx.transpose = T, cell.column = 1
)

h_GSE127465_s <- CreateSeuratObject(counts = h_GSE127465)
cell_IDs <- CellsByIdentities(h_GSE127465_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract tissue type from meta data
Tissue <- as.data.frame(cells[, c(1,3)])
library(tidyverse)
#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="ID")

#assign metadata to cells in seurat object
h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = Tissue,
  col.name = 'Tissue.idents'
)
head(x = h_GSE127465_s[[]])

#repeat for cell_type
cell_type <-as.data.frame(cells[, c(1,9)])
cell_type<-cell_type %>% remove_rownames %>% column_to_rownames(var="ID")

h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = cell_type,
  col.name = 'cell.types'
)
head(x = h_GSE127465_s[[]])

#repeat for cell_subtype
cell_subtype <-as.data.frame(cells[, c(1,11)])
cell_subtype<-cell_subtype %>% remove_rownames %>% column_to_rownames(var="ID")

h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = cell_subtype,
  col.name = 'cell.subtypes'
)
head(x = h_GSE127465_s[[]])

#Downstream analysis for clustering
h_GSE127465_s <- FindVariableFeatures(h_GSE127465_s, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(h_GSE127465_s), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(h_GSE127465_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(h_GSE127465_s)
#h_GSE127465_s <- ScaleData(h_GSE127465_s, features = all.genes) #too long
h_GSE127465_s <- ScaleData(h_GSE127465_s)

#dim reduction
h_GSE127465_s <- RunPCA(h_GSE127465_s, features = VariableFeatures(object = h_GSE127465_s))
DimPlot(h_GSE127465_s, reduction = "pca", group.by = "cell.types")

#selecting PCs
#h_GSE127465_s <- JackStraw(h_GSE127465_s, num.replicate = 100)
##h_GSE127465_s <- ScoreJackStraw(h_GSE127465_s, dims = 1:20)
###JackStrawPlot(h_GSE127465_s, dims = 1:20)

ElbowPlot(h_GSE127465_s,ndims=25)
#choose first 20 PCAs

#clustering
h_GSE127465_s <- FindNeighbors(h_GSE127465_s, dims = 1:20)
h_GSE127465_s <- FindClusters(h_GSE127465_s, resolution = 0.5)
h_GSE127465_s <- RunUMAP(h_GSE127465_s, dims = 1:20)
DimPlot(h_GSE127465_s, reduction = "umap")
DimPlot(h_GSE127465_s, reduction = "umap", group.by = "cell.types")
DimPlot(h_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")
saveRDS(h_GSE127465_s, file = "h_GSE127465_s.rds")
#set active identity to cell type
h_GSE127465_s <- SetIdent(h_GSE127465_s, value = h_GSE127465_s@meta.data$cell.types)
#subset clusters that resemble neutrophils
hNeut_GSE127465_s <- subset(x = h_GSE127465_s, idents = "Neutrophils")

#recluster neutrophils
hNeut_GSE127465_s <- FindVariableFeatures(hNeut_GSE127465_s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(hNeut_GSE127465_s), 10)
Neut_plot1 <- VariableFeaturePlot(hNeut_GSE127465_s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(hNeut_GSE127465_s)
hNeut_GSE127465_s <- ScaleData(hNeut_GSE127465_s, features = neut.all.genes)
hNeut_GSE127465_s <- RunPCA(hNeut_GSE127465_s, features = VariableFeatures(object = hNeut_GSE127465_s))
DimPlot(hNeut_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")
ElbowPlot(hNeut_GSE127465_s)

#choose firts 15 PCAs

hNeut_GSE127465_s <- FindNeighbors(hNeut_GSE127465_s, dims = 1:15)
hNeut_GSE127465_s <- FindClusters(hNeut_GSE127465_s, resolution = 0.5)
hNeut_GSE127465_s <- RunUMAP(hNeut_GSE127465_s, dims = 1:15)
DimPlot(hNeut_GSE127465_s, reduction = "umap", label = T)

DimPlot(hNeut_GSE127465_s, reduction = "umap", label = T, group.by = "Tissue.idents")

hNeut.markers <- FindAllMarkers(hNeut_GSE127465_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hNeut.markers_Top10<- hNeut.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(hNeut.markers_Top10, "hNeut.markers_Top10.csv")

#NETosis genes,
FeaturePlot(hNeut_GSE127465_s, features = c("G0S2","MME","BST1","MPO","CREB5","ENTPD4","CYP4F3","FPR2","HPSE","SELP","CEACAM3","IL17A","VNN3"), min.cutoff = "q10", max.cutoff = "q90")

library(ggplot2)
N_subset_markers <-hNeut.markers_Top10$gene
DoHeatmap(hNeut_GSE127465_s, features = N_subset_markers, size = 3)+ scale_fill_gradientn(colors = c("white", "grey", "red"))

# store the current identities in a new column of meta.data called Cluster.idents
hNeut_GSE127465_s$Cluster.idents <- Idents(hNeut_GSE127465_s)
head(x = hNeut_GSE127465_s[[]])
Idents(hNeut_GSE127465_s) <- "cell.subtypes"
#extract top 10 markers for the subtypes in the paper
hNeut.subtypes.markers <- FindAllMarkers(hNeut_GSE127465_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hNeut.subtypes.markers<- hNeut.subtypes.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(hNeut.subtypes.markers, "h_GSE127465_hNeut.subtypes.markers.csv")

Idents(hNeut_GSE127465_s) <- "Tissue.idents"

###datavisualisation
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)
library(org.Mm.eg.db)
library(clusterProfiler)

plot_stat(hNeut_GSE127465_s, plot_type = "group_count", group_by = "cell.subtypes")

markers <- FindAllMarkers(hNeut_GSE127465_s, logfc.threshold = 0.1, min.pct = 0, only.pos = T)
# Pre-filter features whose detection percentages across the two groups are similar (within
# 0.25)
markers<-FindMarkers(hNeut_GSE127465_s, ident.1 = "tumor", min.diff.pct = 0.25)
markers <- tibble::rownames_to_column(markers, "gene")

plot_heatmap(dataset = hNeut_GSE127465_s, 
             markers = markers$gene[1:40],
             sort_var = c("Tissue.idents","Cluster.idents"),
             anno_var = c("Cluster.idents","Tissue.idents"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("cyan","salmon","yellow","purple"), # color vector
                                "Reds",
                                c("blue","white","red"),                            # Three-color gradient
                                "Greens"))

plot_stat(hNeut_GSE127465_s, plot_type = "group_count", group_by="Cluster.idents")
plot_stat(hNeut_GSE127465_s, plot_type = "prop_fill", group_by = "Tissue.idents" )

#Go Analysis
hNeuts_Markers<-FindAllMarkers(hNeut_GSE127465_s)
plot_cluster_go(hNeuts_Markers, cluster_name = "tumor", org = "human", ont = "CC")
plot_cluster_go(hNeuts_Markers, cluster_name = "blood", org = "human", ont = "CC")


Idents(hNeut_GSE127465_s) <- "Cluster.idents"
plot_cluster_go(hNeut.markers, cluster_name = "2", org = "human", ont = "CC")
plot_cluster_go(hNeut.markers, cluster_name = "0", org = "human", ont = "CC")
plot_cluster_go(hNeut.markers, cluster_name = "1", org = "human", ont = "CC")
plot_cluster_go(hNeut.markers, cluster_name = "4", org = "human", ont = "CC")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

de <- find_diff_genes(dataset = hNeut_GSE127465_s, 
                      clusters = as.character(0:10),
                      comparison = c("Cluster.idents","blood","tumor"),
                      logfc.threshold = 0,   # threshold of 0 is used for GSEA
                      min.cells.group = 1)   # To include clusters with only 1 cell

gsea_res <- test_GSEA(de, 
                      pathway = pathways.hallmark)
