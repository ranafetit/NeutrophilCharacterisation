library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(tidyverse)

expression_matrix <- ReadMtx(
  mtx = "matrix.mtx.gz", features = "genes.tsv.gz",
  cells = "barcodes.tsv.gz"
)
hCRC_s <- CreateSeuratObject(counts = expression_matrix, names.field = 2)

VlnPlot(hCRC_s, features = c("nFeature_RNA", "nCount_RNA"))

#run seurat workflow
hCRC_s <- NormalizeData(hCRC_s)

hCRC_s <- FindVariableFeatures(hCRC_s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRC_s), 10)
plot1 <- VariableFeaturePlot(hCRC_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hCRC_s <- ScaleData(hCRC_s)
hCRC_s <- RunPCA(hCRC_s, features = VariableFeatures(object = hCRC_s))
VizDimLoadings(hCRC_s, dims = 1:2, reduction = "pca")
DimPlot(hCRC_s, reduction = "pca")
ElbowPlot(hCRC_s)

hCRC_s <- FindNeighbors(hCRC_s, dims = 1:20)#SELECT FIRST 20 PCS
hCRC_s <- FindClusters(hCRC_s, resolution = 0.5)
hCRC_s <- RunUMAP(hCRC_s, dims = 1:20)
DimPlot(hCRC_s, reduction = "umap", label = T)

saveRDS(hCRC_s, file = "hCRC_s.rds")
FeaturePlot(hCRC_s, features = c("FCGR3B", "LYZ"), raster = FALSE)

#subsetting only neutrophils based on FCGR3b expression reported in paper
hCRC_neut_s <-subset(x = hCRC_s, idents = c("13","20"))
DimPlot(hCRC_neut_s, reduction = "umap", label = T)

#recluster to analyse

hCRC_neut_s <- FindVariableFeatures(hCRC_neut_s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRC_neut_s), 10)
plot1 <- VariableFeaturePlot(hCRC_neut_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hCRC_neut_s <- ScaleData(hCRC_neut_s)
hCRC_neut_s <- RunPCA(hCRC_neut_s, features = VariableFeatures(object = hCRC_neut_s))
VizDimLoadings(hCRC_neut_s, dims = 1:2, reduction = "pca")
DimPlot(hCRC_neut_s, reduction = "pca")
ElbowPlot(hCRC_neut_s)

hCRC_neut_s <- FindNeighbors(hCRC_neut_s, dims = 1:20)#SELECT FIRST 20 PCS
hCRC_neut_s <- FindClusters(hCRC_neut_s, resolution = 0.5)
hCRC_neut_s <- RunUMAP(hCRC_neut_s, dims = 1:20)
DimPlot(hCRC_neut_s, reduction = "umap")

hCRC.cluster.markers <- FindAllMarkers(hCRC_neut_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(tidyverse)
hCRC.cluster.markers.top10<- hCRC.cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(hCRC.cluster.markers.top10, file="hCRC.neut.top10.markers.csv")

#scoring neutrophil signature in human neutrophil dataset of CRCLM

library(tidyverse)
library(RColorBrewer)

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","GNGT2","THBS1","PTMA")
#Use this list of 20 genes to score cells using the AddModuleScore function:
hCRC_neut_s <- AddModuleScore(hCRC_neut_s,
                                    features = list(t_enriched),
                                    name="T_enriched")
FeaturePlot(hCRC_neut_s,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hCRC_neut_s <- AddModuleScore(hCRC_neut_s,
                                    features = list(h_enriched),
                                    name="H_enriched")
FeaturePlot(hCRC_neut_s,
            features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(hCRC_neut_s,features = c("T_enriched1")) 
VlnPlot(hCRC_neut_s,features = c("H_enriched1")) 

#plotting proportions of cells in each module
hCRC_neut_s@meta.data %>%
  group_by(seurat_clusters,T_enriched1) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=T_enriched1)) +
  geom_col() +
  ggtitle("Percentage  per cluster")

DotPlot(object = hCRC_neut_s, features = c("T_enriched1","H_enriched1"))
write_rds(hCRC_neut_s, file="hCRCLM_neut.rds")
#slingshot analysis

library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(Seurat)
library(patchwork)
library(tidyverse)

hCRC_neut_s.sce <- as.SingleCellExperiment(hCRC_neut_s)
#using slingshot
sce <- slingshot(hCRC_neut_s.sce, clusterLabels = hCRC_neut_s.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')

head(sce$slingPseudotime_1)
pseudo.paths <- slingPseudotime(sce)
head(pseudo.paths)

slingLineages(sce)

library(viridis)
library(scales)
library(RColorBrewer)
#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(hCRC_neut_s$seurat_clusters, hue_pal())


plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')

library(tradeSeq)

counts <- as.matrix(hCRC_neut_s@assays$RNA@counts)
dimred <- hCRC_neut_s@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = hCRC_neut_s.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

dim(counts)

filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = hCRC_neut_s.sce$seurat_clusters, models = sce)

# Define function to plot
library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = hCRC_neut_s.sce$seurat_cluster, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

pseudotime_association <- associationTest(sce)
head(pseudotime_association)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#start vs end: diff expressed genes between two points

startRes <- startVsEndTest(sce)
#We can visualize the estimated smoothers for the third most significant gene.

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[25]]
plotSmoothers(sce, counts, gene = sigGeneStart)

#color the cells in UMAP space with that gene’s expression.
plotGeneCount(lineages, counts, gene = sigGeneStart)

FeaturePlot(hCRC_neut_s, features = c("SPP1","IRF1"))
FeaturePlot(hCRC_neut_s, features = c("CXCL1","CXCL2"))
FeaturePlot(hCRC_neut_s, features = c("CXCR1","CXCR2"))
FeaturePlot(hCRC_neut_s, features = c("CXCR1","MME"), blend = T, cols=c("black","red","green"))

#NETosis genes,
FeaturePlot(hCRC_neut_s, features = c("G0S2","MME","BST1","MPO","CREB5","ENTPD4","CYP4F3","FPR2","HPSE","SELP","CEACAM3","VNN3"), min.cutoff = "q10", max.cutoff = "q90")

DotPlot(hCRC_neut_s, features = c("CXCL1","CXCL2"))
DotPlot(hCRC_neut_s, features = c("CXCR1","CXCR2"))

# Visualize co-expression of two features simultaneously
FeaturePlot(hCRC_neut_s, features = c("H_enriched1", "CXCR1"), blend = TRUE)

FeaturePlot(hCRC_neut_s, features = c("T_enriched1", "CXCL1"), blend = TRUE)
FeaturePlot(hCRC_neut_s, features = c("H_enriched1", "CXCL1"), blend = TRUE)
FeaturePlot(hCRC_neut_s, features = c("T_enriched1", "CXCL2"), blend = TRUE)
FeaturePlot(hCRC_neut_s, features = c("H_enriched1", "T_enriched1"), blend = TRUE, cols = c("black","green","red"))

FeaturePlot(hCRC_neut_s, features = c("IL1B", "IL1R2"), blend = TRUE, cols = c("black","green","red"), max.cutoff = 5)
FeaturePlot(hCRC_neut_s, features = c("CXCL8", "CXCR2"), blend = TRUE, cols = c("black","green","red"))
FeaturePlot(hCRC_neut_s, features = c("CXCL2", "CXCR2"), blend = TRUE, cols = c("black","green","red"))
FeaturePlot(hCRC_neut_s, features = c("CXCL1", "CXCR2"), blend = TRUE, cols = c("black","green","red"))

FeaturePlot(hCRC_neut_s, features = c("ICAM1", "ITGB2"), blend = TRUE, cols = c("black","green","red"), max.cutoff = 3)
FeaturePlot(hCRC_neut_s, features = c("ICAM1", "ITGAX"), blend = TRUE, cols = c("black","green","red"), max.cutoff = 3)

FeaturePlot(hCRC_neut_s, features = c("ADGRE5", "CD55"), blend = TRUE, cols = c("black","green","red"), max.cutoff = 3)

#check diff expressed genes at bifurcation point
earlyDERes <- earlyDETest(sce, knots = c(2, 3))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][231])
plotGeneCount(curves, counts, gene = rownames(earlyDERes)[oEarly][231])

#check diff genes at end
endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[70]]
plotSmoothers(sce, counts, sigGene)
plotGeneCount(curves, counts, gene = sigGene)

#BiocManager::install("dittoSeq")
library(dittoSeq)
dittoBarPlot(
  object = hCRC_neut_s,
  var = "T_enriched1", group.by = "seurat_clusters")

#data visualisation
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

library(org.Mm.eg.db)
library(clusterProfiler)
plot_stat(hCRC_neut_s, plot_type = "group_count", group_by="seurat_clusters")

#Go analysis
plot_cluster_go(hCRC.cluster.markers, cluster_name = "0", org = "human", ont = "CC")
plot_cluster_go(hCRC.cluster.markers, cluster_name = "5", org = "human", ont = "CC")
plot_cluster_go(hCRC.cluster.markers, cluster_name = "1", org = "human", ont = "CC")

#plot_all_cluster_go(hCRC.cluster.markers, org = "human", ont = "CC")



hCRC.cluster.markers.top5<- hCRC.cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

Features<- hCRC.cluster.markers.top5$gene
DoHeatmap(subset(hCRC_neut_s), features = Features, size = 3) + scale_fill_viridis()

DimPlot()

#plot volcanoplot of upregulated genes in cluster 0 (M-specific neutrophil population)

# Select Rows by column value
Cluster0_markers <- hCRC.cluster.markers[hCRC.cluster.markers$cluster == '0',]

library(EnhancedVolcano)
EnhancedVolcano(Cluster0_markers,
                lab = rownames(Cluster0_markers),
                x = 'avg_log2FC',
                y = 'p_val')

#run DEG compared to all clusters

DimPlot(hCRC_neut_s, reduction = "umap")
Idents(hCRC_neut_s)
Idents(hCRC_neut_s)<- hCRC_neut_s@meta.data$seurat_clusters

Cluster0_markers_comp <- FindMarkers(hCRC_neut_s, ident.1 = "0")
EnhancedVolcano(Cluster0_markers_comp,
                lab = rownames(Cluster0_markers_comp),
                x = 'avg_log2FC',
                y = 'p_val',col=c('black', 'black', 'black', 'red3'),FCcutoff = 1.5,
                colAlpha = 1,labSize = 4,boxedLabels = T, drawConnectors = T,selectLab = c('TXNIP','RIPOR2','FKBP5','CEBPD','STK17B','CXCR2','JAML','MME','CXCL8','CXCL2','MIF','CCL3','PI3','C15orf48' ),
                labCol = c('red','red','red','red','red','red','red','red','blue','blue','blue','blue','blue','blue' ), title = "Differentially Expressed Genes in Metastasis-specific Cluster")
