#R18_No.KPN is the seurat object from CRC_other_counts on Zenodo

library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(Seurat)
library(patchwork)
library(tidyverse)

DimPlot(R18_No.KPN, reduction = "umap", label = T, group.by = "seurat_clusters")

#using slingshot
R18_No.KPN.sce <- as.SingleCellExperiment(R18_No.KPN)

sce <- slingshot(R18_No.KPN.sce, clusterLabels = R18_No.KPN.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')

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
cell_colors <- cell_pal(R18_No.KPN$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(R18_No.KPN$seurat_clusters, hue_pal())


plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')

library(tradeSeq)

counts <- as.matrix(R18_No.KPN@assays$RNA@counts)
dimred <- R18_No.KPN@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = R18_No.KPN.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

dim(counts)
#sce <- fitGAM(counts = as.matrix(counts), sds = curves)#too long

filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = R18_No.KPN.sce$seurat_clusters, models = sce)

# Define function to plot
library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters =R18_No.KPN.sce$seurat_cluster, models = sce) + ggplot2::theme(legend.position = "none"), 
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
write.csv(startRes, file="startRes_R18.NO.KPN.csv")
#We can visualize the estimated smoothers for the third most significant gene.

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[273]]
plotSmoothers(sce, counts, gene = sigGeneStart)

#color the cells in UMAP space with that gene’s expression.
plotGeneCount(lineages, counts, gene = sigGeneStart)
