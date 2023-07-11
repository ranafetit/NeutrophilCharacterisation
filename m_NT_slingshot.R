library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(Seurat)
library(patchwork)

m_GSE165276_NT_s<-subset(x = m_GSE165276_s, idents = c("3", "8","13","6","10","11","12","2","7","4"), invert = TRUE)
DimPlot(m_GSE165276_NT_s, reduction = "umap", label = T)

#run slingshot on NT dataset starting with two clusters
m_GSE165276_NT.sce <- as.SingleCellExperiment(m_GSE165276_NT_s)

saveRDS(m_GSE165276_NT_s, file = "m_GSE165276_NT_s.rds")

sce <- slingshot(m_GSE165276_NT.sce, clusterLabels = m_GSE165276_NT.sce$seurat_clusters, omega=TRUE , start.clus=5,reducedDim = 'UMAP')
pseudo.paths <- slingPseudotime(sce)
head(pseudo.paths)

cell_colors <- cell_pal(m_GSE165276_NT_s$Tissue.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(m_GSE165276_NT_s$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
slingLineages(sce)

library(tradeSeq)
#diff expressiona cross pseudotime

counts <- as.matrix(m_GSE165276_NT_s@assays$RNA@counts[m_GSE165276_NT_s@assays$RNA@var.features, ])


dimred <- m_GSE165276_NT_s@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = m_GSE165276_NT.sce$seurat_clusters, dimred, start.clus= 5)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = m_GSE165276_NT.sce$seurat_clusters, models = sce)

# Define function to plot
library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = m_GSE165276_NT.sce$seurat_cluster, models = sce) + ggplot2::theme(legend.position = "none"), 
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
sigGeneStart <- names(sce)[oStart[40]]
plotSmoothers(sce, counts, gene = sigGeneStart)

#color the cells in UMAP space with that gene’s expression.
plotGeneCount(lineages, counts, gene = sigGeneStart)

plotSmoothers(sce, counts, gene = "Spp1")

