library(Seurat)


DimPlot(mNeut_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")
Idents(mNeut_GSE127465_s) <- "Tissue.idents"


DimPlot(m_GSE165276_NT_s, reduction = "umap", group.by = "Tissue.ident")
Idents(m_GSE165276_NT_s) <- "Tissue.ident"

#integrate the datasets
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = c(mNeut_GSE127465_s,m_GSE165276_NT_s))
anchors <- FindIntegrationAnchors(object.list = c(mNeut_GSE127465_s,m_GSE165276_NT_s), anchor.features = features)

NT_NSLC_int <- IntegrateData(anchorset = anchors)
DefaultAssay(NT_NSLC_int) <- "integrated"

NT_NSLC_int@meta.data$Tissue <- NT_NSLC_int@active.ident
Idents(NT_NSLC_int)

# Run the standard workflow for visualization and clustering
NT_NSLC_int <- ScaleData(NT_NSLC_int, verbose = FALSE)
NT_NSLC_int <- RunPCA(NT_NSLC_int, npcs = 30, verbose = FALSE)
NT_NSLC_int <- RunUMAP(NT_NSLC_int, reduction = "pca", dims = 1:30)
NT_NSLC_int <- FindNeighbors(NT_NSLC_int, reduction = "pca", dims = 1:30)
NT_NSLC_int <- FindClusters(NT_NSLC_int, resolution = 0.5)

DimPlot(NT_NSLC_int, reduction = "umap", label=T,group.by = "Tissue",)
DimPlot(NT_NSLC_int, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(NT_NSLC_int, reduction = "umap", label=T,split.by = "Tissue",)


NT_NSLC_int@meta.data$Tissue <- factor(NT_NSLC_int@meta.data$Tissue,
                            levels=c("BM","SP","BL","h","t"))

saveRDS(NT_NSLC_int, file="NT_NSLC_int.rds")

#finding markers per cluster

Cluster.markers <- FindAllMarkers(NT_NSLC_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(tidyverse)
Cluster.markers.top5<- Cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


write.csv(Cluster.markers, "int.cluster.markers_Top10.csv")

#finding markers per tissue type
NT_NSLC_int <- SetIdent(NT_NSLC_int, value = NT_NSLC_int@meta.data$Tissue)
Tissue.markers <- FindAllMarkers(NT_NSLC_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tissue.markers.top5<- Tissue.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(Tissue.markers, "int.tissue.markers_Top10.csv")

#data visualisation
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

library(org.Mm.eg.db)
library(clusterProfiler)

FeaturePlot(NT_NSLC_int, features = c("Ifit3","Ifit1"), min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(NT_NSLC_int, features = c("Pcna","Ier3"), min.cutoff = "q10", max.cutoff = "q90")

plot_stat(NT_NSLC_int, plot_type = "group_count", group_by="Tissue")
plot_stat(NT_NSLC_int, plot_type = "prop_fill", group_by = "Tissue" )
plot_stat(NT_NSLC_int, plot_type = "prop_multi", group_by = "Tissue")

plot_heatmap(dataset = NT_NSLC_int, 
             markers = Tissue.markers.top5,
             sort_var = c("Tissue","seurat_clusters"),
             anno_var = c("seurat_clusters","Tissue"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("cyan","salmon","blue","purple","green")))

plot_heatmap(dataset = NT_NSLC_int, 
             markers = Cluster.markers,
             sort_var = c("seurat_clusters","Tissue"),
             anno_var = c("seurat_clusters", "Tissue"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("cyan","salmon","blue","purple","green")))

DimPlot(NT_NSLC_int, features = c("Lyz1","Ptma", "Pcna", "Cxcl2","Ier3", "Bhlhe40"))

        
#Go analysis
plot_cluster_go(Tissue.markers, cluster_name = "BM", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "BL", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "SP", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "h", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "t", org = "mouse", ont = "CC")

plot_measure(dataset = NT_NSLC_int, 
             measures = c("Spp1"), 
             group_by = "Tissue")+theme_bw()
#GSEA analysis
Idents(NT_NSLC_int) <- "seurat_clusters"
de <- find_diff_genes(dataset = NT_NSLC_int, 
                      clusters = as.character(0:14),
                      comparison = c("Tissue", "h", "t"),
                      logfc.threshold = 0,   # threshold of 0 is used for GSEA
                      min.cells.group = 1)   # To include clusters with only 1 cell

gsea_res <- test_GSEA(de, 
                      pathway = pathways.hallmark)

#slingshot
library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(Seurat)
library(patchwork)

NT_NSLC_int.sce <- as.SingleCellExperiment(NT_NSLC_int)
sce <- slingshot(NT_NSLC_int.sce, clusterLabels = NT_NSLC_int.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')
head(sce$slingPseudotime_1)
pseudo.paths <- slingPseudotime(sce)
head(pseudo.paths)
slingLineages(sce)

#differential exp across pseudotime
library(tradeSeq)
counts <- as.matrix(NT_NSLC_int@assays$RNA@counts)
dimred <- NT_NSLC_int@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = NT_NSLC_int.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
dim(counts)
filt_counts <- counts[rowSums(counts > 10) > ncol(counts)/100, ]
dim(filt_counts)

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves, nknots=10)
plotGeneCount(curves, filt_counts, clusters = NT_NSLC_int$seurat_clusters, models = sce)
# Define function to plot

library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = NT_NSLC_int$seurat_clusters, models = sce) + ggplot2::theme(legend.position = "none"), 
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
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)

#color the cells in UMAP space with that gene’s expression.
plotGeneCount(lineages, counts, gene = sigGeneStart)

endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[3]]
plotSmoothers(sce, counts, sigGene)
plotGeneCount(curves, counts, gene = sigGene)

patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])
plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][1])

plotSmoothers(sce, counts, gene = "spp1")
