library(Seurat)

DimPlot(NT_NSLC_int, reduction = "umap", group.by = "Tissue")
Idents(NT_NSLC_int) <- "Tissue"

DimPlot(R18_KPN, reduction = "umap", group.by = "seurat_clusters")


#integrate the datasets
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = c(R18_KPN,NT_NSLC_int))
anchors <- FindIntegrationAnchors(object.list = c(R18_KPN,NT_NSLC_int), anchor.features = features)

KPN_NT_NSLC_int <- IntegrateData(anchorset = anchors)
DefaultAssay(KPN_NT_NSLC_int) <- "integrated"

Idents(KPN_NT_NSLC_int)
KPN_NT_NSLC_int <- RenameIdents(object = KPN_NT_NSLC_int,  'h' = 'L', 't' = 'L_AC')

KPN_NT_NSLC_int@active.ident <- factor(KPN_NT_NSLC_int@active.ident,
                                       levels=c("BM","SP","BL","L","L_AC","KPN"))

KPN_NT_NSLC_int@meta.data$Tissue.ident <- KPN_NT_NSLC_int@active.ident


# Run the standard workflow for visualization and clustering
KPN_NT_NSLC_int <- ScaleData(KPN_NT_NSLC_int, verbose = FALSE)
KPN_NT_NSLC_int <- RunPCA(KPN_NT_NSLC_int, npcs = 30, verbose = FALSE)
KPN_NT_NSLC_int <- RunUMAP(KPN_NT_NSLC_int, reduction = "pca", dims = 1:30)
KPN_NT_NSLC_int <- FindNeighbors(KPN_NT_NSLC_int, reduction = "pca", dims = 1:30)
KPN_NT_NSLC_int <- FindClusters(KPN_NT_NSLC_int, resolution = 0.5)

DimPlot(KPN_NT_NSLC_int, reduction = "umap", label=T, group.by = "Tissue.ident")


DimPlot(KPN_NT_NSLC_int, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(KPN_NT_NSLC_int, reduction = "umap", label=T,split.by = "Tissue",)

saveRDS(KPN_NT_NSLC_int, file="KPN_NT_NSLC_int.rds")

#finding markers per cluster

Cluster.markers <- FindAllMarkers(KPN_NT_NSLC_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(tidyverse)
Cluster.markers.top5<- Cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(Cluster.markers.top5, file="Cluster.markers.top5.int.csv")

#data visualisation
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

library(org.Mm.eg.db)
library(clusterProfiler)

plot_stat(KPN_NT_NSLC_int, plot_type = "group_count", group_by="Tissue.ident")
plot_stat(KPN_NT_NSLC_int, plot_type = "prop_fill", group_by = "Tissue.ident" )
plot_stat(KPN_NT_NSLC_int, plot_type = "prop_multi", group_by = "Tissue.ident")

plot_heatmap(dataset = KPN_NT_NSLC_int, 
             markers = Cluster.markers.top5,
             sort_var = c("Tissue.ident","seurat_clusters"),
             anno_var = c("seurat_clusters","Tissue.ident"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("cyan1","salmon","blue","purple","green","yellow")))


#finding markers per tissue type
KPN_NT_NSLC_int <- SetIdent(KPN_NT_NSLC_int, value = KPN_NT_NSLC_int@meta.data$Tissue.ident)
Tissue.markers <- FindAllMarkers(KPN_NT_NSLC_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tissue.markers.top5<- Tissue.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(Tissue.markers, "int.tissue.markers_Top5.csv")

plot_heatmap(dataset = KPN_NT_NSLC_int, 
             markers = Tissue.markers.top5,
             sort_var = c("Tissue.ident","seurat_clusters"),
             anno_var = c("seurat_clusters","Tissue.ident"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("cyan1","salmon","blue","purple","green","yellow")))

#Go analysis
plot_cluster_go(Tissue.markers, cluster_name = "BM", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "BL", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "SP", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "L", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "L_AC", org = "mouse", ont = "CC")
plot_cluster_go(Tissue.markers, cluster_name = "KPN", org = "mouse", ont = "CC")

KPN_NT_NSLC_int <- SetIdent(KPN_NT_NSLC_int, value = KPN_NT_NSLC_int@active.ident)
plot_cluster_go(Tissue.markers, cluster_name = "1", org = "mouse", ont = "CC")

#slingshot pseudotime analysis

library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(Seurat)
library(patchwork)

KPN_NT_NSLC_int.sce <- as.SingleCellExperiment(KPN_NT_NSLC_int)
sce <- slingshot(KPN_NT_NSLC_int.sce, clusterLabels = KPN_NT_NSLC_int.sce$seurat_clusters, start.clus = 5, stretch = 0 ,reducedDim = 'UMAP')
head(sce$slingPseudotime_1)
pseudo.paths <- slingPseudotime(sce)
head(pseudo.paths)
slingLineages(sce)

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
library(viridis)
library(scales)

cell_colors <- cell_pal(KPN_NT_NSLC_int$Tissue.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(KPN_NT_NSLC_int$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

#differential gene expression across pseudotime
library(tradeSeq)
dimred <- KPN_NT_NSLC_int@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = KPN_NT_NSLC_int.sce$seurat_clusters, dimred, start.clus= 5)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
counts <- as.matrix(KPN_NT_NSLC_int@assays$RNA@counts)
dim(counts)

sce <- fitGAM(counts = as.matrix(counts), sds = curves)#too long

filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves, nknots=10)
plotGeneCount(curves, filt_counts, clusters = KPN_NT_NSLC_int$seurat_clusters, models = sce)

# Define function to plot

library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = KPN_NT_NSLC_int$seurat_clusters, models = sce) + ggplot2::theme(legend.position = "none"), 
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
sigGeneStart <- names(sce)[oStart[4]]
plotSmoothers(sce, counts, gene = sigGeneStart)

#color the cells in UMAP space with that gene’s expression.
plotGeneCount(lineages, counts, gene = sigGeneStart)

endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[3]]
plotSmoothers(sce, counts, sigGene)
plotGeneCount(curves, counts, gene = sigGene)

#regenerative markers
reg.markers <- c("Anxa1","Lgr5", "Emp1")
DimPlot(KPN_NT_NSLC_int, features = "Anxa1", group.by = Tissue.ident)
VlnPlot(KPN_NT_NSLC_int, features = "Emp1")
DotPlot(KPN_NT_NSLC_int, features = c("Emp1"))

#DGE and GESA
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(DESeq2)

KPN_NT_NSLC_int <- SetIdent(KPN_NT_NSLC_int, value = KPN_NT_NSLC_int@meta.data$Tissue.ident)

KPN.markers <- FindMarkers(KPN_NT_NSLC_int, ident.1 = "KPN", ident.2 = NULL)
KPN.markers <- tibble::rownames_to_column(KPN.markers, "Gene") # Apply rownames_to_column

KPN.gene.list <- KPN.markers$avg_log2FC
names(KPN.gene.list) <- KPN.markers$Gene
KPN.gene.list<-na.omit(KPN.gene.list)

KPN.gene.list = sort(KPN.gene.list, decreasing = TRUE)

gse <- gseGO(geneList=KPN.gene.list, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 1, 
             maxGSSize = 20, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=30)

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

d <- godata('org.Mm.eg.db', ont="BP")
ego2 <- pairwise_termsim(gse, method="Wang", semData = d)
install.packages("ggnewscale")
emapplot(ego2, showCategory = 15)

d_CC <- godata('org.Mm.eg.db', ont="CC")
ego2_CC <- pairwise_termsim(gse, method="Wang", semData = d_CC)
emapplot(ego2_CC, showCategory = 15)

emapplot_cluster(ego2)
emapplot_cluster(ego2_CC)

cnetplot(gse, categorySize="pvalue", foldChange=KPN.gene.list)
gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 1)

#KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(KPN.gene.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Mm.eg.db")
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.KPN <- KPN.markers
KEGG.KPN$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
KEGG.KPN.gene.list <- KEGG.KPN$avg_log2FC

# Name vector with ENTREZ ids
names(KEGG.KPN.gene.list) <- KEGG.KPN$Y

# omit any NA values 
KEGG.KPN.gene.list<-na.omit(KEGG.KPN.gene.list)

# sort the list in decreasing order (required for clusterProfiler)
KEGG.KPN.gene.list = sort(KEGG.KPN.gene.list, decreasing = TRUE)

#Create gseKEGG object
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = KEGG.KPN.gene.list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 20,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, title = "Enriched Pathways")

library(pathview)
install.packages("pathview")
devtools::install_github("javadnoorb/pathview", force=TRUE)
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=KEGG.KPN.gene.list, pathway.id="mmu04620", species = "mmu")

knitr::include_graphics("mmu04620.pathview.png")
dme

#visualise specific markers
Epithelial <- c("Cdx1","Cdh1","Tff3","Hnf4a","Gata2","Cftr")
Mesenchymal <- c("Twist","Snai2","aSma","Vim","Pcolce","Fgf2")
Basal <- c("ITGA5","Ly6d","Krt5","Cav1","Cdkn1c","Krt13","Favp5","Lcn2")
Regenrative <- c("Ly6a","Anxa1","Lgr5","Clu","Emp1")
  
VlnPlot(KPN_NT_NSLC_int, features = Epithelial)
VlnPlot(KPN_NT_NSLC_int, features = Mesenchymal)
VlnPlot(KPN_NT_NSLC_int, features = Basal)

DefaultAssay(KPN_NT_NSLC_int) <- "RNA"
DotPlot(KPN_NT_NSLC_int, features = Regenrative) + RotatedAxis()

