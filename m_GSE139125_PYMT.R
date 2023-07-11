library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

#get exp matrix

WT <- read_tsv("GSM4131336_Wild_Type_Expression_Matrix.txt.gz", col_names = T)
PYMT <- read_tsv("GSM4131337_PYMT_Expression_Matrix.txt.gz", col_names = T)


WT <- as.data.frame(WT)
PYMT <- as.data.frame(PYMT)

#change first column to rownames
WT1 <- WT[,-1]
rownames(WT1) <- WT[,1]

PYMT1 <- PYMT[,-1]
rownames(PYMT1) <- PYMT[,1]

WT_s <- CreateSeuratObject(counts = WT1, project = "WT")
PYMT_s <- CreateSeuratObject(counts = PYMT1, project = "PYMT")

m_GSE139125s <- merge(WT_s, y = PYMT_s, add.cell.ids = c("WT", "PYMT"), project = "GSE139125")


head(colnames(m_GSE139125s))
table(m_GSE139125s$orig.ident)
VlnPlot(m_GSE139125s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#normalise data
m_GSE139125s <- NormalizeData(m_GSE139125s)

#Downstream analysis for clustering

m_GSE139125s<- FindVariableFeatures(m_GSE139125s, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(m_GSE139125s), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m_GSE139125s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(m_GSE139125s)
m_GSE139125s <- ScaleData(m_GSE139125s, features = all.genes) #too long

#dim reduction
m_GSE139125s <- RunPCA(m_GSE139125s, features = VariableFeatures(object = m_GSE139125s))
DimPlot(m_GSE139125s, reduction = "pca")

#select PCs

ElbowPlot(m_GSE139125s,ndims=25)
#select first 15 ones

#clustering
m_GSE139125s <- FindNeighbors(m_GSE139125s, dims = 1:15)
m_GSE139125s <- FindClusters(m_GSE139125s, resolution = 0.3)
m_GSE139125s <- RunUMAP(m_GSE139125s, dims = 1:15)

DimPlot(m_GSE139125s, reduction = "umap", label = T)
DimPlot(m_GSE139125s, reduction = "umap", group.by = "orig.ident")
saveRDS(m_GSE139125s, file = "m_GSE139125s.rds")

#identify neutrophil clusters

FeaturePlot(m_GSE139125s, features = c("Ly6g","Cxcr2"))

#subset neutrophil clusters

mNeut_GSE139125s <-subset(x = m_GSE139125s, idents = c("3","0","1","7"))
DimPlot(mNeut_GSE139125s, reduction = "umap", label = T)
DimPlot(mNeut_GSE139125s, reduction = "umap", label = T, group.by = "orig.ident")


#recluster neutrophils
mNeut_GSE139125s <- FindVariableFeatures(mNeut_GSE139125s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(mNeut_GSE139125s), 10)
Neut_plot1 <- VariableFeaturePlot(mNeut_GSE139125s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(mNeut_GSE139125s)
mNeut_GSE139125s <- ScaleData(mNeut_GSE139125s, features = neut.all.genes)
mNeut_GSE139125s <- RunPCA(mNeut_GSE139125s, features = VariableFeatures(object = mNeut_GSE139125s))

DimPlot(mNeut_GSE139125s, reduction = "pca", group.by = "orig.ident")
ElbowPlot(mNeut_GSE139125s)
#choose firts 15 PCAs

mNeut_GSE139125s <- FindNeighbors(mNeut_GSE139125s, dims = 1:15)
mNeut_GSE139125s <- FindClusters(mNeut_GSE139125s, resolution = 0.5)
mNeut_GSE139125s <- RunUMAP(mNeut_GSE139125s, dims = 1:15)
DimPlot(mNeut_GSE139125s, reduction = "umap", label = T)
DimPlot(mNeut_GSE139125s, reduction = "umap", label = T, group.by = "orig.ident")


tcommon_enriched <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Ifitm6","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1")

#adding some genes
tcommon_enriched <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Ifitm6","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1")

t_enriched <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Cxcr2","Rsad2")
t_enriched_a <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Ptma")

#Use this list of 20 genes to score cells using the AddModuleScore function:
mNeut_GSE139125s <- AddModuleScore(mNeut_GSE139125s,
                                  features = list(t_enriched_a),
                                  name="T_enriched")

# Plot scores
library(RColorBrewer)
FeaturePlot(mNeut_GSE139125s,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



saveRDS(mNeut_GSE139125s, file = "mNeut_GSE139125s.rds")

FeaturePlot(mNeut_GSE139125s, features = "Ptma")

#healthy neut signature
h_enriched <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3","G0s2","Fpr2")

#Use this list of 20 genes to score cells using the AddModuleScore function:
mNeut_GSE139125s <- AddModuleScore(mNeut_GSE139125s,
                                   features = list(h_enriched),
                                   name="h_enriched")

# Plot scores
library(RColorBrewer)
FeaturePlot(mNeut_GSE139125s,
            features = "h_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#plotting proportions of cells in each module
mNeut_GSE139125s@meta.data %>%
  group_by(orig.ident,T_enriched1) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent, fill=T_enriched1)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

mNeut_GSE139125s@meta.data %>%
  group_by(seurat_clusters,T_enriched1) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=T_enriched1)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

VlnPlot(mNeut_GSE139125s,features = c("T_enriched1"), group.by = "orig.ident") 
VlnPlot(mNeut_GSE139125s,features = c("h_enriched1"), group.by = "orig.ident") 

#cellcycle scoring
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

mNeut_GSE139125s <- CellCycleScoring(mNeut_GSE139125s, s.features = s.genes, g2m.features = g2m.genes)
table(mNeut_GSE139125s[[]]$Phase)

VlnPlot(mNeut_GSE139125s,features = c("S.Score","G2M.Score"))
VlnPlot(mNeut_GSE139125s,features = c("S.Score"), group.by = "orig.ident")
VlnPlot(mNeut_GSE139125s,features = c("G2M.Score"), group.by = "orig.ident")
