library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

install.packages('Seurat')
install.Rtools(check = TRUE, check_r_update = TRUE, GUI = TRUE, ...)

#get Blood datasets
BL_1 <- read_tsv("GSM5029335_BL_dataset1.txt.gz", col_names = T)
BL_2 <- read_tsv("GSM5029338_BL_dataset2.txt.gz", col_names = T)

BL_1 <- as.data.frame(BL_1)
BL_2 <- as.data.frame(BL_2)
#change firstcolumn values to row name

BL1 <- BL_1[,-1]
rownames(BL1) <- BL_1[,1]

BL2 <- BL_2[,-1]
rownames(BL2) <- BL_2[,1]

#colnames(BL2) = as.character(c(1:ncol(BL2)))
#rownames(BL2) = as.character(c(1:nrow(BL2)))

#colnames(BL2) = seq(1, ncol(BL2))
#rownames(BL2) = seq(1, nrow(BL2))

#create a combined blood dataset from 2 runs

BL1_s <- CreateSeuratObject(counts = BL1, project = "BL1")
BL1_s

BL2<- as.matrix(BL2)
BL2_s <- CreateSeuratObject(counts = BL2, project = "BL2")

m_GSE165276_Blood_s <- merge(BL1_s, y = BL2_s, add.cell.ids = c("BL1", "BL2"), project = "Blood")

#get BM datasets
BM_1 <- read_tsv("GSM5029336_BM_dataset1.txt.gz", col_names = T)
BM_2 <- read_tsv("GSM5029339_BM_dataset2.txt.gz", col_names = T, show_col_types = FALSE)

BM_1 <- as.data.frame(BM_1)
BM_2 <- as.matrix(BM_2)

BM1 <- BM_1[,-1]
rownames(BM1) <- BM_1[,1]

BM2 <- BM_2[,-1]
rownames(BM2) <- BM_2[,1]

#create a combined blood dataset from 2 runs

BM1_s <- CreateSeuratObject(counts = BM1, project = "BM1")
BM1_s

BM2<- as.matrix(BM2)
BM2_s <- CreateSeuratObject(counts = BM2, project = "BM2")

m_GSE165276_BM_s <- merge(BM1_s, y = BM2_s, add.cell.ids = c("BM1", "BM2"), project = "BM")

#get Spleen datasets
SP_1 <- read_tsv("GSM5029337_SP_dataset1.txt.gz", col_names = T)
SP_2 <- read_tsv("GSM5029340_SP_dataset2.txt.gz", col_names = T, show_col_types = FALSE)

SP_1 <- as.data.frame(SP_1)
SP_2 <- as.data.frame(SP_2)

SP1 <- SP_1[,-1]
rownames(SP1) <- SP_1[,1]

SP2 <- SP_2[,-1]
rownames(SP2) <- SP_2[,1]

SP1_s <- CreateSeuratObject(counts = SP1, project = "SP1")
SP1_s

SP2_s <- CreateSeuratObject(counts = SP2, project = "SP2")
SP2_s

#create a combined dataset1 from blood, bm and sp
m_GSE165276_ds1_s <- merge(BM1_s, y = c(SP1_s, BL1_s), add.cell.ids = c("BM1", "SP1","BL1"), project = "NT_DS1", merge.data=T)

m_GSE165276_ds1_s

GetAssayData(m_GSE165276_ds1_s)[1:10, 1:15]
head(colnames(m_GSE165276_ds1_s))
table(m_GSE165276_ds1_s$orig.ident)


VlnPlot(m_GSE165276_ds1_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#normalise data
m_GSE165276_ds1_s <- NormalizeData(m_GSE165276_ds1_s)

#Downstream analysis for clustering

m_GSE165276_ds1_s <- FindVariableFeatures(m_GSE165276_ds1_s, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(m_GSE165276_ds1_s), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m_GSE165276_ds1_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(m_GSE165276_ds1_s)
m_GSE165276_ds1_s <- ScaleData(m_GSE165276_ds1_s, features = all.genes) #too long

#dim reduction
m_GSE165276_ds1_s <- RunPCA(m_GSE165276_ds1_s, features = VariableFeatures(object = m_GSE165276_ds1_s))
DimPlot(m_GSE165276_ds1_s, reduction = "pca")

#select PCs

ElbowPlot(m_GSE165276_ds1_s,ndims=25)
#select first 15 ones

#clustering
m_GSE165276_ds1_s <- FindNeighbors(m_GSE165276_ds1_s, dims = 1:15)
m_GSE165276_ds1_s <- FindClusters(m_GSE165276_ds1_s, resolution = 0.3)
m_GSE165276_ds1_s <- RunUMAP(m_GSE165276_ds1_s, dims = 1:15)

DimPlot(m_GSE165276_ds1_s, reduction = "umap")
DimPlot(m_GSE165276_ds1_s, reduction = "umap", group.by = "orig.ident")
DimPlot(h_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")
saveRDS(h_GSE127465_s, file = "h_GSE127465_s.rds")