
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

#get exp matrix

#hGSE114727<-fread("GSE114725_rna_imputed.csv.gz")

hGSE114727<- read.delim(file = "GSE114725_rna_imputed.csv.gz", header = TRUE, sep = ",")
hGSE114727[1:4,1:4]

#subset data for patient BC1 and BC5 (containing neutrophils)
hGSE114727_subset<-hGSE114727[hGSE114727$patient %in% c('BC1','BC5'),]
hGSE114727_subset[1:4,1:4]


#create exp matrix
hGSE114727_exp<- hGSE114727_subset[,-1:-4]
metadata <- hGSE114727_subset[,1:5]

hGSE114727_exp <- t(hGSE114727_exp)
hGSE114727_exp[1:5,1:5]

colnames(hGSE114727_exp) <- hGSE114727_exp[1,]
hGSE114727_exp <- hGSE114727_exp[-1, ]

hGSE114727_s <- CreateSeuratObject(counts = hGSE114727_exp, project = "hGSE114727")
write.csv(hGSE114727_exp, file="hGSE114727_exp.csv")
#assign metadata to cells

#extract tissue type from meta data
Tissue <- as.data.frame(metadata[, c(2,5)])
library(tidyverse)
#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_s <- AddMetaData(
  object = hGSE114727_s,
  metadata = Tissue,
  col.name = 'Tissue.ident'
)
head(x = hGSE114727_s[[]])

#extract patient data from meta data
Patient <- as.data.frame(metadata[, c(1,5)])
Patient<-Patient %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_s <- AddMetaData(
  object = hGSE114727_s,
  metadata = Patient,
  col.name = 'Patient'
)
head(x = hGSE114727_s[[]])

#extract cluster annotations by authrors from meta data
Author.cluster <- as.data.frame(metadata[, c(4,5)])
Author.cluster<-Author.cluster %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_s <- AddMetaData(
  object = hGSE114727_s,
  metadata = Author.cluster,
  col.name = 'Author.cluster'
)
head(x = hGSE114727_s[[]])

###run seurat analysis
VlnPlot(hGSE114727_s, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#normalise data
hGSE114727_s <- NormalizeData(hGSE114727_s)
hGSE114727_s <- FindVariableFeatures(hGSE114727_s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hGSE114727_s), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hGSE114727_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hGSE114727_s <- ScaleData(hGSE114727_s)

hGSE114727_s <- RunPCA(hGSE114727_s, features = VariableFeatures(object = hGSE114727_s))
VizDimLoadings(hGSE114727_s, dims = 1:2, reduction = "pca")
DimPlot(hGSE114727_s, reduction = "pca")
#select PCs

ElbowPlot(hGSE114727_s,ndims=25)
#select first 15 ones

#clustering
hGSE114727_s <- FindNeighbors(hGSE114727_s, dims = 1:15)
hGSE114727_s <- FindClusters(hGSE114727_s, resolution = 0.3)
hGSE114727_s <- RunUMAP(hGSE114727_s, dims = 1:15)

DimPlot(hGSE114727_s, reduction = "umap", label = T)
DimPlot(hGSE114727_s, reduction = "umap", group.by = "Tissue.ident")
DimPlot(hGSE114727_s, reduction = "umap", group.by = "Patient")
DimPlot(hGSE114727_s, reduction = "umap", group.by = "Author.cluster", label = T)
saveRDS(hGSE114727_s, file = "hGSE114727_s.rds")

#subset clusters that resemble neutrophils/monocytes according to manuscript
Idents(object = hGSE114727_s) <- "Author.cluster"
hNeut_hGSE114727_s <- subset(x = hGSE114727_s, idents = c("29","36","40","48","68","84","86","89","94"))

#recluster neutrophils
hNeut_hGSE114727_s <- FindVariableFeatures(hNeut_hGSE114727_s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(hNeut_hGSE114727_s), 10)
Neut_plot1 <- VariableFeaturePlot(hNeut_hGSE114727_s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(hNeut_hGSE114727_s)
hNeut_hGSE114727_s <- ScaleData(hNeut_hGSE114727_s, features = neut.all.genes)
hNeut_hGSE114727_s <- RunPCA(hNeut_hGSE114727_s, features = VariableFeatures(object = hNeut_hGSE114727_s))
DimPlot(hNeut_hGSE114727_s, reduction = "pca", group.by = "Tissue.ident")
ElbowPlot(hNeut_hGSE114727_s)

#choose firts 15 PCAs

hNeut_hGSE114727_s <- FindNeighbors(hNeut_hGSE114727_s, dims = 1:15)
hNeut_hGSE114727_s <- FindClusters(hNeut_hGSE114727_s, resolution = 0.5)
hNeut_hGSE114727_s <- RunUMAP(hNeut_hGSE114727_s, dims = 1:15)
DimPlot(hNeut_hGSE114727_s, reduction = "umap", label = T)

DimPlot(hNeut_hGSE114727_s, reduction = "umap", label = F, group.by = "Tissue.ident")
DimPlot(hNeut_hGSE114727_s, reduction = "umap", label = F, group.by = "Patient")

#adding neutrophil score


library(tidyverse)
library(RColorBrewer)
library(Seurat)

Idents(hNeut_hGSE114727_s)

#assign new identities to differentiate by tissue

hNeut_hGSE114727_s <- SetIdent(hNeut_hGSE114727_s, value = hNeut_hGSE114727_s@meta.data$Tissue.ident)
DimPlot(hNeut_hGSE114727_s)

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","GNGT2","THBS1","PTMA")

#Use this list of 20 genes to score cells using the AddModuleScore function:
hNeut_hGSE114727_s <- AddModuleScore(hNeut_hGSE114727_s,
                                    features = list(t_enriched),
                                    name="T_enriched", nbin = 10)

# Plot scores
library(RColorBrewer)
FeaturePlot(hNeut_hGSE114727_s,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(hNeut_hGSE114727_s,features = c("T_enriched1"), group.by = "Tissue.ident") 

h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hNeut_hGSE114727_s <- AddModuleScore(hNeut_hGSE114727_s,
                                    features = list(h_enriched),
                                    name="H_enriched", nbin = 14)

# Plot scores
library(RColorBrewer)
FeaturePlot(hNeut_hGSE114727_s,
            features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

