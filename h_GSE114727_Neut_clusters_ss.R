library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

#get exp matrix

hGSE114727<- read.delim(file = "GSE114725_rna_imputed.csv.gz", header = TRUE, sep = ",")
hGSE114727[1:4,1:4]

#subset data for clusters containing neutrophils
hGSE114727_clsubset<-hGSE114727[hGSE114727$cluster %in% c('89','29'),]
hGSE114727_clsubset[1:4,1:4]

#create exp matrix
hGSE114727_clexp<- hGSE114727_clsubset[,-1:-4]
metadata <- hGSE114727_clsubset[,1:5]

hGSE114727_clexp <- t(hGSE114727_clexp)
hGSE114727_clexp[1:5,1:5]

colnames(hGSE114727_clexp) <- hGSE114727_clexp[1,]
hGSE114727_clexp <- hGSE114727_clexp[-1, ]

hGSE114727_clexp[is.na(hGSE114727_clexp)] <- 0

hGSE114727_cls <- CreateSeuratObject(counts = hGSE114727_clexp, project = "hGSE114727")

#assign metadata to cells
#extract tissue type from meta data
Tissue <- as.data.frame(metadata[, c(2,5)])
library(tidyverse)
#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_cls <- AddMetaData(
  object = hGSE114727_cls,
  metadata = Tissue,
  col.name = 'Tissue.ident'
)
head(x = hGSE114727_cls[[]])

#extract patient data from meta data
Patient <- as.data.frame(metadata[, c(1,5)])
Patient<-Patient %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_cls <- AddMetaData(
  object = hGSE114727_cls,
  metadata = Patient,
  col.name = 'Patient'
)
head(x = hGSE114727_cls[[]])

#extract cluster annotations by authrors from meta data
Author.cluster <- as.data.frame(metadata[, c(4,5)])
Author.cluster<-Author.cluster %>% remove_rownames %>% column_to_rownames(var="cellid")

#assign metadata to cells in seurat object
hGSE114727_cls <- AddMetaData(
  object = hGSE114727_cls,
  metadata = Author.cluster,
  col.name = 'Author.cluster'
)
head(x = hGSE114727_cls[[]])
###run seurat analysis
VlnPlot(hGSE114727_cls, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


##
###
####Try SCTransform
install.packages("sctransform")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)




#normalise data
hGSE114727_cls <- NormalizeData(hGSE114727_cls)
hGSE114727_cls <- FindVariableFeatures(hGSE114727_cls, selection.method = "vst", nfeatures = 2000)

head(hGSE114727_cls[["RNA"]]@var.features)

top10 <- head(VariableFeatures(hGSE114727_cls), 10)

head(VariableFeatures(hGSE114727_cls))
VariableFeatures(hGSE114727_cls)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hGSE114727_cls)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

hGSE114727_cls <- ScaleData(hGSE114727_cls)

hGSE114727_cls <- RunPCA(hGSE114727_cls, features = VariableFeatures(object = hGSE114727_cls))
VizDimLoadings(hGSE114727_cls, dims = 1:2, reduction = "pca")
DimPlot(hGSE114727_cls, reduction = "pca", group.by = "Author.cluster")
#select PCs

ElbowPlot(hGSE114727_cls,ndims=25)
#select first 20 ones

#clustering
hGSE114727_cls <- FindNeighbors(hGSE114727_cls, dims = 1:20)
hGSE114727_cls <- FindClusters(hGSE114727_cls, resolution = 0.3)
hGSE114727_cls <- RunUMAP(hGSE114727_cls, dims = 1:20)

DimPlot(hGSE114727_cls, reduction = "umap", label = T)
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Tissue.ident")
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Patient")
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Author.cluster", label = T)
saveRDS(hGSE114727_cls, file = "hGSE114727_cls.rds")

#scoring for neutrophil signature 

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","GNGT2","THBS1","PTMA")

t_enriched <- c("CDKN1A","PPIA","IFITM1","IL1B","ISG15","CCL4","CD14","IER3","CCL3","IFIT3","PTMA","THBS1","IFIT1","TAGLN2","GNGT2","THBS1")


#Use this list of 20 genes to score cells using the AddModuleScore function:
hGSE114727_cls <- AddModuleScore(hGSE114727_cls,
                                     features = list(t_enriched),
                                     name="T_enriched", nbin = 5)
##

setdiff(t_enriched[[1]], rownames( hGSE114727_cls[['RNA']] ) )

#The above command can return the genes not in the RNA assay.

# Plot scores
library(RColorBrewer)
FeaturePlot(hGSE114727_cls,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(hGSE114727_cls,features = c("T_enriched1"), group.by = "Tissue.ident") 
VlnPlot(hGSE114727_cls,features = c("T_enriched1"), group.by = "seurat_clusters") 
VlnPlot(hGSE114727_cls,features = c("H_enriched1"), group.by = "seurat_clusters") 

DotPlot(hGSE114727_cls,features = c("T_enriched1"), group.by = "Tissue.ident") 


h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hGSE114727_cls <- AddModuleScore(hGSE114727_cls,
                                     features = list(h_enriched),
                                     name="H_enriched", nbin = 5)

# Plot scores
library(RColorBrewer)
FeaturePlot(hGSE114727_cls,
            features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Idents(hGSE114727_cls)
hGSE114727_cls <- SetIdent(hGSE114727_cls, value = hGSE114727_cls@meta.data$Tissue.ident)
DimPlot(hGSE114727_cls)
VlnPlot(hGSE114727_cls,features = c("H_enriched1"), group.by = "Tissue.ident") 

DotPlot(hGSE114727_cls,features = c("H_enriched1"), group.by = "Tissue.ident") 

#
##
###
####
#data visualisation
library(Scillus)
plot_stat(hGSE114727_cls, plot_type = "group_count", group_by="Tissue.ident")
plot_stat(hGSE114727_cls, plot_type = "prop_fill", group_by = "Tissue.ident" )
plot_stat(hGSE114727_cls, plot_type = "prop_multi", group_by = "Tissue.ident")

##work on normal and tumor derived neutrophils only
hGSE114727_t_h<-subset(x = hGSE114727_cls, idents = c("TUMOR", "NORMAL"))
DimPlot(hGSE114727_t_h)

hGSE114727_t_h <- FindVariableFeatures(hGSE114727_t_h, selection.method = "vst", nfeatures = 2000)
hGSE114727_t_h <- ScaleData(hGSE114727_t_h)

hGSE114727_t_h <- RunPCA(hGSE114727_t_h, features = VariableFeatures(object = hGSE114727_t_h))


ElbowPlot(hGSE114727_t_h,ndims=25)
#select first 20 ones

#clustering
hGSE114727_t_h <- FindNeighbors(hGSE114727_t_h, dims = 1:20)
hGSE114727_t_h <- FindClusters(hGSE114727_t_h, resolution = 0.3)
hGSE114727_t_h <- RunUMAP(hGSE114727_t_h, dims = 1:20)

DimPlot(hGSE114727_t_h, reduction = "umap", label = T)
DimPlot(hGSE114727_t_h, reduction = "umap", group.by = "Tissue.ident")
DimPlot(hGSE114727_t_h, reduction = "umap", group.by = "Patient")
DimPlot(hGSE114727_t_h, reduction = "umap", group.by = "Author.cluster", label = T)


#scoring for neutrophil signature

t_enriched <- c("CDKN1A","PPIA","IFITM1","IL1B","ISG15","CCL4","CD14","IER3","CCL3","IFIT3","PTMA","THBS1","IFIT1","TAGLN2","GNGT2","RPS27L")
t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","GNGT2","THBS1","PTMA")
																																							


#Use this list of 20 genes to score cells using the AddModuleScore function:
hGSE114727_t_h <- AddModuleScore(hGSE114727_t_h,
                                 features = list(t_enriched),
                                 name="T_enriched", nbin = 5)

FeaturePlot(hGSE114727_t_h,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(hGSE114727_t_h,features = c("T_enriched1"), group.by = "Tissue.ident") 
VlnPlot(hGSE114727_t_h,features = c("H_enriched1"), group.by = "Tissue.ident") 

VlnPlot(hGSE114727_t_h,features = c("T_enriched1"), group.by = "seurat_clusters") 
VlnPlot(hGSE114727_t_h,features = c("H_enriched1"), group.by = "seurat_clusters") 

DotPlot(hGSE114727_t_h,features = c("H_enriched1", "T_enriched1"), group.by = "Tissue.ident") 


h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hGSE114727_t_h <- AddModuleScore(hGSE114727_t_h,
                                 features = list(h_enriched),
                                 name="H_enriched", nbin = 5)

FeaturePlot(hGSE114727_t_h,
            features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

plot_stat(hGSE114727_t_h, plot_type = "group_count", group_by="Tissue.ident")
plot_stat(hGSE114727_t_h, plot_type = "prop_fill", group_by = "Tissue.ident" )
plot_stat(hGSE114727_t_h, plot_type = "prop_multi", group_by = "Tissue.ident")


hGSE114727_t<-subset(x = hGSE114727_cls, idents = c("TUMOR"))
DimPlot(hGSE114727_t)

write_rds(hGSE114727_t, file="hBCPT_neut.rds")
