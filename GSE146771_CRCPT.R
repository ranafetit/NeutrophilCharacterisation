library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)

GSE146771_metadata <- read_tsv("GSE146771_CRC.Leukocyte.10x.Metadata (1).txt.gz", col_names = T)

#Sys.setenv(VROOM_CONNECTION_SIZE=5000000)
GSE146771 <- read.table("GSE146771_CRC.Leukocyte.10x.TPM.txt.gz",sep="")

head(GSE146771[1:5,1:5])

h_GSE146771_s <- CreateSeuratObject(counts = GSE146771)

cell_IDs <- CellsByIdentities(h_GSE146771_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract Global_Cluster = cell type from meta data
Cell_type <- as.data.frame(GSE146771_metadata[, c(2,15)])

library(tidyverse)
#change column value to row names to match cell IDs in seurat
Cell_type<-Cell_type %>% remove_rownames %>% column_to_rownames(var="CellName")

#assign metadata to cells in seurat object
h_GSE146771_s <- AddMetaData(
  object = h_GSE146771_s,
  metadata = Cell_type,
  col.name = 'Cell.type'
)
head(x = h_GSE146771_s[[]])

# Set identity classes to an existing column in meta data
Idents(object = h_GSE146771_s) <- "Cell.type"
Idents(h_GSE146771_s)

#subset tcells only

hCRCPT_TC<-subset(x = h_GSE146771_s, idents = c("CD4 T cell", "CD8 T cell"))

saveRDS(hCRCPT_TC, file = "hCRCPT_TC.rds")


#Downstream analysis for clustering
hCRCPT_TC <- FindVariableFeatures(hCRCPT_TC, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hCRCPT_TC), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hCRCPT_TC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(hCRCPT_TC)
hCRCPT_TC <- ScaleData(hCRCPT_TC, features = all.genes)

#dim reduction
hCRCPT_TC <- RunPCA(hCRCPT_TC, features = VariableFeatures(object = hCRCPT_TC))
DimPlot(hCRCPT_TC, reduction = "pca", group.by = "Cell.type")

#selecting PCs
ElbowPlot(hCRCPT_TC)

#clustering

hCRCPT_TC <- FindNeighbors(hCRCPT_TC, dims = 1:20)
hCRCPT_TC <- FindClusters(hCRCPT_TC, resolution = 0.5)
hCRCPT_TC <- RunUMAP(hCRCPT_TC, dims = 1:20)
DimPlot(hCRCPT_TC, reduction = "umap")
DimPlot(hCRCPT_TC, reduction = "umap", group.by = "Cell.type")
saveRDS(hCRCPT_TC, file = "hCRCPT_TC.rds")

table(Idents(hCRCPT_TC))
