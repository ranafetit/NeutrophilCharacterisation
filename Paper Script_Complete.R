#Author: Rana Fetit
#This is a structured file that entails the entire computational flow for the analysis undergone in this manuscript from loading the required packages, downloading the required public datasets, their processing and is intended to enable researchers to reproduce the entire analysis in the manuscript.

#The steps to run this file successfully are:
##Download the files as described below from Zenodo, GEOdatabase and National Omics Encyclopaedia
##Follow the processing steps described below for each section 
##Save the processed datasets under the same file names used here in each section

#The following checkpoints will enable users to resume downstream analysis after having saved their work.
##checkpoint 1 at (line 194): The user will have loaded, processed and saved the publicly available dataset GSE127465 (lungs, healthy tumor-bearing_mouse) as "mNeut_GSE127465_s.rds"
##checkpoint 2 at (line 281): The user will have loaded, processed and saved the publicly available dataset GSE165276 (Neutrotime_mouse) as "m_GSE165276_s.rds"
##checkpoint 3 at (line 381): The user will have processed and integrated the datasets: Neutrotime (m_GSE165276_NT_s), mouse lung (mNeut_GSE127465_s) and R18_KPN (CRC_KPN on Zendodo) and saved it as "KPN_NT_NSLC_int.rds"
##checkpoint 4 at (line 511): The user will have loaded, processed and saved the publicly available dataset: GSE139125 (PyMT_mouse)"mNeut_GSE139125s.rds"
##checkpoint 5 at (line 586): The user will have loaded, processed and saved the publicly available dataset: GSE127465 (human; non-small cell lung cancer) as: "h_GSE127465_s.rds"
##checkpoint 6 at (line 689): The user will have loaded, processed and saved the publicly available dataset: (GSE114727, human, breast cancer) as: "hGSE114727_cls.rds"
##checkpoint 7 at (line 878): The user will have scored Healthy_enriched and Tumour_enriched gene signatures in the mouse PyMT dataset (mNeut_GSE139125s), R18_No.KPN (CRC_other_counts on Zenodo), hGSE114727_t_h (breast cancer), as well as  Metastasis-enriched (M_enriched) gene signatures in hNeut_GSE127465_s (human NSCLC) as well as loaded, processed and saved the publicly available dataset: OEP001756 (human, Colorectal cancer liver metastasis) as "hCRC_s.rds"
##checkpoint 8 at (line 1300):The user will have subsetted neutrophils and scored for neutrophil signatures in the public dataset: OEP001756 (human, Colorectal cancer liver metastasis) and saved it as "hCRC_neut_s", performed slingshot analysis for the integrated mouse dataset: KPN_NT_NSLC_int, mNeut_GSE139125s, R18_No.KPN, hNeut_GSE127465_s,and hCRC_neut_s (Colorectal cancer liver metastasis).
##checkpoint 9 at (line 1495):The user will have subsetted neutrophils from Primary tumours (BC and NSCLC) and Metastatic CRC datasets, saved them as "hBCPT_neut.rds","hNSCLCPT_neut.rds" and "hCRCLM_neut.rds" respectively, integrated them and saved them as "T_int.rds", run GO and KEGG analysis, subsetted isolating Tcells from hCRC_s dataset, and saved it as "hCRCLM_TC.rds" 
##checkpoint 10 at (line 1604):The user will have processed publicly available dataset: GSE146771 (Primary tumor, colorectal cancer, Human) and saved it as "hCRCPT_TC.rds", integrated TCell datasets from both primary and metastatic CRC, saved it as "hCRC_TC_PTvLM1.5.rds"
##checkpoint 11 at (line 1934):The user will have run GSE and KEGG analysis on Tcells, merged Neutrophils, Tcells and macrophages (Mph) in human metastatic CRC dataset (CRCLM), saved it as "CRCLM_N_TC_Mph.rds" to run CellChat analysis.

#The expected output is: 
##Panels A-Z in Figure1
##Panels A-H, L and O-R in Figure2
##Panels A-M in Figure3
##Panels A-N in Figure4

##This should match all main Figures 1-5 and all Supplementary Figures S1-S6 in the published manuscript.
#This script has been successfully run as independent files that form the total flow which are also uploaded separately and individually described on the Github page (https://github.com/ranafetit/NeutrophilCharacterisation).


#Load required Libraries
library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(TSCAN)
library(slingshot)
library(scater)
library(ggplot2)
library(patchwork)
library(tibble)
library(Scillus)
library(magrittr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)
library(scales)
library(tradeSeq)
library(RColorBrewer)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(enrichplot)
library(DESeq2)
library(GOSemSim)
library(DOSE)
library(ReactomeGSA)
library(pathview)
library(enrichR)
library(CellChat)
library(NMF)
library(ggalluvial)

#Processing independently generated dataset: R18 (neutrophils from transplant and genetically engineered mouse models of CRC)
##Data provided from R18 research group at the beatson saved as: R18_P_N_Celltypist_Neutrophil_subset.rds
###Metadata, counts and normalised counts are deposited on Zenodo

DimPlot(R18_P_N_Celltypist_Neutrophil_subset, reduction = "umap", group.by = "Genotype")
#subsetting for KPN and other CRC models, deposited on Zenodo as: CRC_KPN and CRC_Other (AKPT, BP, BPN and KP), respectively.

R18_KPN <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "KPN")
saveRDS(R18_KPN, file="R18_KPN.rds") #deposited on Zenodo as: CRC_KPN

R18_No.KPN <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "KPN", invert=T)
saveRDS(R18_No.KPN, file="R18_No.KPN.rds") #deposited on Zenodo as: CRC_Other

#Processing publicly available dataset: GSE127465 (lungs, healthy tumor-bearing_mouse)
##Download count matrix and metadata from GEOdatabase
### save them as: GSE127465_gene_names_mouse_28205.tsv.gz and GSE127465_mouse_cell_metadata_15939x12.tsv.gz respectively to start processing

#create seurat objects from downloaded files
genes <- read_tsv("GSE127465_gene_names_mouse_28205.tsv.gz", col_names = FALSE)
cells <- read_tsv("GSE127465_mouse_cell_metadata_15939x12.tsv.gz", col_names = T)
cells<-as.data.frame(cells)

cells <- tibble::rownames_to_column(cells, "ID")

m_GSE127465 <- ReadMtx(
  mtx = "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz", features = "GSE127465_gene_names_mouse_28205.tsv.gz",
  cells = "GSE127465_mouse_cell_metadata_15939x12.tsv.gz", feature.column = 1, skip.cell = 1, mtx.transpose = T, cell.column = 1
)

m_GSE127465_s <- CreateSeuratObject(counts = m_GSE127465)

cell_IDs <- CellsByIdentities(m_GSE127465_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract tissue type from meta data
Tissue <- as.data.frame(cells[, c(1,2)])

#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="ID")

#assign metadata to cells in seurat object
m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = Tissue,
  col.name = 'Tissue.idents'
)
head(x = m_GSE127465_s[[]])

#repeat for cell_type
cell_type <-as.data.frame(cells[, c(1,10)])
cell_type<-cell_type %>% remove_rownames %>% column_to_rownames(var="ID")

m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = cell_type,
  col.name = 'cell.types'
)
head(x = m_GSE127465_s[[]])

#repeat for cell_subtype
cell_subtype <-as.data.frame(cells[, c(1,11)])
cell_subtype<-cell_subtype %>% remove_rownames %>% column_to_rownames(var="ID")

m_GSE127465_s <- AddMetaData(
  object = m_GSE127465_s,
  metadata = cell_subtype,
  col.name = 'cell.subtypes'
)
head(x = m_GSE127465_s[[]])

#Downstream analysis for clustering

m_GSE127465_s <- FindVariableFeatures(m_GSE127465_s, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(m_GSE127465_s), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m_GSE127465_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scale data
all.genes <- rownames(m_GSE127465_s)
m_GSE127465_s <- ScaleData(m_GSE127465_s, features = all.genes)
#dim reduction
m_GSE127465_s <- RunPCA(m_GSE127465_s, features = VariableFeatures(object = m_GSE127465_s))
DimPlot(m_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")
#selecting PCs
m_GSE127465_s <- JackStraw(m_GSE127465_s, num.replicate = 100)
m_GSE127465_s <- ScoreJackStraw(m_GSE127465_s, dims = 1:20)
JackStrawPlot(m_GSE127465_s, dims = 1:20)
#clustering
m_GSE127465_s <- FindNeighbors(m_GSE127465_s, dims = 1:10)
m_GSE127465_s <- FindClusters(m_GSE127465_s, resolution = 0.5)
m_GSE127465_s <- RunUMAP(m_GSE127465_s, dims = 1:20)
DimPlot(m_GSE127465_s, reduction = "umap")
DimPlot(m_GSE127465_s, reduction = "umap", group.by = "cell.types")

saveRDS(m_GSE127465_s, file = "m_GSE127465_s.rds")

DimPlot(m_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")

#set active identity to cell type
m_GSE127465_s <- SetIdent(m_GSE127465_s, value = m_GSE127465_s@meta.data$cell.types)

#subset clusters that are neutrophils according to original publication
mNeut_GSE127465_s <- subset(x = m_GSE127465_s, idents = "Neutrophils")

#recluster neutrophils
mNeut_GSE127465_s <- FindVariableFeatures(mNeut_GSE127465_s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(mNeut_GSE127465_s), 10)
Neut_plot1 <- VariableFeaturePlot(mNeut_GSE127465_s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(mNeut_GSE127465_s)
mNeut_GSE127465_s <- ScaleData(mNeut_GSE127465_s, features = neut.all.genes)
mNeut_GSE127465_s <- RunPCA(mNeut_GSE127465_s, features = VariableFeatures(object = mNeut_GSE127465_s))
DimPlot(mNeut_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")
ElbowPlot(mNeut_GSE127465_s)
#choose firts 15 PCAs
mNeut_GSE127465_s <- FindNeighbors(mNeut_GSE127465_s, dims = 1:15)
mNeut_GSE127465_s <- FindClusters(mNeut_GSE127465_s, resolution = 0.5)
mNeut_GSE127465_s <- RunUMAP(mNeut_GSE127465_s, dims = 1:15)

DimPlot(mNeut_GSE127465_s, reduction = "umap", label = T)
DimPlot(mNeut_GSE127465_s, reduction = "umap", label = T, group.by = "cell.subtypes")

# store the current identities in a new column of meta.data called Cluster.idents
mNeut_GSE127465_s$Cluster.idents <- Idents(mNeut_GSE127465_s)
head(x = mNeut_GSE127465_s[[]])

saveRDS(mNeut_GSE127465_s, file = "mNeut_GSE127465_s.rds")

#Processing publicly available dataset: GSE165276 (Neutrotime_mouse)
##Download the datasets from GEOdatabase
###save them as: GSM5029335_BL_dataset1.txt.gz, GSM5029338_BL_dataset2.txt.gz, GSM5029336_BM_dataset1.txt.gz, GSM5029336_BM_dataset2.txt.gz, GSM5029337_SP_dataset1.txt.gz, GSM5029337_SP_dataset2.txt.gz)

#create seurat objects from downloaded files
#Bonemarrow (BM)
BM_1 = read.table("GSM5029339_BM_dataset1.txt.gz")
BM1_s <- CreateSeuratObject(counts = BM_1, project = "BM1")

BM_2 = read.table("GSM5029339_BM_dataset2.txt.gz")
BM2_s <- CreateSeuratObject(counts = BM_2, project = "BM2")

#Spleen (SP)
SP_1 = read.table("GSM5029340_SP_dataset1.txt.gz")
SP1_s <- CreateSeuratObject(counts = SP_1, project = "SP1")

SP_2 = read.table("GSM5029340_SP_dataset2.txt.gz")
SP2_s <- CreateSeuratObject(counts = SP_2, project = "SP2")

#Blood (BL)
BL_1 = read.table("GSM5029338_BL_dataset1.txt.gz")
BL1_s <- CreateSeuratObject(counts = BL_1, project = "BL1")

BL_2 = read.table("GSM5029338_BL_dataset2.txt.gz")
BL2_s <- CreateSeuratObject(counts = BL_2, project = "BL2")

#Merge datasets from the two runs per tissue type
m_GSE165276_BL_s <- merge(BL1_s, y = BL2_s, add.cell.ids = c("BL1","BL2"), project = "Blood", merge.data=T)

m_GSE165276_BM_s <- merge(BM1_s, y = BM2_s, add.cell.ids = c("BM1","BM2"), project = "BM", merge.data=T)

m_GSE165276_SP_s <-  merge(SP1_s, y = SP2_s, add.cell.ids = c("SP1","SP2"), project = "SP", merge.data=T)

#create a complete neutrotime dataset and save it as (m_GSE165276_s)
m_GSE165276_s <- merge(m_GSE165276_BL_s, y = c(m_GSE165276_BM_s,m_GSE165276_SP_s), add.cell.ids = c("BL","BM","SP"), project = "NT", merge.data=T)

# Rename cells in a Seurat object according to tissue type
head(x = colnames(x = m_GSE165276_s))
Idents(object = m_GSE165276_s) 

m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `BL1` = "BL")
m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `BL2` = "BL")
m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `BM1` = "BM")
m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `BM2` = "BM")
m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `SP1` = "SP")
m_GSE165276_s <- RenameIdents(object = m_GSE165276_s, `SP2` = "SP")

Idents(object = m_GSE165276_s) 
m_GSE165276_s[["Tissue.ident"]] <- Idents(object = m_GSE165276_s)

##Process Dataset as described in original publication.
VlnPlot(m_GSE165276_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#normalise data
m_GSE165276_s <- NormalizeData(m_GSE165276_s)
#Downstream analysis for clustering
m_GSE165276_s <- FindVariableFeatures(m_GSE165276_s, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(m_GSE165276_s), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m_GSE165276_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scale data
all.genes <- rownames(m_GSE165276_s)
m_GSE165276_s <- ScaleData(m_GSE165276_s, features = all.genes)
#dim reduction
m_GSE165276_s <- RunPCA(m_GSE165276_s, features = VariableFeatures(object = m_GSE165276_s))
DimPlot(m_GSE165276_s, reduction = "pca")
#select PCs
ElbowPlot(m_GSE165276_s,ndims=25)
#select first 15 ones for clustering

#clustering
m_GSE165276_s <- FindNeighbors(m_GSE165276_s, dims = 1:15)
m_GSE165276_s <- FindClusters(m_GSE165276_s, resolution = 0.5)
m_GSE165276_s <- RunUMAP(m_GSE165276_s, dims = 1:15)

DimPlot(m_GSE165276_s, reduction = "umap")
DimPlot(m_GSE165276_s, reduction = "umap", group.by = "orig.ident")
DimPlot(m_GSE165276_s, reduction = "umap", group.by = "Tissue.ident")

m_GSE165276_s[["Cluster.ident"]] <- Idents(object = m_GSE165276_s)

DimPlot(m_GSE165276_s, reduction = "umap", group.by = "Cluster.ident", label=T)

saveRDS(m_GSE165276_s, file = "m_GSE165276_s.rds")

#Subset for the major BM, BL and SP clusters 
m_GSE165276_NT_s<-subset(x = m_GSE165276_s, idents = c("3", "8","13","6","10","11","12","2","7","4"), invert = TRUE)
DimPlot(m_GSE165276_NT_s, reduction = "umap", label = T)

saveRDS(m_GSE165276_NT_s, file = "m_GSE165276_NT_s.rds")#This is the dataset that will be used going forward.

###Integrating Neutrotime (m_GSE165276_NT_s) and mouse lung (mNeut_GSE127465_s) datasets

DimPlot(mNeut_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")
Idents(mNeut_GSE127465_s) <- "Tissue.idents"

DimPlot(m_GSE165276_NT_s, reduction = "umap", group.by = "Tissue.ident")
Idents(m_GSE165276_NT_s) <- "Tissue.ident"

#select features that are repeatedly variable across datasets for integration
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
                                       levels=c("BM","SP","BL","h","t")) #h= healthy lung tissue, t= tunour lung tissue
saveRDS(NT_NSLC_int, file="NT_NSLC_int.rds")

###Subset only for healthy neutrophils in both datasets
NT_NSLC_int_healthy <- subset(x = NT_NSLC_int, idents = "t", invert=T)
DimPlot(NT_NSLC_int_healthy)

#recluster for neutrophils
NT_NSLC_int_healthy <- FindVariableFeatures(NT_NSLC_int_healthy, selection.method = "vst", nfeatures = 2000)
NT_NSLC_int_healthy.all.genes <- rownames(NT_NSLC_int_healthy)
NT_NSLC_int_healthy <- ScaleData(NT_NSLC_int_healthy, features = NT_NSLC_int_healthy.all.genes)
NT_NSLC_int_healthy <- RunPCA(NT_NSLC_int_healthy, features = VariableFeatures(object = NT_NSLC_int_healthy))

DimPlot(NT_NSLC_int_healthy, reduction = "pca")
ElbowPlot(NT_NSLC_int_healthy)
#choose firts 15 PCAs

NT_NSLC_int_healthy <- FindNeighbors(NT_NSLC_int_healthy, dims = 1:15)
NT_NSLC_int_healthy <- FindClusters(NT_NSLC_int_healthy, resolution = 0.5)
NT_NSLC_int_healthy <- RunUMAP(NT_NSLC_int_healthy, dims = 1:15)

DimPlot(NT_NSLC_int_healthy, reduction = "umap", label = T)
DimPlot(NT_NSLC_int_healthy, reduction = "umap", label = T, group.by = "Tissue")

Idents(NT_NSLC_int_healthy)<-NT_NSLC_int_healthy@meta.data$Tissue
Idents(NT_NSLC_int_healthy)
NT_NSLC_int_healthy <- RenameIdents(object = NT_NSLC_int_healthy,  'h' = 'L')

Figure_1A<- DimPlot(NT_NSLC_int_healthy, reduction = "umap", label = T)
Figure_1A+ scale_x_reverse()

###Integrate the datasets: NT_NSLC_int and R18_KPN
DimPlot(NT_NSLC_int, reduction = "umap", group.by = "Tissue")
Idents(NT_NSLC_int) <- "Tissue"

DimPlot(R18_KPN, reduction = "umap", group.by = "seurat_clusters")

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = c(R18_KPN,NT_NSLC_int))
anchors <- FindIntegrationAnchors(object.list = c(R18_KPN,NT_NSLC_int), anchor.features = features)

KPN_NT_NSLC_int <- IntegrateData(anchorset = anchors)
DefaultAssay(KPN_NT_NSLC_int) <- "integrated"

Idents(KPN_NT_NSLC_int)
KPN_NT_NSLC_int <- RenameIdents(object = KPN_NT_NSLC_int,  'h' = 'L', 't' = 'L_AC') #L=Lung tissue, L_AC=Lung Adenocarcinoma

KPN_NT_NSLC_int@active.ident <- factor(KPN_NT_NSLC_int@active.ident,
                                       levels=c("BM","SP","BL","L","L_AC","KPN"))

KPN_NT_NSLC_int@meta.data$Tissue.ident <- KPN_NT_NSLC_int@active.ident

# Run the standard workflow for visualization and clustering
KPN_NT_NSLC_int <- ScaleData(KPN_NT_NSLC_int, verbose = FALSE)
KPN_NT_NSLC_int <- RunPCA(KPN_NT_NSLC_int, npcs = 30, verbose = FALSE)
KPN_NT_NSLC_int <- RunUMAP(KPN_NT_NSLC_int, reduction = "pca", dims = 1:30)
KPN_NT_NSLC_int <- FindNeighbors(KPN_NT_NSLC_int, reduction = "pca", dims = 1:30)
KPN_NT_NSLC_int <- FindClusters(KPN_NT_NSLC_int, resolution = 0.5)

Figure_1B_1 <- DimPlot(KPN_NT_NSLC_int, reduction = "umap", label=T, group.by = "Tissue.ident")
Figure_1B_2 <- DimPlot(KPN_NT_NSLC_int, reduction = "umap", label = TRUE, repel = TRUE)

saveRDS(KPN_NT_NSLC_int, file="KPN_NT_NSLC_int.rds")

#finding cluster markers
Cluster.markers <- FindAllMarkers(KPN_NT_NSLC_int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster.markers.top5<- Cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(Cluster.markers.top5, file="Cluster.markers.top5.int.csv")

#data visualization
supplementary_Figure_S1A <- plot_stat(KPN_NT_NSLC_int, plot_type = "group_count", group_by="Tissue.ident")
supplementary_Figure_S1B <-plot_stat(KPN_NT_NSLC_int, plot_type = "prop_multi", group_by = "Tissue.ident")

supplementary_Figure_S2A <- plot_heatmap(dataset = KPN_NT_NSLC_int, 
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

#assign new identities to differentiate between healthy neutrophils and tumour neutrophils.
head(x=KPN_NT_NSLC_int@meta.data[["Tissue.ident"]])
KPN_NT_NSLC_int@meta.data$T.h.ident <- KPN_NT_NSLC_int@meta.data[["Tissue.ident"]]
KPN_NT_NSLC_int <- SetIdent(KPN_NT_NSLC_int, value = KPN_NT_NSLC_int@meta.data$T.h.ident)
KPN_NT_NSLC_int <- RenameIdents(object = KPN_NT_NSLC_int,  'L_AC' = 'Tumor', 'KPN' = 'Tumor')

KPN_NT_NSLC_int@active.ident <- factor(KPN_NT_NSLC_int@active.ident,
                                       levels=c("BM","SP","BL","L","Tumor"))
DimPlot(KPN_NT_NSLC_int) #now we have one cluster for all tumor-derived neutrophils
DefaultAssay(KPN_NT_NSLC_int) <- "RNA"

#Get top 20 genes enriched in tumor-derived Neutrophils to be used in identifying neutrophil signatures
t_enriched <- FindMarkers(KPN_NT_NSLC_int, ident.1 = "Tumor", verbose = FALSE) %>%
  arrange(-avg_log2FC) %>%
  rownames_to_column(var = "gene") %>%
  pull(gene) %>% 
  .[1:20]

#Processing publicly available dataset: GSE139125 (PyMT_mouse)
##Download the datasets from GEOdatabase
###save them as: GSM4131336_Wild_Type_Expression_Matrix.txt.gz and GSM4131337_PYMT_Expression_Matrix.txt.gz

#create Seurat objects from downloaded files
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

##merge the two datasets into a complete one and save it as (m_GSE139125s)
m_GSE139125s <- merge(WT_s, y = PYMT_s, add.cell.ids = c("WT", "PYMT"), project = "GSE139125")

#Downstream analysis for clustering
head(colnames(m_GSE139125s))
table(m_GSE139125s$orig.ident)
VlnPlot(m_GSE139125s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
m_GSE139125s <- NormalizeData(m_GSE139125s)
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
#dim reduction and PC selection
m_GSE139125s <- RunPCA(m_GSE139125s, features = VariableFeatures(object = m_GSE139125s))
DimPlot(m_GSE139125s, reduction = "pca")
ElbowPlot(m_GSE139125s,ndims=25) #select first 15 ones

#clustering
m_GSE139125s <- FindNeighbors(m_GSE139125s, dims = 1:15)
m_GSE139125s <- FindClusters(m_GSE139125s, resolution = 0.3)
m_GSE139125s <- RunUMAP(m_GSE139125s, dims = 1:15)

DimPlot(m_GSE139125s, reduction = "umap", label = T)
DimPlot(m_GSE139125s, reduction = "umap", group.by = "orig.ident")
saveRDS(m_GSE139125s, file = "m_GSE139125s.rds")

#identifying and subsetting neutrophil clusters as described in original publications
FeaturePlot(m_GSE139125s, features = c("Ly6g","Cxcr2"))
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
ElbowPlot(mNeut_GSE139125s) #choose firts 15 PCAs

mNeut_GSE139125s <- FindNeighbors(mNeut_GSE139125s, dims = 1:15)
mNeut_GSE139125s <- FindClusters(mNeut_GSE139125s, resolution = 0.5)
mNeut_GSE139125s <- RunUMAP(mNeut_GSE139125s, dims = 1:15)
DimPlot(mNeut_GSE139125s, reduction = "umap", label = T)

Figure_1C <- DimPlot(mNeut_GSE139125s, reduction = "umap", label = T, group.by = "orig.ident")
saveRDS(mNeut_GSE139125s, file = "mNeut_GSE139125s.rds") #This will be used onwards for analysis of mouse PyMT dataset.

#Processing publicly available dataset: GSE127465 (human); non-small cell lung tumor and blood 7 patients
##Download the datasets from GEOdatabase
###save them as: GSE127465_gene_names_human_41861.tsv.gz and GSE127465_human_cell_metadata_54773x25.tsv.gz
genes <- read_tsv("GSE127465_gene_names_human_41861.tsv.gz", col_names = FALSE)
cells <- read_tsv("GSE127465_human_cell_metadata_54773x25.tsv.gz", col_names = T)
cells<-as.data.frame(cells)
cells <- tibble::rownames_to_column(cells, "ID")

h_GSE127465 <- ReadMtx(
  mtx = "GSE127465_human_counts_normalized_54773x41861.mtx.gz", features = "GSE127465_gene_names_human_41861.tsv.gz",
  cells = "GSE127465_human_cell_metadata_54773x25.tsv.gz", feature.column = 1, skip.cell = 1, mtx.transpose = T, cell.column = 1
)

h_GSE127465_s <- CreateSeuratObject(counts = h_GSE127465)
cell_IDs <- CellsByIdentities(h_GSE127465_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract tissue type from meta data
Tissue <- as.data.frame(cells[, c(1,3)])
#change column value to row names to match cell IDs in seurat
Tissue<-Tissue %>% remove_rownames %>% column_to_rownames(var="ID")
#assign metadata to cells in seurat object
h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = Tissue,
  col.name = 'Tissue.idents'
)
head(x = h_GSE127465_s[[]])

#repeat for cell_type
cell_type <-as.data.frame(cells[, c(1,9)])
cell_type<-cell_type %>% remove_rownames %>% column_to_rownames(var="ID")

h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = cell_type,
  col.name = 'cell.types'
)
head(x = h_GSE127465_s[[]])

#repeat for cell_subtype
cell_subtype <-as.data.frame(cells[, c(1,11)])
cell_subtype<-cell_subtype %>% remove_rownames %>% column_to_rownames(var="ID")

h_GSE127465_s <- AddMetaData(
  object = h_GSE127465_s,
  metadata = cell_subtype,
  col.name = 'cell.subtypes'
)
head(x = h_GSE127465_s[[]])

#Downstream analysis for clustering
h_GSE127465_s <- FindVariableFeatures(h_GSE127465_s, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(h_GSE127465_s), 10)
plot1 <- VariableFeaturePlot(h_GSE127465_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scale data and perform dim reduction
all.genes <- rownames(h_GSE127465_s)
h_GSE127465_s <- ScaleData(h_GSE127465_s)
h_GSE127465_s <- RunPCA(h_GSE127465_s, features = VariableFeatures(object = h_GSE127465_s))
DimPlot(h_GSE127465_s, reduction = "pca", group.by = "cell.types")
ElbowPlot(h_GSE127465_s,ndims=25)#choose first 20 PCAs
#clustering
h_GSE127465_s <- FindNeighbors(h_GSE127465_s, dims = 1:20)
h_GSE127465_s <- FindClusters(h_GSE127465_s, resolution = 0.5)
h_GSE127465_s <- RunUMAP(h_GSE127465_s, dims = 1:20)
DimPlot(h_GSE127465_s, reduction = "umap")
DimPlot(h_GSE127465_s, reduction = "umap", group.by = "cell.types")
DimPlot(h_GSE127465_s, reduction = "umap", group.by = "Tissue.idents")
saveRDS(h_GSE127465_s, file = "h_GSE127465_s.rds")

#set active identity to cell type
h_GSE127465_s <- SetIdent(h_GSE127465_s, value = h_GSE127465_s@meta.data$cell.types)

#subset clusters that are neutrophils as identified in original publication
hNeut_GSE127465_s <- subset(x = h_GSE127465_s, idents = "Neutrophils")

#recluster neutrophils
hNeut_GSE127465_s <- FindVariableFeatures(hNeut_GSE127465_s, selection.method = "vst", nfeatures = 2000)
Neut_top10 <- head(VariableFeatures(hNeut_GSE127465_s), 10)
Neut_plot1 <- VariableFeaturePlot(hNeut_GSE127465_s)
Neut_plot2 <- LabelPoints(plot = Neut_plot1, points = Neut_top10, repel = TRUE)
neut.all.genes <- rownames(hNeut_GSE127465_s)
hNeut_GSE127465_s <- ScaleData(hNeut_GSE127465_s, features = neut.all.genes)
hNeut_GSE127465_s <- RunPCA(hNeut_GSE127465_s, features = VariableFeatures(object = hNeut_GSE127465_s))
DimPlot(hNeut_GSE127465_s, reduction = "pca", group.by = "Tissue.idents")
ElbowPlot(hNeut_GSE127465_s)#choose firts 15 PCAs
hNeut_GSE127465_s <- FindNeighbors(hNeut_GSE127465_s, dims = 1:15)
hNeut_GSE127465_s <- FindClusters(hNeut_GSE127465_s, resolution = 0.5)
hNeut_GSE127465_s <- RunUMAP(hNeut_GSE127465_s, dims = 1:15)
DimPlot(hNeut_GSE127465_s, reduction = "umap", label = T)

Figure_1I <- DimPlot(hNeut_GSE127465_s, reduction = "umap", label = T, group.by = "Tissue.idents")

#Processing publicly available dataset: (GSE114727, human, breast cancer) 
##Download the datasets from GEOdatabase
###save expression matrix as: GSE114725_rna_imputed.csv.gz
hGSE114727<- read.delim(file = "GSE114725_rna_imputed.csv.gz", header = TRUE, sep = ",")
hGSE114727[1:4,1:4]

#subset data for clusters containing neutrophils as identified in original publication
hGSE114727_clsubset<-hGSE114727[hGSE114727$cluster %in% c('89','29'),]
hGSE114727_clsubset[1:4,1:4]
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

#run standard seurat analysis
VlnPlot(hGSE114727_cls, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hGSE114727_cls <- NormalizeData(hGSE114727_cls)
hGSE114727_cls <- FindVariableFeatures(hGSE114727_cls, selection.method = "vst", nfeatures = 2000)
head(hGSE114727_cls[["RNA"]]@var.features)
top10 <- head(VariableFeatures(hGSE114727_cls), 10)
head(VariableFeatures(hGSE114727_cls))
VariableFeatures(hGSE114727_cls)
plot1 <- VariableFeaturePlot(hGSE114727_cls)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

hGSE114727_cls <- ScaleData(hGSE114727_cls)
hGSE114727_cls <- RunPCA(hGSE114727_cls, features = VariableFeatures(object = hGSE114727_cls))
VizDimLoadings(hGSE114727_cls, dims = 1:2, reduction = "pca")
DimPlot(hGSE114727_cls, reduction = "pca", group.by = "Author.cluster")
#select PCs
ElbowPlot(hGSE114727_cls,ndims=25) #select first 20 ones
#clustering
hGSE114727_cls <- FindNeighbors(hGSE114727_cls, dims = 1:20)
hGSE114727_cls <- FindClusters(hGSE114727_cls, resolution = 0.3)
hGSE114727_cls <- RunUMAP(hGSE114727_cls, dims = 1:20)
DimPlot(hGSE114727_cls, reduction = "umap", label = T)
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Tissue.ident")
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Patient")
DimPlot(hGSE114727_cls, reduction = "umap", group.by = "Author.cluster", label = T)
saveRDS(hGSE114727_cls, file = "hGSE114727_cls.rds")

#hGSE114727_cls contains the 2 neutrophil clusters (29 and 89) from patients BC1-8 in Normal, tumour tissue as well as lymphnode and blood.
##We are interested in normal and tumor derived neutrophils only
hGSE114727_t_h<-subset(x = hGSE114727_cls, idents = c("TUMOR", "NORMAL"))
DimPlot(hGSE114727_t_h)
hGSE114727_t_h <- FindVariableFeatures(hGSE114727_t_h, selection.method = "vst", nfeatures = 2000)
hGSE114727_t_h <- ScaleData(hGSE114727_t_h)
hGSE114727_t_h <- RunPCA(hGSE114727_t_h, features = VariableFeatures(object = hGSE114727_t_h))
ElbowPlot(hGSE114727_t_h,ndims=25)#select first 20 ones
#clustering
hGSE114727_t_h <- FindNeighbors(hGSE114727_t_h, dims = 1:20)
hGSE114727_t_h <- FindClusters(hGSE114727_t_h, resolution = 0.3)
hGSE114727_t_h <- RunUMAP(hGSE114727_t_h, dims = 1:20)

Figure_1L <- DimPlot(hGSE114727_t_h, reduction = "umap", group.by = "Patient")

#scoring neutrophil signatures derived from KPN and L_AC on other models
##Scoring for Healthy_enriched and Tumour_enriched gene signatures in the mouse PyMT dataset
###Use the list of identified genes to score cells using the AddModuleScore function
h_enriched <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3","G0s2","Fpr2")
mNeut_GSE139125s <- AddModuleScore(mNeut_GSE139125s,
                                   features = list(h_enriched),
                                   name="h_enriched")

Figure_1D <- FeaturePlot(mNeut_GSE139125s,
            features = "h_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

t_enriched_a <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Ptma")
mNeut_GSE139125s <- AddModuleScore(mNeut_GSE139125s,
                                   features = list(t_enriched_a),
                                   name="T_enriched")

Figure_1E <- FeaturePlot(mNeut_GSE139125s,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Scoring for Healthy_enriched and Tumour_enriched gene signatures in the R18_No.KPN dataset
##R18_No.KPN is the seurat object from CRC_other_counts on Zenodo
###Use the list of identified genes to score cells using the AddModuleScore function
Figure_1F <- DimPlot(R18_No.KPN)

h_enriched_a <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3","G0s2","Fpr2")
R18_No.KPN <- AddModuleScore(R18_No.KPN,
                             features = list(h_enriched_a),
                             name="H_enriched", ctrl=50)
Figure_1G <- FeaturePlot(R18_No.KPN,
                        features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

t_enriched_a <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Ptma")
R18_No.KPN <- AddModuleScore(R18_No.KPN,
                             features = list(t_enriched_a),
                             name="T_enriched", ctrl=50)
Figure_1H <- FeaturePlot(R18_No.KPN,
            features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Scoring for Healthy_enriched and Tumour_enriched gene signatures in other mouse CRC models
R18_AKPT <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "AKPT")
R18_KP <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "KP")
R18_BPN <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "BPN")
R18_BP <- subset(x = R18_P_N_Celltypist_Neutrophil_subset, idents = "BP")

h_enriched_a <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3","G0s2","Fpr2")
t_enriched_a <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Ptma")

R18_AKPT <- AddModuleScore(R18_AKPT,
                           features = list(h_enriched_a),
                           name="H_enriched", ctrl=50)
R18_KP <- AddModuleScore(R18_KP,
                         features = list(h_enriched_a),
                         name="H_enriched", ctrl=50)
R18_BP <- AddModuleScore(R18_BP,
                         features = list(h_enriched_a),
                         name="H_enriched", ctrl=50)
R18_BPN <- AddModuleScore(R18_BPN,
                          features = list(h_enriched_a),
                          name="H_enriched", ctrl=50)

Supplementary_Figure_S2B1 <-FeaturePlot(R18_AKPT,
                                        features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2B2 <-FeaturePlot(R18_KP,
                                        features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2B3 <-FeaturePlot(R18_BP,
                                        features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2B4 <-FeaturePlot(R18_BPN,
                                        features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
R18_AKPT <- AddModuleScore(R18_AKPT,
                           features = list(t_enriched_a),
                           name="T_enriched", ctrl=50)
R18_KP <- AddModuleScore(R18_KP,
                         features = list(t_enriched_a),
                         name="T_enriched", ctrl=50)
R18_BP <- AddModuleScore(R18_BP,
                         features = list(t_enriched_a),
                         name="T_enriched", ctrl=50)
R18_BPN <- AddModuleScore(R18_BPN,
                          features = list(t_enriched_a),
                          name="T_enriched", ctrl=50)

Supplementary_Figure_S2C1 <- FeaturePlot(R18_AKPT,
                                         features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2C2 <- FeaturePlot(R18_KP,
                                         features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2C3 <- FeaturePlot(R18_BP,
                                         features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Supplementary_Figure_S2C4 <- FeaturePlot(R18_BPN,
                                         features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Scoring for Healthy_enriched, Tumour_enriched and Metastasis-enriched (M_enriched) gene signatures in hNeut_GSE127465_s (human NSCLC)
##Use the list of identified genes to score cells using the AddModuleScore function

#assign new identities to differentiate between neutrophils from blood and tumor.
hNeut_GSE127465_s <- SetIdent(hNeut_GSE127465_s, value = hNeut_GSE127465_s@meta.data$Tissue.idents)

h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hNeut_GSE127465_s <- AddModuleScore(hNeut_GSE127465_s,
                                    features = list(h_enriched),
                                    name="h_enriched")
Figure_1J <- FeaturePlot(hNeut_GSE127465_s,
                         features = "h_enriched1", max.cutoff = 10) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","THBS1","PTMA")
hNeut_GSE127465_s <- AddModuleScore(hNeut_GSE127465_s,
                                    features = list(t_enriched),
                                    name="T_enriched")
Figure_1K <- FeaturePlot(hNeut_GSE127465_s,
                         features = "T_enriched1",max.cutoff = 15) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


m_enriched <- c("CD74","TXNIP","STK17B","CXCR2","FKBP5","CXCR1","CTSS","CEBPD","JAML","IRF1")
hNeut_GSE127465_s <- AddModuleScore(hNeut_GSE127465_s,
                                    features = list(m_enriched),
                                    name="M_enriched")
Supplementary_Figure_S3A <- FeaturePlot(hNeut_GSE127465_s,
                                        features = "M_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Scoring for neutrophil signatures in the human dataset: hGSE114727_t_h (breast cancer)
t_enriched <- c("CDKN1A","PPIA","IFITM1","IL1B","ISG15","CCL4","CD14","IER3","CCL3","IFIT3","PTMA","THBS1","IFIT1","TAGLN2","GNGT2","RPS27L")
hGSE114727_t_h <- AddModuleScore(hGSE114727_t_h,
                                 features = list(t_enriched),
                                 name="T_enriched", nbin = 5)
Figure_1N <- FeaturePlot(hGSE114727_t_h,
            features = "T_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
hGSE114727_t_h <- AddModuleScore(hGSE114727_t_h,
                                 features = list(h_enriched),
                                 name="H_enriched", nbin = 5)

Figure_1M <- FeaturePlot(hGSE114727_t_h,
            features = "H_enriched1") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Processing publicly available dataset: OEP001756 (human, Colorectal cancer liver metastasis) 
##download matrix, features and barcodes from National Omics Encyclopaedia
###Save datasets as: matrix.mtx.gz, genes.tsv.gz and barcodes.tsv.gz

#create expression matrix and run standard seurat workflow
expression_matrix <- ReadMtx(
  mtx = "matrix.mtx.gz", features = "genes.tsv.gz",
  cells = "barcodes.tsv.gz"
)
hCRC_s <- CreateSeuratObject(counts = expression_matrix, names.field = 2)

VlnPlot(hCRC_s, features = c("nFeature_RNA", "nCount_RNA"))
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

#subsetting only neutrophils based on the marker expression reported in original publication
FeaturePlot(hCRC_s, features = c("FCGR3B", "LYZ"), raster = FALSE)
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

Figure_2A <- DimPlot(hCRC_neut_s, reduction = "umap")

hCRC.cluster.markers <- FindAllMarkers(hCRC_neut_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hCRC.cluster.markers.top10<- hCRC.cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(hCRC.cluster.markers.top10, file="hCRC.neut.top10.markers.csv")

#Scoring for neutrophil signatures in the human dataset: OEP001756 (Colorectal cancer liver metastasis)
t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","GNGT2","THBS1","PTMA")
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

Figure_2B <- FeaturePlot(hCRC_neut_s, features = c("H_enriched1", "T_enriched1"), blend = TRUE, cols = c("black","green","red"))

saveRDS(hCRC_neut_s, file="hCRC_neut_s.rds")

#plot volcanoplot of upregulated genes in cluster 0 (M-specific neutrophil population)
##run DEG compared to all clusters
DimPlot(hCRC_neut_s, reduction = "umap")
Idents(hCRC_neut_s)
Idents(hCRC_neut_s)<- hCRC_neut_s@meta.data$seurat_clusters

Cluster0_markers_comp <- FindMarkers(hCRC_neut_s, ident.1 = "0")

Figure_2L <- EnhancedVolcano(Cluster0_markers_comp,
                lab = rownames(Cluster0_markers_comp),
                x = 'avg_log2FC',
                y = 'p_val',col=c('black', 'black', 'black', 'red3'),FCcutoff = 1.5,
                colAlpha = 1,labSize = 4,boxedLabels = T, drawConnectors = T,selectLab = c('TXNIP','RIPOR2','FKBP5','CEBPD','STK17B','CXCR2','JAML','MME','CXCL8','CXCL2','MIF','CCL3','PI3','C15orf48' ),
                labCol = c('red','red','red','red','red','red','red','red','blue','blue','blue','blue','blue','blue' ), title = "Differentially Expressed Genes in Metastasis-specific Cluster")

#slingshot pseudotime analysis
##slingshot analysis for the integrated mouse dataset: KPN_NT_NSLC_int
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

cell_colors <- cell_pal(KPN_NT_NSLC_int$Tissue.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(KPN_NT_NSLC_int$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

#differential gene expression across pseudotime
dimred <- KPN_NT_NSLC_int@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = KPN_NT_NSLC_int.sce$seurat_clusters, dimred, start.clus= 5)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
counts <- as.matrix(KPN_NT_NSLC_int@assays$RNA@counts)
dim(counts)
sce <- fitGAM(counts = as.matrix(counts), sds = curves)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves, nknots=10)
Figure_1O <- plotGeneCount(curves, filt_counts, clusters = KPN_NT_NSLC_int$seurat_clusters, models = sce)

# Define function to plot
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

#start and end tests:
startRes <- startVsEndTest(sce)
#Visualize the estimated smoothers for the Nth most significant gene.
##color the cells in UMAP space with that gene’s expression
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]

Figure_1S <- plotGeneCount(lineages, counts, gene = sigGeneStart)
Figure_1W <-plotSmoothers(sce, counts, gene = sigGeneStart)

##slingshot analysis for the dataset: mNeut_GSE139125s
PYMT.markers <- FindAllMarkers(mNeut_GSE139125s)
PYMT.markers_Top10<- PYMT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

mNeut_GSE139125s.sce <- as.SingleCellExperiment(mNeut_GSE139125s)

#using slingshot
sce <- slingshot(mNeut_GSE139125s.sce, clusterLabels = mNeut_GSE139125s.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')
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
cell_colors <- cell_pal(mNeut_GSE139125s$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(mNeut_GSE139125s$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')

counts <- as.matrix(mNeut_GSE139125s@assays$RNA@counts[mNeut_GSE139125s@assays$RNA@var.features, ])
dimred <- mNeut_GSE139125s@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = mNeut_GSE139125s.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
dim(counts)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

Figure_1P <- plotGeneCount(curves, filt_counts, clusters = mNeut_GSE139125s.sce$seurat_clusters, models = sce)

# Define function to plot
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = mNeut_GSE139125s.sce$seurat_cluster, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

pseudotime_association <- associationTest(sce)

head(pseudotime_association)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)
feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#start vs end Test
startRes <- startVsEndTest(sce)
write.csv(startRes, file="startRes_PYMT.csv")
#visualize the estimated smoothers for selected most significant gene.
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[22]]

Figure_1T <- plotGeneCount(lineages, counts, gene = sigGeneStart)
Figure_1X <- plotSmoothers(sce, counts, gene = sigGeneStart)

##slingshot analysis for the dataset: R18_No.KPN
#R18_No.KPN is the seurat object from CRC_other_counts on Zenodo

R18_No.KPN.sce <- as.SingleCellExperiment(R18_No.KPN)
sce <- slingshot(R18_No.KPN.sce, clusterLabels = R18_No.KPN.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')
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
cell_colors <- cell_pal(R18_No.KPN$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(R18_No.KPN$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')

counts <- as.matrix(R18_No.KPN@assays$RNA@counts)
dimred <- R18_No.KPN@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = R18_No.KPN.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
dim(counts)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

Figure_1Q <- plotGeneCount(curves, filt_counts, clusters = R18_No.KPN.sce$seurat_clusters, models = sce)

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

#start vs end test
startRes <- startVsEndTest(sce)
write.csv(startRes, file="startRes_R18.NO.KPN.csv")
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[559]]

Figure_1U <- plotGeneCount(lineages, counts, gene = sigGeneStart)
Figure_1Y <- plotSmoothers(sce, counts, gene = sigGeneStart)

#slingshot analysis for the human dataset: hNeut_GSE127465_s

hNeut_GSE127465_s.sce <- as.SingleCellExperiment(hNeut_GSE127465_s)
sce <- slingshot(hNeut_GSE127465_s.sce, clusterLabels = hNeut_GSE127465_s.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')
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
cell_colors <- cell_pal(hNeut_GSE127465_s$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(hNeut_GSE127465_s$seurat_clusters, hue_pal())

plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
counts <- as.matrix(hNeut_GSE127465_s@assays$RNA@counts)
dimred <- hNeut_GSE127465_s@reductions$umap@cell.embeddings
lineages <- getLineages(data = dimred, clusterLabels = hNeut_GSE127465_s.sce$seurat_clusters, dimred)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
dim(counts)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

Figure_1R <- plotGeneCount(curves, filt_counts, clusters = hNeut_GSE127465_s.sce$seurat_clusters, models = sce)

plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = mNeut_GSE139125s.sce$seurat_cluster, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

pseudotime_association <- associationTest(sce)
head(pseudotime_association)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)
feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#start vs end test
startRes <- startVsEndTest(sce)
write.csv(startRes, file="startRes_h_NSCLC.csv")
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[316]]

Figure_1V <- plotGeneCount(lineages, counts, gene = sigGeneStart)
Figure_1Z <- plotSmoothers(sce, counts, gene = sigGeneStart)

#Slingshot Pseudotime analysis in the human dataset: hCRC_neut_s (Colorectal cancer liver metastasis)
hCRC_neut_s.sce <- as.SingleCellExperiment(hCRC_neut_s)
sce <- slingshot(hCRC_neut_s.sce, clusterLabels = hCRC_neut_s.sce$seurat_clusters,stretch = 0 ,reducedDim = 'UMAP')
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

cell_colors_clust <- cell_pal(hCRC_neut_s$seurat_clusters, hue_pal())
plot(reducedDims(sce)$UMAP, col = cell_colors_clust, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')

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

#End Test to examine differentially expressed genes at lineage ends
endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[70]]

Figure_2C <- plotGeneCount(curves, counts, gene = sigGene)
Figure_2D <- plotSmoothers(sce, counts, sigGene)

#start vs end Test
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[25]]

Figure_2F <- plotGeneCount(lineages, counts, gene = sigGeneStart)
Figure_2G <- plotSmoothers(sce, counts, gene = sigGeneStart)
Figure_2E <- FeaturePlot(hCRC_neut_s, features = c("TXNIP", "CXCR2"), blend = TRUE, cols = c("GREY","RED","GREEN"))
Figure_2H <- FeaturePlot(hCRC_neut_s, features = c("CXCL8", "IL1B"), blend = TRUE, cols = c("GREY","RED","GREEN"))


#subset neutrophils from Primary tumours (NSCLC and BC) and Metastatic CRC datasets processed above

hGSE114727_t<-subset(x = hGSE114727_cls, idents = c("TUMOR"))
DimPlot(hGSE114727_t)

write_rds(hGSE114727_t, file="hBCPT_neut.rds") #This dataset contains neutrophils from primary tumour, human Breast cancer dataset

hNSCLCPT_neut<-subset(x = hNeut_GSE127465_s, idents = "tumor")
write_rds(hNSCLCPT_neut, file="hNSCLCPT_neut.rds")#This dataset contains neutrophils from primary tumour, human NSCLC dataset

write_rds(hCRC_neut_s, file="hCRCLM_neut.rds")#This dataset contains neutrophils from metastatic human CRC dataset.

#load and view the datasets

DimPlot(hBCPT_neut)
DimPlot(hNSCLCPT_neut)
DimPlot(hCRCLM_neut)

#integrate the datasets
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = c(hBCPT_neut,hCRCLM_neut,hNSCLCPT_neut))
anchors <- FindIntegrationAnchors(object.list = c(hBCPT_neut,hCRCLM_neut,hNSCLCPT_neut), anchor.features = features)

T_int.list =  list(hBCPT_neut,hCRCLM_neut,hNSCLCPT_neut)
all_genes <- lapply(T_int.list, row.names) %>% Reduce(intersect, .) 

T_int <- IntegrateData(anchorset = anchors, features.to.integrate = all_genes)
DefaultAssay(T_int) <- "integrated"
Idents(T_int)

T_int <- RenameIdents(object = T_int,  'tumor' = 'PT_NSCLC', 'TUMOR' = 'PT_BC')
T_int <- RenameIdents(object = T_int,  '12' = 'M_CRC')

T_int@active.ident <- factor(T_int@active.ident,
                             levels=c("PT_NSCLC","PT_BC","M_CRC"))
T_int@meta.data$Tumor.type <- T_int@active.ident

# Run the standard workflow for visualization and clustering
T_int <- ScaleData(T_int, verbose = FALSE)
T_int <- RunPCA(T_int, npcs = 30, verbose = FALSE)
T_int <- RunUMAP(T_int, reduction = "pca", dims = 1:30)
T_int <- FindNeighbors(T_int, reduction = "pca", dims = 1:30)
T_int <- FindClusters(T_int, resolution = 0.5)

Figure_2O <- DimPlot(T_int, reduction = "umap", label=T, repel = T)
Figure_2P1 <- FeaturePlot(T_int, features = c("IL1B", "CXCR2"), blend = TRUE, cols = c("GREY","RED","GREEN"), max.cutoff = 10)
Figure_2P2 <- FeaturePlot(T_int, features = c("TXNIP", "CXCR2"), blend = TRUE, cols = c("GREY","RED","GREEN"), max.cutoff = 10)

# Find differentially expressed features between M_CRC (metastatic) and all other neutrophils derived from Primary tumours
M.CRC.markers <- FindMarkers(T_int, ident.1 = "M_CRC")
write.csv(M.CRC.markers, file="MCRC_markers.csv")

Figure_2Q <- EnhancedVolcano(M.CRC.markers,
                lab = M.CRC.markers$Gene,
                x = 'avg_log2FC',
                y = 'p_val')

#GO and KEGG analysis and exploring implicated pathways
M.CRC.markers <- tibble::rownames_to_column(M.CRC.markers, "Gene") # Apply rownames_to_column
M.CRC.gene.list <- M.CRC.markers$avg_log2FC
names(M.CRC.gene.list) <- M.CRC.markers$Gene
M.CRC.gene.list<-na.omit(M.CRC.gene.list)
M.CRC.gene.list = sort(M.CRC.gene.list, decreasing = TRUE)
gse <- gseGO(geneList=M.CRC.gene.list, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 1, 
             maxGSSize = 37, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=30)
d <- godata('org.Hs.eg.db', ont="BP")
ego2 <- pairwise_termsim(gse, method="Wang", semData = d)
install.packages("ggnewscale")
emapplot(ego2, showCategory = 15)
d_CC <- godata('org.Hs.eg.db', ont="CC")
ego2_CC <- pairwise_termsim(gse, method="Wang", semData = d_CC)

emapplot(ego2_CC, showCategory = 15)
emapplot_cluster(ego2)
emapplot_cluster(ego2_CC)

cnetplot(gse, categorySize="pvalue", foldChange=M.CRC.gene.list)
gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 1)

#KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(M.CRC.gene.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CRC <- M.CRC.markers
KEGG.CRC <- KEGG.CRC[-24,]
KEGG.CRC$Y = dedup_ids$ENTREZID
# Create a vector 
KEGG.CRC.gene.list <- KEGG.CRC$avg_log2FC
# Name vector with ENTREZ ids
names(KEGG.CRC.gene.list) <- KEGG.CRC$Y
# omit any NA values 
KEGG.CRC.gene.list<-na.omit(KEGG.CRC.gene.list)
# sort the list in decreasing order (required for clusterProfiler)
KEGG.CRC.gene.list = sort(KEGG.CRC.gene.list, decreasing = TRUE)
#Create gseKEGG object
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = KEGG.CRC.gene.list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 30,
               pvalueCutoff = 0.01,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, title = "Enriched Pathways")
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=KEGG.KPN.gene.list, pathway.id="mmu04620", species = "mmu")
knitr::include_graphics("mmu04620.pathview.png")
dme

#run reactome GSA
gsva_result <- analyse_sc_clusters(T_int, verbose = TRUE)
gsva_result
# Combines and returns the pathways of all clusters
pathway_expression <- pathways(gsva_result)
# Simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
# Find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
# Calculate the maximum expression difference for each pathway
max_difference$diff <- max_difference$max - max_difference$min
# Sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
# Check the first few pathways and their maximum expression differences
head(max_difference)
# Plot the expression of a single pathway across the different cell types
plot_obj <- plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-70171")
print(plot_obj)

max.diff.pathways <-max_difference
max.diff.pathways<- tibble::rownames_to_column(max_difference, "PathwayID") # Apply rownames_to_column
pathwayIDs<- max.diff.pathways$PathwayID
plot_gsva_heatmap(gsva_result, pathway_ids = pathwayIDs[1:40])
plot_gsva_pca(gsva_result)

Figure_2R <- DEenrichRPlot(object = T_int, ident.1 = "M_CRC",ident.2 = "PT_NSCLC",enrich.database = "GO_Biological_Process_2018",max.genes = 50)
Figure_2S <-DEenrichRPlot(object = T_int, ident.1 = "M_CRC",ident.2 = "PT_NSCLC",enrich.database = "KEGG_2021_Human",max.genes = 50)

#isolating Tcells from hCRC_s dataset (metastatic CRC dataset, human processed above) 
#Use markers in original publication using a threshhold >1,5 not >3 to include all cells.
#CD4 for CD4 T cells
##FOXP3 for Tregs
###CD8A for CD8 T cells
####SLC4A10 for MAIT cells, 

FeaturePlot(hCRC_s, features = c("CD4","CD8A","FOXP3","SLC4A10"))
#subsetting Tcells based on marker expression reported in paper
hCRC_CD4_TC_s <-subset(x = hCRC_s, subset = CD4 > 1.5)
hCRC_CD8_TC_s <-subset(x = hCRC_s, subset = CD8A > 1.5)
hCRC_FOXP3_TC_s <-subset(x = hCRC_s, subset = FOXP3 > 1.5)
hCRC_mait_TC_s <-subset(x = hCRC_s, subset = SLC4A10 > 1.5)
#merge objects
hCRCLM_TC <- merge(x = hCRC_CD4_TC_s, y = list(hCRC_CD8_TC_s, hCRC_FOXP3_TC_s,hCRC_mait_TC_s))
#recluster to analyse
hCRCLM_TC <- FindVariableFeatures(hCRCLM_TC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRCLM_TC), 10)
plot1 <- VariableFeaturePlot(hCRCLM_TC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

hCRCLM_TC <- ScaleData(hCRCLM_TC)
hCRCLM_TC <- RunPCA(hCRCLM_TC, features = VariableFeatures(object = hCRCLM_TC))
VizDimLoadings(hCRCLM_TC, dims = 1:2, reduction = "pca")
DimPlot(hCRCLM_TC, reduction = "pca")
ElbowPlot(hCRCLM_TC)
hCRCLM_TC <- FindNeighbors(hCRCLM_TC, dims = 1:20)#SELECT FIRST 20 PCS
hCRCLM_TC <- FindClusters(hCRCLM_TC, resolution = 0.5)
hCRCLM_TC <- RunUMAP(hCRCLM_TC, dims = 1:20)
DimPlot(hCRCLM_TC, reduction = "umap", label = T)

FeaturePlot(hCRCLM_TC, features = c("CD4","CD8A","FOXP3","SLC4A10"))

# Rename identity classes
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `13` = "CD4 T cell")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `17` = "Treg")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `7` = "CD8 T cell")
Idents(hCRCLM_TC)
hCRCLM_TC$Cell.type <- Idents(hCRCLM_TC)
DimPlot(hCRCLM_TC, reduction = "umap", group.by = "Cell.type", label = T, label.size = 5)
saveRDS(hCRCLM_TC, file="hCRCLM_TC.rds")

#Processing publicly available dataset: GSE146771 (Primary tumor, colorectal cancer, Human)
##Download metadata and TPM from GEO database
### Save them as: GSE146771_CRC.Leukocyte.10x.Metadata (1).txt.gz and GSE146771_CRC.Leukocyte.10x.TPM.txt.gz
GSE146771_metadata <- read_tsv("GSE146771_CRC.Leukocyte.10x.Metadata (1).txt.gz", col_names = T)
GSE146771 <- read.table("GSE146771_CRC.Leukocyte.10x.TPM.txt.gz",sep="")

head(GSE146771[1:5,1:5])
h_GSE146771_s <- CreateSeuratObject(counts = GSE146771)
cell_IDs <- CellsByIdentities(h_GSE146771_s)
cell_IDs <- as.data.frame(cell_IDs)
names(cell_IDs) <- "cell_ID"
cells$ID <- cell_IDs$cell_ID

#extract Global_Cluster = cell type from meta data
Cell_type <- as.data.frame(GSE146771_metadata[, c(2,15)])
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

#subset tcells as annotated in metadata from original publication
hCRCPT_TC<-subset(x = h_GSE146771_s, idents = c("CD4 T cell", "CD8 T cell"))
saveRDS(hCRCPT_TC, file = "hCRCPT_TC.rds")

#Downstream analysis for clustering
hCRCPT_TC <- FindVariableFeatures(hCRCPT_TC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hCRCPT_TC), 10)
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

#Integrate TCell datasets from primary and metastatic CRC

DimPlot(hCRCLM_TC, reduction = "umap", group.by = "Cell.type")
DimPlot(hCRCPT_TC, reduction = "umap", group.by = "Cell.type")

# Rename identity classes
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `CD4 T cell` = "LM_CD4")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `CD8 T cell` = "LM_CD8")
hCRCLM_TC <- RenameIdents(object = hCRCLM_TC, `Treg` = "LM_Treg")

Idents(hCRCPT_TC) <- hCRCPT_TC$Cell.type
hCRCPT_TC <- RenameIdents(object = hCRCPT_TC, `CD4 T cell` = "PT_CD4")
hCRCPT_TC <- RenameIdents(object = hCRCPT_TC, `CD8 T cell` = "PT_CD8")

hCRCPT_TC$Cell.type <- Idents(hCRCPT_TC)
hCRCLM_TC$Cell.type <- Idents(hCRCLM_TC)

Idents(hCRCLM_TC) <-"Cell.type"
Idents(hCRCPT_TC) <-"Cell.type"

#Integrate datasets
features <- SelectIntegrationFeatures(object.list = c(hCRCPT_TC,hCRCLM_TC))
anchors <- FindIntegrationAnchors(object.list = c(hCRCPT_TC,hCRCLM_TC), anchor.features = features)
hCRC_TC_PTvLM1.5 <- IntegrateData(anchorset = anchors)
DefaultAssay(hCRC_TC_PTvLM1.5) <- "integrated"
Idents(hCRC_TC_PTvLM1.5)

hCRC_TC_PTvLM1.5@active.ident <- factor(hCRC_TC_PTvLM1.5@active.ident,
                                        levels=c("PT_CD4","PT_CD8","LM_CD4","LM_CD8","LM_Treg"))

hCRC_TC_PTvLM1.5@meta.data$Cell.type <- hCRC_TC_PTvLM1.5@active.ident
Idents(hCRC_TC_PTvLM1.5) <- hCRC_TC_PTvLM1.5$Cell.type
hCRC_TC_PTvLM1.5 <- RenameIdents(object = hCRC_TC_PTvLM1.5, `PT_CD4` = "PT")
hCRC_TC_PTvLM1.5 <- RenameIdents(object = hCRC_TC_PTvLM1.5, `PT_CD8` = "PT")
hCRC_TC_PTvLM1.5 <- RenameIdents(object = hCRC_TC_PTvLM1.5, `LM_CD4` = "LM")
hCRC_TC_PTvLM1.5 <- RenameIdents(object = hCRC_TC_PTvLM1.5, `LM_CD8` = "LM")
hCRC_TC_PTvLM1.5 <- RenameIdents(object = hCRC_TC_PTvLM1.5, `LM_Treg` = "LM")

Idents(hCRC_TC_PTvLM1.5)
hCRC_TC_PTvLM1.5@meta.data$Tumor.type <- hCRC_TC_PTvLM1.5@active.ident

# Run the standard workflow for visualization and clustering
hCRC_TC_PTvLM1.5 <- ScaleData(hCRC_TC_PTvLM1.5, verbose = FALSE)
hCRC_TC_PTvLM1.5 <- RunPCA(hCRC_TC_PTvLM1.5, npcs = 30, verbose = FALSE)
hCRC_TC_PTvLM1.5 <- RunUMAP(hCRC_TC_PTvLM1.5, reduction = "pca", dims = 1:30)
hCRC_TC_PTvLM1.5 <- FindNeighbors(hCRC_TC_PTvLM1.5, reduction = "pca", dims = 1:30)
hCRC_TC_PTvLM1.5 <- FindClusters(hCRC_TC_PTvLM1.5, resolution = 0.5)

Figure_3A <- DimPlot(hCRC_TC_PTvLM1.5, reduction = "umap", label=T, group.by = "Cell.type")
Figure_3B <- DimPlot(hCRC_TC_PTvLM1.5, reduction = "umap", label=T, group.by = "Tumor.type")

saveRDS(hCRC_TC_PTvLM1.5, file="hCRC_TC_PTvLM1.5.rds")

#finding markers per  tumor type
hCRC_TC_PTvLM1.5 <- SetIdent(hCRC_TC_PTvLM1.5, value = hCRC_TC_PTvLM1.5@meta.data$Cell.type)
PTCD8.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "PT_CD8")
LMCD8.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "LM_CD8")

PTCD4.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "PT_CD4")
LMCD4.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "LM_CD4")

CD4.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "LM_CD4", ident.2 = "PT_CD4")
CD8.markers1.5 <- FindMarkers(hCRC_TC_PTvLM1.5, ident.1 = "LM_CD8", ident.2 = "PT_CD8")

Figure_3C <- EnhancedVolcano(CD4.markers1.5,
                             lab = CD4.markers1.5$Gene,
                             x = 'avg_log2FC',
                             y = 'p_val', FCcutoff = 0.5, title = "CD4 Tcells LM v PT")

CD4.markers1.5 <- tibble::rownames_to_column(CD4.markers1.5, "Gene") # Apply rownames_to_column

CD4.gene.list1.5 <- CD4.markers1.5$avg_log2FC
names(CD4.gene.list1.5) <- CD4.markers1.5$Gene
CD4.gene.list1.5<-na.omit(CD4.gene.list1.5)
CD4.gene.list1.5 = sort(CD4.gene.list1.5, decreasing = TRUE)

gse1.5 <- gseGO(geneList=CD4.gene.list1.5, 
                ont ="ALL", 
                keyType = "SYMBOL",
                minGSSize = 1, 
                maxGSSize = 200, 
                pvalueCutoff = 0.01, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "none")

Figure_3D <- dotplot(gse1.5, showCategory=15)

d <- godata('org.Hs.eg.db', ont="CC")
CD41.5 <- pairwise_termsim(gse1.5, method="Wang", semData = d)
install.packages("ggnewscale")
emapplot(CD41.5, showCategory = 15)
emapplot_cluster(CD41.5)
cnetplot(gse1.5, categorySize="pvalue", foldChange=CD4.gene.list1.5)

#KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
CD4ids1.5<-bitr(names(CD4.gene.list1.5), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# remove duplicate IDS 
CD4dedup_ids1.5 = CD4ids1.5[!duplicated(CD4ids1.5[c("SYMBOL")]),]
# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CD41.5 <- CD4.markers1.5
KEGG.CD41.5$Y = CD4dedup_ids1.5$ENTREZID
CD4dedup_ids1.5 <- CD4dedup_ids1.5 %>% 
  rename(
    Gene = SYMBOL)

#identify unmatching rows
anti_join(KEGG.CD41.5, CD4dedup_ids1.5, by="Gene")
#Remove rows where Gene is not equal to genes not present in CD4dedup_ids
KEGG.CD4.filtered1.5 <- subset(KEGG.CD41.5,Gene !='MARCH1' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='H1FX' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='MT-ND3' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='MT-ND1' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='H3F3A' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='HIST1H4C' )
KEGG.CD4.filtered1.5 <- subset(KEGG.CD4.filtered1.5,Gene !='H2AFZ' )

# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CD4.filtered1.5$Y = CD4dedup_ids1.5$ENTREZID
## Create a vector of the gene unuiverse
KEGG.CD4.gene.list1.5 <- KEGG.CD4.filtered1.5$avg_log2FC
# Name vector with ENTREZ ids
names(KEGG.CD4.gene.list1.5) <- KEGG.CD4.filtered1.5$Y
# omit any NA values 
KEGG.CD4.gene.list1.5<-na.omit(KEGG.CD4.gene.list1.5)
# sort the list in decreasing order (required for clusterProfiler)
KEGG.CD4.gene.list1.5 = sort(KEGG.CD4.gene.list1.5, decreasing = TRUE)
#Create gseKEGG object
kegg_organism = "hsa"
kk2.CD41.5 <- gseKEGG(geneList     = KEGG.CD4.gene.list1.5,
                      organism     = kegg_organism,
                      nPerm        = 10000,
                      minGSSize    = 3,
                      maxGSSize    = 200,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "none",
                      keyType       = "ncbi-geneid")

Figure_3E <- dotplot(kk2.CD41.5, title = "CD4 T cells Enriched Pathways")

#CellChat Analysis in neutrophils and TCells of human Metastatic CRC dataset

DimPlot(hCRC_neut_s, reduction = "umap", label = T)

# Rename identity classes
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `12` = "Neutrophils")
Idents(hCRC_neut_s)
hCRC_neut_s$Cell.type <- Idents(hCRC_neut_s)

#Merge Neutrophil and Tcells in CRCLM to run CellChat
CRCLM_N_TC <- merge(hCRC_neut_s, y = hCRCLM_TC, add.cell.ids = c("Neutrophils", "Tcells"), project = "CRCLM")
table(CRCLM_N_TC$Cell.type)
DimPlot(CRCLM_N_TC, reduction = "umap", label = T)

#recluster to analyse
CRCLM_N_TC <- FindVariableFeatures(CRCLM_N_TC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CRCLM_N_TC), 10)
plot1 <- VariableFeaturePlot(CRCLM_N_TC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

CRCLM_N_TC <- ScaleData(CRCLM_N_TC)
CRCLM_N_TC <- RunPCA(CRCLM_N_TC, features = VariableFeatures(object = CRCLM_N_TC))
VizDimLoadings(CRCLM_N_TC, dims = 1:2, reduction = "pca")
DimPlot(CRCLM_N_TC, reduction = "pca")
ElbowPlot(CRCLM_N_TC)

CRCLM_N_TC <- FindNeighbors(CRCLM_N_TC, dims = 1:20)#SELECT FIRST 20 PCS
CRCLM_N_TC <- FindClusters(CRCLM_N_TC, resolution = 0.5)
CRCLM_N_TC <- RunUMAP(CRCLM_N_TC, dims = 1:20)
DimPlot(CRCLM_N_TC, reduction = "umap", label = T)

DimPlot(CRCLM_N_TC, reduction = "umap", label = T, group.by = "Cell.type")
# Create cellchat object
cellchat <- createCellChat(object = CRCLM_N_TC, group.by = "Cell.type", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

Figure_3F <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
Figure_3G <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Supplementary_Figure_S3B
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}

#visualisation of the different pathways & interactions
#pathways.show <- c("MHC-I", "MHC-II","MIF","CD45" ,"GALECTIN", "CLEC","CXCL" , "ANNEXIN","ADGRE5","CCL","CD99", "ITGB2","LCK","IL1","TNF","VCAM","CD86","ICAM", "IFN-II","PECAM1") 
pathways.show <- c("ICAM")
par(mfrow=c(1,1))

Supplementary_Figure_S4 <- netVisual_aggregate(cellchat, signaling = pathways.show , layout = "circle") #repeat for the other 19 pathways listed above
vertex.receiver = seq(1,2)

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
## visualize cell-cell communication mediated by a single ligand-receptor pair
###extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.

pathways.show <- c("CXCL")
Figure_3H <- netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show first ligand-receptor pair
Figure_3k1 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[2,] # show second ligand-receptor pair
Figure_3k2 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[3,] # show third ligand-receptor pair
Figure_3k3 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[4,] # show fourth ligand-receptor pair
Figure_3k4 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[5,] # show fifth ligand-receptor pair
Figure_3k5 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[6,] # show sixth ligand-receptor pair
Figure_3k6 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#Repeat for IL1 pathway
pathways.show <- c("IL1")
Figure_3I <- netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.IL1 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[1,] # show one ligand-receptor pair
Figure_3L <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#Repeat for TNF pathway
pathways.show <- c("TNF")
Figure_3J <- netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.TNF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TNF[1,] # show one ligand-receptor pair
Figure_3M2 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.TNF[2,] # show one ligand-receptor pair
Figure_3M1 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
Supplementary_Figure_S5A <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4), remove.isolate = FALSE)
Supplementary_Figure_S5B <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:4), remove.isolate = FALSE)
Supplementary_Figure_S5C <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4), remove.isolate = FALSE)
Supplementary_Figure_S5D <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3), remove.isolate = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
Figure_4C2 <- netAnalysis_signalingRole_scatter(cellchat)

#which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#how multiple cell groups and signaling pathways coordinate to function?
#run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")

#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 7
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
Figure_4A <- netAnalysis_river(cellchat, pattern = "outgoing")

#incoming patterns
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
Figure_4B <- netAnalysis_river(cellchat, pattern = "incoming")
saveRDS(cellchat, file = "cellchat_hCRCLM_N_Tc.rds")

#CellChat Analysis in human Primary CRC dataset
## Create cellchat object
cellchat <- createCellChat(object = hCRCPT_TC, group.by = "Cell.type", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
Figure_4C1 <- netAnalysis_signalingRole_scatter(cellchat)

##Analysis and characterisation of Neutrophil subtypes in human metastatic CRC dataset
DimPlot(hCRC_neut_s, reduction = "umap", group.by = "seurat_clusters", label = T)
Idents(object = hCRC_neut_s) <- "seurat_clusters"

##rename clusters accoring to subtype
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `0` = "TXNIP+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `1` = "Inf_reg_Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `2` = "COX+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `3` = "HSP+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `4` = "IFN+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `5` = "Pro-apoptotic_Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `6` = "PLPP3+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `7` = "RPS+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `8` = "ARG1+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `9` = "MT+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `10` = "Activated_Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `11` = "HLA-DR+Neut")
hCRC_neut_s <- RenameIdents(object = hCRC_neut_s, `12` = "Inf_reg_Neut")

Figure_4D <- DimPlot(hCRC_neut_s, reduction = "umap", label = T, label.size = 5)

Idents(hCRC_neut_s)
hCRC_neut_s$Cell.type <- Idents(hCRC_neut_s)


#Merge Neutrophil, Tcells and macrophages (Mph) in human metastatic CRC dataset (CRCLM) to run CellChat

CRCLM_N_TC_Mph <- merge(hCRC_neut_s, y = list(hCRCLM_TC,hCRC_Mph_s1), add.cell.ids = c("Neutrophils", "Tcells","Macrophages"), project = "CRCLM")
table(CRCLM_N_TC_Mph$Cell.type)

#recluster to analyse
CRCLM_N_TC_Mph <- FindVariableFeatures(CRCLM_N_TC_Mph, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CRCLM_N_TC_Mph), 10)
plot1 <- VariableFeaturePlot(CRCLM_N_TC_Mph)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

CRCLM_N_TC_Mph <- ScaleData(CRCLM_N_TC_Mph)
CRCLM_N_TC_Mph <- RunPCA(CRCLM_N_TC_Mph, features = VariableFeatures(object = CRCLM_N_TC_Mph))
VizDimLoadings(CRCLM_N_TC_Mph, dims = 1:2, reduction = "pca")
DimPlot(CRCLM_N_TC_Mph, reduction = "pca")
ElbowPlot(CRCLM_N_TC_Mph)
CRCLM_N_TC_Mph <- FindNeighbors(CRCLM_N_TC_Mph, dims = 1:20)#SELECT FIRST 20 PCS
CRCLM_N_TC_Mph <- FindClusters(CRCLM_N_TC_Mph, resolution = 0.5)
CRCLM_N_TC_Mph <- RunUMAP(CRCLM_N_TC_Mph, dims = 1:20)

DimPlot(CRCLM_N_TC_Mph, reduction = "umap", label = T, group.by = "Cell.type", repel = T, label.box = T)
DimPlot(CRCLM_N_TC_Mph, reduction = "umap", label = T, split.by  = "active.ident", repel = T, label.box = T)

saveRDS(CRCLM_N_TC_Mph, file="CRCLM_N_TC_Mph.rds")

# Create cellchat object and run cellchat analysis
cellchat <- createCellChat(object = CRCLM_N_TC_Mph, group.by = "Cell.type", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB #use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

Figure_4E <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}


Figure_4F <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, idents.use = c("TXNIP+Neut"), arrow.size = 1)
Figure_4G <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, idents.use = c("SPP1+Mph"), arrow.size = 1)
Figure_4H <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, sources.use = c("M1-like-Mph"), arrow.size = 1)
Figure_4I <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, targets.use = c("M2-like-Mph"), arrow.size = 1)

#visualisation of the different pathways & interactions
cellchat@netP[["pathways"]]
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show , layout = "circle")

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
#and visualize cell-cell communication mediated by a single ligand-receptor pair
Figure_4J <- netAnalysis_contribution(cellchat, signaling = pathways.show)

#visualize the cell-cell communication mediated by a single ligand-receptor pair.
#extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL<- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[6,] 
Figure_4M <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[8,] 
Figure_4K <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

LR.show <- pairLR.CXCL[9,] 
Figure_4L <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
Supplementary_Figure_S6A <- netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1:9,11:19), remove.isolate = FALSE)
Supplementary_Figure_S6B <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10,12:19), remove.isolate = FALSE)
Supplementary_Figure_S6C <- netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:11,13:19), remove.isolate = FALSE)
Supplementary_Figure_S6D <- netVisual_bubble(cellchat, sources.use = 17, targets.use = c(1:16,18:19), remove.isolate = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

#which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

Figure_4N <- ht1 + ht2

saveRDS(cellchat, file = "cellchat_hCRCLM_N_Tc_Mph.rds")

