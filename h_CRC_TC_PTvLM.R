library (Seurat)

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

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = c(hCRCPT_TC,hCRCLM_TC))
anchors <- FindIntegrationAnchors(object.list = c(hCRCPT_TC,hCRCLM_TC), anchor.features = features)

hCRC_TC_PTvLM <- IntegrateData(anchorset = anchors)
DefaultAssay(hCRC_TC_PTvLM) <- "integrated"

Idents(hCRC_TC_PTvLM)

hCRC_TC_PTvLM@active.ident <- factor(hCRC_TC_PTvLM@active.ident,
                                       levels=c("PT_CD4","PT_CD8","LM_CD4","LM_CD8","LM_Treg"))

hCRC_TC_PTvLM@meta.data$Cell.type <- hCRC_TC_PTvLM@active.ident
Idents(hCRC_TC_PTvLM) <- hCRC_TC_PTvLM$Cell.type
hCRC_TC_PTvLM <- RenameIdents(object = hCRC_TC_PTvLM, `PT_CD4` = "PT")
hCRC_TC_PTvLM <- RenameIdents(object = hCRC_TC_PTvLM, `PT_CD8` = "PT")
hCRC_TC_PTvLM <- RenameIdents(object = hCRC_TC_PTvLM, `LM_CD4` = "LM")
hCRC_TC_PTvLM <- RenameIdents(object = hCRC_TC_PTvLM, `LM_CD8` = "LM")
hCRC_TC_PTvLM <- RenameIdents(object = hCRC_TC_PTvLM, `LM_Treg` = "LM")

Idents(hCRC_TC_PTvLM)
hCRC_TC_PTvLM@meta.data$Tumor.type <- hCRC_TC_PTvLM@active.ident

# Run the standard workflow for visualization and clustering
hCRC_TC_PTvLM <- ScaleData(hCRC_TC_PTvLM, verbose = FALSE)
hCRC_TC_PTvLM <- RunPCA(hCRC_TC_PTvLM, npcs = 30, verbose = FALSE)
hCRC_TC_PTvLM <- RunUMAP(hCRC_TC_PTvLM, reduction = "pca", dims = 1:30)
hCRC_TC_PTvLM <- FindNeighbors(hCRC_TC_PTvLM, reduction = "pca", dims = 1:30)
hCRC_TC_PTvLM <- FindClusters(hCRC_TC_PTvLM, resolution = 0.5)

DimPlot(hCRC_TC_PTvLM, reduction = "umap", label=T, group.by = "Cell.type")
DimPlot(hCRC_TC_PTvLM, reduction = "umap", label=T)
DimPlot(hCRC_TC_PTvLM, reduction = "umap", label=T, group.by = "Tumor.type")

saveRDS(hCRC_TC_PTvLM, file="hCRC_TC_PTvLM.rds")
#data visualisation
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

library(org.Hs.eg.db)
library(clusterProfiler)

plot_stat(hCRC_TC_PTvLM, plot_type = "group_count", group_by="Tumor.type")
plot_stat(hCRC_TC_PTvLM, plot_type = "prop_fill", group_by = "Tumor.type" )
plot_stat(hCRC_TC_PTvLM, plot_type = "prop_multi", group_by = "Tumor.type")

#finding markers per  tumor type
hCRC_TC_PTvLM <- SetIdent(hCRC_TC_PTvLM, value = hCRC_TC_PTvLM@meta.data$Cell.type)
PTCD8.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "PT_CD8")
LMCD8.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "LM_CD8")

PTCD4.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "PT_CD4")
LMCD4.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "LM_CD4")

CD4.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "LM_CD4", ident.2 = "PT_CD4")
CD8.markers <- FindMarkers(hCRC_TC_PTvLM, ident.1 = "LM_CD8", ident.2 = "PT_CD8")


library(EnhancedVolcano)

EnhancedVolcano(CD4.markers,
                lab = rownames(CD4.markers),
                x = 'avg_log2FC',
                y = 'p_val', FCcutoff = 0.5, title = "CD4 Tcells LM v PT")


EnhancedVolcano(CD8.markers,
                lab = rownames(CD8.markers),
                x = 'avg_log2FC',
                y = 'p_val',  FCcutoff = 0.3,title = "CD8 Tcells LM v PT")

#goanalysis
library(Scillus)
library(tidyverse)
library(Seurat)
library(magrittr)

library(org.Hs.eg.db)
library(clusterProfiler)


#DGE and GESA
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DESeq2)

CD4.markers <- tibble::rownames_to_column(CD4.markers, "Gene") # Apply rownames_to_column

CD4.gene.list <- CD4.markers$avg_log2FC
names(CD4.gene.list) <- CD4.markers$Gene
CD4.gene.list<-na.omit(CD4.gene.list)

CD4.gene.list = sort(CD4.gene.list, decreasing = TRUE)

gse <- gseGO(geneList=CD4.gene.list, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 1, 
             maxGSSize = 200, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

CD8.markers <- tibble::rownames_to_column(CD8.markers, "Gene") # Apply rownames_to_column

CD8.gene.list <- CD8.markers$avg_log2FC
names(CD8.gene.list) <- CD8.markers$Gene
CD8.gene.list<-na.omit(CD8.gene.list)

CD8.gene.list = sort(CD8.gene.list, decreasing = TRUE)

gseCD8 <- gseGO(geneList=CD8.gene.list, 
             ont ="ALL", 
             keyType = "SYMBOL",
             minGSSize = 1, 
             maxGSSize = 29, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=15)
dotplot(gseCD8, showCategory=15)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

d <- godata('org.Hs.eg.db', ont="CC")
CD4 <- pairwise_termsim(gse, method="Wang", semData = d)
install.packages("ggnewscale")
emapplot(CD4, showCategory = 15)

CD8 <- pairwise_termsim(gseCD8, method="Wang", semData = d)
emapplot(CD8, showCategory = 15)

emapplot_cluster(CD4)
emapplot_cluster(CD8)

cnetplot(gse, categorySize="pvalue", foldChange=CD4.gene.list)

cnetplot(gseCD8, categorySize="pvalue", foldChange=CD8.gene.list)

#KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
CD4ids<-bitr(names(CD4.gene.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# remove duplicate IDS 
CD4dedup_ids = CD4ids[!duplicated(CD4ids[c("SYMBOL")]),]

# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CD4 <- CD4.markers
KEGG.CD4$Y = CD4dedup_ids$ENTREZID


CD4dedup_ids <- CD4dedup_ids %>% 
  rename(
   Gene = SYMBOL)

#identify unmatching rows
library(dplyr)
anti_join(KEGG.CD4, CD4dedup_ids, by="Gene")

#Remove rows where Gene is not equal to genes not present in CD4dedup_ids
KEGG.CD4.filtered <- subset(KEGG.CD4,Gene !='MARCH1' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='H1FX' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='MT-ND3' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='MT-ND1' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='H3F3A' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='HIST1H4C' )
KEGG.CD4.filtered <- subset(KEGG.CD4.filtered,Gene !='H2AFZ' )

# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CD4.filtered$Y = CD4dedup_ids$ENTREZID


## Create a vector of the gene unuiverse
KEGG.CD4.gene.list <- KEGG.CD4.filtered$avg_log2FC

# Name vector with ENTREZ ids
names(KEGG.CD4.gene.list) <- KEGG.CD4.filtered$Y

# omit any NA values 
KEGG.CD4.gene.list<-na.omit(KEGG.CD4.gene.list)

# sort the list in decreasing order (required for clusterProfiler)
KEGG.CD4.gene.list = sort(KEGG.CD4.gene.list, decreasing = TRUE)

#Create gseKEGG object
kegg_organism = "hsa"
kk2.CD4 <- gseKEGG(geneList     = KEGG.CD4.gene.list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2.CD4, title = "CD4 T cells Enriched Pathways")

library(pathview)
install.packages("pathview")
devtools::install_github("javadnoorb/pathview", force=TRUE)
# Produce the native KEGG plot (PNG)
dme.cd4 <- pathview(gene.data=KEGG.CD4.gene.list, pathway.id="hsa04722", species = "hsa")

knitr::include_graphics("hsa04722.pathview.png")
dme.cd4

###Repeat KEGG analysis for CD8
#KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
CD8ids<-bitr(names(CD8.gene.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# remove duplicate IDS 
CD8dedup_ids = CD8ids[!duplicated(CD8ids[c("SYMBOL")]),]

# Create a new column in df2 with the corresponding ENTREZ IDs
KEGG.CD8 <- CD8.markers
KEGG.CD8$Y = CD8dedup_ids$ENTREZID


CD8dedup_ids <- CD8dedup_ids %>% 
  rename(
    Gene = SYMBOL)

#identify unmatching rows
library(dplyr)
anti_join(KEGG.CD8, CD8dedup_ids, by="Gene")

## Create a vector of the gene unuiverse
KEGG.CD8.gene.list <- KEGG.CD8$avg_log2FC

# Name vector with ENTREZ ids
names(KEGG.CD8.gene.list) <- KEGG.CD8$Y

# omit any NA values 
KEGG.CD8.gene.list<-na.omit(KEGG.CD8.gene.list)

# sort the list in decreasing order (required for clusterProfiler)
KEGG.CD8.gene.list = sort(KEGG.CD8.gene.list, decreasing = TRUE)

#Create gseKEGG object
kegg_organism = "hsa"
kk2.CD8 <- gseKEGG(geneList     = KEGG.CD8.gene.list,
                   organism     = kegg_organism,
                   nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 100,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "none",
                   keyType       = "ncbi-geneid")


