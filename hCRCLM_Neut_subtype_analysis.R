##Neutrophil subtypes in CRCLM

DimPlot(hCRC_neut_s, reduction = "umap", group.by = "seurat_clusters", label = T)

Idents(object = hCRC_neut_s) <- "seurat_clusters"
DimPlot(hCRC_neut_s, reduction = "umap", label = T)

#data visualisation
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)
plot_stat(hCRC_neut_s, plot_type = "group_count", group_by="seurat_clusters")
#plot_stat(hCRC_neut_s, plot_type = "prop_fill",, group_by="seurat_clusters" )
#plot_stat(hCRC_neut_s, plot_type = "prop_multi", group_by = "seurat_clusters")


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

DimPlot(hCRC_neut_s, reduction = "umap", label = T, label.size = 5)
Idents(hCRC_neut_s)
hCRC_neut_s$Cell.type <- Idents(hCRC_neut_s)

plot_stat(hCRC_neut_s, plot_type = "group_count", group_by="Cell.type")

FeaturePlot(hCRCLM_TC, features = "SPP1", max.cutoff = 5)
DimPlot(hCRCLM_TC, reduction = "umap", label = T, label.size = 5)
