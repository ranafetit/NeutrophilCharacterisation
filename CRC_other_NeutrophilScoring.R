library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyverse)
library(ggplot2)

#scoring tumor-specific neut signature derived from KPN and L_AC on other CRC models
#R18_No.KPN is the seurat object from CRC_other_counts on Zenodo

# list genes enriched in tumor-derived Neuts
t_enriched_a <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1","Ifit1","Il1b","Ccl3","Ptma")

#Use this list of 20 genes to score cells using the AddModuleScore function:
R18_No.KPN <- AddModuleScore(R18_No.KPN,
                             features = list(t_enriched_a),
                             name="T_enriched", ctrl=50)
library(RColorBrewer)
# Plot scores
FeaturePlot(R18_No.KPN,
            features = "T_enriched1", ) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#healthy neut signature
h_enriched_a <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3","G0s2","Fpr2")

#Use this list of 20 genes to score cells using the AddModuleScore function:
R18_No.KPN <- AddModuleScore(R18_No.KPN,
                             features = list(h_enriched_a),
                             name="H_enriched", ctrl=50)

# Plot scores
library(RColorBrewer)
FeaturePlot(R18_No.KPN,
            features = "H_enriched1", split.by = "Genotype") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


VlnPlot(R18_No.KPN,features = c("T_enriched1"), group.by = "Genotype") 
VlnPlot(R18_No.KPN,features = c("h_enriched1"), group.by = "Genotype") 

#cellcycle scoring
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

R18_No.KPN <- CellCycleScoring(R18_No.KPN, s.features = s.genes, g2m.features = g2m.genes)
table(R18_No.KPN[[]]$Phase)

VlnPlot(R18_No.KPN,features = c("S.Score","G2M.Score"))
VlnPlot(R18_No.KPN,features = c("S.Score"), group.by = "Genotype")
VlnPlot(R18_No.KPN,features = c("G2M.Score"), group.by = "Genotype")


