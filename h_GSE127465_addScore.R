#scoring neutrophil signature in human neutrophil dataset of lung tumor

library(tidyverse)
library(RColorBrewer)
library(Seurat)

Idents(hNeut_GSE127465_s)

#assign new identities to differentiate bw neuts from blood and tumor.

hNeut_GSE127465_s <- SetIdent(hNeut_GSE127465_s, value = hNeut_GSE127465_s@meta.data$Tissue.idents)
DimPlot(hNeut_GSE127465_s)


#scoring for established gene signature in Tumor specific 
t_enriched <- c("Wfdc17","Cdkn1a","Ppia","Ifitm1","Tagln2","Ifit3","Slfn4","Gngt2","Ifitm6","Isg15","Cxcl2","Ccl4","Cd14","Rps27l","Thbs1")

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","GNGT2","THBS1","PTMA")

t_enriched <- c("CDKN1A","PPIA","IFITM1","TAGLN2","ISG15","CXCL8","CCL4","CD14","RPS27L","IER3","CCL3","IFIT3","IFIT1","IL1B","WFDC1","THBS1","PTMA")


# list genes enriched in tumor-derived Neuts
tcommon_enriched <- c("Cxcl2","Thbs1","Ccl4","Ccl3","Cd14","Cdkn1a","Ppia","Gngt2","Ier3","Rps27l","Ptma")



FeaturePlot(hNeut_GSE127465_s, features = "SIGLECF")

#Use this list of 20 genes to score cells using the AddModuleScore function:
hNeut_GSE127465_s <- AddModuleScore(hNeut_GSE127465_s,
                                   features = list(t_enriched),
                                   name="T_enriched")

# Plot scores
library(RColorBrewer)
FeaturePlot(hNeut_GSE127465_s,
            features = "T_enriched1",max.cutoff = 15) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



#healthy neut signature
h_enriched <- c("Mmp8","Ifitm6","S100a6","Lyz1","Lyz2","Ctla2a","Chil3")
h_enriched <- c("MMP8","IFITM2","IFITM3","S100A6","LYZ","CTLA4","CHI3L1","G0S2","FPR2")
#NO IFITM6 IN HUMANS
#NO CTLA2A IN HUMNAS, ALTERNATIVELY CTLA4
#add some netosis genes
#By contrast, rodents lack a direct homologue of IL-8, but the chemokines CXCL1/KC, CXCL2/MIP-2, and CXCL5-6/LIX are regarded as functional homologues of IL-8 

#Use this list of 20 genes to score cells using the AddModuleScore function:
hNeut_GSE127465_s <- AddModuleScore(hNeut_GSE127465_s,
                                   features = list(h_enriched),
                                   name="h_enriched")

# Plot scores
library(RColorBrewer)
FeaturePlot(hNeut_GSE127465_s,
            features = "h_enriched1", max.cutoff = 10) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(hNeut_GSE127465_s,features = c("T_enriched1"), group.by = "Tissue.idents") 
VlnPlot(hNeut_GSE127465_s,features = c("h_enriched1"), group.by = "Tissue.idents") 
table(hNeut_GSE127465_s$T_enriched1)

#cellcycle scoring
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

hNeut_GSE127465_s <- CellCycleScoring(hNeut_GSE127465_s, s.features = s.genes, g2m.features = g2m.genes)
table(hNeut_GSE127465_s[[]]$Phase)

VlnPlot(hNeut_GSE127465_s,features = c("S.Score","G2M.Score"))
VlnPlot(hNeut_GSE127465_s,features = c("S.Score"), group.by = "Tissue.idents")
VlnPlot(hNeut_GSE127465_s,features = c("G2M.Score"), group.by = "Tissue.idents")

hNSCLCPT_neut<-subset(x = hNeut_GSE127465_s, idents = "tumor")
write_rds(hNSCLCPT_neut, file="hNSCLCPT_neut.rds")
