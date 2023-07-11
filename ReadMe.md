# Characterising neutrophil subtypes in cancer using human and murine single-cell RNA sequencing datasets

## Project description 

This repository contains all the Rscripts generated in analysing both publicly available and independently generated scRNAsequencing datasets of healthy neutrophils as well as neutrophils in lung, breast and colorectal cancer (CRC) to identify neutrophil transcriptomic subtypes and their developmental lineages in health and malignancy. 

Public datasets were retrieved from the GEOdatabase and National Omics Encyclopaedia

Accession numbers:

- GSE127465
- GSE165276
- GSE139125
- GSE114727
- OEP001756
- GSE146771

Independent datasets include neutrophils from transplant and genetically engineered mouse models of CRC recapitulating primary tumor setting. Metadata, counts and normalised counts are deposited on Zenodo. 

- CRC_KPN
- CRC_Other (AKPT, BP, BPN and KP)

## Workflow
Datasets were retrieved from the GEOdatabase and National Omics Encyclopaedia and processed using the package Seurat (version 4.3.0) on R (versions 3.17 and 4.1.1). Datasets were integrated by RPCA using the IntegrateData function then scaled and normalised. Dimension reduction was performed using PCA followed by clustering using the FindNeighbours and FindClusters functions. Marker genes for individual clusters were determined using the FindAllMarkers function and neutrophils were isolated by authors using the cluster identities and markers assigned in the original publications. Datasets were integrated to establish and test neutrophil gene-signatures using the AddModuleScore function. Pseudo-time analysis was performed using the Rpackage Slingshot (version 2.8.0) to identify neutrophil lineages. Gene expression along the different trajectories was performed using the Rpackage TradeSeq (version 1.14.0). Gene Set Enrichment (GSE), Gene Ontology (GO) and KEGG analyses were performed using the Rpackages ClusterProfiler (version 4.8.1) and EnrichR (version 3.2). Ligand-receptor (L-R) interactions and signalling pathways between neutrophils and other immune cell populations in primary and metastatic sites were investigated using the Rpackage CellChat (version 1.6.1). 

## Links for software processing pipelines on Github

| Software | Link |
| ------ | ------ |
| Seurat | [https://github.com/satijalab/seurat] |
| Slingshot | [https://github.com/kstreet13/slingshot]|
| TradeSeq | [https://github.com/statOmics/tradeSeq] |
| ClusterProfiler | [https://github.com/YuLab-SMU/clusterProfiler]|
| EnrichR | [https://github.com/wjawaid/enrichR] |
| CellChat | [https://github.com/sqjin/CellChat]|

## Rscripts in this repository
- [KPN_NT_NSCLC_int] - Integrating neutrophil datasets from CRC KPN, Neutrotime (NT) and non-small cell lung cancer (NSCLC). All mouse.
- [m_GSE139125_PYMT_slingshot] - Pseudotime analysis using slingshot on PYMT breast cancer dataset (mouse)
- [CRC_other_slingshot] - Pseudotime analysis using slingshot on other CRC genotypes.
- [CRC_other_NeutrophilScoring] - Scoring neutrophil signatures in other CRC genotypes.
- [CellChat_hCRC_Tc_PT] - Ligand-Receptor interactions and cell communication analysis using CellChat in T-cells (Tc) of primary tumor (PT) datasets of human CRC
- [CellCHat_N_TC_CRCLM] - Ligand-Receptor interactions and cell communication analysis using CellChat between T-cells (Tc) and Neutrophils (N) in human liver metastatic CRC tissue (CRCLM)
- [GSE146771_CRCPT] - Analysis of Neutrophils in primary tumor CRC tissue.
- [h_CRC_TC_PTvLM] - Comparison of Tcells (TC) between human CRC primary tumor and liver metastasis (LM)
- [h_CRCLM_analysis] - Analysis of Neutrophils in metastatic CRC dataset.
- [h_GSE114727_analysis] - Analysis of Neutrophils in primary tumor breast cancer GSE114727 dataset.
- [h_GSE114727_Neut_clusters_ss] - Analysis of Neutrophils in primary tumor breast cancer GSE114727 dataset subsetted by neutrophil clusters identified in original paper.
- [h_GSE127465_addScore] - Scoring neutrophil signatures in human Primary tumor lung cancer GSE127465 dataset.
- [h_GSE127465_analysis] - Analysis of Neutrophils in human Primary tumor lung cancer GSE127465 dataset to identify and subset neutrophils.
- [h_GSE127465_pseudotime] - Pseudotime analysis using slingshot of Neutrophils in human Primary tumor lung cancer GSE127465 dataset
- [hCRCLM_Mph] - Analysis of Macrophages in metastatic CRC dataset.
- [hCRCLM_Neut_subtype_analysis] - Characterising neutrophil subtypes in metastatic CRC dataset.
- [m_GSE139125_PYMT] - Analysis of mouse PYMT breast cancer dataset to identify and subset neutrophils.
- [m_GSE139125_PYMT_slingshot] - Pseudotime analysis using slingshot of neutrophils in mouse PYMT breast cancer dataset
- [m_GSE165276_NT] - Analysis of neutrophils from healthy Blood, bone marrow and spleen in Neutrotime dataset
- [m_NT_NSLC_integrated] - Integrating neutrophils from Neutrotime (NT) dataset and non-small cell lung cancer (NSCLC) dataset; both mouse.
- [m_NT_slingshot] - Pseudotime analysis using slingshot of healthy neutrophils in mouse neutrotime dataset.

## Credits

This work has been done by Rana Fetit, postdoctoral researcher in the Steele Lab (M48) at the Beatson Institute for Cancer Research, Glasgow.

## Corresponding authors
Rana Fetit, r.fetit@beatson.gla.ac.uk; rana.fetit@gmail.com
Colin W. Steele, colin.steele@glasgow.ac.uk

## License

MIT


