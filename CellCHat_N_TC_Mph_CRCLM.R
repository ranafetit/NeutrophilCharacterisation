#Merge Neutrophil, Tcells and Mph in CRCLM to run CellChat

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
DimPlot(CRCLM_N_TC_Mph, reduction = "umap", label = T)

DimPlot(CRCLM_N_TC_Mph, reduction = "umap", label = T, group.by = "Cell.type", repel = T, label.box = T)
DimPlot(CRCLM_N_TC_Mph, reduction = "umap", label = T, split.by  = "active.ident", repel = T, label.box = T)

saveRDS(CRCLM_N_TC_Mph, file="CRCLM_N_TC_Mph.rds")

devtools::install_github("sqjin/CellChat")

library(CellChat)

# Create cellchat object
cellchat <- createCellChat(object = CRCLM_N_TC_Mph, group.by = "Cell.type", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, idents.use = c("SPP1+Mph"), arrow.size = 1)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, sources.use = c("Activated_Neut"), arrow.size = 1)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, targets.use = c("CD8 T cell"), arrow.size = 1)

#visualisation of the different pathways & interactions
cellchat@netP[["pathways"]]

#pathways.show <- c("MHC-I", "MHC-II","MIF","CD45" ,"GALECTIN", "CLEC","CXCL" , "ANNEXIN","ADGRE5","CCL","CD99", "ITGB2","LCK","IL1","TNF","VCAM","CD86","ICAM", "IFN-II","PECAM1") 
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show , layout = "circle")

vertex.receiver = seq(1,2)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
#and visualize cell-cell communication mediated by a single ligand-receptor pair

netAnalysis_contribution(cellchat, signaling = pathways.show)

#visualize the cell-cell communication mediated by a single ligand-receptor pair.
#extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.

pathways.show <- c("MHC-II") 
pairLR.MHC2<- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MHC2[1:5,] # show one ligand-receptor pair

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

pathways.show <- c("ADGRE5") 
pairLR.ADGRE5<- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ADGRE5[1,] # show one ligand-receptor pair

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:19), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 19, targets.use = c(1:18), remove.isolate = FALSE)


plotGeneExpression(cellchat, signaling = "MIF")
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("TGFb") 
plotGeneExpression(cellchat, signaling = "TGFb")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb"))
gg2

#which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL","IL1","TNF","IFN-II"),pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL","IL1","TNF","IFN-II"),pattern = "incoming")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL","IL1","TNF"),pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL","IL1","TNF"),pattern = "incoming")


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-1","CLEC","CD99","LCK","MIF"),pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-1","CLEC","CD99","LCK","MIF"),pattern = "incoming")

ht1+ht2

#how multiple cell groups and signaling pathways coordinate to function?
library(NMF)
library(ggalluvial)

#un selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")

#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 7

nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#incoming patterns
selectK(cellchat, pattern = "incoming")
nPatterns = 8
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")

#Identify signaling groups based on their functional similarity
reticulate::py_install(envname="Renv", packages ='umap-learn')
install.packages("uwot")
library(uwot)

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#reticulate::py_install(packages = 'umap-learn')
#reticulate::use_python("/Users/USERS/opt/miniconda3/bin/python", required=T)


cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
#saveRDS(cellchat, file = "cellchat_hCRCLM_N_Tc.rds")


DotPlot(object = CRCLM_N_TC_Mph, features = "SPP1", group.by = "Cell.type")
Idents(CRCLM_N_TC_Mph)
