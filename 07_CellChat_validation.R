#Lucile Massenet-Regad 
#Created 2023-04-23


rm(list=ls())
setwd(dir = work.dir)

output.dir="analyses/3_Downstream_analysis/"
dir.create(path = paste0(work.dir,output.dir,'8_CellChat'), recursive = TRUE, showWarnings = FALSE)

#libraries
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)


# Load seurat object, retrieve matrix, average expression for each cluster
seurat=readRDS("data/ccRCC_nCount1000_integrated_Harmony_res1_PC50_2022-02-08.rds")
Idents(seurat)=seurat$Celltype_Harmony2

seurat.tum = subset(seurat, cells=which(seurat$Tissue=="Tum"))
data=seurat.tum@assays$RNA@data
meta.data= seurat.tum@meta.data

# Create CellChat object
cellchat <- createCellChat(object = data, meta = meta.data, group.by = "Celltype_Harmony2")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 

# CellChat DB
CellChatDB <- CellChatDB.human # see how add specific interactions

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- CellChat::subsetData(cellchat, features = NULL) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# compute the communication probability and infer cellular communication networ
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat) #df at the level of LR . precise slot.name= "netP" to access at level of signaling pathway.


## Systems analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat2)



netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", cluster.rows = T,  width = 20 , height = 20)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",  width = 20 , height = 20)

interac_oi_OUT = data.frame("interaction_name"=c("EDN1_EDNRA", "EDN1_EDNRB", "CD70_CD27", "SPP1_CD44", "SPP1_ITGA4_ITGB1", "SPP1_ITGA5_ITGB1", 
                                            "SPP1_ITGA8_ITGB1", "SPP1_ITGA9_ITGB1", "SPP1_ITGAV_ITGB1", "SPP1_ITGAV_ITGB5", "CCL28_CCR10", 
                                            "RARRES2_CMKLR1", "RARRES2_GPR1", "VEGFA_VEGFR1", "VEGFA_VEGFR2","VEGFA_VEGFR1R2","VEGFA_VEGFR3",  "TGFA_EGFR_ERBB2", "TGFA_EGFR"))
interac_oi_IN =  data.frame("interaction_name"=c("TGFA_EGFR",
                                            "AREG_EGFR", "EREG_EGFR", "HBEGF_EGFR", "HGF_MET"))
netVisual_bubble(cellchat_10,  pairLR.use = interac_oi_OUT, 
                 remove.isolate = FALSE, sort.by.source = T, sort.by.source.priority = F, thresh = 0.01)
netVisual_bubble(cellchat_10,  pairLR.use = interac_oi_IN, 
                 remove.isolate = FALSE, sort.by.target = F, sort.by.source.priority = T, thresh = 0.01)

# Specific communication
saveRDS(cellchat, file = "analyses/3_Downstream_analysis/8_CellChat/cellchat_LM_data_0.10trim.rds")

