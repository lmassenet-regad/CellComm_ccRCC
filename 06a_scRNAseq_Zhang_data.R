#Lucile Massenet-Regad - scRNAseq ccRCC PhD project
#Created 2021-09-17, last modified: 2021-11-23
#Integration - Comparison cohort to LM022 LML027 and Young dataset

rm(list=ls())
setwd(dir = "~/Documents/PhD_TUMOR_ccRCC/BIOINFO/Zhang et al, 2021, PNAS/")

library(clustree)
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)

# # files <- list.files("GSE159115_csv/")
# # meta.all=readxl::read_excel("Metadata.xlsx")
# # 
# # Seurat=list()
# # for(i in 1:11){
# #   #i=1
# #   id=substr(files[i],12, 19)
# #   data=t(read.csv(paste0("GSE159115_csv/", files[i]), header=T))
# #   colnames(data)=data[1,]
# #   data=data[,-1]
# #   seurat <- CreateSeuratObject(counts = data , project = id )
# #   seurat <- NormalizeData(seurat)
# #   seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
# #   seurat
# # }             
# 
# #saveRDS(Seurat, "Seurat_list_object.rds")

#FUNCTION NEEDED FOR PREPROCESSING  
seurat.preprocess <- function(seuratOb=seuratOb, the.seed=1337L, n.dim=50){
  seurat=seuratOb
  # QC
  pdf(paste0("analyses/QC_",as.character(unique(seurat$orig.ident)), "_Preprocessing_individual_seurat_objects_Zhang.pdf"))
  seurat[["percent.mt"]]=PercentageFeatureSet(seurat, pattern="^MT.")
  pa <- VlnPlot(seurat, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) 
  plot(pa)
  pb <- FeatureScatter(seurat, feature1="nCount_RNA", feature2="nFeature_RNA")
  plot(pb)
  metadata <- as.data.frame(seurat@meta.data)
  p1 <- metadata %>%  ggplot(aes(x=nCount_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    xlab("Cell UMIs") +
    geom_vline(xintercept = 1000)+
    annotate("text", x=500, y=1, label= paste0(round((sum(metadata$nCount<1000)/dim(metadata)[1])*100,0), ' %') ) +
    ggplot2::ggtitle(as.character(unique(seurat$orig.ident)))
  plot(p1)
  
  p2 <- metadata %>%  ggplot(aes(x=nFeature_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    xlab("Cell Features") +
    geom_vline(xintercept = 200)+
    annotate("text", x=500, y=1, label= paste0(round((sum(metadata$nFeature_RNA<200)/dim(metadata)[1])*100,0), ' %') ) +
    ggplot2::ggtitle(as.character(unique(seurat$orig.ident)))
  plot(p2)
  
  # subset according to QC
  seurat=subset(seurat, subset=nCount_RNA >1000 & nFeature_RNA>200 & percent.mt <20)
  
  # cell cycle QC
  s.genes <- Seurat::cc.genes$s.genes
  s.genes <- s.genes[s.genes %in% rownames(seurat)] # genes in dataset
  
  g2m.genes <- Seurat::cc.genes$g2m.genes
  g2m.genes <- g2m.genes[g2m.genes %in% rownames(seurat)] # genes in dataset
  
  # compute the scoring
  seurat<- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, assay = "RNA", nbin = 10, seed = the.seed)
  seurat$Seurat.SmG2M.Score <- seurat$S.Score - seurat$G2M.Score
  
  
  # doublet prediction
  seurat$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(seurat)))$scDblFinder.class
  seurat$scDblFinder.class <- unname(seurat$scDblFinder.class == "doublet")
  print('scDblFinder doublets :')
  print(table(seurat$scDblFinder.class))
  
  # scds
  seurat$scds.score <- scds::cxds_bcds_hybrid(Seurat::as.SingleCellExperiment(seurat))$hybrid_score
  seurat$scds.class <- unname(seurat$scds.score > 1)
  print('scds-hybrid doublets :')
  print(table(seurat$scds.class))
  
  # union of scDblFinder & scds & Manual incompatibility
  seurat$doublets_consensus.class <- seurat$scDblFinder.class | seurat$scds.class #| seurat$ManDblt
  print('Consensus doublets :')
  print(table(seurat$doublets_consensus.class))
  
  df_S.cycle_doublets <- data.frame(Phase = seurat$Phase, doublets = seurat$doublets_consensus.class)
  print(table(df_S.cycle_doublets))
  print(paste0(round((nrow(df_S.cycle_doublets[df_S.cycle_doublets$Phase == "G2M" &
                                                 df_S.cycle_doublets$doublets == TRUE,])/
                        nrow(df_S.cycle_doublets[df_S.cycle_doublets$doublets == TRUE,])*100),2),
               "% of doublets are in G2M phase by Seurat."))
  
  # plot PCA ## Careful, not good way to normalise +++``
  seurat.tmp = seurat #  temporary object only to visualise how doublet and cell cycling are.
  seurat.tmp <- NormalizeData(seurat.tmp, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.tmp <- FindVariableFeatures(seurat.tmp, selection.method = "vst", nfeatures = 2000)
  seurat.tmp <- ScaleData(seurat.tmp)
  seurat.tmp <- RunPCA(seurat.tmp, npcs = 100, verbose = FALSE)
  ElbowPlot(seurat.tmp, ndims = 100) 
  seurat.tmp <- RunUMAP(seurat.tmp, dims = 1:n.dim, seed.use = the.seed)
  
  plot(DimPlot(seurat.tmp, reduction = "umap", group.by = "Phase") + ggplot2::ggtitle("Cell Phase (Seurat)"))
  plot(DimPlot(seurat.tmp, group.by = "scDblFinder.class") + ggplot2::ggtitle("Cell doublets (scDblFinder)"))
  plot(DimPlot(seurat.tmp, group.by = "scds.class") + ggplot2::ggtitle("Cell doublets (scds)")) 
  plot(DimPlot(seurat.tmp, group.by = "doublets_consensus.class") + ggplot2::ggtitle("Cell doublets (union)")) 
  dev.off()
  
  # Filtering cell doublets
  seurat <- seurat[, !seurat$doublets_consensus.class]
  seurat= NormalizeData(seurat)
  seurat=FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  
  return(seurat)
} 


## 1 :  QC and processing individual datasets

Seurat=readRDS("Seurat_list_object.rds")

for(i in 1:length(Seurat)){
  seurat2= seurat.preprocess(Seurat[[i]])
  Seurat[[i]]=seurat2
}


## 2 - Integration using Seurat CCA

kidney.anchors <- FindIntegrationAnchors(object.list = Seurat, dims=1:50)
kidney.combined <- IntegrateData(anchorset = kidney.anchors, dims = 1:50)
saveRDS(kidney.combined, "Seurat_Integrated_20220217.rds")

rm(kidney.anchors)
rm(Seurat)
#rm(data)

# Apply Seurat pipeline
#===================================

kidney.combined = readRDS( "Seurat_Integrated_20220217.rds")
kidney.combined <- ScaleData(kidney.combined, verbose = FALSE)
kidney.combined <- RunPCA(kidney.combined, npcs = 100, verbose = FALSE)

DefaultAssay(kidney.combined) <- "integrated"

pdf(paste0("analyses/clustering/Clustering_CCA_Zhang_",Sys.Date(),".pdf"),useDingbats = F)
VizDimLoadings(object = kidney.combined, dims = 1:2)
DimHeatmap(object = kidney.combined, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(object =kidney.combined,ndims = 100)  #50 enough to explain much variance 
dev.off()

kidney.combined <- FindNeighbors(object = kidney.combined, dims = 1:50)


# Clustree - 2021-09-20
# ==============================================================================
#Silhouette
seurat.test=kidney.combined

for (res in seq(0.1,1.4,0.1)){
  #res=0.1
  seurat.test <- FindClusters(object = seurat.test, resolution = res)

# Clusttree
pdf(paste0("analyses/clustering/Clustering_CCA_Clustree_Zhang_PC50_",Sys.Date(),".pdf"),useDingbats = F,height = 15,width=10)
exprs <- "data"
prefix = "integrated_snn_res."
args <- list()
gene_names <- rownames(seurat.test@assays$RNA@data)
for (node_aes in c("node_colour", "node_size", "node_alpha")) {
  if (node_aes %in% names(args)) {
    node_aes_value <- args[[node_aes]]
    if (node_aes_value %in% gene_names) {
      aes_name <- paste0(exprs, "_", node_aes_value)
      seurat.test@meta.data[aes_name] <-
        slot(seurat.test, exprs)[node_aes_value, ]
      args[[node_aes]] <- aes_name
    }
  }
}
args$x <- seurat.test@meta.data
args$prefix <- prefix
do.call(clustree, args)
dev.off()
}

# 4- UMAP and clustering 
#====================================
DefaultAssay(kidney.combined)="integrated" #important for umap and findcluster if integration performed
kidney.combined <- RunUMAP(kidney.combined, dims = 1:50)
kidney.combined <- FindClusters(object = kidney.combined, resolution = 0.5)
DimPlot(kidney.combined, reduction = "umap", label=T)


kidney.combined$Tissue=NA
kidney.combined$Tissue[which(kidney.combined$orig.ident %in% c("SI_18854", "SI_18855", "SI_19703", "SI_22368", 
                                                               "SI_22604", "SI_23459", "SI_23843"))]="Tum"
kidney.combined$Tissue[which(kidney.combined$orig.ident %in% c("SI_18856", "SI_19704", "SI_22369", "SI_22605"))]="Hty"

kidney.combined$Patient=NA
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_18854"))]="p2005"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_18855","SI_18856"))]="p2006"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_19703","SI_19704"))]="p2007"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_22368", "SI_22369"))]="p2017"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_22604", "SI_22605"))]="p2022"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c("SI_23459"))]="p2023"
kidney.combined$Patient[which(kidney.combined$orig.ident %in% c( "SI_23843"))]="p2026"


# 5 - celltype annotation
#====================================
DefaultAssay(kidney.combined) <- "RNA"

seurat.markers<- FindAllMarkers(kidney.combined,only.pos = TRUE, min.diff.pct = 0.25, logfc.threshold = 0.5)
metascape <- seurat.markers %>% dplyr::filter(p_val_adj<=0.05 & avg_log2FC > 1)  %>% as.data.frame 
write.csv(metascape, paste0("analyses/clustering/Integration_CCA_Markers__Zhang_res0.6_",Sys.Date(),".csv"))

kidney.combined$Celltype_CCA = NA

kidney.combined$Celltype_CCA[Idents(kidney.combined)=="0"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="1"] ="MMAC_1"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="2"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="3"] = "Endoth_PLVAP"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="4"] = "Fibro"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="5"] ="T_cells"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="6"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="7"] = "Endoth_ACKR1"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="8"] = "DC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="9"] = "Prolif"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="10"] =  "MMAC_2"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="11"] = "Endoth_PDGFRB"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="12"] = "PlasmaC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="13"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="14"] = "NK"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="15"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="16"] = "TumC"
kidney.combined$Celltype_CCA[Idents(kidney.combined)=="17"] = "Mast"

# UMAP and clustering graph + SAVE DATA 

pdf(file = paste0("analyses/clustering/CCA_UMAP_","res0.5","_PC50_",Sys.Date(),".pdf"), useDingbats = F, width = 10, height = 10 )
DimPlot(kidney.combined, reduction = "umap", label = T, label.size = 4) + labs(title= "res 0.5")
DimPlot(kidney.combined, reduction = "umap", group.by = "orig.ident", label = T, label.size = 3) 
DimPlot(kidney.combined, reduction = "umap", group.by = "Tissue", label = F, label.size = 3) 
DimPlot(kidney.combined, reduction = "umap", group.by = "Patient", label = F, label.size = 3) 

pt=as.data.frame(table(kidney.combined$Celltype_CCA, kidney.combined$Tissue))
ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) 

pt=as.data.frame(table(kidney.combined$Celltype_CCA, kidney.combined$Patient))
ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) 

DimPlot(kidney.combined, reduction = "umap", label=T)
DimPlot(kidney.combined , reduction = "umap", group.by = "integrated_snn_res.0.5")
DimPlot(kidney.combined , reduction = "umap", group.by = "Celltype_CCA", label=T)
dev.off()

saveRDS(kidney.combined, "Seurat_Integrated_20220705.rds")


