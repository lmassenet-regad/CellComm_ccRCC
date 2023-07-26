#2021-06-04 Young et al. - last modified 2/2/22
# Analysis - preprocessing, integration, clustering

rm(list=ls())

# Load libraries
#===================================
library(Seurat)
library(clustree)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

work.dir= "~/Documents/PhD_TUMOR_ccRCC/BIOINFO/Young et al_2018_Kidney/"
dir.create(paste0(work.dir,"/analyses/1_QC&Integration_Harmony/"), recursive = TRUE, showWarnings = FALSE)
output.dir=paste0(work.dir, "analyses/1_QC&Integration_Harmony/")
input.dir=paste0(work.dir, "data/")
setwd(dir = "~/Documents/PhD_TUMOR_ccRCC/BIOINFO/Young et al_2018_Kidney/")

#Parameters 
the.seed<-1337L
min.counts=500
min.features=200
n.dims=30
pct.mito=20
min.cells=2

##==========================================================
#### Change metadata files from Seurat object - redo everything
#==========================================================

#Load supplementary data
raw.data=Read10X(data.dir =paste0(input.dir,"NIHMS78997-supplement-Supplementary_data")) #Warning message indicates that the last line of the file doesn't end with an End Of Line (EOL) character 
raw.data@Dimnames[[2]]=sapply(strsplit(as.character(raw.data@Dimnames[[2]]), "___", fixed=T), "[[", 1)
# Load metadata from supp tables + cluster info (supp table)
metadata1=as.data.frame(readxl::read_excel(paste0(input.dir,"SUPP_TABLES_aat1699-Young-TablesS1-S12-revision1.xlsx"), sheet="TableS6 - Sample manifest"))
colnames(metadata1)=metadata1[1,]
metadata1=metadata1[-1,]
metadata2=as.data.frame(readxl::read_excel(paste0(input.dir,"SUPP_TABLES_aat1699-Young-TablesS1-S12-revision1.xlsx"),sheet="TableS11 - Cell manifest"))
colnames(metadata2)=metadata2[1,]
metadata2=metadata2[-1,]
metadata=dplyr::right_join(metadata2, metadata1, by="SangerID")
cluster_info=as.data.frame(readxl::read_excel(paste0(input.dir,"SUPP_TABLES_aat1699-Young-TablesS1-S12-revision1.xlsx"), sheet="TableS2 - Cluster info"))
colnames(cluster_info)=cluster_info[1,]
cluster_info=cluster_info[-1,]
colnames(cluster_info)[1]="ClusterID"
metadata=dplyr::right_join(cluster_info[,c(1,8,9)], metadata, by="ClusterID")

# #check that order is the same / make unique barcode in the metadata and in Seurat object
all(metadata$barcode%in%raw.data@Dimnames[[2]])
pos=match(raw.data@Dimnames[[2]],metadata$barcode)
metadata=metadata[pos,]
all(metadata$barcode==raw.data@Dimnames[[2]])

raw.data@Dimnames[[2]]=make.unique(metadata$DropletID) #to have unique

QC.passed=metadata$DropletID[which(metadata$Experiment%in%c("RCC1","RCC2", "VHL_RCC") & 
                                     metadata$Organ=="Kidney")] #Take only cell of interest

data=raw.data[,QC.passed]
metadata=dplyr::filter(metadata, metadata$DropletID%in%QC.passed)
rownames(metadata)=make.unique(metadata$DropletID)


##==========================================================
#### Create Seurat Object - QC
#==========================================================

seurat=CreateSeuratObject(data, assay = "RNA", min.cells = min.cells, min.features = 0, meta.data = metadata)
table(Idents(seurat))

##### Visualize QC metrics as a violin plot
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[LS][[:digit:]]") #ribosomal genes

pdf(paste0(output.dir, "01_QC_Simple_metrics_beforeFiltering_Young.pdf"), useDingbats = F)
VlnPlot(seurat, features = c("nCount_RNA"), ncol = 1, pt.size=0.00) + NoLegend() + scale_y_log10() +  geom_hline(yintercept = min.counts)
VlnPlot(seurat, features = c("nFeature_RNA"), ncol = 1, pt.size=0.00) + NoLegend() + scale_y_log10() +  geom_hline(yintercept = min.features)
VlnPlot(seurat, features = c("percent.mt"), ncol = 1, pt.size=0.00) + NoLegend() +  geom_hline(yintercept = pct.mito)
VlnPlot(seurat, features = c("percent.rb"), ncol = 1, pt.size=0.00) + NoLegend() 
FeatureScatter(seurat, feature1 = "percent.rb", feature2 = "percent.mt") + NoLegend()

FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QCpass") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)

metadata3 <- as.data.frame(seurat@meta.data)
metadata3 %>%  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = min.counts)+
  annotate("text", x=100, y=1, label= paste0("<min.counts =", round((sum(metadata3$nCount_RNA<min.counts)/dim(metadata3)[1])*100,0), ' %') )
dev.off()

seurat <- subset(seurat, subset = nFeature_RNA >min.features & nCount_RNA > min.counts & percent.mt < pct.mito)  

pdf(paste0(output.dir, "01_QC_Simple_metrics_afterFiltering_Young.pdf"), useDingbats = F)
VlnPlot(seurat, features = c("nCount_RNA"), ncol = 1, pt.size=0.00) + NoLegend() + scale_y_log10() +  geom_hline(yintercept = min.counts)
VlnPlot(seurat, features = c("nFeature_RNA"), ncol = 1, pt.size=0.00) + NoLegend() + scale_y_log10() +  geom_hline(yintercept = min.features)
VlnPlot(seurat, features = c("percent.mt"), ncol = 1, pt.size=0.00) + NoLegend() +  geom_hline(yintercept = pct.mito)
VlnPlot(seurat, features = c("percent.rb"), ncol = 1, pt.size=0.00) + NoLegend() 
FeatureScatter(seurat, feature1 = "percent.rb", feature2 = "percent.mt") + NoLegend()

FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "QCpass") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)

metadata3 <- as.data.frame(seurat@meta.data)
metadata3 %>%  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = min.counts)+
  annotate("text", x=100, y=1, label= paste0("<min.counts =", round((sum(metadata3$nCount_RNA<min.counts)/dim(metadata3)[1])*100,0), ' %') )
dev.off()

###### Cell cycle prediction  # Alternative method to Seurat for cell cycle prediction -> Cyclone
#------------------------------------------------------------------------
s.genes <- Seurat::cc.genes$s.genes
s.genes <- s.genes[s.genes %in% rownames(seurat)] # genes in dataset

g2m.genes <- Seurat::cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% rownames(seurat)] # genes in dataset

# compute the scoring
seurat<- CellCycleScoring(object = seurat, s.features = s.genes, g2m.features = g2m.genes, assay = "RNA", nbin = 10, seed = the.seed)
seurat$Seurat.SmG2M.Score <- seurat$S.Score - seurat$G2M.Score


###### Doublets predicion + check with cell cycle
#------------------------------------------------------------------------
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
seurat$doublets_consensus.class <- seurat$scDblFinder.class | seurat$scds.class
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
seurat.tmp <- RunUMAP(seurat.tmp, dims = 1:n.dims, seed.use = the.seed)

pdf(paste0(output.dir, "02_QC_CellCycle_Doublets_PCA_temp_Young_.pdf"), useDingbats = F)
DimPlot(seurat.tmp, reduction = "umap", group.by = "Phase") + ggplot2::ggtitle("Cell Phase (Seurat)")
DimPlot(seurat.tmp, group.by = "scDblFinder.class") + ggplot2::ggtitle("Cell doublets (scDblFinder)") 
DimPlot(seurat.tmp, group.by = "scds.class") + ggplot2::ggtitle("Cell doublets (scds)") 
DimPlot(seurat.tmp, group.by = "doublets_consensus.class") + ggplot2::ggtitle("Cell doublets (union)") 
DimPlot(seurat.tmp, group.by = "TissueDiseaseState") + ggplot2::ggtitle("Tissue state") 
DimPlot(seurat.tmp, group.by = "Experiment") + ggplot2::ggtitle("Patient - need to normalise") 

FeatureScatter(seurat.tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "doublets_consensus.class") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)+
  ggplot2::ggtitle("QC -metrics - Cell doublets (union)") 

FeatureScatter(seurat.tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "TissueDiseaseState") + 
  scale_y_log10() + scale_x_log10()+ geom_vline(xintercept = min.counts) + geom_hline(yintercept = min.features)+
  ggplot2::ggtitle("QC -metrics - Cell doublets (union)") 
dev.off()

# Filte ring cell doublets
seurat <- seurat[, !seurat$doublets_consensus.class]

##==========================================================
#### Normalize, Integration Harmony nd Dim reduction
# ==============================================================================

seurat$batch=paste0(seurat$Experiment,"_", seurat$TissueDiseaseState)

# Normalisation WITH HARMONY
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)

seurat <- RunPCA(seurat, npcs = n.dims+20, verbose = FALSE)
ElbowPlot(seurat, ndims = 100) #n=30 fine
seurat <- RunUMAP(seurat, dims = 1:n.dims, seed.use = the.seed)


p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = "Experiment", pt.size = .1) +
                    NoLegend() +
                    ggtitle("Before harmony"))
p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", pt.size = .1, group.by = "Experiment") + NoLegend()) + ggplot2::ggtitle("Before Integration")
p2 <- AugmentPlot(VlnPlot(object = seurat, features = "PC_1", group.by = "Experiment", pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))

seurat <- suppressWarnings(harmony::RunHarmony(object =seurat,  group.by.vars="Experiment",
                                               reduction = "pca", assay.use="RNA", dims.use=1:n.dims+20, plot_convergence = T, verbose = T))

p3 <- AugmentPlot(DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by ="Experiment") + NoLegend())+ ggplot2::ggtitle("After Integration")
p4 <- AugmentPlot(VlnPlot(object = seurat, features = "harmony_1", group.by = "Experiment", pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))


seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:n.dims, reduction.name = "harmony_umap")
p5 <- AugmentPlot(DimPlot(seurat, reduction = "harmony_umap", group.by = "Experiment", pt.size = .1) +
                    NoLegend() +
                    ggtitle("After harmony"))
pdf(paste0(output.dir, "03_Integration_Harmony_PC", n.dims+20,"_Young.pdf"), useDingbats = F, width = 10, height = 10)
ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4))
plot_grid(p0,p5)
dev.off()

# Apply Seurat pipeline
#===================================
pdf(paste0(output.dir,"03_Clustering_Young_",Sys.Date(),".pdf"),useDingbats = F)
VizDimLoadings(object = seurat, dims = 1:2,reduction = "harmony" )
DimHeatmap(object = seurat, dims = 1:15, cells = 500, balanced = TRUE, reduction = "harmony")
ElbowPlot(object =seurat, ndims = n.dims+20, reduction = "harmony")  
dev.off()

seurat <- FindNeighbors(object = seurat, dims = 1:n.dims, reduction = "harmony")

# 4- UMAP and clustering  - chose resolution
#====================================
res=0.7
seurat <- FindClusters(object = seurat, resolution = res)

#Plot and save UMAP
pdf(paste0(output.dir,"04_Clustering_UMAP_", res, "_",Sys.Date(),".pdf"),useDingbats = F, height = 10,width=10)
DimPlot(seurat, reduction = "harmony_umap", label=T) + ggplot2::ggtitle(paste0("resolution = ", res)) + NoLegend()
DimPlot(seurat, reduction = "harmony_umap", group.by="SangerID",label=F) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
DimPlot(seurat, reduction = "harmony_umap", group.by="Experiment",label=F) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
DimPlot(seurat, reduction = "harmony_umap", group.by="TissueDiseaseState",label=T, pt.size = 0.01 ) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
DimPlot(seurat, reduction = "harmony_umap", split.by="TissueDiseaseState",  group.by="RNA_snn_res.0.7",label=T, pt.size = 0.01 ) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
dev.off()

# 5 - DEG - celltype annotation 
#====================================
DefaultAssay(seurat) <- "RNA" # SCT or RNA
if (DefaultAssay(seurat)== "SCT"){
  seurat <- PrepSCTFindMarkers(seurat) # correct counts when multiple SCT slots
  seurat.markers<- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, assay = "SCT")
}else if (DefaultAssay(seurat)== "RNA"){
  seurat.markers<- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, assay = "RNA")
} else stop ("DefaultAssay(seurat) should be RNA or SCT")

metascape <- seurat.markers %>% dplyr::filter(p_val_adj<=0.05 & avg_log2FC > 0.5)  %>% as.data.frame 
write.csv(metascape, paste0(output.dir, "05_DEG_CCA_Markers_res",res,"_RNAassay",Sys.Date(),".csv"))

seurat.markers %>% dplyr::filter(p_val_adj<=0.05 & avg_log2FC > 1) %>% arrange(cluster) %>% group_by(cluster) %>% top_n(20, avg_log2FC)  %>% View()



genes.notIMM=unique(c("PTPRC", "CA9","NDUFA4L2", "SLC17A3","LDHA", "NNMT", "ALDOB", "SLC34A1", "PDZK1IP1", "ANPEP",
                      "SLC22A8", "SLC17A3", "SLC16A9", "SLC7A13", "SLC13A3","ACTA2", "RGS5", "NOTCH3", "MCAM", "PDGFRB", 
                      "TAGLN", "PODXL","PTPRO","PTPRQ", "PTGDS", "CLIC5", "F3", "DCN", "MME", 
                      "RNASE1","RAMP2", "ENG", "PECAM1", "PTPRB",  "PTPRB",  "VCAM1", "PLVAP", "SLC14A1", "SOST", "MEG3", "EHD3", "TGFBR2", 
                      "KNG1", "UMOD", "CLDN16", "PCP4", "DEFB1","DUSP9", "SLC12A1", "ATP6V1G3","ATP6V0D2","EPCAM","LGALS3","TMEM213", "SLC4A1","SLC26A4", 
                      "CLCNKB", "FXYD4", "AQP2", "AQP3", "CD24", "EPCAM", "CD9","HSD11B2","CALB1","MAL","CLDN8", "KCNJ1"))
genes.IMM= unique(c("PTPRC","CD3D", "CD8A","CD8B", "NCAM1", "GNLY", "GZMB", "KLRD1","XCL1", "PRF1", "TRAC", "IL7R", "FOXP3", "TIGIT", 
                    "CD14","LYZ", "S100A12", "FCGR3A", "FCGR3B", "FCGR1A", "CD68", "FCGR2A", "APOE", "ITGAM", 
                    "CST3", "ITGAX", "IRF7", "IRF8", "CLEC4C", "TCF4", "LILRA4", "MS4A1", "CD79A", "CD19", "IGHM",
                    "IGHG2", "JCHAIN", "KIT", "MS4A2", "CPA3", "TPSAB1", "S100A8", "S100A9", "FCGR3B", "CSF3R", "MKI67", "TOP2A", "PCNA"))

pdf(paste0(output.dir, "05_DotPlot_CellAssignment_",res,"_PC",n.dims, ".pdf"), width = 15, height = 9,useDingbats = F)
DotPlot(seurat, features=as.vector(genes.notIMM), assay = "RNA") + RotatedAxis() 
DotPlot(seurat, features=as.vector(genes.IMM), assay = "RNA") + RotatedAxis() 
dev.off()

#Cluster assignment - res 0.7 - Harmony batch correction - by Experiment
seurat$Celltype_Harmony = NA
seurat$Celltype_Harmony[Idents(seurat)=="0"] = "CD8T"
seurat$Celltype_Harmony[Idents(seurat)=="1"] = "PT_GPX3"
seurat$Celltype_Harmony[Idents(seurat)=="2"] = "CD4T"
seurat$Celltype_Harmony[Idents(seurat)=="3"] = "NK"
seurat$Celltype_Harmony[Idents(seurat)=="4"] = "PT_MT1G"
seurat$Celltype_Harmony[Idents(seurat)=="5"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="6"] = "Unassign"
seurat$Celltype_Harmony[Idents(seurat)=="7"] = "Mac_C1QA"
seurat$Celltype_Harmony[Idents(seurat)=="8"] = "Mono_c"
seurat$Celltype_Harmony[Idents(seurat)=="9"] = "TumC_1"
seurat$Celltype_Harmony[Idents(seurat)=="10"] = "Treg"
seurat$Celltype_Harmony[Idents(seurat)=="11"] = "TumC_2"
seurat$Celltype_Harmony[Idents(seurat)=="12"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="13"] = "NK_CD160"
seurat$Celltype_Harmony[Idents(seurat)=="14"] = "Mono_nc"
seurat$Celltype_Harmony[Idents(seurat)=="15"] = "Epith"
seurat$Celltype_Harmony[Idents(seurat)=="16"] = "Fibro"
seurat$Celltype_Harmony[Idents(seurat)=="17"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="18"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="19"] = "PlasmaC"
seurat$Celltype_Harmony[Idents(seurat)=="20"] = "Epith"
seurat$Celltype_Harmony[Idents(seurat)=="21"] =  "RBC"
seurat$Celltype_Harmony[Idents(seurat)=="22"] = "Epith"
seurat$Celltype_Harmony[Idents(seurat)=="23"] = "Prolif"
seurat$Celltype_Harmony[Idents(seurat)=="24"] = "Mast"
seurat$Celltype_Harmony[Idents(seurat)=="25"] = "Unassign"


seurat$Tissue=seurat$TissueDiseaseState
seurat$Patient=seurat$Experiment

pdf(file = paste0(output.dir, "06_Final_plots_Integrated_Young_res",res, "_PC",n.dims,"_",Sys.Date(),".pdf"), useDingbats = F, width = 8, height = 8 )
DimPlot(seurat, reduction = "harmony_umap", label = T, label.size = 4) + labs(title= "res 0.6")
DimPlot(seurat, reduction = "harmony_umap", group.by = "Experiment", label = F, label.size = 4)
DimPlot(seurat, reduction = "harmony_umap", group.by = "TissueDiseaseState", label = F, label.size = 4)

DimPlot(seurat, reduction = "harmony_umap", group.by = "Celltype_CCA", label = F, label.size = 4)
FeaturePlot(seurat, reduction = "harmony_umap", features = "nCount_RNA")
FeaturePlot(seurat, reduction = "harmony_umap", features = "nFeature_RNA")

pt=as.data.frame(table(Idents(seurat), seurat$batch, seurat$Celltype_CCA, seurat$Tissue, 
                       seurat$Experiment))


ggplot(pt, aes(x = Freq , y = Var3 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))

ggplot(pt, aes(x = Freq , y = Var3 , fill = Var4)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))

ggplot(pt, aes(x = Freq , y = Var3 , fill = Var5)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))

dev.off()


pdf(file = paste0(output.dir, "06_Final_plots_SplitTissue_Integration_Young_UMAP_res",res, "_PC",n.dims,"_",Sys.Date(),".pdf"), useDingbats = F, width = 14, height = 8 )
DimPlot(seurat, reduction = "harmony_umap", split.by = "Tissue")
DimPlot(seurat, reduction = "harmony_umap", group.by = "Celltype_CCA", split.by = "Tissue")
DimPlot(seurat, reduction = "harmony_umap", group.by = "Patient", split.by = "Tissue")
dev.off()

pdf(file = paste0(output.dir, "06_Final_plots_QC_",Sys.Date(),".pdf"), useDingbats = F, width = 14, height = 8)
VlnPlot(seurat, features = "nFeature_RNA" , pt.size=0.01, split.by =  "Tissue", group.by="Celltype_CCA")
VlnPlot(seurat, features = "nCount_RNA" , pt.size=0.0, split.by =  "Tissue", group.by="Celltype_CCA") + scale_y_log10()
VlnPlot(seurat, features = "percent.mt", pt.size=0.01, split.by =  "Tissue", group.by="Celltype_CCA")
VlnPlot(seurat, features = "percent.rb" , pt.size=0.01, split.by =  "Tissue", group.by="Celltype_CCA")
dev.off()

pdf(file = paste0(output.dir, "06_Final_cells_QC_metrics",Sys.Date(),".pdf"), useDingbats = F, width = 20, height = 25)
metadata <- as.data.frame(seurat@meta.data)
metadata %>%  ggplot(aes(color=batch, x=nCount_RNA, fill= batch)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800)+
  annotate("text", x=1000, y=1, label= paste0(round((sum(metadata$nCount_RNA<800)/dim(metadata)[1])*100,0), ' %') )

metadata %>%  ggplot(aes(color=batch, x=nCount_RNA, fill=batch)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)+
  facet_wrap(~ Celltype_CCA, ncol=3, scales = "free")


metadata %>%  ggplot(aes(color=batch, x=nFeature_RNA, fill= batch)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  annotate("text", x=100, y=1, label= paste0(round((sum(metadata$nFeature_RNA<200)/dim(metadata)[1])*100,0), ' %') )
dev.off()

write.csv(table( seurat$Celltype_CCA, seurat$batch), file= paste0(output.dir, "06_Table_CCA_clustering",res,"_RNAassay_",Sys.Date(),".csv"))

# #SAVE DATA 
saveRDS(seurat, paste0("data/ccRCC_Young_integrated_Harmony_2_res0.7_PC50_",Sys.Date(),".rds"))


