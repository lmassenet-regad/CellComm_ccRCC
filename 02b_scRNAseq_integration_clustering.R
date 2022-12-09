#Lucile Massenet-Regad - scRNAseq ccRCC PhD project
#Created 2021-05-21, last modified: 2022-02-1
#Data analysis - samples integration 

rm(list=ls())
work.dir= work.dir
output.dir =output.dir

# Load libraries
#===================================
library(dplyr)
library(ggplot2)
library(clustree)
library(Seurat)
library(readxl)
library(gridExtra)
library(ComplexHeatmap)

the.seed<-1337L
n.dims=30

#==========================================================================================#
######  A:  Integration tumor / normal tissue - using data after preprocessing script ######
#==========================================================================================#

# Integration with Harmony
#===================================
seurat.tum1=readRDS("data/LM022_Tum_res0.2_PC50_data_2022-02-04.rds") # Hty : Juxtatumoral, Tum: Tumoral
seurat.hty1=readRDS("data/LM022_Hty_res0.2_PC50_data_2022-02-04.rds")
seurat.tum2=readRDS("data/LM027_Tum_res0.2_PC50_data_2022-02-04.rds")
seurat.hty2=readRDS("data/LM027_Hty_res0.2_PC50_data_2022-02-04.rds")
seurat.tum3=readRDS("data/LM029_Tum_res0.2_PC50_data_2022-02-04.rds")
seurat.hty3=readRDS("data/LM029_Hty_res0.2_PC50_data_2022-02-04.rds")

seurat=merge(seurat.tum1, c(seurat.tum2, seurat.tum3, seurat.hty1, seurat.hty2, seurat.hty3), project = "combined")
seurat=subset(seurat, subset= nCount_RNA>1000)

table(seurat$orig.ident)
DefaultAssay(seurat)

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)

seurat <- RunPCA(seurat, npcs = n.dims+20, verbose = FALSE)
ElbowPlot(seurat, ndims = 100) #n=30 fine
seurat <- RunUMAP(seurat, dims = 1:n.dims, seed.use = the.seed)


p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = "orig.ident", pt.size = .1) +
                    NoLegend() +
                    ggtitle("Before harmony"))
p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", pt.size = .1, group.by = "orig.ident") + NoLegend()) + ggplot2::ggtitle("Before Integration")
p2 <- AugmentPlot(VlnPlot(object = seurat, features = "PC_1", group.by = "orig.ident", pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))

seurat <- suppressWarnings(harmony::RunHarmony(object =seurat,  group.by.vars="orig.ident",
                                               reduction = "pca", assay.use="RNA", dims.use=1:n.dims+20, plot_convergence = T, verbose = T))

p3 <- AugmentPlot(DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by ="orig.ident") + NoLegend())+ ggplot2::ggtitle("After Integration")
p4 <- AugmentPlot(VlnPlot(object = seurat, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))


seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:n.dims, reduction.name = "harmony_umap")
p5 <- AugmentPlot(DimPlot(seurat, reduction = "harmony_umap", group.by = "orig.ident", pt.size = .1) +
                    NoLegend() +
                    ggtitle("After harmony"))
pdf(paste0(output.dir, "03_Integration_Harmony_PC", n.dims+20,".pdf"), useDingbats = F, width = 13, height = 10)
ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4))
plot_grid(p0,p5)
dev.off()

# Apply Seurat pipeline
#===================================
pdf(paste0(output.dir,"03_Clustering_",Sys.Date(),".pdf"),useDingbats = F)
VizDimLoadings(object = seurat, dims = 1:2,reduction = "harmony" )
DimHeatmap(object = seurat, dims = 1:15, cells = 500, balanced = TRUE, reduction = "harmony")
ElbowPlot(object =seurat, ndims = n.dims+20, reduction = "harmony")  
dev.off()

seurat <- FindNeighbors(object = seurat, dims = 1:n.dims, reduction = "harmony")

# Clustree & Silhouette
seurat.test=seurat
pdf(paste0(output.dir,"03_Clustering_Silhouette_PC",n.dims,"_",Sys.Date(),".pdf"),useDingbats = F)
for (res2 in seq(0.1,1.4,0.1)){
  print(res2)
  #res2=0.1
  seurat.test <- FindClusters(object = seurat.test, resolution = res2)
  si <- cluster::silhouette(x= as.numeric(Idents(seurat.test)),
                            dist = dist(FetchData(seurat.test,
                                                  vars=c(as.vector(paste0("PC_", c(1:n.dims))) ))) )

  #plot(si)
  plot(si, cex.names = 0.5)
}
dev.off()

# Clusttree
pdf(paste0(output.dir,"03_Clustering_Clustree_PC",n.dims,"Young_",Sys.Date(),".pdf"),useDingbats = F, height = 15,width=10)
exprs <- "data"
prefix = "RNA_snn_res."
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

#DimPlot(seurat.test, group.by="RNA_snn_res.0.5", label = T)

DimPlot(seurat.test, reduction="harmony_umap", group.by="RNA_snn_res.1", label = T)

# 4- UMAP and clustering  - chose resolution
#====================================
res=1
seurat <- FindClusters(object = seurat, resolution = res)

seurat$Tissue=stringi::stri_sub(seurat$orig.ident,-3,-1)
seurat$Patient=stringi::stri_sub(seurat$orig.ident,1,5)

#Plot and save UMAP
pdf(paste0(output.dir,"04_Clustering_UMAP_", res, "_",Sys.Date(),".pdf"),useDingbats = F, height = 10,width=10)
DimPlot(seurat, reduction = "harmony_umap", label=T) + ggplot2::ggtitle(paste0("resolution = ", res)) + NoLegend()
DimPlot(seurat, reduction = "harmony_umap", group.by="orig.ident",label=F) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
DimPlot(seurat, reduction = "harmony_umap", group.by="Tissue",label=F) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
DimPlot(seurat, reduction = "harmony_umap", group.by="Patient",label=T, pt.size = 0.01 ) + ggplot2::ggtitle(paste0("Harmony integration - resolution = ", res))
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
write.csv(metascape, paste0(output.dir, "05_DEG_Harmony_Markers_res",res,"_RNAassay",Sys.Date(),".csv"))

seurat.markers %>% dplyr::filter(p_val_adj<=0.05 & avg_log2FC > 1) %>% arrange(cluster) %>% group_by(cluster) %>% top_n(20, avg_log2FC)  %>% View()

#Manual gene list for annotation
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


#Cluster assignment - res 1
seurat$Celltype_Harmony = NA
seurat$Celltype_Harmony[Idents(seurat)=="0"] = "CD4T"
seurat$Celltype_Harmony[Idents(seurat)=="1"] = "CD8T"
seurat$Celltype_Harmony[Idents(seurat)=="2"] = "PT_MT1G"
seurat$Celltype_Harmony[Idents(seurat)=="3"] = "PT_GPX3"
seurat$Celltype_Harmony[Idents(seurat)=="4"] = "NK_T"
seurat$Celltype_Harmony[Idents(seurat)=="5"] = "NK"
seurat$Celltype_Harmony[Idents(seurat)=="6"] = "MMAC_VCAN"
seurat$Celltype_Harmony[Idents(seurat)=="7"] = "MMAC_FCGR3A"
seurat$Celltype_Harmony[Idents(seurat)=="8"] = "cDC2"
seurat$Celltype_Harmony[Idents(seurat)=="9"] = "TumC"
seurat$Celltype_Harmony[Idents(seurat)=="10"] = "NK_CD160"
seurat$Celltype_Harmony[Idents(seurat)=="11"] = "B_cell"
seurat$Celltype_Harmony[Idents(seurat)=="12"] = "Prolif"
seurat$Celltype_Harmony[Idents(seurat)=="13"] = "Mast"
seurat$Celltype_Harmony[Idents(seurat)=="14"] = "Neutrop"
seurat$Celltype_Harmony[Idents(seurat)=="15"] = "MMAC_C1QA"
seurat$Celltype_Harmony[Idents(seurat)=="16"] = "Treg"
seurat$Celltype_Harmony[Idents(seurat)=="17"] = "Epith"
seurat$Celltype_Harmony[Idents(seurat)=="18"] = "pDC"
seurat$Celltype_Harmony[Idents(seurat)=="19"] = "cDC1"
seurat$Celltype_Harmony[Idents(seurat)=="20"] = "PlasmaC"
seurat$Celltype_Harmony[Idents(seurat)=="21"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="22"] = "Fibro"

pdf(file = paste0(output.dir, "04_Final_plots_Integration_HARMONY_UMAP_res",res, "_PC",n.dims,"_",Sys.Date(),".pdf"), useDingbats = F, width = 8, height = 8 )
DimPlot(seurat, reduction = "harmony_umap", label = T, label.size = 4) + labs(title= "res 0.7")
DimPlot(seurat, reduction = "harmony_umap", group.by = "orig.ident", label = F, label.size = 4) 
DimPlot(seurat, reduction = "harmony_umap", group.by = "Tissue", label = F, label.size = 4) 
DimPlot(seurat, reduction = "harmony_umap", group.by = "Patient", label = F, label.size = 4) 
DimPlot(seurat, reduction = "harmony_umap", group.by = "Celltype_Harmony", label = F, label.size = 4) 
FeaturePlot(seurat, reduction = "harmony_umap", features = "nCount_RNA")
FeaturePlot(seurat, reduction = "harmony_umap", features = "nFeature_RNA")

pt=as.data.frame(table(Idents(seurat), seurat$orig.ident))
pt$Var2=factor(pt$Var2, levels =c("LM022_Tum","LM027_Tum", "LM029_Tum","LM022_Hty","LM027_Hty","LM029_Hty"))
pt$Var3=substr(pt$Var2, 7,9)

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var3)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) 


pt=as.data.frame(table(seurat$Celltype_Harmony, seurat$orig.ident))
pt$Var2=factor(pt$Var2, levels =c("LM022_Tum","LM027_Tum","LM022_Hty","LM027_Hty", "LM029_Tum","LM029_Hty"))
pt$Var3=substr(pt$Var2, 7,9)

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var3)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) 

dev.off()


pdf(file = paste0(output.dir, "04_Final_plots_SplitTissue_Integration_LM022_27_29_Harmony_UMAP_res",res, "_PC",n.dims,"_",Sys.Date(),".pdf"), useDingbats = F, width = 14, height = 8 )
DimPlot(seurat, reduction = "harmony_umap", split.by = "Tissue") 
DimPlot(seurat, reduction = "harmony_umap", group.by = "Celltype_Harmony", split.by = "Tissue")
DimPlot(seurat, reduction = "harmony_umap", group.by = "Patient", split.by = "Tissue")
dev.off()

pdf(file = paste0(output.dir, "04_Final_plots_QC_",Sys.Date(),".pdf"), useDingbats = F, width = 14, height = 8)
VlnPlot(seurat, features = "nFeature_RNA" , pt.size=0.01, split.by =  "Tissue", group.by="Celltype_Harmony")
VlnPlot(seurat, features = "nCount_RNA" , pt.size=0.0, split.by =  "Tissue", group.by="Celltype_Harmony") + scale_y_log10() 
VlnPlot(seurat, features = "percent.mt", pt.size=0.01, split.by =  "Tissue", group.by="Celltype_Harmony")
VlnPlot(seurat, features = "percent.rb" , pt.size=0.01, split.by =  "Tissue", group.by="Celltype_Harmony")
dev.off()

pdf(file = paste0(output.dir, "04_Final_cells_QC_metrics",Sys.Date(),".pdf"), useDingbats = F, width = 14, height = 15)
metadata <- as.data.frame(seurat@meta.data)
metadata %>%  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800)+
  annotate("text", x=1000, y=1, label= paste0(round((sum(metadata$nCount_RNA<800)/dim(metadata)[1])*100,0), ' %') )

metadata %>%  ggplot(aes(color=orig.ident, x=nCount_RNA, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)+
  facet_wrap(~ Celltype_Harmony, ncol=3)


metadata %>%  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  annotate("text", x=100, y=1, label= paste0(round((sum(metadata$nFeature_RNA<200)/dim(metadata)[1])*100,0), ' %') )
dev.off()

write.csv(table( seurat$Celltype_Harmony, seurat$orig.ident), file= paste0(output.dir, "04_Table_CCA_clustering",res,"_RNAassay_",Sys.Date(),".csv"))

#SAVE DATA
saveRDS(seurat, paste0("data/ccRCC_nCount1000_integrated_Harmony_res1_PC50_",Sys.Date(),".rds"))


#==========================================================================================#
######                         B: Paper figures                                       ######
#==========================================================================================#

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpattern)

DimPlot(seurat, reduction="harmony_umap", label=T, group.by="RNA_snn_res.1") + NoLegend() # Figure 1C

#Supp Figure 2A - BEFORE AND AFTER HARMONY
display.brewer.pal(n = 6, name = 'RdBu', )
brewer.pal(n = 6, name = 'RdBu')
DimPlot(seurat, reduction = "harmony_umap", group.by = 'orig.ident', pt.size = 0.1, label=FALSE,
        cols = c("LM022_Hty"="#2166AC","LM022_Tum"="#B2182B", "LM027_Hty"= "#67A9CF", 
                 "LM027_Tum"="#EF8A62", "LM029_Hty"=  "#D1E5F0"  , "LM029_Tum"="#FDDBC7" ))


# Supp Figure 2B
Idents(seurat)=seurat$RNA_snn_res.1
markers=read.csv("05_DEG_Harmony_Markers_res1_RNAassay2022-02-08.csv", header = T)
rownames(markers)=markers$X
markers=markers[,-1]
markers %>% filter(p_val_adj<=0.05 & avg_log2FC > 1) %>% group_by(cluster) %>% top_n(3, avg_log2FC) -> filter.markers
seurat.test.averages = AverageExpression(seurat, return.seurat = TRUE)
seurat.test.averages = ScaleData(seurat.test.averages, features = filter.markers$gene, assay = "RNA")
DoHeatmap(seurat.test.averages, features = filter.markers$gene, assay = "RNA", raster = FALSE, draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

#Supp Figure 2C - Repartition per sample / per cluster
pt=as.data.frame(table(seurat$Celltype_Harmony, seurat$Tissue, seurat$Patient))
ggplot(pt , aes(y = Var1 , x = ifelse(test = Var2 == "Tum", yes = Freq, no = -Freq), fill = Var3, pattern=Var2)) +
  geom_col_pattern(color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(Tum = "stripe", Hty = "none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_bw(base_size = 15) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))


ggplot(pt , aes(y = Var1 , x = Freq, fill = Var3, pattern=Var2)) +
  geom_col_pattern(color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(Tum = "stripe", Hty = "none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_bw(base_size = 15) +
  ylab("Cell type") +
  xlab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1))


head(pt)
pt$Var4=paste0(pt$Var2, "_", pt$Var3)

head(pt2)

ggplot(pt , aes(x = Var1 , y = Freq , fill = Var4)) +
  geom_col(position = "fill", width = 0.5)+
  theme_bw(base_size = 15) +
  xlab("Cell type") +
  ylab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=45, hjust = 1))+
  scale_fill_brewer(palette = "RdBu", direction = -1)
# scale_color_manual(values = brewer.pal(n = 6, name = "RdBu"))


# Figure 1D

cytokines=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cytokines")
receptors=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cyt_receptors")
checkpoints=readxl::read_excel("../Liste_molecules.xlsx", sheet="Checkpoints")
chemokines=readxl::read_excel("../Liste_molecules.xlsx", sheet="Chemokines")
GOI=unique(c(unique(receptors$Symbol), unique(cytokines$Symbol), unique(chemokines$Symbol),unique(checkpoints$Symbol)))

meta.data=FetchData(object = seurat, slot="data", assay="RNA",  vars = c("Tissue", "Celltype_Harmony"))
meta.data$Sample.ID=rownames(meta.data)
data=as.data.frame(seurat[["RNA"]]@data)
data=data %>% filter(rownames(data) %in% GOI)
head(rownames(data))

family="receptors"
if (family=="checkpoints"){
  list=checkpoints$Symbol
}else if (family=="chemokines"){
  list=chemokines$Symbol
}else if (family=="cytokines"){
  list=cytokines$Symbol
}else if (family=="receptors"){
  list=receptors$Symbol
}else {print("not correct family name")}

data.oi=data %>% filter(rownames(data) %in% list ) # 63 chemokines detected
metagene.chemo=data.frame("Sample.ID"=colnames(data.oi) ,"score"=colSums(data.oi)) 
metagene.chemo=left_join(metagene.chemo, meta.data)

#heatmap
median.matrix = metagene.chemo %>%
  group_by(Celltype_Harmony, Tissue, .add=T) %>%
  summarise_if(is.numeric, median, na.rm=TRUE)
test=reshape2::dcast(median.matrix, formula = Celltype_Harmony ~ 
                       Tissue , value.var = "score", drop = T) 
rownames(test)= test$Celltype_Harmony
test$ratio=test$Tum/test$Hty
test$log2FC=log2(test$ratio)
test
write.csv(test, paste0("analyses/3_Downstream_analysis/0_Scores_results_",family,".csv"))

col_fun = circlize::colorRamp2(c(round(min(test[,2:3])), round(min(test[,2:3]) + (max(test[,2:3]) -min(test[,2:3]))/2), round(max(test[,2:3]))), c("blue","white", "red"))
ht1=Heatmap(as.matrix(test[,2:3]), cluster_rows = T, cluster_columns = F, name=family,
            show_column_dend = TRUE , show_row_dend = T, show_column_names = T,  show_row_names = T, col = col_fun)
ht1

#add pvalue info + save file
anno_df=ggpubr::compare_means(formula=score~Tissue, data=metagene.chemo,  group.by="Celltype_Harmony",  p.adjust.method="BH" )  %>%
  mutate(Family = family)
write.csv(anno_df, paste0("analyses/3_Downstream_analysis/0_comparison_scores_pvalue_",family,".csv"))

