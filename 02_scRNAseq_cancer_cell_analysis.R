#Lucile Massenet-Regad 
#Data analysis post integration 

rm(list=ls())
work.dir=work.dir
setwd(dir = work.dir)
output.dir="analyses/3_Downstream_analysis/"
input.dir="data/"

######  Load libraries ###### 
#===================================#
library(dplyr)
library(Seurat)
library(readxl)
library(gridExtra)
library(ComplexHeatmap)
library(EnhancedVolcano)

library(ggplot2)
library(ggrepel)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)

###### List of molecules of interest ###### 
#===================================#
all=unique(c(unique(receptors$Symbol), unique(cytokines$Symbol), unique(chemokines$Symbol),unique(checkpoints$Symbol)))
db=as.data.frame(readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20230412.xlsx"))
mol=unique(c(unique(db$`Ligand 1`), unique(db$`Ligand 2`), unique(db$`Receptor 1`),unique(db$`Receptor 2`),unique(db$`Receptor 3`)))
GOI=unique(c(mol,all)) 

###### Load data ###### 
#======================#
seurat=readRDS(paste0("data/ccRCC_nCount1000_integrated_Harmony_res1_PC50_2022-02-08.rds"))
assay="RNA"


## A- Cancer cells reclustering -  2 populations - only tumor sample ##### 
#============================================================================#
Idents(seurat)=seurat$Celltype_Harmony
sub2=subset(seurat, idents = "ccRCC", cells=which(seurat$Tissue=="Tum"))
sub2 <- RunPCA(sub2, dims = 1:50)
sub2 <- harmony::RunHarmony(sub2, dims = 1:50, group.by.vars="orig.ident",
                            reduction = "pca", assay.use="RNA", dims.use=1:50,
                            plot_convergence = T, verbose = T, reduction.save = "harmony2")

sub2 <- RunUMAP(sub2, dims = 1:50, reduction = "harmony2", reduction.name = "harmony2_umap")

sub2 <- FindNeighbors(sub2, dims = 1:50, reduction = "harmony2")
sub2  <- FindClusters(sub2, resolution = 0.1)

DimPlot(sub2, reduction = "harmony2_umap", label=T) # Figure 2A

FeatureScatter(sub2, feature1 = "nFeature_RNA", feature2 = "percent.rb") # Figure 2B

VlnPlot(sub2, features = c("CA9"), group.by = "Celltype_Harmony2",
        pt.size = 0.1, 
        y.max = 3.5, # add the y-axis maximum value - otherwise p-value hidden
) + ggpubr::stat_compare_means(comparisons = list(c("ccRCC1", "ccRCC2")), label = "p.format") + #"p.signif" or "p.format"
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)) # Figure 2C - left

comm.mol=intersect(unique(mol), rownames(sub2))
GOIinSeurat=list(rownames(sub2)[which(rownames(sub2)%in%GOI)])
sub2=AddModuleScore(sub2, features = GOIinSeurat, name = "CommMol")
VlnPlot(sub2, features = c("CommMol1"), group.by = "Celltype_Harmony2", pt.size = 0.1, y.max = 0.2) + 
  ggpubr::stat_compare_means(comparisons = list(c("ccRCC1", "ccRCC2")), label = "p.format") + #"p.signif" or "p.format"
  theme(axis.text.x = element_text(angle = 0, hjust=0.5))  + 
  stat_summary(fun.y = median, 
               fun.ymin = median, 
               fun.ymax = median, 
               geom = "crossbar", 
               width = 0.5) # Figure 2C - right

seurat.markers2<- FindAllMarkers(sub2,only.pos = TRUE, min.pct = 0.1, 
                                 logfc.threshold = 0.25, assay = assay,
                                 recorrect_umi = FALSE)
markers <- seurat.markers2 %>% dplyr::filter(p_val_adj<=0.05 )  %>% as.data.frame. # Supplementary Table S3A

#Add ccRCC1 and ccRCC2 label to original datasets
seurat.test$Celltype_Harmony2=seurat.test$Celltype_Harmony
seurat.test$Celltype_Harmony2[WhichCells(sub2, idents = c("0"))]="ccRCC1"
seurat.test$Celltype_Harmony2[WhichCells(sub2, idents = c("1"))]="ccRCC2"

pt=as.data.frame(table(Idents(sub2), sub2$orig.ident))
pt$Var2=factor(pt$Var2, levels =c("LM022_Tum","LM027_Tum", "LM029_Tum","LM022_Hty","LM027_Hty","LM029_Hty"))

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion of each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) # Figure 2A - bottom


## B - DEG ccRCC2  compared to PT cells (from juxtatumors) ##### 
#============================================================================#
# For communication molecules
Idents(seurat.com2)=seurat.com2$Celltype_Harmony2
sub4=subset(seurat.com2,idents=c( "PT_GPX3", "PT_MT1G", "ccRCC2"))

DefaultAssay(sub4)=assay
sub4$Cell=NA
sub4$Cell[which(Idents(sub4)%in%c("ccRCC@"))]="ccRCC"
sub4$Cell[which(Idents(sub4)%in%c("PT_GPX3", "PT_MT1G"))]="PT"
Idents(sub4)=paste0(sub4$Cell, "_", sub4$Tissue)
markers4=FindMarkers(sub4, ident.1 ="ccRCC_Tum", ident.2 = "PT_Hty", only.pos = F, logfc.threshold = 0.25, 
                     min.pct = 0.1, assay = assay, recorrect_umi = FALSE)
markers4$gene=rownames(markers4)
markers4 =markers4 %>% filter(p_val_adj <0.05) # Supplementary Table S3B


## C- Cancer cell-specific communication genes analysis ##### 
#============================================================================#

seurat.tum = subset(seurat.test, cells=which(seurat.test$Tissue=="Tum"))
Idents(seurat.tum)=seurat.tum$Celltype_Harmony2
seurat.tum.com = DietSeurat(seurat.tum, features = GOI, data = TRUE, scale.data =T,  counts = T, assays = assay)
Idents(seurat.tum.com)=seurat.tum.com$Celltype_Harmony2

# Remove tumor cells cluster 
markers=data.frame()
cellOI=unique(Idents(seurat.tum.com))[-c(6,12)] #  to remove tumor cells clusters

for (cell in cellOI ){
  print(cell)
  markers3a=FindMarkers(seurat.tum.com, ident.1 = c("ccRCC2"), ident.2 = cell , 
                        assay = assay, min.pct = 0.05, logfc.threshold = 0.0, only.pos = T)
  markers3a$gene=rownames(markers3a)
  markers3a$type=paste0("VS_",cell)  
  markers=rbind(markers,markers3a)
}
markers$p_val_adj2=p.adjust(p =markers$p_val , method = "bonferroni") 
markers$diff.pct=markers$pct.1-markers$pct.2

markers2=markers%>%filter(avg_log2FC > 0.25  | (avg_log2FC > 0.1 & pct.2 <0.1 & diff.pct>0.05) )
markers2=markers2%>%filter( pct.1 >0.1)

fr.gene=data.frame(table(markers2$gene))
fr.gene %>% filter(Freq==length(cellOI)) %>% pull(Var1)

GenesOI = unique(c( as.vector(fr.gene %>% filter(Freq==length(cellOI)) %>% pull(Var1) )))
GenesOI=GenesOI[-c(13,25)] # not among DEG PT_hty vs cancer cells

#Figure 2D
Idents(seurat.tum) <- factor(x = Idents(seurat.tum), levels = unique(Idents(seurat.tum))[order(levels(Idents(seurat.tum)), decreasing = T)])
DotPlot(seurat.tum,  features = unique(GenesOI), cols = c("blue", "red"),  dot.min=0.1, assay = "RNA", idents = Idents(seurat.tum), scale.by="size", ) + RotatedAxis() 

