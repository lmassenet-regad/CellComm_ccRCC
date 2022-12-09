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
cytokines=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cytokines")
receptors=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cyt_receptors")
checkpoints=readxl::read_excel("../Liste_molecules.xlsx", sheet="Checkpoints")
chemokines=readxl::read_excel("../Liste_molecules.xlsx", sheet="Chemokines")
all=unique(c(unique(receptors$Symbol), unique(cytokines$Symbol), unique(chemokines$Symbol),unique(checkpoints$Symbol)))
db=as.data.frame(readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20220705.xlsx"))
mol=unique(c(unique(db$`Ligand 1`), unique(db$`Ligand 2`), unique(db$`Receptor 1`),unique(db$`Receptor 2`),unique(db$`Receptor 3`)))
GOI=unique(c(mol,all)) 

###### Load data ###### 
#======================#
seurat.test=readRDS(paste0("data/ccRCC_nCount1000_integrated_Harmony_res1_PC50_2022-02-08.rds"))
assay="RNA"

#Supp Figure 2D
Idents(seurat.test)=seurat.test$RNA_snn_res.1
VlnPlot(seurat.test, features = c("CA9", "NDUFA4L2", "NNMT","EGLN3"), pt.size = 0, ncol = 2)


## A- TUMOR CELLS reclustering -  2 populations - only tumor sample ##### 
#============================================================================#
Idents(seurat.test)=seurat.test$Celltype_Harmony
sub2=subset(seurat.test, idents = "TumC", cells=which(seurat.test$Tissue=="Tum"))
sub2 <- RunPCA(sub2, dims = 1:50)
sub2 <- harmony::RunHarmony(sub2, dims = 1:50, group.by.vars="orig.ident",
                            reduction = "pca", assay.use="RNA", dims.use=1:50,
                            plot_convergence = T, verbose = T, reduction.save = "harmony2")

sub2 <- RunUMAP(sub2, dims = 1:50, reduction = "harmony2", reduction.name = "harmony2_umap")

sub2 <- FindNeighbors(sub2, dims = 1:50, reduction = "harmony2")
sub2  <- FindClusters(sub2, resolution = 0.1)
DimPlot(sub2, reduction = "harmony2_umap", label=T) # Figure 2A

FeatureScatter(sub2, feature1 = "nFeature_RNA", feature2 = "percent.rb") # Figure 2B

GOIinSeurat=list(rownames(sub2)[which(rownames(sub2)%in%GOI)])
sub2=AddModuleScore(sub2, features = GOIinSeurat, name = "CommMol")
VlnPlot(sub2, features="CommMol1", assay = assay, pt.size = 0.5) +stat_summary(fun= median, geom='point', size = 25, colour = "black", shape = 95) + 
  ggpubr::stat_compare_means(label = "p.format") + NoLegend() # Figure 2C

DefaultAssay(sub2)=assay 
seurat.markers2<- FindAllMarkers(sub2,only.pos = TRUE, min.pct = 0.1, 
                                 logfc.threshold = 0.25, assay = assay,
                                 recorrect_umi = FALSE) #Differential expression marker list TumC1 vs TumC2

markers <- seurat.markers2 %>% dplyr::filter(p_val_adj<=0.05 )  %>% as.data.frame
write.csv(markers , "/DEG_TumC_subcluster_TumSample_only_0.25logFC_pvalue0.05.csv") #Supp Table 2D

#Add TumC1 and TumC2 label to original dataset + plots ######
seurat.test$Celltype_Harmony2=seurat.test$Celltype_Harmony
seurat.test$Celltype_Harmony2[WhichCells(sub2, idents = c("0"))]="TumC1"
seurat.test$Celltype_Harmony2[WhichCells(sub2, idents = c("1"))]="TumC2"

pt=as.data.frame(table(Idents(sub2), sub2$orig.ident))
pt$Var2=factor(pt$Var2, levels =c("LM022_Tum","LM027_Tum", "LM029_Tum","LM022_Hty","LM027_Hty","LM029_Hty"))

ggplot(pt, aes(x = Freq , y = Var1 , fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Cell type") +
  xlab("Proportion of each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=0, hjust = 1)) # Figure 2A - suite


## B- DEG TumC (from tumoral tissue) compared to PT cells (from healthy tissue) ##### 
#============================================================================#
# For communication molecules
Idents(seurat.com2)=seurat.com2$Celltype_Harmony2
sub4=subset(seurat.com2,idents=c( "PT_GPX3", "PT_MT1G", "TumC2"))

DefaultAssay(sub4)=assay
sub4$Cell=NA
sub4$Cell[which(Idents(sub4)%in%c("TumC2"))]="TumC"
sub4$Cell[which(Idents(sub4)%in%c("PT_GPX3", "PT_MT1G"))]="PT"
Idents(sub4)=paste0(sub4$Cell, "_", sub4$Tissue)
markers4=FindMarkers(sub4, ident.1 ="TumC_Tum", ident.2 = "PT_Hty", only.pos = F, logfc.threshold = 0.25, 
                     min.pct = 0.1, assay = assay, recorrect_umi = FALSE)
markers4$gene=rownames(markers4)
markers4 =markers4 %>% filter(p_val_adj <0.05)
write.csv(markers4 %>% filter(p_val_adj <0.05), file = paste0(output.dir,"DEG_TumC_",assay,"/DEG_TumC2_CA9_Tum_vs_PT_Hty_inDB_logFC0.25_minpct0.1_padj0.05.csv"))


###### VISUALISATION DEG  Volcano Plot ###### 
# for Volcano Plot, need to redo FinMarkers without logFC threshold. 
mat=FindMarkers(sub4, ident.1 ="TumC_Tum", ident.2 = "PT_Hty", only.pos = F, logfc.threshold = 0, 
                min.pct = 0.1, assay = assay, recorrect_umi = FALSE)
mat$gene=rownames(mat)
mat=mat%>% filter(gene %in% GOI)  # DEG filtered from all the DEG TumC2 vs PT - change only 2 genes in comparison to take directly only the GOI for DEG analysis
mat2=mat %>% filter(abs(avg_log2FC)>0.25 & p_val_adj <0.05)
EnhancedVolcano::EnhancedVolcano(mat,
                                 lab = mat$gene,
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj', 
                                 FCcutoff = 0.25, legendPosition = "none", labSize = 3 , pointSize = 2, pCutoff = 0.05,   
                                 drawConnectors = TRUE, widthConnectors = 0.25,  axisLabSize = 10) #Figure 2D


###### VISUALISATION DEG Chart Plot - DEG communication list ###### 
df=as.data.frame(readxl::read_excel("analyses/0_BILAN_molecules_v2_2022.xlsx", sheet="DEG TumC vs PT")) #Supp Table
df$Family[is.na(df$Family)] <- "Other"
df2=as.data.frame(table(df$Family))

#get the positions for the labels
df3 <- df2 %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

pie <- ggplot(df2, aes(x="", y=Freq, fill=Var1))+
  geom_col(color = "black") + coord_polar("y", start=0)  +
  scale_y_continuous(
    breaks=df3$pos,    
    labels=paste0(df2$Freq)) + 
  theme_minimal() + theme(axis.text.x=element_text(color='black'), legend.position = "none") # the labels
pie # Figure 2E

###### PERFORM GENE ENRICHMENT ANALYSIS ###### 
df=read.csv("analyses/3_Downstream_analysis/DEG_TumC_RNA/DEG_TumC2_CA9_Tum_vs_PT_Hty_logFC0.25_minpct0.1_padj0.05.csv")
dim(df)

df=df %>% filter(avg_log2FC >0) # UP or DOWN
df=df %>% filter(gene %in% GOI)

H = msigdbr(species = "Homo sapiens", category = "H")
H.symbol= select(H, gs_name, gene_symbol)
#run enrichment
enrich.H=enricher(gene=df$gene, TERM2GENE = H.symbol)

#extract result
head(enrich.H@result)
enrich.H.df=enrich.H@result %>% separate (BgRatio, into=c("size.term", "size.cat"), sep="/") %>%
  separate (GeneRatio, into=c("size.overlap.term", "size.overlap.cat"), sep="/") %>%
  #convert character to numeric
  mutate_at(vars("size.term", "size.cat","size.overlap.term", "size.overlap.cat" ), as.numeric) %>%
  #calculate informative ratio size.overlap.term /size.term
  mutate("ratio"=size.overlap.term/size.term) %>%
  mutate(Description=gsub("HALLMARK_", "", Description), 
         Description=gsub("_", " ", Description))

enrich.H.df2=enrich.H.df %>% filter(p.adjust < 0.05) 
print(dim(enrich.H.df2)[[1]])

ggplot(enrich.H.df2, aes(x=reorder(Description, ratio), #reorder by ratio : most enriched term
                         y=ratio)) + geom_col() + theme_classic() + coord_flip() +
  labs(x="Gene set", y="Significant gene in set / total number of gene in set", 
       title = "Enrichment in DEG TumC vs PT ")

enrich.H.df2["HALLMARK_APOPTOSIS", "geneID"]

write.csv(enrich.H.df2, file = paste0(output.dir,"DEG_TumC_RNA", "/DEG_TumC2_CA9_Tum_vs_PT_Hty_ALL_DOWN_HALLMARK_enrichment_annotations.csv"))

## Plot most enriched term (ncount)
enrich.H.df3= enrich.H.df %>% filter(Count >8)
ggplot(enrich.H.df3, aes(x=reorder(Description, Count), #reorder by ratio : most enriched term
                         y=Count)) + geom_col() + theme_classic() + coord_flip() +
  labs(x="Gene set", y="number of gene in set", 
       title = "Enrichment in DEG TumC2 vs PT ") # Figure 2F



## C- ANALYSE SPECIFIC GENES TO TUMOR CELLS IN SEURAT.TUM OBJECT ##### 
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
  markers3a=FindMarkers(seurat.tum.com, ident.1 = c("TumC2"), ident.2 = cell , 
                        assay = assay, min.pct = 0.05, logfc.threshold = 0.0, only.pos = T)
  #recorrect_umi = FALSE) # TumC2_CA9 vs 
  markers3a$gene=rownames(markers3a)
  markers3a$type=paste0("VS_",cell)  
  markers=rbind(markers,markers3a)
}
markers$p_val_adj2=p.adjust(p =markers$p_val , method = "bonferroni") #adviced by Lilith
markers$diff.pct=markers$pct.1-markers$pct.2

markers2=markers%>%filter(avg_log2FC > 0.25  | (avg_log2FC > 0.1 & pct.2 <0.1 & diff.pct>0.05) )
markers2=markers2%>%filter( pct.1 >0.1)

fr.gene=data.frame(table(markers2$gene))
fr.gene %>% filter(Freq==length(cellOI)) %>% pull(Var1)

GenesOI = unique(c( as.vector(fr.gene %>% filter(Freq==length(cellOI)) %>% pull(Var1) )))
GenesOI=GenesOI[-c(13,25)] # not among DEG PT_hty vs tumor cells

#DotPlot - Figure 2G
Idents(seurat.tum) <- factor(x = Idents(seurat.tum), levels = unique(Idents(seurat.tum))[order(levels(Idents(seurat.tum)), decreasing = T)])
DotPlot(seurat.tum,  features = unique(GenesOI), cols = c("blue", "red"),  dot.min=0.1, assay = "RNA", idents = Idents(seurat.tum), scale.by="size", ) + RotatedAxis() 

