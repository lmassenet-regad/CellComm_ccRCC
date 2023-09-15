#Lucile Massenet-Regad 
#Created 2021-05-21, last modified: 2023-09-15
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
seurat.markers<- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, assay = "RNA")

seurat.markers %>% dplyr::filter(p_val_adj<=0.05 & avg_log2FC > 0.5)  %>% as.data.frame  #Table S1


#Manual gene list for annotation
genes.assignment=c("CD3D", "CD4", "IL7R", "FOXP3", "CD8A", "GZMK",  "MKI67", "TOP2A", "GNLY", "KLRD1", "NCAM1", 
              "MS4A1", "CD79A",  "IGHM", "JCHAIN", 
              "C1QA", "MRC1", "CD14", "LYZ","FCGR3A", "HLA-DRA",
              "CD1C", "IDO1", "CLEC9A", "IRF7", "CLEC4C", "LILRA4",
              "KIT", "TPSAB1", "S100A8", "S100A9", "FCGR3B", 
              "ALDOB", "SLC34A1", "PDZK1IP1", "ANPEP", 
              "CA9","NDUFA4L2","NNMT", "EGLN3", 
              "EPCAM", "DEFB1", "UMOD", 
              "ENG", "PECAM1", "PTPRB",
              "PDGFRB", "ACTA2", "RGS5")

DotPlot(seurat, features=as.vector(genes.assignment), assay = "RNA") + coord_flip() # Figure 1C

#Cluster assignment - res 1
seurat$Celltype_Harmony = NA
seurat$Celltype_Harmony[Idents(seurat)=="0"] = "CD4T"
seurat$Celltype_Harmony[Idents(seurat)=="1"] = "CD8T"
seurat$Celltype_Harmony[Idents(seurat)=="2"] = "PT_MT1G"
seurat$Celltype_Harmony[Idents(seurat)=="3"] = "PT_GPX3"
seurat$Celltype_Harmony[Idents(seurat)=="4"] = "NK_T"
seurat$Celltype_Harmony[Idents(seurat)=="5"] = "NK"
seurat$Celltype_Harmony[Idents(seurat)=="6"] = "Mono_c"
seurat$Celltype_Harmony[Idents(seurat)=="7"] = "Mono_nc"
seurat$Celltype_Harmony[Idents(seurat)=="8"] = "cDC2"
seurat$Celltype_Harmony[Idents(seurat)=="9"] = "ccRCC"
seurat$Celltype_Harmony[Idents(seurat)=="10"] = "NK_CD160"
seurat$Celltype_Harmony[Idents(seurat)=="11"] = "B_cell"
seurat$Celltype_Harmony[Idents(seurat)=="12"] = "Prolif"
seurat$Celltype_Harmony[Idents(seurat)=="13"] = "Mast"
seurat$Celltype_Harmony[Idents(seurat)=="14"] = "Neutrop"
seurat$Celltype_Harmony[Idents(seurat)=="15"] = "Mac_C1QA"
seurat$Celltype_Harmony[Idents(seurat)=="16"] = "Treg"
seurat$Celltype_Harmony[Idents(seurat)=="17"] = "Epith"
seurat$Celltype_Harmony[Idents(seurat)=="18"] = "pDC"
seurat$Celltype_Harmony[Idents(seurat)=="19"] = "cDC1"
seurat$Celltype_Harmony[Idents(seurat)=="20"] = "PlasmaC"
seurat$Celltype_Harmony[Idents(seurat)=="21"] = "Endoth"
seurat$Celltype_Harmony[Idents(seurat)=="22"] = "Fibro"


#SAVE DATA
saveRDS(seurat, paste0("data/ccRCC_nCount1000_integrated_Harmony_res1_PC50_",Sys.Date(),".rds"))



#==========================================================================================#
######  B: Differential analysis for each cluster to update ICELLNET database         ######
#==========================================================================================#

# All cell types except cancer cells - Conduct differential analyses between tissue
db=as.data.frame(readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20210830.xlsx"))
mol=unique(c(unique(db$`Ligand 1`), unique(db$`Ligand 2`), unique(db$`Receptor 1`),unique(db$`Receptor 2`),unique(db$`Receptor 3`)))

seurat.com = DietSeurat(seurat, features = mol, data = TRUE, scale.data =T,  counts = T, assays = c(assay) )
Idents(seurat.com)=seurat.com$Celltype_Harmony
cellOI=unique(Idents(seurat.com))[-c(6, 10)] # Exclude cluster of tumor cells and prolif

markers.deg=data.frame()
for (cell in cellOI){
  seurat.sub=subset(seurat.com, cells=which(seurat.com$Celltype_Harmony==cell))
  Idents(seurat.sub)=seurat.sub$Tissue
  sub.marker=FindMarkers(object = seurat.sub, 
                         ident.1 = "Tum", ident.2 = "Hty", 
                         logfc.threshold = 0.1, min.pct =0.1, 
                         recorrect_umi = FALSE) 
  sub.marker$gene=rownames(sub.marker)
  sub.marker$diff.pct=sub.marker$pct.1 - sub.marker$pct.2
  sub.marker=sub.marker%>% filter(p_val_adj<0.1)
  write.csv(sub.marker, file=paste0(output.dir, "DEG_TUMvsHTY_",assay, "/DEG_",cell,"_TUMvsHTY_inDBicellnet20210830",Sys.Date(),".csv"))
  if (dim(sub.marker)[[1]]>0){
    sub.marker$Cell=cell
    markers.deg=rbind(markers.deg,sub.marker)
  }
  
}
write.csv(markers.deg, file=paste0(output.dir, "DEG_TUMvsHTY_",assay, "/ALL_summaryDEG_TUMvsHTY_inDBicellnet20210830_",Sys.Date(),".csv"))


# INTERSECTION DEG all genes with NATMI file -> to identify putative missing interactions
#============================================================================#
NATMIdb=readxl::read_excel("NATMI_DB_41467_2020_18873_MOESM4_ESM.xlsx", sheet="literature_support") #From NATMI paper, PMID 33024107 - Supplementary Data 1 
NATMIgene=unique(c(NATMIdb$`Ligand gene symbol`), c(NATMIdb$`Receptor gene symbol`))

# intersection of previously DEG analyses (cells from tumors versus juxtatumors, for each cell type) with NATMI DB
ALL_tum_hty_deg=read.csv("ALL_summaryDEG_TUMvsHTY_2022-02-15.csv")
ALL_tum_hty_deg=ALL_tum_hty_deg %>% filter(ALL_tum_hty_deg$p_val_adj<0.05 & abs(ALL_tum_hty_deg$avg_log2FC)>0.25)
intersect(ALL_tum_hty_deg$gene, NATMIgene) 
setdiff(intersect(ALL_tum_hty_deg$gene, NATMIgene), mol)


#intersection of NATMI DB with the list of cancer DEG (ccRCC2 compared to juxtatumoral proximal tubules) 
TumCdeg=read.csv("DEG_TumC2_CA9_Tum_vs_PT_Hty_logFC0.25_minpct0.1_padj0.05.csv") #  Supplementary Table S4
TumCdeg=TumCdeg %>% filter(TumCdeg$p_val_adj<0.05 & abs(TumCdeg$avg_log2FC)>0.25)
intersect(TumCdeg$gene, NATMIgene) # Comm molecules TumC2 DEG 
setdiff(intersect(TumCdeg$gene, NATMIgene), mol)


#==========================================================================================#
######                         C: Paper figures                                       ######
#==========================================================================================#

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpattern)

#### Figure 1B ####
DimPlot(seurat, reduction="harmony_umap", label=T, group.by="RNA_snn_res.1") + NoLegend() # Figure 1B
DimPlot(seurat, reduction="harmony_umap", label=F, split.by ="Tissue") #Figure 1B

#### Figure 1C ####
levels(seurat) <- c("0","16","1","12","4","5","10","11","20",
                    "15", "6","7","8","19", "18", "13","14",
                    "2","3","9","17", "21", "22")

DotPlot(seurat, features=as.vector(genes.final), assay = "RNA") + coord_flip()

#### Figure 1D-1E ####
listes_goi=readxl::read_excel("~/Liste_broad_communication_molecules.csv")

meta.data=FetchData(object = seurat, slot="data", assay="RNA",  vars = c("Tissue", "Celltype_Harmony"))
meta.data$Sample.ID=rownames(meta.data)
seurat.com = DietSeurat(seurat, features = GOI, data = TRUE, scale.data =T,  counts = T, assays ="RNA") # 549 over 604 unique genes detected
data=as.data.frame(seurat.com[["RNA"]]@data)
data.tum=data[,meta.data$Sample.ID[which(meta.data$Tissue=="Tum")]]

names=unique(meta.data$Celltype_Harmony)
names=names[order(names)]

score.mat=matrix(ncol = dim(listes_goi)[2], nrow =  length(names))
colnames(score.mat)=colnames(listes_goi)
rownames(score.mat)=names

for (i in colnames(listes_goi)){
  
  data.oi= data.tum %>% filter(rownames(data) %in% as.character(unlist(listes_goi[i]))) 
  metagene=data.frame("Sample.ID"=colnames(data.oi) ,"score"=colSums(data.oi)) 
  metagene=left_join(metagene, meta.data)
  mean.matrix = metagene %>%
    group_by(Celltype_Harmony, Tissue, .add=T) %>%
    summarise_if(is.numeric, mean, na.rm=TRUE)
  rownames(mean.matrix)=paste0(mean.matrix$Celltype_Harmony)
  mean.matrix=mean.matrix[names,]
  
  score.mat[,i]=mean.matrix$score
}

score.mat= score.mat[complete.cases(score.mat),]
heat.data = apply(score.mat, MARGIN=2, FUN = function(x) (x-mean(x))/sd(x)) # center reduction step

ht1=ComplexHeatmap::Heatmap(t(heat.data), cluster_rows = T, cluster_columns = T ,
                            show_column_dend = TRUE , show_row_dend = T, show_column_names = T,  show_row_names = T, name = "z-score")
ht1 # Figure 1D
ht1 = draw(ht1)
rows=row_dend(ht1) 
cols=column_dend(ht1)
order_row=rownames(t(heat.data))[unlist(rows)]
order_cols=colnames(t(heat.data))[unlist(cols)]

# For Figure 1E
score.mat2=matrix(ncol = dim(listes_goi)[2], nrow =  length(names))
colnames(score.mat2)=colnames(listes_goi)
rownames(score.mat2)=names

score.pvalue=matrix(ncol = dim(listes_goi)[2], nrow =  length(names))
colnames(score.pvalue)=colnames(listes_goi)
rownames(score.pvalue)=names

for (i in colnames(listes_goi)){
  data.oi= data %>% filter(rownames(data) %in% as.character(unlist(listes_goi[i])))
  metagene=data.frame("Sample.ID"=colnames(data.oi) ,"score"=colSums(data.oi)) 
  metagene=left_join(metagene, meta.data)
  mean.matrix = metagene %>%
    group_by(Celltype_Harmony, Tissue, .add=T) %>%
    summarise_if(is.numeric, mean, na.rm=TRUE)
  
  test=reshape2::dcast(mean.matrix, formula = Celltype_Harmony ~ 
                         Tissue , value.var = "score", drop = T) 
  rownames(test)= test$Celltype_Harmony

  test["ccRCC", "Hty"]=mean(x=metagene$score[which(metagene$Celltype_Harmony %in% c("PT_GPX3", "PT_MT1G") & metagene$Tissue =="Hty")])
  metagene.ccRCC=filter(metagene, (metagene$Celltype_Harmony %in% c("PT_GPX3", "PT_MT1G") & metagene$Tissue =="Hty") | (metagene$Celltype_Harmony=="ccRCC"))
  
  test$ratio=test$Tum/test$Hty
  anno_df=ggpubr::compare_means(formula=score~Tissue, data=metagene,  group.by="Celltype_Harmony", p.adjust.method = "none")
  test=left_join(test, anno_df[c("Celltype_Harmony", "p")])
  rownames(test)= test$Celltype_Harmony
  test["ccRCC", "p"]=wilcox.test(formula=score~Tissue, data=metagene.ccRCC)$p.value
  test$p.adjust=p.adjust(test$p, method = "fdr")
  
  write.csv(test, paste0("analyses/3_Downstream_analysis/0_comparison_scores_pvalue_",i,"_", Sys.Date() ,".csv"))
  score.mat2[,i]=test$ratio
  score.pvalue[,i]=test$p.adjust
}

log2score=log2(score.mat2)

ht3=ComplexHeatmap::Heatmap(t(log2score), row_order = order_row, column_order = order_cols, 
                            show_column_dend = TRUE , show_row_dend = T, show_column_names = T,  show_row_names = T, name = "log2FC", cell_fun = function(j, i, x, y, w, h, fill){
                              if (t(score.pvalue)[i, j] < 0.001 & abs(t(log2score)[i,j])>0.25 ) {
                                grid.text("***", x, y)
                              } else if(t(score.pvalue)[i, j] < 0.01 & abs(t(log2score)[i,j])>0.25) {
                                grid.text("**", x, y)
                              }
                            }) #highlight significant and abs(logFC) >0.25
ht3 # Figure 1E



#### Figure S1C #### 
display.brewer.pal(n = 6, name = 'RdBu', )
brewer.pal(n = 6, name = 'RdBu')
DimPlot(seurat, reduction = "harmony_umap", group.by = 'orig.ident', pt.size = 0.1, label=FALSE,
        cols = c("LM022_Hty"="#2166AC","LM022_Tum"="#B2182B", "LM027_Hty"= "#67A9CF", 
                 "LM027_Tum"="#EF8A62", "LM029_Hty"=  "#D1E5F0"  , "LM029_Tum"="#FDDBC7" ))

#### Figure S1D #### 
pt=as.data.frame(table(seurat$Celltype_Harmony, seurat$Tissue, seurat$Patient))
pt$Var4=paste0(pt$Var2, "_", pt$Var3)

ggplot(pt , aes(x = Var1 , y = Freq , fill = Var4)) +
  geom_col(position = "fill", width = 0.5)+
  theme_bw(base_size = 15) +
  xlab("Cell type") +
  ylab("Proportion or each sample") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle=45, hjust = 1))+
  scale_fill_brewer(palette = "RdBu", direction = -1)

