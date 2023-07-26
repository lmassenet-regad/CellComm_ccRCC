#Lucile Massenet-Regad -
#Created 2021-05-27, last modified: 2022-07-05


rm(list=ls())
setwd(dir = work.dir)

output.dir="analyses/3_Downstream_analysis/"
assay="RNA" #SCT or RNA

#libraries
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
library(Seurat)
library(ComplexHeatmap)
library(circlize)


dir.create(path = paste0(work.dir,output.dir,'ICELLNET_', assay), 
           recursive = TRUE, showWarnings = FALSE)


db=as.data.frame(readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20220705.xlsx", sheet = 1))
db.name.couple=name.lr.couple(db, type="Family")

# load seurat object, retrieve matrix, average expression for each cluster
seurat=readRDS("data/ccRCC_Young_integrated_Harmony_2_res0.7_PC50_2022-07-05.rds")
table(Idents(seurat))

Idents(seurat)=seurat$Celltype_CCA
seurat.tum = subset(seurat, cells=which(seurat$Tissue=="Tumour"))
rm(seurat)

Idents(seurat.tum)=seurat.tum$Celltype_CCA
table(Idents(seurat.tum))

data <- sc.data.cleaning(object = seurat.tum, db = db, filter.perc =10, save_file = T, force.file = F, path= paste0(output.dir,'ICELLNET_', assay, "/"))
data.icell= as.data.frame(gene.scaling(data, n=1, db=db))

# central cell 1
CC2="TumC_2"
target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), "Class"=colnames(data.icell))
rownames(target)=target$ID

PC=c("PlasmaC", "Treg","CD4T",  "CD8T", "Unassign_T", "NK", "NK_CD160", "Prolif",
     "MMAC_1", "MMAC_2", "MMAC_3", "Mast",
     "Epith","PT_GPX3", "PT_MT1G","TumC_1","TumC_2", "Endoth","Fibro")


#direction communication score
direction="out"

# data selection 
CC2.data=as.data.frame(data.icell[,c(CC2, "Symbol")], row.names = rownames(data.icell))
PC.data=as.data.frame(data.icell[,c(PC, "Symbol")], row.names = rownames(data.icell))

#icellnet score
score.computation.1= icellnet.score(direction=direction, PC.data=PC.data, 
                                    CC.data= CC2.data,  
                                    PC.target = target, PC=PC, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]
ymax=max(score1)

#BARPLOT

my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding", "ECM")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  , "ECM"="slateblue4",
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "slateblue1",  "other" = "#908F90",  "NA"="#908F90")

p1=LR.family.score(lr=lr1, my.family=my.family, db.couple=db.name.couple, plot=T, title=paste0(CC2,"-", direction), family.col=family.col, ymax = ymax) +
  ggplot2::theme(axis.text.x=element_text(angle=90, hjust = 1))
p1
ggsave(paste0(output.dir,"ICELLNET_", assay,"/", CC2,"_", direction, "_barplot_",Sys.Date(),".pdf"))

# B -  HEATMAPS - see differences of specific channels
data.cell=lr1
data.cell=data.cell[complete.cases(data.cell),]
data.cell=data.cell[which(rowSums(data.cell)>0),]
dim(data.cell)
write.csv(data.cell,paste0(work.dir,output.dir,'ICELLNET_', assay, "/ICELLNET_LR_TumC_Scores_",direction,"_",Sys.Date(),".csv"))

#-----------------------------------------------
# Contribution TumC  (for each interaction) compared to other cells - analysis considering all the db (not only interactions >0 for TumC)
#========================
test=matrix(nrow=length(rownames(db)), ncol=length(PC))
rownames(test)=db.name.couple[,1]
colnames(test)=colnames(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], couple=db.name.couple[1,1], db=db, plot=F))
for (i in db.name.couple[,1]){
  if (direction =="out"){
    test[i,]=rowSums(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], couple=i, db=db, plot=F)) # direction OUT
  }
  if (direction =="in"){
    test[i,]=colSums(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], couple=i, db=db, plot=F)) # direction IN
  }
}
write.csv(test,paste0(work.dir,output.dir,'ICELLNET_', assay, "/ICELLNET_sumScores_",direction,"_ALL_CELLS_",Sys.Date(),".csv"))





######## 2021-08-04 FIND CRITERIA TO SELECT specific tumor cell interactions #####
##################################################################################


int.spe <- function(LR.mat, CoI, thresh){  # function to determine specific interactions for a cell type
  colMax= apply(LR.mat, 1, which.max)
  a=which(colnames(LR.mat)==CoI)
  LR.mat2=LR.mat[which(colMax==a),]
  #create dataframe to store maximum value of CoI
  value=data.frame("vmax"=LR.mat2[,a])
  rownames(value)=rownames(LR.mat2)
  value$vmax2=apply(LR.mat2[,-a], 1, function(x) (max(x)+0.1))
  value$vmax2_nb=apply(LR.mat2[,-a], 1, which.max)
  value$vmax2_cell=colnames(LR.mat2[,-a])[value$vmax2_nb]
  value=value%>% dplyr::select(-c("vmax2_nb"))
  value$vmean_pos=apply( LR.mat2[,-a], 1, function(x) mean(x[x!=0]))
  value$ratio=value$vmax/value$vmax2
  value = value %>% filter(ratio >thresh)
  return(value)
}

thresh=1.5
LR.mat=as.data.frame(read.csv(paste0(paste0( output.dir,"/ICELLNET_RNA/ICELLNET_sumScores_",direction,"_ALL_CELLS_2022-10-26.csv"))))
rownames(LR.mat)=LR.mat$X
head(LR.mat)
LR.mat=LR.mat[,-c(1,18)] # do not forget to remove X and TumC1-> to study TumC2 specificity
test=int.spe(LR.mat, CoI="TumC_2", thresh=thresh)
write.csv(test, paste0(work.dir,output.dir,'ICELLNET_', assay, "/LR_spé_tumC_",direction,"_thresh_",thresh,"_ALL_CELLS_", Sys.Date(),".csv"))
test
dim(test)

### 2nd step # ADD info of DEG -> keep only interaction where TumC upregulated compared to PT
DEG=as.data.frame(read.csv("analyses/2_DEG_TumC/DEG_TumC2_CA9_Tum_vs_PT_Hty_inDB_logFC0.1_minpct0.1_padj0.05.csv"))
DEG.gene=dplyr::filter(DEG, DEG$avg_log2FC>0.25 & as.numeric(DEG$p_val_adj)<0.05) %>% pull(X)

if (direction=="out"){db.DEG <- filter(db, db$`Ligand 1` %in% DEG.gene | db$`Ligand 2` %in% DEG.gene)}
if (direction=="in"){db.DEG <- filter(db, db$`Receptor 1` %in% DEG.gene| db$`Receptor 2` %in% DEG.gene | db$`Receptor 3` %in% DEG.gene)}
DEG.group= icellnet::name.lr.couple(db.DEG, type="Family")[,1]

test2=test %>% filter(rownames(test) %in% DEG.group)
test2 
dim(test2)
dplyr::setdiff(rownames(test), rownames(test2))

write.csv(test2, paste0(work.dir,output.dir,'ICELLNET_', assay, "/LR_spé_tumC_DEG",direction,"_thresh_",thresh,"_ALL_CELLS_", Sys.Date(),".csv"))

