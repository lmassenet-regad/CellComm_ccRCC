#Lucile Massenet-Regad - scRNAseq ccRCC PhD project
#Created 2021-05-27, last modified: 2023-09-15
# ICELLNET on Zhang dataset

rm(list=ls())
work.dir="~/Documents/PhD_TUMOR_ccRCC/BIOINFO/Zhang et al, 2021, PNAS/"
setwd(dir = work.dir)

output.dir="analyses/"
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


db=as.data.frame(readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20230412.xlsx", sheet = 1))
db2=db
db.name.couple=name.lr.couple(db2, type="Family")

#### DATASET DEPENDANT STEP
# load seurat object, retrieve matrix, average expression for each cluster
seurat=readRDS("Seurat_Integrated_20230830.rds")

Idents(seurat)=seurat$Celltype_CCA
seurat.tum = subset(seurat, cells=which(seurat$Tissue=="Tum"))
rm(seurat)

table(Idents(seurat.tum))

data <- sc.data.cleaning(object = seurat.tum, db = db2, filter.perc =10, save_file = T, force.file = F, path= paste0(work.dir,output.dir,'ICELLNET_', assay, "/"))
data.icell= as.data.frame(gene.scaling(data, n = 1, db=db2))

## CYTOKINE RECEPTOR DETECTION
data.stat=read.csv(paste0(output.dir,'ICELLNET_', assay, "/scRNAseq_statsInfo_for_ICELLNET.csv"))

# central cell 1
CC1="TumC"
target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), "Class"=colnames(data.icell))
rownames(target)=target$ID

PC=c("PlasmaC", "T_cells",  "NK", "MMAC_1", "MMAC_2", "DC",
     "Mast","Endoth_ACKR1", "Endoth_PLVAP", "Endoth_PDGFRB", "Fibro", "Prolif", "TumC")

#direction communication score
direction="in"

# data selection 
CC1.data=as.data.frame(data.icell[,c(CC1, "Symbol")], row.names = rownames(data.icell))
PC.data=as.data.frame(data.icell[,c(PC, "Symbol")], row.names = rownames(data.icell))

#icellnet score
score.computation.1= icellnet.score(direction=direction, PC.data=PC.data, 
                                    CC.data= CC1.data,  
                                    PC.target = target, PC=PC, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db2)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]
ymax=max(score1)

#BARPLOT

my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding", "ECM")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  , "ECM"="slateblue4",
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "slateblue1",  "other" = "#908F90",  "NA"="#908F90")

p1=LR.family.score(lr=lr1, my.family=my.family, db.couple=db.name.couple, plot=T, title=paste0(CC1,"-", direction), family.col=family.col, ymax = ymax) +
  ggplot2::theme(axis.text.x=element_text(angle=90, hjust = 1))
p1
#ggsave(paste0(output.dir,"ICELLNET_", assay,"/", CC1,"_", direction, "_barplot_",Sys.Date(),".pdf"))

# B -  HEATMAPS - see differences of specific channels
data.cell=lr1
data.cell=data.cell[complete.cases(data.cell),]
data.cell=data.cell[which(rowSums(data.cell)>0),]
dim(data.cell)
write.csv(data.cell,paste0(work.dir,output.dir,'ICELLNET_', assay, "/ICELLNET_TumC_Scores_",direction,"_",Sys.Date(),".csv"))

# -----------------------------------------------
# Contribution TumC  (for each interaction) compared to other cells - analysis considering all the db (not only interactions >0 for TumC)
# ========================
test=matrix(nrow=length(rownames(db)), ncol=length(PC))
rownames(test)=db.name.couple[,1]
colnames(test)=colnames(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], couple="CD70 / CD27", db=db, plot=F))
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
LR.mat=as.data.frame(read.csv(paste0(paste0("analyses/ICELLNET_RNA/ICELLNET_sumScores_",direction,"_ALL_CELLS_2023-08-30.csv")))) 
rownames(LR.mat)=LR.mat$X
head(LR.mat)
LR.mat=LR.mat[,-c(1, 6, 13)] # do not forget to remove X and Prolif, and Epith, that are contaminated by Tumor cells-> to study TumC specificity

test=int.spe(LR.mat, CoI="ccRCC", thresh=thresh)
#write.csv(test, paste0(work.dir,output.dir,'ICELLNET_', assay, "/LR_spé_tumC_",direction,"_thresh_",thresh,"_ALL_CELLS_", Sys.Date(),".csv"))
test
dim(test)

DEG=as.data.frame(read.delim("analyses/DEG_TumC_Tum_vs_Hty_inDB_logFC0.25_minpct0.1_padj0.05.csv", sep = ",", header = T))
DEG.gene=dplyr::filter(DEG, DEG$avg_log2FC>0.25 & as.numeric(DEG$p_val_adj)<0.05) %>% pull(X)

if (direction=="out"){db.DEG <- filter(db, db$`Ligand 1` %in% DEG.gene | db$`Ligand 2` %in% DEG.gene)}
if (direction=="in"){db.DEG <- filter(db, db$`Receptor 1` %in% DEG.gene| db$`Receptor 2` %in% DEG.gene | db$`Receptor 3` %in% DEG.gene)}
DEG.group= icellnet::name.lr.couple(db.DEG, type="Family")[,1]

test2=test %>% filter(rownames(test) %in% DEG.group)
test2 
dim(test2)
setdiff(rownames(test), rownames(test2))

write.csv(test2, paste0(work.dir,output.dir,'ICELLNET_', assay, "/LR_spé_tumC_DEG",direction,"_thresh_",thresh,"_ALL_CELLS_", Sys.Date(),".csv"))
