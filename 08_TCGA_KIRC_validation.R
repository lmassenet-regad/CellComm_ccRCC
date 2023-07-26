# Created by Lucile Massenet-Regad, 16/09/2019


library(dplyr)
library(corrplot)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(FactoMineR)
library(dendextend)
library(gridExtra)
library(DESeq2)
library(reshape2)


# 1 - Create metadata file from downloaded info
#--------------------------------------------

#Read metadata, clinical data and sample data sheet
sample.data=read.csv("~/TCGA/Data/Downloaded_used/gdc_sample_sheet.2019-07-25.tsv", sep="\t")
clinical.csv=read.delim("~/TCGA/Data/Downloaded_used/clinical_data_downloaded_TCGA_portal_json_convertcsv_v2.csv", sep=";", header = T)
metadata.csv=read.delim("~/TCGA/Data/Downloaded_used/metadata_downloaded_TCGA_portal_json_convertcsv.csv", sep="\t", header = T)


# define a unique ID (future colnames of the matrix)
length(unique(as.factor(sample.data$Case.ID))) #9055 --> nb of patients
duplicate.patient=sample.data$Case.ID[which(duplicated(as.factor(sample.data$Case.ID))==T)] # patient with 2 samples (ex : control/tumor, or two parts of the same tumor, etc..) 
length(duplicate.patient) #790 --> OK (9055+790 = 9845 files in total)

length(unique(as.factor(sample.data$Sample.ID)))
duplicate.sampleID=sample.data$Sample.ID[duplicated(as.factor(sample.data$Sample.ID))] # different parts of tumor from same patient ?
sample.data[which(sample.data$Sample.ID==duplicate.sampleID[1]),]
length(duplicate.sampleID) #35 --> OK nb of files in total

# create target file
metadata.csv$File.Name=metadata.csv$file_name
metadata.csv$case_id=metadata.csv$associated_entities.0.case_id
target.int=dplyr::left_join(sample.data, metadata.csv, by="File.Name")
target.all.param=dplyr::left_join(target.int, clinical.csv, by="case_id")
length(unique(as.factor(target.all.param$Sample.ID))) #9809
target.all.param$Sample.ID=as.factor(make.unique(as.character(target.all.param$Sample.ID), sep="-")) # unique sample ID

param=c(colnames(sample.data), c("case_id", "treatments.treatment_type", "treatments.state", "treatments.therapeutic_agents", "treatments.treatment_or_therapy",
                                 "treatments.1.treatment_type", "diagnoses.icd_10_code","diagnoses.days_to_diagnosis","diagnoses.tissue_or_organ_of_origin","diagnoses.progression_or_recurrence",
                                 "diagnoses.prior_malignancy","diagnoses.synchronous_malignancy","diagnoses.site_of_resection_or_biopsy","diagnoses.days_to_last_follow_up","demographic.updated_datetime",
                                 "demographic.created_datetime","gender","state.y","submitter_id.y","year_of_birth", "race","days_to_birth", "ethnicity", "vital_status", "age_at_index", "year_of_death", "diagnoses.tumor_stage", "associated_entities.0.entity_submitter_id"))

target.all.reduced=target.all.param[,colnames(target.all.param)%in%param]

save(target.all.param, file="~/TCGA/Data/target_all_param_20191030.Rdata")
save(target.all.reduced, file="~/TCGA/Data/target_all_reduced_20191030.Rdata")


#  2 - Create data_matrix of the samples
#--------------------------------------------
rm(list=ls())
load("target_all_param_20191030.Rdata")
load( "target_all_reduced_20191030.Rdata")
target=target.all.reduced
rownames(target)=target$Sample.ID
name.fileID=target.all.reduced$File.ID

setwd("~/TCGA_Raw_data/")
file=read.table(paste0("Data_files/",target$File.ID[1],"/",target$File.Name[1]))
colnames(file)=c("EnsemblID",as.character(target$Sample.ID[1]))
data=file
for (i in 2:length(target$Sample.ID)){
  file=read.table(paste0("Data_files/",target$File.ID[i],"/",target$File.Name[i]))
  colnames(file)=c("EnsemblID",as.character(target$Sample.ID[i]))
  data=left_join(data,file, by="EnsemblID")
}

save(data, file=paste0("~/TCGA_Raw_data/TCGA_allsamples_rawdata.Rdata"))

# 3 - Filtering & Normalisation by DESeq2
#--------------------------------------------
rownames(data)=data$EnsemblID
data=data[,-1]

rownames(target)=target$Sample.ID #unique needed! 
data=data[,rownames(target)]

#filtering
filter<-apply(data, 1 , function(x) length(x[x>5])>=0.25*dim(data)[2])
data.filtered = data[filter,]

# Create DESeq2 object
library(DESeq2)

#to check : -->OK
all(rownames(target)%in%colnames(data))
all(rownames(target)==colnames(data))
#if FALSE, you have to reorder the data
#data=data[,rownames(target)]
#all(rownames(target)==colnames(data))


dds <- DESeqDataSetFromMatrix(countData = data.filtered,
                              colData = target,
                              design = ~1) #Rq : design not used to do normalisation steps
dds

#estimate size factors to have normalised counts
dds=estimateSizeFactors(dds)
norm.data=counts(dds, normalized=TRUE)
sum(duplicated(rownames(norm.data))) #no duplicate


# 4 - Add gene names and save file
#--------------------------------------------
# Convert EnsemblID to gene symbols
res.sig=data.frame("ensemblID"=rownames(norm.data))
res.sig$ensembl <- gsub('\\..+$', '', res.sig$ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#to see what you can retrieve --> listAttributes(ensembl)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res.sig$ensembl,
                  mart = ensembl )

idx <- match( res.sig$ensembl, genemap$ensembl_gene_id )
res.sig$entrez <- genemap$entrezgene_id[ idx ]
res.sig$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

all(rownames(norm.data)==res.sig$ensemblID)
all(!duplicated(res.sig$ensemblID))
all(!duplicated(res.sig$hgnc_symbol)) 
sum(!duplicated(res.sig$hgnc_symbol))
sum(is.na(res.sig$hgnc_symbol))

duplicate=res.sig$ensemblID[(duplicated(res.sig$hgnc_symbol) | is.na(res.sig$hgnc_symbol))]
res.sig=res.sig %>% filter(!(rownames(res.sig)%in%duplicate) )
res.sig=res.sig[!is.na(res.sig$hgnc_symbol),]

norm.data=as.data.frame(norm.data) %>% filter(!(rownames(norm.data)%in%duplicate))

all(rownames(norm.data)==res.sig$ensemblID)
rownames(norm.data)=res.sig$hgnc_symbol

# #save normalised data
save(norm.data, file="TCGA_all_samples_normDESeq2.Rdata")


# 5 - TCGA communication gene expression analysis - for Supplementary Figure S3D
#--------------------------------------------
db_list=as.data.frame(readxl::read_excel("~/Desktop/LR_selected_DB.xlsx"))
rownames(db_list)= db_list$Pair
data.scale= t(apply(matrix, MARGIN=1,  FUN = function(x) (x/max(x)))) 

# COMPUTE multiplication for each LR pair
mult=matrix(nrow = dim(db_list)[1],  ncol = dim(data.scale)[2])
colnames(mult)=colnames(data.scale)
rownames(mult)=rownames(db_list)

for (j in 1:dim(db_list)[1]){
  pair=db_list[j,"Pair"]
  ligand=db_list[j,"Ligand"]
  if (is.na(db_list[j,"Rec_2"])){
    receptor=db_list[j,"Rec_1"]
  } else { receptor=c(db_list[j,"Rec_1"],  db_list[j,"Rec_2"])}
  data_meta = data.scale[which(rownames(data.scale) %in% c(ligand, receptor)), ]
  if (is.na(db_list[j,"Rec_2"])){
    mult[pair,]=data_meta[which(rownames(data_meta)== ligand), ]* data_meta[which(rownames(data_meta)== db_list[j,"Rec_1"]),]
  } else {
    mult[pair,]=data_meta[which(rownames(data_meta)== ligand), ]* sqrt(data_meta[which(rownames(data_meta)== db_list[j,"Rec_1"]),] * data_meta[which(rownames(data_meta)== db_list[j,"Rec_2"]),])
  }
}

mult_df=as.data.frame(reshape2::melt(mult))
colnames(mult_df)=c("Pair", "Sample.ID", "value")
mult_df=dplyr::right_join(mult_df, target[,c("Sample.ID", "Sample.Type")])
mult_df$Sample.Type=as.character(mult_df$Sample.Type)

# Compute logFC 
test2= mult_df %>% group_by(Pair, Sample.Type, .add=T) %>% summarise_if(is.numeric, mean, na.rm=TRUE)
test3 = test2 %>% reshape2::dcast(Pair ~ Sample.Type) 
test3$FC = test3$`Primary Tumor` / test3$`Solid Tissue Normal`
test3$log2FC = log2(test3$FC)

#compute pvalue
stat= ggpubr::compare_means(value ~ Sample.Type, data = mult_df, group.by = "Pair", method = "wilcox.test", p.adjust.method = "bonferroni")
all(stat$Pair == test3$Pair) #check
test3$padj=stat$p.adj

#volcano plot
EnhancedVolcano::EnhancedVolcano(test3,lab = test3$Pair,
                                 x = 'log2FC',
                                 y = 'padj', 
                                 FCcutoff = 1, legendPosition = "top", labSize = 5 , pointSize = 4 , pCutoff = 0.001, 
                                 drawConnectors = TRUE, widthConnectors = 0.25,  axisLabSize = 10)
