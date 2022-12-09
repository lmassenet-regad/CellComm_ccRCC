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


# 5 - TCGA communication gene analysis - bulk global communication score
#--------------------------------------------
rm(list=ls())
load("target_all_param_20191030.Rdata")
load( "target_all_reduced_20191030.Rdata")
load("TCGA_all_samples_normDESeq2.Rdata") #norm.data


target=target.all.reduced
rownames(target)=target$Sample.ID
target=filter(target.all.param,target.all.param$Sample.Type%in% c("Solid Tissue Normal","Primary Tumor") 
              & !(Project.ID %in% c("TCGA-DLBC", "TCGA-LGG","TCGA-OV", "TCGA-CESC", "TCGA-GBM", "TCGA-CESC", "TCGA-PAAD", "TCGA-SARC", "TCGA-SKCM"))) #no control of not enough controls (<10...)
rownames(target)=target$Sample.ID

head(norm.data[1:5,1:5])
data2=norm.data[,which(colnames(norm.data)%in%as.character(target$Sample.ID))]

# Extract cytokines and cytokines receptors genes from data matrix
cytokines=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cytokines")
receptors=readxl::read_excel("../Liste_molecules.xlsx", sheet="Cyt_receptors")
checkpoints=readxl::read_excel("../Liste_molecules.xlsx", sheet="Checkpoints")
chemo=readxl::read_excel("../Liste_molecules.xlsx", sheet="Chemokines")

#all cyt or rec 
cyt_oi=cytokines %>%  pull(Symbol)  %>% unique()
rec_oi=receptors %>%  pull(Symbol)  %>% unique()
check_oi= checkpoints %>%  pull(Symbol)  %>% unique()
chemo_oi=chemo %>%  pull(Symbol)  %>% unique()
all=unique(c(cyt_oi, rec_oi, check_oi,  chemo_oi))  # to change according to plot of interest

#Extract molecules genes 
data.cyt=as.data.frame(data2) %>% filter(rownames(data2) %in% all)
metagene.cyt=data.frame("Sample.ID"=colnames(data.cyt) ,"score"=colSums(data.cyt)) 
target$Project.ID= sapply(strsplit(as.character(target$Project.ID), split='-', fixed=TRUE), function(x) (x[2]))
metagene.cyt=dplyr::right_join(metagene.cyt,target[,c("Sample.ID","Project.ID", "Sample.Type")])
str(metagene.cyt)
summary(metagene.cyt)


#Exclude metastasis - keep only Solid tissue normal and Primary tumors
test=metagene.cyt%>%
  group_by(Project.ID, Sample.Type)%>% 
  summarise(Mean=mean(score), Median=median(score))

test2=test%>%
  group_by(Project.ID) %>%
  mutate(deltaMed = Median[Sample.Type == 'Primary Tumor'] - Median[Sample.Type == 'Solid Tissue Normal']) %>%
  mutate(deltaMean = Mean[Sample.Type == 'Primary Tumor'] - Mean[Sample.Type == 'Solid Tissue Normal'])%>%
  filter(Sample.Type=="Primary Tumor")

# Figure 1A
ggplot(test2, aes(x=deltaMed, y=Median, color=Project.ID, label=Project.ID)) + 
  geom_point() + ggrepel::geom_text_repel(size=3) + theme_classic() + theme(legend.position = "none") + 
  labs(x = "Score difference between tumor and healthy", y = "Median of communication score",  title = "Comm mol") 

