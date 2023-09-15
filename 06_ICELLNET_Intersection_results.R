##Comparison of specific interaction + DEG on tumor side, in the 3 datasets

library(ggplot2)
library(VennDiagram)
library(dplyr)

setwd(path)

library(ComplexHeatmap)
col_fun = circlize::colorRamp2(c(0, 100), c("white", "red"))


#### Load results for each dataset ####

# direction = in
zhang.in = as.data.frame(read.csv("../../Zhang et al, 2021, PNAS/analyses/ICELLNET_Harmony_RNA/LR_spé_tumC_DEGin_thresh_1.5_ALL_CELLS_2023-08-30.csv"))
young.in = as.data.frame(read.csv("../../Young et al_2018_Kidney/analyses/3_Downstream_analysis/ICELLNET_RNA/LR_spé_tumC_DEGin_thresh_1.5_ALL_CELLS_2023-04-20.csv"))
soumelis.in =as.data.frame(read.csv("../analyses/3_Downstream_analysis/ICELLNET_RNA/LR_spé_tumC_DEGin_thresh_1.5_ALL_CELLS_2023-06-14.csv"))

zhang.in = zhang.in %>% filter(ratio >1.5) %>% pull(X) %>% as.vector()
young.in= young.in %>% filter( ratio >1.5) %>% pull(X) %>% as.vector()
soumelis.in = soumelis.in %>% filter( ratio >1.5) %>% pull(X) %>% as.vector()

# direction = out
zhang.out = as.data.frame(read.csv("../../Zhang et al, 2021, PNAS/analyses/ICELLNET_Harmony_RNA/LR_spé_tumC_DEGout_thresh_1.5_ALL_CELLS_2023-08-30.csv"))
young.out = as.data.frame(read.csv("../../Young et al_2018_Kidney/analyses/3_Downstream_analysis/ICELLNET_RNA/LR_spé_tumC_DEGout_thresh_1.5_ALL_CELLS_2023-04-20.csv"))
soumelis.out = as.data.frame(read.csv("../analyses/3_Downstream_analysis/ICELLNET_RNA/LR_spé_tumC_DEGout_thresh_1.5_ALL_CELLS_2023-06-14.csv"))

zhang.out = zhang.out %>% filter(ratio >1.5) %>% pull(X) %>% as.vector()
young.out= young.out%>% filter( ratio >1.5) %>% pull(X) %>% as.vector()
soumelis.out = soumelis.out %>% filter( ratio >1.5) %>% pull(X) %>% as.vector()

#### Venn diagramm - IN ####
library(ggVennDiagram)

x = list("Zhang"= zhang.in, "Young"= young.in, "Soumelis"=soumelis.in)
ggVennDiagram(x, label_alpha=0, label="count", set_size = 3) +  scale_fill_gradient(low="white",high = "white") + scale_color_brewer(palette = "Set1") +
  ggtitle("IN") 


overlap=calculate.overlap(x=list(zhang.in, young.in, soumelis.in))
# names overlap : "a5"  "a2" "a4" "a6" "a1" "a3" "a7"
#meaning = c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")
names(overlap) <- c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")
overlap %>% unlist() %>% as.data.frame() %>% write.csv( "IN_DEG_overlap_VD_20230906.csv")

#### Venn diagramm - OUT ####

x2 = list("Zhang"= zhang.out, "Young"= young.out, "Soumelis"=soumelis.out)
ggVennDiagram(x2, label_alpha=0, label="count") +  
  scale_fill_gradient(low="white",high = "white") + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("OUT")


overlap=calculate.overlap(x=list(zhang.out, young.out, soumelis.out))
# names overlap : "a5"  "a2" "a4" "a6" "a1" "a3" "a7"
#meaning = c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")
names(overlap) <- c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")
overlap %>% unlist() %>% as.data.frame() %>% write.csv( "OUT_DEG_overlap_VD_20230906.csv")


##### INTERPRETATION & VISUALISATION #####


#### Specific ICELLNET plot - IN ####  - 2022 07 20 - updated 2022 10 27

LR.in=read.csv("IN_DEG_overlap_VD_20230906.csv")
overlap_in=calculate.overlap(x=list(zhang.in, young.in, soumelis.in))
names(overlap_in) <- c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")


icell.zhang.in=as.data.frame(read.csv("../../Zhang et al, 2021, PNAS/analyses/ICELLNET_Harmony_RNA/ICELLNET_TumC_Scores_in_2023-08-30.csv"))
icell.young.in = as.data.frame(read.csv("../../Young et al_2018_Kidney/analyses/3_Downstream_analysis/ICELLNET_RNA/ICELLNET_LR_TumC_Scores_in_2023-04-20.csv"))
icell.soumelis.in =as.data.frame(read.csv("../analyses/3_Downstream_analysis/ICELLNET_RNA/ICELLNET_LR_TumC_Scores_in_2023-04-18.csv"))
rownames(icell.zhang.in)=icell.zhang.in$X
icell.zhang.in=icell.zhang.in[,-1]
rownames(icell.young.in)=icell.young.in$X
icell.young.in=icell.young.in[,-1]
rownames(icell.soumelis.in)=icell.soumelis.in$X
icell.soumelis.in=icell.soumelis.in[,-1]

# Soumelis dataset
icell.soumelis.in2=icell.soumelis.in[c(overlap_in$ZYS, overlap_in$ZS, overlap_in$YS),]
icell.soumelis.in2=icell.soumelis.in2[which(rowSums(icell.soumelis.in2)>0),]
Heatmap(as.matrix(icell.soumelis.in2), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("soumelis-in"), col = col_fun,
        row_names_gp = gpar(fontsize =10))

# Zhang dataset --careful numbers
icell.zhang.in2=icell.zhang.in[c(overlap_in$ZYS, overlap_in$ZS, overlap_in$ZY),]
Heatmap(as.matrix(icell.zhang.in2), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("zhang-in_shared"), col = col_fun,
        row_names_gp = gpar(fontsize =10))

# Young dataset -- careful numbers
icell.young.in3=icell.young.in[c(overlap_in$ZYS, overlap_in$ZY, overlap_in$YS),]
Heatmap(as.matrix(icell.young.in3), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("young-in_shared"), col = col_fun,
        row_names_gp = gpar(fontsize =10))


#### Specific ICELLNET plot - OUT #### - 2022 10 27
LR.out=read.csv("OUT_DEG_overlap_VD_20230420.csv")
overlap_out=calculate.overlap(x=list(zhang.out, young.out, soumelis.out))
names(overlap_out) <- c("ZYS","ZY", "ZS", "YS", "Z", "Y", "S")

icell.zhang.out=as.data.frame(read.csv("../../Zhang et al, 2021, PNAS/analyses/ICELLNET_Harmony_RNA/ICELLNET_TumC_Scores_out_2023-08-30.csv"))
icell.young.out = as.data.frame(read.csv("../../Young et al_2018_Kidney/analyses/3_Downstream_analysis/ICELLNET_RNA/ICELLNET_LR_TumC_Scores_out_2023-04-20.csv"))
icell.soumelis.out =as.data.frame(read.csv("../analyses/3_Downstream_analysis/ICELLNET_RNA/ICELLNET_LR_TumC_Scores_out_2023-04-18.csv"))
rownames(icell.zhang.out)=icell.zhang.out$X
icell.zhang.out=icell.zhang.out[,-1]
rownames(icell.young.out)=icell.young.out$X
icell.young.out=icell.young.out[,-1]
rownames(icell.soumelis.out)=icell.soumelis.out$X
icell.soumelis.out=icell.soumelis.out[,-1]

icell.soumelis.out2=icell.soumelis.out[c(overlap_out$ZYS, overlap_out$ZS, overlap_out$YS),]
Heatmap(as.matrix(icell.soumelis.out2), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("soumelis-out-shared"), col = col_fun,
        row_names_gp = gpar(fontsize =10))


icell.zhang.out2=icell.zhang.out[c(overlap_out$ZYS, overlap_out$ZS, overlap_out$ZY),]
Heatmap(as.matrix(icell.zhang.out2), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("zhang-out"), col = col_fun,
        row_names_gp = gpar(fontsize =10))


icell.young.out2=icell.young.out[c(overlap_out$ZYS, overlap_out$ZY, overlap_out$YS),]
Heatmap(as.matrix(icell.young.out2), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("young-out"), col = col_fun,
        row_names_gp = gpar(fontsize =10))



icell.young.in3=icell.young.in[rownames(icell.young.in)%in%LR.in[1:8,2],]
icell.young.in3=icell.young.in3[LR.in[1:8,2],]
icell.young.in3=icell.young.in3[which(rowSums(icell.young.in3)>0),]
Heatmap(as.matrix(icell.young.in3), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("young-in_shared"), col = col_fun,
        row_names_gp = gpar(fontsize =10))

icell.young.in3=icell.young.in[rownames(icell.young.in)%in%LR.in[9:12,2],]
icell.young.in3=icell.young.in3[which(rowSums(icell.young.in3)>0),]
Heatmap(as.matrix(icell.young.in3), cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,
        column_title =paste0("young-in_spe"), col = col_fun,
        row_names_gp = gpar(fontsize =10))





