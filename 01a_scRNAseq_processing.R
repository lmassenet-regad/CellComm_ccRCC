#Author: Lucile Massenet-Regad 
#Created 2021-05-17
#inspired from pipeline script at SincellTE2022

rm(list=ls())
work.dir=work.dir
file.name="p27_Tumoral"
sample.name <- "p27_Tum"

# Set working directory + create output directory
setwd(dir = work.dir)
dir.create(path = paste(c(work.dir,'analyses/1_Preprocessing_data/', sample.name), collapse = '/'), 
           recursive = TRUE, showWarnings = FALSE)
input.dir <- paste0("data/Raw/",file.name,"/raw_feature_bc_matrix/")
output.dir <- paste0("analyses/1_Preprocessing_data/", sample.name,"/")

#Load Library
library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)

###### Parameters
#------------------------------------------------------------------------
# Computational Parameters
the.seed <- 1337L

# Analysis Parameters
## Barcode-level QC cell
pcmito.max <- 20
min.features <- 200
min.counts <- 500
## Feature-level QC
min.cells <- 2
#Dimension reduction and normalisation
features.n <- 3000
dim.max <- 50 # to compute PCA
n.dim=30 # for clustering / UMAP
res=0.2 # resolution for clustering
norm="LogNormalize" #logNormalize or SCT

###### Empty drop removal + creation of Seurat object
#------------------------------------------------------------------------
# Loading raw count matrix
scmat <- Seurat::Read10X(input.dir) 

# Removing empty droplets
bc_rank <- DropletUtils::emptyDrops(m = scmat)
scmat_filtered <- scmat[, which(bc_rank$FDR < 1E-03)]

# Kneeplot
## Make the dataframe with all droplets (before filtering) and the number of UMI for each droplet, and if the droplets are filtered or not by emptydrops
nb_umi_by_barcode <- data.frame(nb_umi = Matrix::colSums(scmat), barcodes = colnames(scmat))
nb_umi_by_barcode <- nb_umi_by_barcode %>% arrange(desc(nb_umi)) %>% dplyr::mutate(num_barcode = seq.int(ncol(scmat)))
nb_umi_by_barcode$droplets_state <- "Empty Droplets"
nb_umi_by_barcode[nb_umi_by_barcode$barcodes %in% colnames(scmat_filtered), "droplets_state"] <- "Full Droplets"

nb_umi_by_barcode=nb_umi_by_barcode[seq(1, nrow(nb_umi_by_barcode), 20), ] # only for graph visualisation - to limit nb of points 

# Creation of the Seurat object
# Rq: Filtering features with ultra-low expression (less than 5 cells)
seurat <- Seurat::CreateSeuratObject(counts = scmat_filtered, project = sample.name, min.cells=min.cells)
seurat

# ## Draw kneeplot
pdf(paste0(output.dir, "01_Kneeplot_EmptyDrops_",sample.name, ".pdf"), useDingbats = F)
ggplot2::ggplot(nb_umi_by_barcode, ggplot2::aes(y = nb_umi, x = num_barcode, color = droplets_state)) +
  ggplot2::geom_point() +
  ggplot2::theme(legend.title = ggplot2::element_blank()) +
  ggplot2::scale_y_log10(name = "Number of UMI by droplet (log scale)") +
  ggplot2::scale_x_log10(name = "Droplet rank (log scale)") +
  ggplot2::expand_limits(x = 0, y = 0) + #because log 
  ggplot2::scale_colour_manual(values = c("cyan3","royalblue4"), guide = 'legend') +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c(0, 0), shape = c(16, 16))))+
  ggplot2::ggtitle(paste0("Number of cells = ", dim(seurat)[[2]]))
dev.off()

# Cleaning
rm(bc_rank,nb_umi_by_barcode,scmat,scmat_filtered)
invisible(gc())

###### Basic metrics-QC
#------------------------------------------------------------------------
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-") #mitochondrial genes
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[LS][[:digit:]]") #ribosomal genes

pdf(paste0(output.dir, "01_QC_Simple_metrics_beforeFiltering_", sample.name, ".pdf"), useDingbats = F)
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



# Filter
seurat <- subset(seurat, subset = nFeature_RNA > min.features & nCount_RNA > min.counts  & percent.mt < pcmito.max)  

pdf(paste0(output.dir, "01_QC_Simple_metrics_afterFiltering_", sample.name, ".pdf"), useDingbats = F)
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

# union of scDblFinder & scds
seurat$doublets_consensus.class <- seurat$scDblFinder.class | seurat$scds.class 
print('Consensus doublets :')
print(table(seurat$doublets_consensus.class))

df_S.cycle_doublets <- data.frame(Phase = seurat$Phase, doublets = seurat$doublets_consensus.class)
print(table(df_S.cycle_doublets))
print(paste0(round((nrow(df_S.cycle_doublets[df_S.cycle_doublets$Phase == "G2M" &
                                               df_S.cycle_doublets$doublets == TRUE,])/
                      nrow(df_S.cycle_doublets[df_S.cycle_doublets$doublets == TRUE,])*100),2),
             "% of doublets are in G2M phase by Seurat."))

# plot PCA 
seurat.tmp = seurat #  temporary object only to visualise how doublet and cell cycling are.
seurat.tmp <- NormalizeData(seurat.tmp, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.tmp <- FindVariableFeatures(seurat.tmp, selection.method = "vst", nfeatures = 2000)
seurat.tmp <- ScaleData(seurat.tmp)
seurat.tmp <- RunPCA(seurat.tmp, npcs = 100, verbose = FALSE)
ElbowPlot(seurat.tmp, ndims = 100) 
seurat.tmp <- RunUMAP(seurat.tmp, dims = 1:n.dim, seed.use = the.seed)

pdf(paste0(output.dir, "02_QC_CellCycle_Doublets_PCA_temp_",sample.name, ".pdf"), useDingbats = F)
DimPlot(seurat.tmp, reduction = "umap", group.by = "Phase") + ggplot2::ggtitle("Cell Phase (Seurat)")
DimPlot(seurat.tmp, group.by = "scDblFinder.class") + ggplot2::ggtitle("Cell doublets (scDblFinder)") 
DimPlot(seurat.tmp, group.by = "scds.class") + ggplot2::ggtitle("Cell doublets (scds)") 
DimPlot(seurat.tmp, group.by = "ManDblt") + ggplot2::ggtitle("Manual markers incompatibility -> likely doublets") 
DimPlot(seurat.tmp, group.by = "doublets_consensus.class") + ggplot2::ggtitle("Cell doublets (union)") 
dev.off()

###### Filtering + Normalisation + Dimension reduction + Clustering
#------------------------------------------------------------------------
# Filtering cell doublets
seurat <- seurat[, !seurat$doublets_consensus.class]

# Graph - QC - after Filtering
pdf(paste0(output.dir, "02_QC_After_filtering_Mitochodrial_genes_",sample.name, ".pdf"), useDingbats = F)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb" ), 
        ncol = 4, pt.size=0.01, group.by =  "orig.ident")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mt", feature2 = "percent.rb")
dev.off()

#Normalisation
if (norm=="LogNormalize"){
  seurat <- NormalizeData( seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat <- FindVariableFeatures( seurat, selection.method = "vst", nfeatures = 2000)
  seurat <- ScaleData( seurat)
}else if (norm=="SCT"){
  #Results are saved in a new assay (named SCT by default) with counts being (corrected) counts, data being log1p(counts), scale.data being pearson residuals;
  seurat <- suppressWarnings(Seurat::SCTransform(object = seurat, assay = "RNA", seed.use = the.seed, variable.features.n = features.n, return.only.var.genes = TRUE)) #alternative to the NormalizeData, FindVariableFeatures, ScaleData
}

# Identify the top 20 HVGs
HVG <- Seurat::VariableFeatures(seurat)
top20 <- head(HVG, 20)
top20
# Quantify mitochondrial/ribosomal genes in HVGs:
nb_mito_in_HVG <- length(grep("^MT-", HVG, value = TRUE))
print(paste0("There are ", nb_mito_in_HVG, " mitochondrial genes in the highly variable genes (", round(x = (nb_mito_in_HVG/length(HVG)*100), 2), "%)."))
nb_ribo_in_HVG <- length(grep("^RP[LS][[:digit:]]", HVG, value = TRUE))
print(paste0("There are ", nb_ribo_in_HVG, " ribosomal genes in the highly variable genes (", round(x = (nb_ribo_in_HVG/length(HVG)*100), 2), "%)."))

# PCA
seurat <- Seurat::RunPCA(object = seurat, verbose = FALSE, npcs = dim.max, seed.use = the.seed)
ElbowPlot(seurat, ndims = dim.max)

###### Clustering
#------------------------------------------------------------------------
seurat <- FindNeighbors(seurat, dims = 1:n.dim)

#Silhouette
pdf(paste0(output.dir, "03_Clustering_Silhouette_PC", n.dim,"_",sample.name, ".pdf"), useDingbats = F)
for (res2 in seq(0.1,1.4,0.1)){
  print(paste0("resolution=",res2))
  seurat <- FindClusters(object = seurat, resolution = res2, force.recalc = T)
  seurat[[paste0("ClusterNames_", res)]] <- Idents(object = seurat)
  si <- cluster::silhouette(x= as.numeric(Idents(seurat)),
                            dist = dist(FetchData(seurat, vars=as.vector(paste0("PC_", c(1:n.dim))) )))
  plot(si, cex.names = 0.5)
}
dev.off()

# Clustree
pdf(paste0(output.dir, "03_Clustering_Clustree_PC", n.dim,"_",sample.name, ".pdf"), useDingbats = F, height = 15,width=10)
exprs <- "data"
prefix = "SCT_snn_res."
args <- list()
gene_names <- rownames(seurat@assays$RNA@data)
for (node_aes in c("node_colour", "node_size", "node_alpha")) {
  if (node_aes %in% names(args)) {
    node_aes_value <- args[[node_aes]]
    if (node_aes_value %in% gene_names) {
      aes_name <- paste0(exprs, "_", node_aes_value)
      seurat@meta.data[aes_name] <-
        slot(seurat, exprs)[node_aes_value, ]
      args[[node_aes]] <- aes_name
    }
  }
}
args$x <- seurat@meta.data
args$prefix <- prefix
do.call(clustree, args)
dev.off()

## STOP 1 - check resolution choice ! 
# Clustering : choice of resolutions
seurat <- FindClusters(seurat, resolution = res)
seurat <-RunUMAP(seurat, dims = 1:n.dim, seed.use = the.seed)
DimPlot(seurat, reduction = "umap", label = T)

# Save UMAP
pdf(paste0(output.dir, "04_UMAP_",res,"_PC", n.dim, "_",sample.name, ".pdf"), useDingbats = F)
DimPlot(seurat, reduction = "umap", group.by = paste0("SCT_snn_res.", res), label = T, label.size = 4) + labs(title= paste("resolution", res))
dev.off()

#QC after clustering
pdf(paste0(output.dir, "04_QC_After_clustering_",sample.name, ".pdf"), useDingbats = F)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb" ), 
        ncol = 4, pt.size=0.01, group.by =  "orig.ident")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mt", feature2 = "percent.rb")
DimPlot(seurat, reduction = "umap", group.by = "Phase", label = T, label.size = 4) + labs(title= "QC cell cycle")
FeaturePlot(seurat, reduction = "umap", features = "nCount_RNA") + labs(title= "QC nCount")
FeaturePlot(seurat, reduction = "umap", features = "nFeature_RNA") + labs(title= "QC nFeature")
dev.off()

###### Save seurat object
#------------------------------------------------------------------------
saveRDS(seurat, file = paste0("data/",sample.name,"_res", res,"_PC", dim.max, "_data_",Sys.Date(),".rds"))
