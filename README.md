# CellComm_ccRCC

This repository host the R scripts used in the paper "Systems analysis of cell communication identifies angiogenin-dependent modulation of tumor proliferation and inflammation in clear cell renal cell carcinoma".

## Index:
ICELLNET ligand/receptor interaction database used in the scripts can be found on ICELLNET github.
List of cytokines, receptors, checkpoints and chemokines is provided as xlsx file.

### 01_TCGA_analysis: 
Script from downloaded TCGA data. Steps include matrix creation, data normalisation with DESeq2, and analysis until Figure 1A generation.

### 02a_scRNAseq_processing:
Input: cellranger output files - exemple with one of the sample
Output: processed individual scRNAseq data. 

### 02b_scRNAseq_integration_clustering:
Input: processed individual scRNAseq data. 
Output: one merged scRNAseq with integrated samples and consensus clustering 
This script also contains code for the Figure 1C-1D and Supp. Figure 2A-C

### 03_scRNAseq_Tumor_cell_analysis:
Input: merged scRNAseq with integrated samples and consensus clustering 
- Clustering of tumor cells (TumC1 and TumC2) + differential expression analysis (Figure 2A-C)
- Analysis of differential expression between TumC2 and proximal tubules (PT) from juxtatumoral samples (Figure 2D-F)
- Identification of tumor cell specific genes compared to all other clusters in the dataset, in tumoral samples (Figure 2G)

### 04_scRNAseq_ICELLNET:
Input: merged scRNAseq with integrated samples and consensus clustering +  ICELLNET ligand/receptor interaction database
Cell cell communication analysis using [ICELLNET package](https://github.com/soumelis-lab/ICELLNET)


### 05a_scRNAseq_Young_data:
Reanalysis of [Young et al. data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6104812/)

### 05b_ICELLNET_Young_data
Cell cell communication analysis to identify tumor cell specific interactions.


### 06a_scRNAseq_Zhang_data:
Reanalysis of [Zhang et al. data](https://www.pnas.org/doi/10.1073/pnas.2103240118)

### 06b_ICELLNET_Zhang_data
Cell cell communication analysis to identify tumor cell specific interactions.

### 07_ICELLNET_Intersection_results:
Provide code used to generate Figure 3 and Supp. Figures 3.