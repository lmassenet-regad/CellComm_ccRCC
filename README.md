# CellComm_ccRCC

This repository hosts the R scripts used in the paper "Large-scale analysis of cell-cell communication reveals angiogenin-dependent tumor progression in clear cell renal cell carcinoma".

## Index:
ICELLNET ligand/receptor interaction database used in the scripts can be found on [ICELLNET github](https://github.com/soumelis-lab/ICELLNET).

### 01a_scRNAseq_processing:
Input: cellranger output files - exemple with one of the sample. Data can be retrived [here from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222703).



Output: processed individual scRNAseq data. 

### 01b_scRNAseq_integration_clustering:
Input: processed individual scRNAseq data 


Output: one merged scRNAseq with integrated samples and consensus clustering.
This script also contains code for the Figure 1B-1D and Supp. Figure S1.

### 02_scRNAseq_cancer_cell_analysis:
Input: merged scRNAseq with integrated samples and consensus clustering 
- A: Clustering of cancer cells (ccRCC1 and ccRCC2) + differential expression analysis (Figure 2A-C, Supplementary Table S3A)
- B:  Analysis of differential expression between ccRCC2 and proximal tubules (PT) from juxtatumoral samples (Supplementary Table S3B)
- C: Identification of cancer cell specific genes compared to all other clusters in the dataset, in tumoral samples (Figure 2D)

### 03_scRNAseq_ICELLNET:
Input: merged scRNAseq with integrated samples and consensus clustering +  ICELLNET ligand/receptor interaction database


Cell cell communication analysis using [ICELLNET package](https://github.com/soumelis-lab/ICELLNET)

### 04a_scRNAseq_Young_data:
Reanalysis of [Young et al. data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6104812/)

### 04b_ICELLNET_Young_data:
Cell cell communication analysis to identify tumor cell specific interactions.


### 05a_scRNAseq_Zhang_data:
Reanalysis of [Zhang et al. data](https://www.pnas.org/doi/10.1073/pnas.2103240118)

### 05b_ICELLNET_Zhang_data:
Cell cell communication analysis to identify tumor cell specific interactions.

### 06_ICELLNET_Intersection_results:
Provide code used to generate Figure 3 and Supp. Figures 3.

### 07_CellChat_validation:
Cell-cell communication analysis performed on original data, using [CellChat](https://github.com/sqjin/CellChat) inference method 

### 08_TCGA_KIRC_validation: 
Script from downloaded TCGA data. Steps include matrix creation, data normalisation with DESeq2, and analysis for Supplementary Figure S3D generation.
