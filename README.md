# Experimental isolation of neural stem cells: A reference benchmark that summarizes the impact of the isolation method on the heterogeneity of the isolated cell subpopulations based on scRNA-seq data analysis.

## Project Overview

The goal of this project was to integrate scRNA-seq datasets from various Neural Stem Cell (NSCs) studies in order to investigate how the experimental isolation method of choice impacts the isolated NSC subpopulations. NSCs experiments are usually performed in mouse lines that have been uniquely labeled with gene promoters. Under the valid assumption that NSCs are considerably heterogeneous, the gene promoter choice could affect the isolated NSC subpopulations.


## Reference Materials

The publications that can be found in this repository and were used for the integrated analysis are the following:

- **Kalamakis et al., 2019**
- **Mizrak et al. 2020**
- **Xie et al., 2014**
- **Dulken et al., 2017**
- **Hamed et al., 2022**


## Dataset Description

- **Description**: All datasets can be found through the online material of the publications mentioned in the Reference Materials.



## Methodology

- **Data Collection**: scRNA-seq datasets from published studies focusing on neural stem cells and neurogenesis.
- **Preprocessing**: Removed, filtered, and transformed data in order to ensure compatibility across datasets for integration based on the scope of the research. Only cells with at least 200 and less that 7500 genes were kept.
- **Integration**: Utilized the `Seurat` package in R to perform data integration and analysis, including:
  - Identification of highly variable genes across cells.
  - Identification of the different cell subpopulations based on gene expression variability and known biological markers from published literature.
  - Dimensionality reduction using PCA and UMAP.
  - Gene expression visualizations of the identified cell subpopulations (Violin plots, Volcano plots).


## Script Details

- `integration_script.R`: All preprocessing and analysis steps of the scRNA-seq datasets are found in this script. 

