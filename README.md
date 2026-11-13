# ACTA2 Identifies a Macrophage-Depleted, Stromal-Excluded Leiomyosarcoma Subtype with Compartmentalized Immune Checkpoints and Therapeutic Stratification

This repository contains the complete R scripts used to generate the main and supplementary figures and supplementary tables for our manuscript. 

## Study Design & Core Analyses
To ensure robust and reproducible findings, this study integrated **dual bulk transcriptomic cohorts** (TCGA-SARC and GSE21050) with **single-cell RNA sequencing** (GSE212527). The core analytical pipeline includes:
- **Dual-algorithm Deconvolution**: Cross-cohort validation using both ESTIMATE and xCell algorithms.
- **Single-cell Analysis**: Utilization of the standard Seurat pipeline (PCA for dimensionality reduction, Harmony for batch effect correction, and graph-based clustering), followed by inferCNV, AddModuleScore, and CellChat to identify cellular origins and ligand-receptor communication.
- **Partial Correlation Analysis**: Rigorous statistical adjustment to separate compositional (cell abundance) effects from cell-intrinsic regulation.

## Data Sources
- **Bulk RNA-seq**: TCGA-SARC and GEO (GSE21050). Processed data can be downloaded from GDC and GEO.
- **Single-cell RNA-seq**: GEO (GSE212527).

## Scripts and Outputs Overview

| Script | Main Analysis / Generated Figures & Tables |
| :--- | :--- |
| `gitLMS-TCGA-SARC.R` | **TCGA Bulk Cohort (Discovery)**:<br>• Spearman correlation heatmap (Fig. 1a)<br>• GSVA Hallmark pathways (Fig. 1b)<br>• xCell immune infiltration heatmap (Fig. 1c)<br>• ESTIMATE and tumor purity analysis (Fig. 2a-d)<br>• Ligand-receptor gene correlation with ACTA2 (Fig. 5k)<br>• Single-gene checkpoint analysis (Supplementary Table S3) |
| `gitLMS-GEO.R` | **GEO Bulk Cohort (Validation)**:<br>• Cross-cohort validation of ESTIMATE scores (Fig. 3a-d)<br>• Cross-cohort validation of xCell cell types (Fig. 4a-f)<br>• Cross-cohort single-gene checkpoint analysis (Fig. 6a)<br>• Partial correlation adjustment for macrophage abundance/purity (Fig. 6b)<br>• Compositional adjudication data (Supplementary Table S4)<br>• Partial correlation of ligand-receptor genes (Supplementary Table S6) |
| `gitLMS-scRNA.R` | **Single-cell RNA-seq Analysis**:<br>• QC metrics (Supplementary Fig. S1) and cell type annotation (Supplementary Fig. S2)<br>• UMAP visualization, CNV heatmap, ACTA2 expression distribution (Fig. 5a-5e)<br>• Macrophage content vs. differentiation score (Fig. 5f)<br>• CellChat forward and reverse ligand-receptor dot plots and circle plots (Fig. 5g-5j)<br>• Single-cell checkpoint expression dot plot (Fig. 6c)<br>• Ligand-receptor interactions (Supplementary Tables S1, S2)<br>• Checkpoint expression statistics across cell types (Supplementary Table S5) |

## Requirements
- R version 4.5.3 (or newer)
- Key packages: `Seurat`, `CellChat`, `inferCNV`, `ppcor`, `xCell`, `ESTIMATE`, `GSVA`, `hgu133plus2.db`

```R
# CRAN packages
install.packages(c("dplyr", "ggplot2", "ggpubr", "ggrepel", "patchwork", 
                   "tidyverse", "pheatmap", "ppcor", "cowplot"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor and specific packages
BiocManager::install(c("Seurat", "inferCNV", "CellChat", "GSVA", 
                       "ESTIMATE", "xCell", "hgu133plus2.db"))
```
## Usage & Reproducibility Notes

1. **Paths**: Please note that the scripts expect a specific project directory structure. The submitted scripts may contain legacy Windows absolute paths and must be edited (e.g., updating `setwd()` and input file paths) before execution on another computer.
2. **Intermediate Objects**: Large intermediate objects (e.g., `Seurat .rds` files, processed bulk matrices) are not included in this repository due to size limits. They can be reconstructed using the scripts provided or obtained from the corresponding author upon reasonable request.
3. **Bulk Validation**: Some cross-cohort plots (e.g., combining TCGA and GEO for Fig 6a and 6b) may require the generation of specific data frames in their respective scripts prior to plotting.
4. **Interpretation**: The code supports the computational associations reported in the manuscript.
5. **Execution Order**: We strongly recommend running the scripts in the following sequence: gitLMS-TCGA-SARC.R (Discovery) → gitLMS-GEO.R (Validation) → gitLMS-scRNA.R (Single-cell). This ensures that all cross-cohort data frames are correctly generated prior to plotting.

## Citation

If you use these scripts, please cite the associated manuscript and the original datasets used in the analysis.
