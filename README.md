# Differential Gene Expression Analysis: Lactation Project (GSE60450)

This repository contains an R script (`Actividad_2.R`) designed for the statistical processing and analysis of RNA-seq data. The study compares transcriptomic profiles in mouse mammary cells across different physiological states (Lactation vs. Pregnancy).

## 📊 Pipeline Overview

The script implements a standard Bioinformatics workflow using the Bioconductor suite, specifically the **edgeR** package:

1.  **Data Loading & Cleaning:** Importing raw count matrices and sample metadata.
2.  **Exploratory Data Analysis (EDA):** Calculating maximum values per sample and the percentage of expressed genes.
3.  **Genomic Annotation:** Mapping Entrez IDs to official Gene Symbols using `org.Mm.eg.db`.
4.  **Preprocessing:**
    * Creation of the `DGEList` object.
    * Filtering low-expressed genes using `filterByExpr`.
    * Library size normalization using the TMM method (`calcNormFactors`).
5.  **Quality Control Visualization:**
    * LogCPM boxplots (pre and post-normalization).
    * Multi-dimensional Scaling (MDS) plots to assess sample clustering.
6.  **Statistical Modeling:**
    * Estimation of common, trended, and tagwise dispersion.
    * Fitting Quasi-Likelihood (QL) Generalized Linear Models (GLM).
7.  **Hypothesis Testing:**
    * Defining specific contrasts (e.g., Basal Lactate vs. Basal Pregnant).
    * Identifying Differentially Expressed Genes (DEGs) based on FDR and Log Fold Change.

## 🛠 Requirements

To run this script, you need to install the following R packages:

```r
install.packages(c("ggplot2", "statmod"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "org.Mm.eg.db"))
