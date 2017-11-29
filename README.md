# *GDCRNATools* - An R package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC


## Introduction

The [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) maintains standardized genomic, clinical, and biospecimen data from National Cancer Institute (NCI) programs including [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/) and [Therapeutically Applicable Research To Generate Effective Treatments (TARGET)](https://ocg.cancer.gov/programs/target), It also accepts high quality datasets from non-NCI supported cancer research programs, such as genomic data from the [Foundation Medicine](https://www.foundationmedicine.com/).

**GDCRNATools** is an R/Bioconductor package which provides a standard, easy-to-use and comprehensive pipeline for downloading, organizing, and integrative analyzing RNA expression data in the GDC portal with an emphasis on deciphering the lncRNA-mRNA related ceRNA regulatory network in cancer.

Many analyses can be perfomed using GDCRNATools, including differential gene expression analysis ([limma](http://bioconductor.org/packages/release/bioc/html/limma.html)[ref](1), [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)[ref](2), and [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)[ref](3)), univariate survival analysis (CoxPH and KM), competing endogenous RNA network analysis (hypergeometric test, Pearson correlation analysis, regulation similarity analysis, sensitivity Pearson partial  correlation[ref](4)), and functional enrichment analysis(GO, KEGG, DO). Besides some routine visualization methods such as volcano plot, scatter plot, and bubble plot, etc., three simple shiny apps are developed in GDCRNATools allowing users visualize the results on a local webpage. All the figures are plotted based on [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package unless otherwise specified.

This user-friendly package allows researchers perform the analysis by simply running a few functions and integrate their own pipelines such as molecular subtype classification, [weighted correlation network analysis (WGCNA)](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)[ref](5), and TF-miRNA co-regulatory network analysis, etc. into the workflow easily. This could open a door to accelerate the study of crosstalk among different classes of RNAs and their regulatory relationships in cancer.


## Installation
`GDCRNATools` can be installed from github directly using the `install_github()` function in `devtools` package. So users should have `devtools` installed before installing `GDCRNATools`
devtools::install_github(repo='Jialab-UCR/GDCRNATools')
