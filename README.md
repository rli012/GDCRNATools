# *GDCRNATools* - An R/Bioconductor package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC


# Introduction

The [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) maintains standardized genomic, clinical, and biospecimen data from National Cancer Institute (NCI) programs including [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/) and [Therapeutically Applicable Research To Generate Effective Treatments (TARGET)](https://ocg.cancer.gov/programs/target), It also accepts high quality datasets from non-NCI supported cancer research programs, such as genomic data from the [Foundation Medicine](https://www.foundationmedicine.com/).

`GDCRNATools` is an R package which provides a standard, easy-to-use and comprehensive pipeline for downloading, organizing, and integrative analyzing RNA expression data in the GDC portal with an emphasis on deciphering the lncRNA-mRNA related ceRNA regulatory network in cancer.

Many analyses can be perfomed using `GDCRNATools`, including differential gene expression analysis ([limma](http://bioconductor.org/packages/release/bioc/html/limma.html)(Ritchie et al. 2015), [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)(Robinson, McCarthy, and Smyth 2010), and [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)(Love, Huber, and Anders 2014)), univariate survival analysis (CoxPH and KM), competing endogenous RNA network analysis (hypergeometric test, Pearson correlation analysis, regulation similarity analysis, sensitivity Pearson partial  correlation), and functional enrichment analysis(GO, KEGG, DO). Besides some routine visualization methods such as volcano plot, scatter plot, and bubble plot, etc., three simple shiny apps are developed in GDCRNATools allowing users visualize the results on a local webpage.


This user-friendly package allows researchers perform the analysis by simply running a few functions and integrate their own pipelines such as molecular subtype classification, [weighted correlation network analysis (WGCNA)](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)(Langfelder and Horvath 2008), and TF-miRNA co-regulatory network analysis, etc. into the workflow easily.


# Installation
`GDCRNATools` is now under review in Bioconductor. Users can install the package locally.

## On Windows system
* Download the package [GDCRNATools_0.99.5.tar.gz](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_0.99.5.tar.gz)
* Make sure you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed
* Add R and Rtools to the Path Variable on the Environment Variables panel, including

  c:\program files\Rtools\bin

  c:\program files\Rtools\gcc-4.6.3\bin

  c:\program files\R\R.3.x.x\bin\i386

  c:\program files\R\R.3.x.x\bin\x64 

* Open a command prompt. Type R CMD INSTALL GDCRNATools_0.99.5.tar.gz

## On Linux and Mac systems
Run the following command in R
```R
install.packages('GDCRNATools_0.99.5.tar.gz', repos = NULL, type='source')
```

If `GDCRNATools` cannot be installed due to the lack of dependencies, please run the following code ahead to install those pacakges simutaneously or separately:
```R
source("https://bioconductor.org/biocLite.R")

### install packages simutaneously ###
biocLite(c('limma', 'edgeR', 'DESeq2', 'clusterProfiler', 'DOSE', 'org.Hs.eg.db', 'biomaRt', 'BiocParallel','GenomicDataCommons'))
install.packages(c('shiny', 'jsonlite', 'rjson', 'survival', 'survminer', 'ggplot2', 'gplots', 'Hmisc', 'DT'))

### install packages seperately ###
biocLite('limma')
biocLite('edgeR')
biocLite('DESeq2')
biocLite('clusterProfiler')
biocLite('DOSE')
biocLite('org.Hs.eg.db')
biocLite('biomaRt')
biocLite('BiocParallel')
biocLite('GenomicDataCommons')

install.packages('shiny')
install.packages('jsonlite')
install.packages('rjson')
install.packages('survival')
install.packages('survminer')
install.packages('ggplot2')
install.packages('gplots')
install.packages('Hmisc')
install.packages('DT')
```


# Manual
The comprehensive manual of `GDCRNATools` is available here: [GDCRNATools Manual](http://htmlpreview.github.io/?https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_manual.html)
