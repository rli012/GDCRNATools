# *GDCRNATools* - an R/Bioconductor package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC


## 1. Introduction

The [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) maintains standardized genomic, clinical, and biospecimen data from National Cancer Institute (NCI) programs including [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/) and [Therapeutically Applicable Research To Generate Effective Treatments (TARGET)](https://ocg.cancer.gov/programs/target), It also accepts high quality datasets from non-NCI supported cancer research programs, such as genomic data from the [Foundation Medicine](https://www.foundationmedicine.com/).

`GDCRNATools` is an R package which provides a standard, easy-to-use and comprehensive pipeline for downloading, organizing, and integrative analyzing RNA expression data in the GDC portal with an emphasis on deciphering the lncRNA-mRNA related ceRNA regulatory network in cancer.



## 2. Manual and R script
The comprehensive manual of `GDCRNATools` is available here: [GDCRNATools Manual](http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html)

R code of the workflow is available here: [GDCRNATools Workflow](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R)



## 3. Installation
`GDCRNATools` has been accepted by [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/GDCRNATools.html) as the development version. However, before the next release of R/Bioconductor in mid-April, it is highly recommended to install `GDCRNATools` by downloading and installing locally, rather than via Bioconductor.

Please download the compressed package here: [GDCRNATools_0.99.16.tar.gz](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_0.99.16.tar.gz)



### On Windows system
* Make sure that your R is installed in 'c:\program files'  
* Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in 'c:\program files'  
* Add R and Rtools to the Path Variable on the Environment Variables panel, including  

  c:\program files\Rtools\bin

  c:\program files\Rtools\gcc-4.6.3\bin

  c:\program files\R\R.3.x.x\bin\i386

  c:\program files\R\R.3.x.x\bin\x64 

* Run the following code in R
```R
install.packages('GDCRNATools_0.99.16.tar.gz', repos = NULL, type='source')
```

### On Linux and Mac systems
Just run the following code in R
```R
install.packages('GDCRNATools_0.99.16.tar.gz', repos = NULL, type='source')
```

### Note
If `GDCRNATools` cannot be installed due to the lack of dependencies, please run the following code ahead to install those pacakges either simutaneously or separately:
```R
source("https://bioconductor.org/biocLite.R")

### install packages simutaneously ###
biocLite(c('limma', 'edgeR', 'DESeq2', 'clusterProfiler', 'DOSE', 'org.Hs.eg.db', 'biomaRt', 'BiocParallel', 'GenomicDataCommons'))
install.packages(c('shiny', 'jsonlite', 'rjson', 'survival', 'survminer', 'ggplot2', 'gplots', 'Hmisc', 'DT', 'matrixStats', 'xml2'))

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
install.packages('matrixStats')
install.packages('xml2')
```
