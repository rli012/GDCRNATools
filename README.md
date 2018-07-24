# *GDCRNATools* - an R/Bioconductor package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC


## 1. Introduction

The [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) maintains standardized genomic, clinical, and biospecimen data from National Cancer Institute (NCI) programs including [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/) and [Therapeutically Applicable Research To Generate Effective Treatments (TARGET)](https://ocg.cancer.gov/programs/target), It also accepts high quality datasets from non-NCI supported cancer research programs, such as genomic data from the [Foundation Medicine](https://www.foundationmedicine.com/).

`GDCRNATools` is an R/Bioconductor package which provides a standard, easy-to-use and comprehensive pipeline for downloading, organizing, and integrative analyzing RNA expression data in the GDC portal with an emphasis on deciphering the lncRNA-mRNA related ceRNA regulatory network in cancer.



## 2. Manual and R script
The comprehensive manual of `GDCRNATools` is available here: [GDCRNATools Manual](http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html)

R code of the workflow is available here: [GDCRNATools Workflow](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R)



## 3. Installation
### 3.1 Installation via Bioconductor
*[`GDCRNATools`](http://bioconductor.org/packages/devel/bioc/html/GDCRNATools.html) requires R(>=3.5.0) and Bioconductor(>=3.7)*. To install this package, start R and enter:

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GDCRNATools")
```

### 3.2 Installation locally
Please download the compressed package here: [GDCRNATools_1.1.3.tar.gz](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_1.1.3.tar.gz)


#### 3.2.1 On Windows system
* Make sure that your R is installed in 'c:\program files'  
* Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in 'c:\program files'  
* Add R and Rtools to the Path Variable on the Environment Variables panel, including  

  c:\program files\Rtools\bin

  c:\program files\Rtools\gcc-4.6.3\bin

  c:\program files\R\R.3.x.x\bin\i386

  c:\program files\R\R.3.x.x\bin\x64 

* Run the following code in R
```R
install.packages('GDCRNATools_1.1.3.tar.gz', repos = NULL, type='source')
```

#### 3.2.2 On Linux and Mac systems
Just run the following code in R
```R
install.packages('GDCRNATools_1.1.3.tar.gz', repos = NULL, type='source')
```

### 3.3 Note
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
