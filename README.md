# *GDCRNATools* - an R/Bioconductor package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC


* **The [GDCRNATools Manual](http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html) and R code of the [GDCRNATools Workflow](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R) has been updated in 10-30-2018.**

* **If you use `GDCRNATools` in your published research, please cite:**  
Li, R., Qu, H., Wang, S., Wei, J., Zhang, L., Ma, R., Lu, J., Zhu, J., Zhong, W., and Jia, Z. (2018). GDCRNATools: an R/Bioconductor package for integrative analysis of lncRNA, miRNA and mRNA data in GDC. Bioinformatics 34, 2515-2517. https://doi.org/10.1093/bioinformatics/bty124.

* **Please add my WeChat: li-rui-dong or email to rli012@ucr.edu if you have further questions.**

***

## 1. Introduction

The [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) maintains standardized genomic, clinical, and biospecimen data from National Cancer Institute (NCI) programs including [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov/) and [Therapeutically Applicable Research To Generate Effective Treatments (TARGET)](https://ocg.cancer.gov/programs/target), It also accepts high quality datasets from non-NCI supported cancer research programs, such as genomic data from the [Foundation Medicine](https://www.foundationmedicine.com/).

`GDCRNATools` is an R/Bioconductor package which provides a standard, easy-to-use and comprehensive pipeline for downloading, organizing, and integrative analyzing RNA expression data in the GDC portal with an emphasis on deciphering the lncRNA-mRNA related ceRNA regulatory network in cancer.



## 2. Manual and R script (updated in 10-30-2018)
The comprehensive manual of `GDCRNATools` is available here: [GDCRNATools Manual](http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html)

R code of the workflow is available here: [GDCRNATools Workflow](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R)



## 3. Installation
### 3.1 Installation via Bioconductor
* The stable release version of [`GDCRNATools`](https://bioconductor.org/packages/release/bioc/html/GDCRNATools.html) requires R(>=3.5.0) and Bioconductor(>=3.8). Please start R and enter:

```R
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GDCRNATools")
```

* To install the development version of [`GDCRNATools`](https://bioconductor.org/packages/devel/bioc/html/GDCRNATools.html), please update your R and Biocondutor to the latest version and run:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GDCRNATools", version = "devel")
```

### 3.2 Installation locally
Please download the compressed package here: [GDCRNATools_1.1.5.tar.gz](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools_1.1.5.tar.gz)


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
install.packages('GDCRNATools_1.1.5.tar.gz', repos = NULL, type='source')

```

#### 3.2.2 On Linux and Mac systems
Just run the following code in R
```R
install.packages('GDCRNATools_1.1.5.tar.gz', repos = NULL, type='source')

```

### 3.3 Note
If `GDCRNATools` cannot be installed due to the lack of dependencies, please run the following code ahead to install those pacakges either simutaneously or separately:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

### install packages simutaneously ###
BiocManager::install(c('limma', 'edgeR', 'DESeq2', 'clusterProfiler', 'DOSE', 'org.Hs.eg.db', 'biomaRt', 'BiocParallel', 'GenomicDataCommons'))
install.packages(c('shiny', 'jsonlite', 'rjson', 'survival', 'survminer', 'ggplot2', 'gplots', 'Hmisc', 'DT', 'matrixStats', 'xml2'))

### install packages seperately ###
BiocManager::install('limma')
BiocManager::install('edgeR')
BiocManager::install('DESeq2')
BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
BiocManager::install('org.Hs.eg.db')
BiocManager::install('biomaRt')
BiocManager::install('BiocParallel')
BiocManager::install('GenomicDataCommons')

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


## 4. Frequently Asked Questions.
**Q1: gdcRNADownload() function doesn't work with the following error:**  
Error in FUN(X[[i]], ...):  
    unused arguments(desination_dir=directory, overwrite=TRUE)  
**A1:** This error occurs when the default API method for downloading fails. Please add   `method='gdc-client'` to the gdcRNADownload() function.

```R
####### Download RNAseq data #######
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client', ### use 'gdc-client' to download data
               directory      = rnadir)
```

**Q2: gdcRNAMerge() doesn't work with the following error:**  
Error in open.connection(file, 'rt'): cannot open the connection  
In addition: Warning message:  
In open.connection(file, 'rt'):  
    cannot open compressed file 'TCGA-XXXX/RNAseq/xxx-xxx-xxx-xxx.htseq.counts.gz', probable reason 'No such file or directory'.  
**A2:** This is usually because the data for different samples are downloaded in separate   folders. Please add `organized=FALSE` to the gdcRNAMerge() function.

```R
####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir,
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')
```
