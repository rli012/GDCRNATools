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
* ADD R and Rtools to the Path Variable on the Environment Variables panel, including

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
A simply manual of `GDCRNATools` is available below. Users are also highly recommended to download the comprhensive manual in _.html_ format and view on local computer [GDCRNATools Manual](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.html)





<a href="https://www.codecogs.com/eqnedit.php?latex=p=1-\sum_{k=0}^m&space;\frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p=1-\sum_{k=0}^m&space;\frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}" title="p=1-\sum_{k=0}^m \frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}" /></a>


here <a href="https://www.codecogs.com/eqnedit.php?latex=m" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m" title="m" /></a> is the number of shared miRNAs, <a href="https://www.codecogs.com/eqnedit.php?latex=N" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N" title="N" /></a> is the total number of miRNAs in the database, <a href="https://www.codecogs.com/eqnedit.php?latex=n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n" title="n" /></a> is the number of miRNAs targeting the lncRNA, <a href="https://www.codecogs.com/eqnedit.php?latex=K" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K" title="K" /></a> is the number of miRNAs targeting the protein coding gene.



<a href="https://www.codecogs.com/eqnedit.php?latex=$$Regulation\&space;similarity\&space;score&space;=&space;1-\frac{1}{M}&space;\sum_{k=1}^M&space;[{\frac{|corr(m_k,l)-corr(m_k,g)|}{|corr(m_k,l)|&plus;|corr(m_k,g)|}}]^M$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$Regulation\&space;similarity\&space;score&space;=&space;1-\frac{1}{M}&space;\sum_{k=1}^M&space;[{\frac{|corr(m_k,l)-corr(m_k,g)|}{|corr(m_k,l)|&plus;|corr(m_k,g)|}}]^M$$" title="$$Regulation\ similarity\ score = 1-\frac{1}{M} \sum_{k=1}^M [{\frac{|corr(m_k,l)-corr(m_k,g)|}{|corr(m_k,l)|+|corr(m_k,g)|}}]^M$$" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=m" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m" title="m" /></a> is the total number of shared miRNAs, <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a> is the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>th shared miRNAs, <a href="https://www.codecogs.com/eqnedit.php?latex=corr(m_k,&space;l)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?corr(m_k,&space;l)" title="corr(m_k, l)" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=corr(m_k,&space;g)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?corr(m_k,&space;g)" title="corr(m_k, g)" /></a> represents the Pearson correlation between the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>th miRNA and lncRNA, the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>th miRNA and mRNA, respectively




<a href="https://www.codecogs.com/eqnedit.php?latex=Sensitivity\&space;correlation&space;=&space;corr(l,g)-\frac{1}{M}\sum_{k=1}^M&space;{\frac{corr(l,g)-corr(m_k,l)corr(m_k,g)}{\sqrt{1-corr(m_k,l)^2}\sqrt{1-corr(m_k,g)^2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Sensitivity\&space;correlation&space;=&space;corr(l,g)-\frac{1}{M}\sum_{k=1}^M&space;{\frac{corr(l,g)-corr(m_k,l)corr(m_k,g)}{\sqrt{1-corr(m_k,l)^2}\sqrt{1-corr(m_k,g)^2}}}" title="Sensitivity\ correlation = corr(l,g)-\frac{1}{M}\sum_{k=1}^M {\frac{corr(l,g)-corr(m_k,l)corr(m_k,g)}{\sqrt{1-corr(m_k,l)^2}\sqrt{1-corr(m_k,g)^2}}}" /></a>


where <a href="https://www.codecogs.com/eqnedit.php?latex=M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?M" title="M" /></a> is the total number of shared miRNAs, <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a> is the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>th shared miRNAs, <a href="https://www.codecogs.com/eqnedit.php?latex=corr(l,&space;g)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?corr(l,&space;g)" title="corr(l, g)" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=corr(m_k,&space;l)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?corr(m_k,&space;l)" title="corr(m_k, l)" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=corr(m_k,&space;g)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?corr(m_k,&space;g)" title="corr(m_k, g)" /></a> represents the Pearson correlation between the long non-coding RNA and the protein coding gene, the kth miRNA and lncRNA, the kth miRNA and mRNA, respectively



![](vignettes/figures/TCGA-PRAD.shinyCorPlot.gif)



![](vignettes/figures/network.png)


![](vignettes/figures/TCGA-PRAD.shinyKMPlot.gif)


![](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/figures/gdcEnrichPlot.GO.bar.png)

![](https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/figures/gdcEnrichPlot.GO.bubble.png)

![](vignettes/figures/TCGA-PRAD.shinyPathview.gif)


## sessionInfo
```{r sessionInfo}
sessionInfo()
```

## References
Chou, Chih-Hung, Sirjana Shrestha, Chi-Dung Yang, Nai-Wen Chang, Yu-Ling Lin, Kuang-Wen Liao, Wei-Chi Huang, et al. 2017. “MiRTarBase Update 2018: A Resource for Experimentally Validated MicroRNA-Target Interactions.” Nucleic Acids Research, November, gkx1067–gkx1067. doi:10.1093/nar/gkx1067.

Colaprico, Antonio, Tiago C. Silva, Catharina Olsen, Luciano Garofano, Claudia Cava, Davide Garolini, Thais S. Sabedot, et al. 2016. “TCGAbiolinks: An R/Bioconductor Package for Integrative Analysis of TCGA Data.” Nucleic Acids Research 44 (8): e71. doi:10.1093/nar/gkv1507.

Furi’o-Tar’i, Pedro, Sonia Tarazona, Toni Gabald’on, Anton J. Enright, and Ana Conesa. 2016. “SpongeScan: A Web for Detecting MicroRNA Binding Elements in LncRNA Sequences.” Nucleic Acids Research 44 (Web Server issue): W176–W180. doi:10.1093/nar/gkw443.

Jeggari, Ashwini, Debora S Marks, and Erik Larsson. 2012. “MiRcode: A Map of Putative MicroRNA Target Sites in the Long Non-Coding Transcriptome.” Bioinformatics 28 (15): 2062–3. doi:10.1093/bioinformatics/bts344.

Langfelder, Peter, and Steve Horvath. 2008. “WGCNA: An R Package for Weighted Correlation Network Analysis.” BMC Bioinformatics 9 (December): 559. doi:10.1186/1471-2105-9-559.

Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. 2014. “Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read Counts.” Genome Biology 15 (February): R29. doi:10.1186/gb-2014-15-2-r29.

Li, Jun-Hao, Shun Liu, Hui Zhou, Liang-Hu Qu, and Jian-Hua Yang. 2014. “StarBase V2.0: Decoding MiRNA-CeRNA, MiRNA-NcRNA and Protein–RNA Interaction Networks from Large-Scale CLIP-Seq Data.” Nucleic Acids Research 42 (Database issue): D92–D97. doi:10.1093/nar/gkt1248.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (December): 550. doi:10.1186/s13059-014-0550-8.

Luo, Weijun, and Cory Brouwer. 2013. “Pathview: An R/Bioconductor Package for Pathway-Based Data Integration and Visualization.” Bioinformatics 29 (14): 1830–1. doi:10.1093/bioinformatics/btt285.

Paci, Paola, Teresa Colombo, and Lorenzo Farina. 2014. “Computational Analysis Identifies a Sponge Interaction Network Between Long Non-Coding RNAs and Messenger RNAs in Human Breast Cancer.” BMC Systems Biology 8 (July): 83. doi:10.1186/1752-0509-8-83.

Ritchie, Matthew E., Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, and Gordon K. Smyth. 2015. “Limma Powers Differential Expression Analyses for RNA-Sequencing and Microarray Studies.” Nucleic Acids Research 43 (7): e47. doi:10.1093/nar/gkv007.

Robinson, Mark D., and Alicia Oshlack. 2010. “A Scaling Normalization Method for Differential Expression Analysis of RNA-Seq Data.” Genome Biology 11 (March): R25. doi:10.1186/gb-2010-11-3-r25.

Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. 2010. “EdgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” Bioinformatics 26 (1): 139–40. doi:10.1093/bioinformatics/btp616.

Yu, Guangchuang, Li-Gen Wang, Yanyan Han, and Qing-Yu He. 2012. “ClusterProfiler: An R Package for Comparing Biological Themes Among Gene Clusters.” OMICS : A Journal of Integrative Biology 16 (5): 284–87. doi:10.1089/omi.2011.0118.

Yu, Guangchuang, Li-Gen Wang, Guang-Rong Yan, and Qing-Yu He. 2015. “DOSE: An R/Bioconductor Package for Disease Ontology Semantic and Enrichment Analysis.” Bioinformatics 31 (4): 608–9. doi:10.1093/bioinformatics/btu684.

Zhu, Yitan, Peng Qiu, and Yuan Ji. 2014. “TCGA-Assembler: An Open-Source Pipeline for TCGA Data Downloading, Assembling, and Processing.” Nature Methods 11 (6): 599–600. doi:10.1038/nmeth.2956.
