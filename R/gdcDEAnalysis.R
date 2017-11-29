##' @title Differential gene expression analysis
##' @description Performs differential gene expression analysis by \pkg{limma}, \pkg{edgeR}, and \pkg{DESeq2}
##' @param counts a dataframe or numeric matrix of raw counts data generated from \code{\link{gdcRNAMerge}}
##' @param group a vector giving the group that each sample belongs to
##' @param comparison a character string specifying the two groups being compared. \cr
##'   Example: \code{comparison='PrimaryTumor-SolidTissueNormal'}
##' @param method one of \code{'limma'}, \code{'edgeR'}, and \code{'DESeq2'}. Default is \code{'limma'} \cr
##'   Note: It may takes long time for \code{method='DESeq2'} with a single core
##' @param n.cores a numeric value of cores to be used for \code{method='DESeq2'} to accelate the analysis process. Default is \code{NULL}
##' @param filter logical, whether to filter out low expression genes. If \code{TRUE}, only genes 
##'   with \code{cpm > 1} in more than half of the samples will be kept. Default is \code{TRUE}
##' @import edgeR
##' @importFrom limma makeContrasts
##' @importFrom limma voom
##' @importFrom limma lmFit
##' @importFrom limma contrasts.fit
##' @importFrom limma eBayes
##' @importFrom limma topTable
##' @importFrom DESeq2 DESeqDataSetFromMatrix
##' @importFrom DESeq2 DESeq
##' @importFrom DESeq2 results
##' @importFrom BiocParallel MulticoreParam
##' @importFrom BiocParallel register
##' @return A dataframe containing Ensembl gene ids/miRBase v21 mature miRNA ids, gene symbols, biotypes, fold change on the log2 scale,
##'   p value, and FDR etc. of all genes/miRNAs of analysis.
##' @note It may takes long time for \code{method='DESeq2'} with a single core. Please use multiple cores if possible
##' @references
##'   Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010 Jan 1;26(1):139-40. \cr
##'   Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research. 2015 Jan 20;43(7):e47-e47. \cr
##'   Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology. 2014 Dec 5;15(12):550.
##' @export
##' @author Ruidong Li and Han Qu
gdcDEAnalysis <- function(counts, group, comparison, method='limma', n.cores=NULL, filter=TRUE) {
  dge = DGEList(counts = counts)
  keep <- rowSums(cpm(dge) > 1) >= 0.5*length(group)
  
  if (method == 'DESeq2') {
    
    cat ('DE analysis using DESeq2 may take long time with a single core')
    
    coldata <- data.frame(group)
    dds <- DESeqDataSetFromMatrix(countData = readCounts,
                                  colData = coldata,
                                  design = ~ group)
    
    dds$group <- factor(dds$group, levels = rev(strsplit(comparison, '-', fixed=T)[[1]]))
    dds <- dds[keep, ]
    
    if (! is.null(n.cores)) {
      register(MulticoreParam(n.cores))
      dds <- DESeq(dds, parallel=TRUE)
    } else {
      dds <- DESeq(dds)
    }
    
    res <- results(dds)
    #res <- lfcShrink(dds, coef=2, res=res)
    
    DEGAll <- data.frame(res)
    colnames(DEGAll) <- c('baseMean', 'logFC', 'lfcSE', 'stat', 'PValue', 'FDR')
    
  } else if (method %in% c('edgeR', 'limma')) {
    group <- factor(group)
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    contrast.matrix <- makeContrasts(contrasts=comparison, levels=design)
    
    dge <- dge[keep,,keep.lib.sizes = TRUE]
    dge <- calcNormFactors(dge)
    
    if (method == 'edgeR') {
      dge <- estimateDisp(dge, design)
      fit <- glmFit(dge, design)
      lrt <- glmLRT(fit, contrast=contrast.matrix)
      
      DEGAll <- lrt$table
      
    } else if (method == 'limma') {
      v <- voom(dge, design=design, plot = FALSE)
      
      fit <- lmFit(v, design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      DEGAll <- topTable(fit2, coef=1, n = Inf)
      colnames(DEGAll) <- c('logFC', 'AveExpr', 't', 'PValue', 'FDR', 'B')
    }
  }
  
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = 'fdr')
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  
  if (startsWith(rownames(counts)[1], 'ENSG')) {
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, group=degList$group, DEGAll)
    
    keep <- which(! is.na(degOutput$symbol))
    degOutput <- degOutput[keep,]
    return(degOutput)
  } else {
    return (DEGAll)
  }
}


###
deAnalysislimma <- function(v, design, contrast.matrix, type='RNAseq') {
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  DEGAll <- topTable(fit2, coef=1, n = Inf)
  
  o <- order(DEGAll$adj.P.Val)
  DEGAll <- DEGAll[o,]
  
  if (type=='miRNAs') {
    return (DEGAll)
  } else {
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, group=degList$group, DEGAll)
    
    keep <- which(! is.na(degOutput$symbol))
    
    degOutput <- degOutput[keep,]
    return(degOutput)
  }
}


##' @title Report differentially expressed genes/miRNAs
##' @description Report genes/miRNAs that are differentially expressed satisfying a given threshold
##' @param deg A dataframe of DE analysis result from \code{\link{gdcDEAnalysis}}
##' @param gene.type one of \code{'all'}, \code{'long_non_coding'}, \code{'protein_coding'}, and \code{'miRNAs'}. Default is \code{'all'}
##' @param fc a numeric value specifying the threshold of fold change
##' @param pval a nuemric value specifying the threshold of p value
##' @return A dataframe or numeric matrix of differentially expressed genes/miRNAs
##' @export
##' @author Ruidong Li and Han Qu
gdcDEReport <- function(deg, gene.type='all', fc=2, pval=0.01) {
  sig <- abs(deg$logFC)>log(fc,2) & deg$FDR<pval
  
  degFinal <- deg[sig,]
  
  if (gene.type=='all' | gene.type=='miRNAs') {
    return (degFinal)
  } else {
    keep <- degFinal$group == gene.type
    degFinal <- degFinal[keep,]
    return (degFinal)
  }
}