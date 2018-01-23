##' @title TMM normalization and voom transformation
##' @description 
##'     Normalize raw counts data by TMM implemented in \pkg{edgeR} 
##'     and then transform it by \code{\link[limma]{voom}} in \pkg{limma}
##' @param counts raw counts of RNA/miRNA expression data
##' @param filter logical, whether to filter out low-expression genes. 
##'     If \code{TRUE}, only genes with \code{cpm > 1} in more than half 
##'     of the samples will be kept. Default is \code{TRUE}
##' @return A dataframe or numeric matrix of TMM normalized and 
##'     \code{\link[limma]{voom}} transformed expression values on 
##'     the log2 scale
##' @references
##'     Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package 
##'     for differential expression analysis of digital gene expression data. 
##'     Bioinformatics. 2010 Jan 1;26(1):139-40. \cr
##'     Law CW, Chen Y, Shi W, Smyth GK. Voom: precision weights unlock 
##'     linear model analysis tools for RNA-seq read counts. Genome biology. 
##'     2014 Feb 3;15(2):R29.
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### Normalization #######
##' rnaMatrix <- matrix(sample(1:100,100), 4, 25)
##' rnaExpr <- gdcVoomNormalization(counts=rnaMatrix, filter=FALSE)

gdcVoomNormalization <- function(counts, filter=TRUE) {
    expr = DGEList(counts = counts)
    expr = calcNormFactors(expr)
    if (filter==TRUE) {
        ## filter out low expression genes
        keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(counts)
        nGenes <- as.numeric(summary(keepALL)[2]) + 
            as.numeric(summary(keepALL)[3])
        nKeep <- summary(keepALL)[3]
        
        cat (paste('Total Number of genes: ', nGenes, '\n', sep=''))
        cat (paste('Number of genes for downstream analysis: ', nKeep, 
            '\n', sep=''))
        
        exprALL <- expr[keepALL,,keep.lib.sizes = TRUE]
        v <- voom(exprALL, design=NULL, plot = FALSE)$E
        
    } else if (filter==FALSE) {
        v <- voom(expr, design=NULL, plot = FALSE)$E
    }
    return (v)
}
