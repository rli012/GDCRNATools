##' @title Match samples in metadata and expression matrix
##' @description Check if samples in the metadata and expression data match
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @return A logical value. If \code{TRUE}, all the samples matched
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' genes <- c('ENSG00000000938','ENSG00000000971','ENSG00000001036',
##'         'ENSG00000001084','ENSG00000001167','ENSG00000001460')
##' 
##' samples <- c('TCGA-2F-A9KO-01', 'TCGA-2F-A9KP-01',
##'             'TCGA-2F-A9KQ-01', 'TCGA-2F-A9KR-01',
##'             'TCGA-2F-A9KT-01', 'TCGA-2F-A9KW-01')
##'              
##' metaMatrix <- data.frame(sample_type=rep('PrimaryTumor',6),
##'                         sample=samples,
##'                         days_to_death=seq(100,600,100),
##'                         days_to_last_follow_up=rep(NA,6))
##' rnaExpr <- matrix(c(2.7,7.0,4.9,6.9,4.6,2.5,
##'                     0.5,2.5,5.7,6.5,4.9,3.8,
##'                     2.1,2.9,5.9,5.7,4.5,3.5,
##'                     2.7,5.9,4.5,5.8,5.2,3.0,
##'                     2.5,2.2,5.3,4.4,4.4,2.9,
##'                     2.4,3.8,6.2,3.8,3.8,4.2),6,6)
##' rownames(rnaExpr) <- genes
##' colnames(rnaExpr) <- samples
##' gdcMatchSamples(metadata=metaMatrix, rna.expr=rnaExpr)
gdcMatchSamples <- function(metadata, rna.expr) {
    correctNum <- sum(metadata$sample==colnames(rna.expr))

    if (correctNum == nrow(metadata) & correctNum == ncol(rna.expr)) {
        ('Samples matched')
    }
}
